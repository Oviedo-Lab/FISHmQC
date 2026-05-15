
// Rcpp
// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <random>
#include <nlopt.hpp>
#include <cmath>
using namespace Rcpp;
using namespace Eigen;

struct Codebook {
  int N_bits;
  std::vector<uint64_t> barcodes;
  std::vector<std::string> species;
  std::vector<int> blanks;
};

struct ST_data {
  std::vector<int> bc_counts; // Vector of length N_barcodes, giving total counts for each barcode across all cells
  Codebook cb;
  std::unordered_map<int, std::vector<uint64_t>> correction_table_inverted;
  Eigen::Matrix<uint64_t, Dynamic, Dynamic> transform_flips_all;
};

struct FlipRates {
  std::vector<double> rate10; // P(1 -> 0) for each bit, so length = N_bits
  std::vector<double> rate01; // P(0 -> 1) for each bit, so length = N_bits
  MatrixXd corr; // N_bits x N_bits matrix of luminance noise correlations between bits
}; 

struct SpotSim {
  std::vector<uint64_t> barcodes_true;
  std::vector<uint64_t> barcodes_read; 
  std::vector<uint64_t> barcodes_corrected;
  std::vector<int> labels_true;
  std::vector<int> labels_read;
  std::vector<int> labels_corrected;
  std::vector<int> cell_ids;
  std::vector<VectorXd> lum;
  std::vector<double> rate01_pre;
  std::vector<double> rate10_pre;
  std::vector<double> rate01_post;
  std::vector<double> rate10_post;
  std::vector<int> read_counts;
  std::vector<int> corrected_counts;
  std::vector<int> true_counts;
  std::vector<double> CR;
  std::vector<double> PPV;
};

MatrixXd rmvnorm(
    int n,
    const VectorXd& mu,
    const MatrixXd& Sigma
  ) {
    int d = mu.size();
    
    // Cholesky decomposition
    Eigen::LLT<MatrixXd> llt(Sigma);
    MatrixXd L = llt.matrixL();
    
    // RNG
    std::mt19937 rng(std::random_device{}());
    std::normal_distribution<double> norm(0.0, 1.0);
    
    // Standard normal samples
    MatrixXd Z(n, d);
    
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < d; ++j) {
        Z(i, j) = norm(rng);
      }
    }
    
    // Correlate + shift mean
    MatrixXd X =
      (Z * L.transpose()).rowwise() +
      mu.transpose();
    
    return X;
  }

std::vector<int> grep_idx(
    Rcpp::CharacterVector x,
    std::string s
  ) {
    std::vector<int> idx;
    for (int i = 0; i < x.size(); ++i) {
      std::string xi = Rcpp::as<std::string>(x[i]);
      if (xi.find(s) != std::string::npos) {idx.push_back(i);}
    }
    return idx;
  } 


/*
 * Functions to analytically compute expected count read from flip rates and true counts
 */

// Hold and flip rates by bit
double H(int bit, double rate10, double rate01) {return 1.0 - rate01 + bit * (rate01 - rate10);}
double F(int bit, double rate10, double rate01) {return rate01 + bit * (rate10 - rate01);}
double TR(
    uint64_t bc,
    uint64_t transform_flips,
    std::vector<double> rate10,
    std::vector<double> rate01
  ) {
    double transform_rate = 1.0;
    for (int b = 0; b < rate10.size(); ++b) {
      int bit = (bc >> b) & 1ULL;
      int flip = (transform_flips >> b) & 1ULL;
      double h = H(bit, rate10[b], rate01[b]);
      transform_rate *= h + flip * (F(bit, rate10[b], rate01[b]) - h);
    }
    return transform_rate; 
  }

Eigen::Matrix<uint64_t, Dynamic, Dynamic> find_transform_flips_all(
    std::vector<uint64_t> barcodes,
    int N_bits
  ) {
    int N_barcodes = barcodes.size();
    // Initialize matrix to hold transforming flips
    // ... T(i,j) = bit flips needed to get from barcode i to barcode j, and T(i,j) = T(j,i)
    Eigen::Matrix<uint64_t, Dynamic, Dynamic> T(N_barcodes, N_barcodes);
    for (int j = 0; j < N_barcodes; ++j) {
      for (int i = j; i < N_barcodes; ++i) {
        T(i, j) = (barcodes[i] ^ barcodes[j]) & ((1ULL << N_bits) - 1);
      }
    }
    return T;
  }

// Number of B correctly read
double expected_preserved_spots(
    uint64_t bc, // True barcode, as a packed uint64_t
    const std::vector<double>& rate10, // bit-flip rates, P(1 -> 0) for each bit, so length = N_bits
    const std::vector<double>& rate01, // bit-flip rates, P(0 -> 1) for each bit, so length = N_bits
    int bc_count_true // Number of spots correctly decoded as bc
  ) {
    double holding_rate = 1.0;
    for (int b = 0; b < rate10.size(); ++b) {
      int bit = (bc >> b) & 1ULL;
      // Assume bit-flips are independent, so multiply hold/flip rates across bits
      holding_rate *= H(bit, rate10[b], rate01[b]);
    }
    return holding_rate * (double)bc_count_true;
  }

// Number of B' incorrectly read as B
double expected_misreads(
    int k, // Index of true barcode
    const std::vector<double>& rate10, // bit-flip rates, P(1 -> 0) for each bit, so length = N_bits
    const std::vector<double>& rate01, // bit-flip rates, P(0 -> 1) for each bit, so length = N_bits
    const std::vector<int>& bc_counts_true, // Vector of same length as true_barcodes, giving number of spots with each true barcode
    const std::vector<uint64_t>& true_barcodes, // Vector giving all possible true spot barcodes
    const Eigen::Matrix<uint64_t, Dynamic, Dynamic>& transform_flips_all
  ) {
    double misread_count = 0.0;
    for (int i = 0; i < true_barcodes.size(); ++i) {
      // Skip true barcode, since we want misreads
      if (i == k) {continue;} 
      // Get bar code of interest
      uint64_t bct = true_barcodes[i];
      // Get flips needed to transform bct into bc 
      uint64_t flips = transform_flips_all(i, k);
      // Find transform rate, assuming independent bit flips, so multiply across bits
      double transform_rate = TR(bct, flips, rate10, rate01);
      misread_count += transform_rate * (double)bc_counts_true[i];
    }
    return misread_count;
  }


std::vector<uint64_t> find_transform_flips(
    std::vector<uint64_t> barcodes, // Possible true barcodes (e.g., from codebook)
    uint64_t observed_barcode, // A single observed barcode (e.g., from a spot)
    int N_bits
  ) {
    int N_barcodes = barcodes.size();
    std::vector<uint64_t>  T;
    T.reserve(N_barcodes);
    for (int i = 0; i < N_barcodes; ++i) {
      uint64_t flips = (barcodes[i] ^ observed_barcode) & ((1ULL << N_bits) - 1);
      T.push_back(flips);
    }
    return T; // For each possible true barcode, the misread flips which would turn it into observed_barcode
  }
  

// Number of B' incorrectly read as something corrected to B
double expected_miscorrections(
    int k,                                        // Index of true barcode
    const std::vector<double>& rate10,            // bit-flip rates, P(1 -> 0) for each bit, so length = N_bits
    const std::vector<double>& rate01,            // bit-flip rates, P(0 -> 1) for each bit, so length = N_bits
    const std::vector<int>& bc_counts_true,       // Vector of same length as true_barcodes, giving number of spots with each true barcode
    const std::vector<uint64_t>& true_barcodes,   // Vector giving all possible true spot barcodes
    const std::vector<uint64_t>& corrected_to_k   // vector of barcodes that would be corrected to barcode indexed by k
  ) {
    int N_bits = rate10.size();
    double misread_count = 0.0;
    // For each possible barcode misread that would be corrected to the barcode indexed by k ...
    for (uint64_t bcr : corrected_to_k) {
      // Get bit-flips required to transform each true barcode into this misread barcode
      std::vector<uint64_t> transform_flips = find_transform_flips(true_barcodes, bcr, N_bits);
      for (int i = 0; i < true_barcodes.size(); ++i) {
        if (i == k) {continue;} // Skip true barcode, since we want misreads
        // Get bar code of interest
        uint64_t bct = true_barcodes[i];
        // Get flips needed to transform bct into bcr 
        uint64_t flips = transform_flips[i];
        // Find transform rate, assuming independent bit flips, so multiply across bits
        double transform_rate = TR(bct, flips, rate10, rate01);
        // Find expected count of spot misreads corrected to barcode k, from true barcode i, and add to total misread count
        misread_count += transform_rate * (double)bc_counts_true[i];
      }
    }
    return misread_count;
  }

// Estimate expected barcode counts after correction, assuming no flip-rate correlations or over-dispersion, as a function of flip rates
std::vector<double> expected_corrected_counts(
    const std::vector<double>& rate10,                                               // bit-flip rates, P(1 -> 0) for each bit, so length = N_bits
    const std::vector<double>& rate01,                                               // bit-flip rates, P(0 -> 1) for each bit, so length = N_bits
    const std::vector<int>& bc_counts_true,                                          // Vector of same length as true_barcodes, giving number of spots with each true barcode
    const std::vector<uint64_t>& true_barcodes,                                      // Vector giving all possible true spot barcodes
    const std::unordered_map<int, std::vector<uint64_t>>& correction_table_inverted, // Inverted correction table mapping each corrected barcode index to vector of misread barcodes that would be corrected to it
    const Eigen::Matrix<uint64_t, Dynamic, Dynamic>& transform_flips_all
  ) {
    int N_barcodes = true_barcodes.size();
    std::vector<double> expected_counts(N_barcodes);
    for (int k = 0; k < N_barcodes; ++k) {
      expected_counts[k] = 
        expected_preserved_spots(true_barcodes[k], rate10, rate01, bc_counts_true[k]) + 
        expected_misreads(k, rate10, rate01, bc_counts_true, true_barcodes, transform_flips_all) + 
        expected_miscorrections(k, rate10, rate01, bc_counts_true, true_barcodes, correction_table_inverted.at(k));
    }
    return expected_counts;
  }

// Idea: for each vector of actual observed corrected barcode counts bc_counts, set all blanks to zero to approximate true_barcodes, 
// then optimize rate10 and rate01 to minimize the distance between ecc = expected_corrected_counts(rate10, rate01, true_barcodes, correction_table_inverted, transform_flips_all) and bc_counts.
// Here, "distance" is (negative log) likelihood of bc_counts given ecc as the mean of (an on overdispersed??) Poisson distribution. 
// Ignore over-dispersion for now, and handle that in the DG simulations??
// Use this as initial values for the flip rates, and set correlations of flip rates as zero to start. 

double observed_counts_nll(
    const std::vector<double>& rates, 
    std::vector<double>& grad,
    void* data                    
  ) {
   
    // Extract N_bits 
    int N_bits = rates.size() / 2;
    std::vector<double> rate10(rates.begin(), rates.begin() + N_bits);
    std::vector<double> rate01(rates.begin() + N_bits, rates.end());
    
    // Extract data
    const auto* d = static_cast<ST_data*>(data);
    const auto& bc_counts = d->bc_counts;
    const auto& blanks = d->cb.blanks;
    const auto& true_barcodes = d->cb.barcodes;
    const auto& correction_table_inverted = d->correction_table_inverted;
    const auto& transform_flips_all = d->transform_flips_all;
    
    // Set blanks to zero to approximate true barcodes
    std::vector<int> bc_counts_true = bc_counts; 
    for (int idx : blanks) {
      bc_counts_true[idx] = 0;
    }
    
    // Compute expected corrected counts from these flip rates
    std::vector<double> ecc = expected_corrected_counts(
      rate10, 
      rate01, 
      bc_counts_true, 
      true_barcodes, 
      correction_table_inverted, 
      transform_flips_all
    );
    
    // Compute negative log likelihood of observed barcodes given expected corrected counts as mean of Poisson distribution
    double nll = 0.0;
    for (int i = 0; i < bc_counts.size(); ++i) {
      double lambda = ecc[i];
      lambda = std::max(lambda, 1e-12); // Avoid log(0) 
      int k = bc_counts[i];
      nll += lambda - k * std::log(lambda) + std::lgamma(k + 1); // Poisson negative log likelihood
    }
    
    // ... and return
    return nll;
  }

// Optimize rate10 and rate01 to minimize distance between expected_corrected_counts and observed_barcodes
std::pair<std::vector<double>, std::vector<double>> estimate_flip_rates_initial(
    ST_data& STdata,
    double ctol = 1e-7,
    int max_evals = 1000
  ) {
   
    // Initialize rate10 and rate01 with some reasonable starting values, e.g., 0.01 for all bits
    int N_bits = STdata.cb.N_bits; 
    std::vector<double> rates(N_bits*2, 0.01); // First N_bits are rate10, second N_bits are rate01
    size_t n = rates.size();
    
    // Set up NLopt optimizer
    nlopt::opt opt(nlopt::LN_SBPLX, n); 
    opt.set_min_objective(observed_counts_nll, &STdata);
    opt.set_ftol_rel(ctol);       // stop when iteration changes objective fn value by less than this fraction 
    opt.set_maxeval(max_evals);   // Maximum number of evaluations to try
    // ... add bounds 
    std::vector<double> lb(n, 0.0);
    std::vector<double> ub(n, 1.0);
    opt.set_lower_bounds(lb);
    opt.set_upper_bounds(ub);
    
    // Fit model
    int success_code = 0;
    double min_fx;
    try {
      nlopt::result sc = opt.optimize(rates, min_fx);
      success_code = static_cast<int>(sc);
    } catch (std::exception& e) {  
      if (false) {
        Rcpp::Rcout << "Optimization failed: " << e.what() << std::endl;
      }  
      success_code = 0;
    }  
    
    // Extract rate10 and rate01 from rates vector and return
    std::vector<double> rate10(rates.begin(), rates.begin() + N_bits);
    std::vector<double> rate01(rates.begin() + N_bits, rates.end());
    return {rate10, rate01};
  }


/*
 * Next up: remake the set_true_spot_info and generate_simulated_spots functions; 
 *  ... done. set_true_spot_info is now make_true_bc_counts, and generate_simulated_spots is now make_SpotSim.
 *  ... eliminated use of MASS? 
 * releated: \set_luminance_noise_correlation_matrix
 */

std::vector<int> decode_lum(
    const VectorXd& x
  ) {
    const int n = x.size();
    std::vector<int> y(n);
    for (int i = 0; i < n; ++i) {y[i] = (x[i] > 0.0);}
    return y;
  }

uint64_t pack_decode_lum(
    const VectorXd& lum
  ) {
    uint64_t packed = 0ULL;
    for (int i = 0; i < lum.size(); ++i) {
      packed |= (uint64_t)(lum[i] > 0.0) << i;
    }
    return packed;
  }

std::pair<std::vector<double>, std::vector<double>> compute_flip_rates(
    const std::vector<uint64_t>& barcodes_true,
    const std::vector<uint64_t>& barcodes_read,
    int N_bits
  ) {
    
    std::vector<uint64_t> n01(N_bits, 0); // 0 -> 1
    std::vector<uint64_t> n10(N_bits, 0); // 1 -> 0
    
    std::vector<uint64_t> denom0(N_bits, 0);
    std::vector<uint64_t> denom1(N_bits, 0);
    
    const size_t N_barcodes = barcodes_true.size();
    
    for (size_t i = 0; i < N_barcodes; ++i) {
      
      uint64_t t = barcodes_true[i];
      uint64_t r = barcodes_read[i];
      uint64_t diff = t ^ r;
      
      for (int b = 0; b < N_bits; ++b) {
        uint64_t mask = 1ULL << b;
        bool true_bit = t & mask;
        bool flipped  = diff & mask;
        if (true_bit) {
          ++denom1[b];
          if (flipped) {++n10[b];}
        } else {
          ++denom0[b];
          if (flipped) {++n01[b];}
        }
      }
    }
    
    std::vector<double> rate01(N_bits);
    std::vector<double> rate10(N_bits);
    
    for (int b = 0; b < N_bits; ++b) {
      rate01[b] = denom0[b] > 0 ? double(n01[b]) / denom0[b] : 0.0;
      rate10[b] = denom1[b] > 0 ? double(n10[b]) / denom1[b] : 0.0;
    }
    
    return {rate01, rate10};
  }

MatrixXi make_true_bc_counts(
    const ArrayXd& gene_kernel_rates, // 1 D vector, of length N_genes, from empirical observations
    const ArrayXd& gamma_variance, // 1 D vector, length N_genes, from empirical observations 
    const std::vector<int> gene_cols,
    int N_barcodes,
    int N_cells
  ) {
    
    int N_genes = gene_cols.size();
    
    // Convert mean/variance -> gamma params
    ArrayXd gamma_rt = gene_kernel_rates / gamma_variance;
    ArrayXd gamma_sp = gene_kernel_rates.square() / gamma_variance;
    
    // Output matrix
    MatrixXi barcode_counts_true = MatrixXi::Zero(N_cells, N_barcodes);
    
    // RNG
    std::mt19937 rng(std::random_device{}());
    
    // Precompute gamma distributions
    std::vector<std::gamma_distribution<double>> gamma_dists;
    gamma_dists.reserve(N_genes);
    for (int g = 0; g < N_genes; ++g) {
      gamma_dists.emplace_back(
        gamma_sp[g],          // shape
        1.0 / gamma_rt[g]     // scale = 1/rate
      );
    }
    
    // Simulate counts
    for (int cell = 0; cell < N_cells; ++cell) {
      for (int g = 0; g < N_genes; ++g) {
        // Gamma-Poisson draw
        double lambda = gamma_dists[g](rng);
        std::poisson_distribution<int> pois(lambda);
        barcode_counts_true(cell, gene_cols[g]) = pois(rng);
      }
    }
    
    return barcode_counts_true;
  }

SpotSim make_SpotSim(
    const MatrixXi& bc_counts, // random draw seeded by empirical observation, e.g., a MERFISH run
    const Codebook& cb, 
    const FlipRates& fr, 
    const std::unordered_map<uint64_t, int>& correction_table
  ) {
    
    // Collect basic information
    int total_spots = bc_counts.sum();
    int N_barcodes = bc_counts.cols();
    int N_cells = bc_counts.rows();
    int N_bits = cb.N_bits;
    
    // Initialize spot sim
    SpotSim sim;
    
    // Reserve space for barcodes and labels
    sim.barcodes_true.reserve(total_spots);
    sim.barcodes_read.reserve(total_spots);
    sim.barcodes_corrected.reserve(total_spots);
    sim.labels_true.reserve(total_spots);
    sim.labels_read.reserve(total_spots);
    sim.labels_corrected.reserve(total_spots);
    sim.cell_ids.reserve(total_spots);
    sim.lum.reserve(total_spots);
    
    // Compute inverse of flip rates using normal-distribution quantiles (inverse CDF)
    std::vector<double> rate10_inv(N_bits);
    std::vector<double> rate01_inv(N_bits);
    for (int b = 0; b < N_bits; ++b) {
      rate10_inv[b] = R::qnorm(fr.rate10[b], 0.0, 1.0, 1, 0); // qnorm(p, mean = 0, sd = 1, lower_tail = true, log_p = false)
      rate01_inv[b] = R::qnorm(1 - fr.rate01[b], 0.0, 1.0, 1, 0); 
    }
    
    // Set bit means and noise by barcode
    MatrixXd bit_means(N_barcodes, N_bits); // Make RowMajor?
    MatrixXd bit_noise(N_barcodes, N_bits);
    for (int i = 0; i < N_barcodes; ++i) {
      for (int b = 0; b < N_bits; ++b) {
        // Extract bit
        int bit = (cb.barcodes[i] >> b) & 1ULL;
        // Convert: 0 -> -1, 1 ->  1
        double m = (double)bit * 2.0 - 1.0;
        bit_means(i, b) = m;
        if (bit == 1) {
          bit_noise(i, b) = -m / rate10_inv[b];
        } else {
          bit_noise(i, b) = -m / rate01_inv[b];
        }
      }
    }
    
    // For each barcode ...
    for (int j = 0; j < N_barcodes; j++) {
      // ... make spots for all cells in one pass
      int count = bc_counts.col(j).sum();
      
      if (count >= 1) {
        // Build noise diagonal matrix
        MatrixXd noise = bit_noise.row(j).asDiagonal();
        
        // Build noised covariance matrix for this barcode
        MatrixXd corr_noised = noise * fr.corr * noise; 
        
        // Sample luminance levels from multivariate normal 
        MatrixXd lum = rmvnorm(count, bit_means.row(j).transpose(), corr_noised);
        
        // For each cell ...
        for (int i = 0; i < N_cells; ++i) {
          
          // ... simulate the number of spots for this barcode
          int bc_count = bc_counts(i, j); 
          for (int k = 0; k < bc_count; ++k) {
            
            // Grab spot and decode 
            VectorXd spot = lum.row(i*bc_count + k).transpose();
            uint64_t spot_bc = pack_decode_lum(spot);
            uint64_t spot_bc_corrected = spot_bc;
            int label_corrected;
            int label_read = -1;
            
            // Correct decoding
            auto it = correction_table.find(spot_bc);
            if (it == correction_table.end()) {
              label_corrected = -1; // uncorrectable
            } else {
              label_corrected = it->second;
            }
            
            if (label_corrected >= 0) {
              spot_bc_corrected = cb.barcodes[label_corrected];
              if (spot_bc_corrected == spot_bc) {
                label_read = label_corrected;
              }
            }
            
            // Save cell ID and spot-luminance information
            sim.cell_ids.push_back(i);
            sim.lum.push_back(spot);
            
            // Save labels 
            sim.labels_true.push_back(j);
            sim.labels_read.push_back(label_read);
            sim.labels_corrected.push_back(label_corrected);
            
            // Savd barcodes
            sim.barcodes_true.push_back(cb.barcodes[j]);
            sim.barcodes_read.push_back(spot_bc);
            sim.barcodes_corrected.push_back(spot_bc_corrected);
          }
        }
      }
    }
    
    // Compute flip rates pre- and post-correction
    std::pair<std::vector<double>, std::vector<double>> fr_pre = compute_flip_rates(sim.barcodes_true, sim.barcodes_read, N_bits);
    std::pair<std::vector<double>, std::vector<double>> fr_post = compute_flip_rates(sim.barcodes_true, sim.barcodes_corrected, N_bits);
    sim.rate01_pre = fr_pre.first;
    sim.rate10_pre = fr_pre.second;
    sim.rate01_post = fr_post.first;
    sim.rate10_post = fr_post.second;
    
    // Count barcodes
    std::vector<int> read_counts(N_barcodes, 0);
    std::vector<int> corrected_counts(N_barcodes, 0);
    std::vector<int> true_counts(N_barcodes, 0);
    for (size_t i = 0; i < sim.labels_corrected.size(); ++i) {
      int dc = sim.labels_read[i];
      if (dc >= 0) {++read_counts[dc];}
      int cc = sim.labels_corrected[i];
      if (cc < 0) {continue;}
      ++corrected_counts[cc];
      if (cc == sim.labels_true[i]) {++true_counts[cc];}
    }
    sim.read_counts = read_counts;
    sim.corrected_counts = corrected_counts;
    sim.true_counts = true_counts;
    
    // Compute CR and PPV
    std::vector<double> CR(N_barcodes, 0.0);
    std::vector<double> PPV(N_barcodes, 0.0);
    for (size_t i = 0; i < N_barcodes; ++i) {
      PPV[i] = corrected_counts[i] > 0 ? double(true_counts[i]) / double(corrected_counts[i]) : 0.0;
      CR[i] = corrected_counts[i] > 0 ? double(read_counts[i]) / double(corrected_counts[i]) : 0.0;
    }
    sim.CR = CR;
    sim.PPV = PPV;
    
    return sim;
  }

uint64_t pack(
    std::vector<int> bits
  ) {
    uint64_t packed = 0;
    for (int i = 0; i < bits.size(); ++i) {packed |= (uint64_t(bits[i]) << i);} 
    return packed;
  }

Codebook pack_codebook(
    const IntegerMatrix& codebook
  ) {
    int N_bits = codebook.ncol();
    int N_barcodes = codebook.nrow();
    CharacterVector species = rownames(codebook);
    Codebook cb;
    cb.N_bits = N_bits;
    cb.blanks = grep_idx(species, "Blank");
    cb.barcodes.reserve(N_barcodes);
    cb.species.reserve(N_barcodes);
    for (int i = 0; i < N_barcodes; ++i) {
      std::vector<int> bits(N_bits);
      for (int j = 0; j < N_bits; ++j) {bits[j] = codebook(i, j);}
      cb.barcodes.push_back(pack(bits));
      cb.species.push_back(as<std::string>(species[i]));
    }
    return cb;
  }

void generate_neighbors(
    uint64_t x,
    int n_bits,
    int dist,
    int start_bit,
    std::vector<uint64_t>& out
  ) {
    if (dist == 0) {
      out.push_back(x);
      return;
    }
    for (int b = start_bit; b < n_bits; ++b) {
      uint64_t flipped = x ^ (1ULL << b); // Flip bit b
      generate_neighbors(
        flipped,
        n_bits,
        dist - 1,
        b + 1,
        out
      );
    }
  }

std::vector<uint64_t> neighbors(
    uint64_t x,
    int n_bits,
    int hamming_dist
  ) {
    std::vector<uint64_t> out;
    generate_neighbors(
      x,
      n_bits,
      hamming_dist,
      0,
      out
    );
    return out;
  }

std::unordered_map<uint64_t, int> build_correction_table(
    const Codebook& cb,
    int max_correctable_Hamming_distance
  ) {
    std::unordered_map<uint64_t, int> correction_table;
    for (size_t i = 0; i < cb.barcodes.size(); ++i) {
      correction_table[cb.barcodes[i]] = i; // Exact match
      // Generate all barcodes within max_correctable_Hamming_distance
      for (int d = 1; d <= max_correctable_Hamming_distance; ++d) {
        std::vector<uint64_t> n = neighbors(cb.barcodes[i], cb.N_bits, d);
        for (int j = 0; j < n.size(); ++j) {
          if (correction_table.count(n[j])) {
            // If this neighbor is already mapped to a different barcode, we have a tie
            if (correction_table[n[j]] != i) {
              correction_table[n[j]] = -2; // Mark as ambiguous
            }
          } else {
            correction_table[n[j]] = i; // Map neighbor to original barcode index
          }
        }
      }
    }
    return correction_table;
  }

std::unordered_map<int, std::vector<uint64_t>> invert_lookup_table(
    const std::unordered_map<uint64_t, int>& lut
  ) {
    std::unordered_map<int, std::vector<uint64_t>> inverted;
    for (const auto& kv : lut) {
      uint64_t barcode = kv.first;
      int label = kv.second;
      inverted[label].push_back(barcode);
    }
    return inverted;
  } 

std::vector<int> decode_spots(
    std::vector<uint64_t> spots,
    const std::unordered_map<uint64_t, int>& correction_table
  ) {
    std::vector<int>  labels(spots.size());
    for (int i = 0; i < spots.size(); ++i) {
      auto it = correction_table.find(spots[i]);
      if (it == correction_table.end()) {
        labels[i] = -1; // uncorrectable
      } else {
        labels[i] = it->second;
      }
    }
    return labels;
  }


// Hamming distance
// int dist = __builtin_popcountll(a ^ b);



ST_data load_STdata(
    NumericVector bc_countsR,
    IntegerMatrix codebook,
    int max_correctable_Hamming_distance
  ) {
    
    Codebook cb = pack_codebook(codebook);
    
    int N_barcodes = bc_countsR.size();
    std::vector<int> bc_counts;
    bc_counts.reserve(N_barcodes);
    for (int i = 0; i < N_barcodes; ++i) {
      bc_counts.push_back((int)bc_countsR[i]);
    }
    
    std::unordered_map<uint64_t, int> correction_table = build_correction_table(cb, max_correctable_Hamming_distance);
    std::unordered_map<int, std::vector<uint64_t>> correction_table_inverted = invert_lookup_table(correction_table);
    
    Eigen::Matrix<uint64_t, Dynamic, Dynamic> transform_flips_all = find_transform_flips_all(cb.barcodes, cb.N_bits);
    
    return {bc_counts, cb, correction_table_inverted, transform_flips_all};
    
  }

// [[Rcpp::export]]
List estimate_flip_rates_initial_R( 
    NumericVector bc_countsR,
    IntegerMatrix codebook,
    int max_correctable_Hamming_distance
  ) {
    auto STdata = load_STdata(bc_countsR, codebook, max_correctable_Hamming_distance);
    auto flip_rates = estimate_flip_rates_initial(STdata);
    return List::create(
      _["rate10"] = flip_rates.first,
      _["rate01"] = flip_rates.second
    );
  }

