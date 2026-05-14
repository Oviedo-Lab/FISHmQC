
// Rcpp
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <random>
#include <nlopt.hpp>
#include <boost/math/distributions/normal.hpp>
using namespace Rcpp;
using namespace Eigen;

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

struct Codebook {
  int N_bits;
  std::vector<uint64_t> barcodes;
  std::vector<string> species;
};

struct FlipRates {
  std::vector<double> one_zero; // P(1 -> 0) for each bit, so length = N_bits
  std::vector<double> zero_one; // P(0 -> 1) for each bit, so length = N_bits
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
};

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

// [[Rcpp::export]]
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

// [[Rcpp::export]]
SpotSim make_SpotSim(
    const MatrixXi& bc_counts, // random draw seeded by empirical observation, e.g., a MERFISH run
    const Codebook& cb, 
    const FlipRates& fr, 
    const std::unordered_map<uint64_t, int>& correction_table
  ) {
    
    // Collect basic information
    int total_spots = bc_counts.sum();
    int N_barcodes = bc_counts.ncol();
    int N_cells = bc_counts.nrow();
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
    
    // Compute inverse of flip rates using normal distribution quantiles (inverse CDF)
    std::vector<double> one_zero_inv(N_bits);
    std::vector<double> zero_one_inv(N_bits);
    for (int b = 0; b < N_bits; ++b) {
      one_zero_inv[b] = R::qnorm(fr.one_zero[b], 0.0, 1.0, 1, 0); // qnorm(p, mean=0, sd=1, lower_tail=true, log_p=false)
      zero_one_inv[b] = R::qnorm(1 - fr.zero_one[b], 0.0, 1.0, 1, 0); 
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
          bit_noise(i, b) = -m / one_zero_inv[b];
        } else {
          bit_noise(i, b) = -m / zero_one_inv[b];
        }
      }
    }
    
    for (int j = 0; j < N_barcodes; j++) {
      
      // For each barcode, make spots for all cells in one pass
      int count = bc_counts.col(j).sum();
      
      if (count >= 1) {
        
        // Build noise diagonal matrix
        MatrixXd noise = bit_noise.row(j).asDiagonal();
        
        // Build noise covariance matrix for this barcode
        MatrixXd corr_noised = noise * fr.corr * noise; 
        
        // Sample from multivariate normal 
        MatrixXd lum = rmvnorm(count, bit_means.row(j).transpose(), corr_noised);
        
        for (int i = 0; i < N_cells; ++i) {
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
    CharacterVector species = codebook(m);
    Codebook cb;
    cb.N_bits = N_bits;
    cb.barcodes.reserve(N_barcodes);
    cb.species.reserve(N_barcodes);
    for (int i = 0; i < N_barcodes; ++i) {
      std::vector<int> bits(N_bits);
      for (int j = 0; j < N_bits; ++j) {bits[j] = codebook(i, j);}
      cb.barcodes.push_back(pack(bits));
      cb.species.push_back(species[i]);
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










// [[Rcpp::export]]
NumericVector find_HD_by_row(
    const NumericMatrix& codebook, 
    const NumericVector& spot
  ) {
    int n = codebook.nrow();
    NumericVector HD = NumericVector(n);
    for (int i = 0; i < n; i++) {
      HD[i] = sum(abs(codebook(i, _) - spot));
    }
    return HD;
  }

// [[Rcpp::export]]
List find_spot_labels_fast(
    const IntegerMatrix& codebook, 
    const IntegerMatrix& spots,
    const int& max_correctable_Hamming_distance
  ) {
    int n_spots = spots.nrow();
    int n_codes = codebook.nrow();
    int m = codebook.ncol(); // dimension of each spot/codeword
    IntegerVector spot_labels(n_spots, NA_INTEGER);
    IntegerVector spot_labels_uncorrected(n_spots, NA_INTEGER);
    
    for (int i = 0; i < n_spots; i++) {
      int min_HD = m + 1;  // max possible HD is m
      int min_idx = -1;
      bool tie = false;
      
      for (int j = 0; j < n_codes; j++) {
        int h = 0;
        for (int k = 0; k < m; k++) {
          h += (static_cast<int>(spots(i, k)) != static_cast<int>(codebook(j, k)));
          if (h > min_HD) break;
        }
        if (h < min_HD) {
          min_HD = h;
          min_idx = j;
          tie = false;
        } else if (h == min_HD) {
          tie = true;
        }
        if (min_HD == 0) {
          break;
        }
      }
      
      if (min_HD <= max_correctable_Hamming_distance && !tie) {
        spot_labels[i] = min_idx + 1; // 1-indexed for R
        if (min_HD == 0) {
          spot_labels_uncorrected[i] = min_idx + 1; // 1-indexed for R
        }
      }
      
    }
    
    List spot_labels_list = List::create(
      Named("spot_labels") = spot_labels,
      Named("spot_labels_uncorrected") = spot_labels_uncorrected
    );
    
    return spot_labels_list;
  }


// [[Rcpp::export]]
IntegerVector find_nearby_barcodes_fast(
  const int& bci,
  const IntegerMatrix& codebook,
  const int& max_correctable_Hamming_distance
  ) {
    int n_barcodes = codebook.nrow();
    int n_bits = codebook.ncol(); 
    IntegerVector barcode = codebook(bci, _);
    IntegerVector nearby_barcodes;
    for (int i = 0; i < n_barcodes; i++) {
      IntegerVector codebook_i = codebook(i, _);
      int HD = 0;
      for (int k = 0; k < n_bits; k++) {
        HD += (static_cast<int>(codebook_i[k]) != static_cast<int>(barcode[k]));
        if (HD > max_correctable_Hamming_distance) {
          break; // No need to continue if HD exceeds the limit
        }
      }
      if (HD <= max_correctable_Hamming_distance) {
        nearby_barcodes.push_back(i + 1); // R is 1-indexed
      }
    }
    return(nearby_barcodes);
  }

// Function to compute PPV for each barcode
// ... want to answer the question: for each barcode bc, if a spot is decoded as bc, what is the probability this is correct?
// [[Rcpp::export]]
NumericVector compute_PPV_fast(
    const List& barcodes_by_cell_true_list,
    const List& spots_bc_labels_list,
    const int& N_barcodes
  ) {
    
    // Initialize vectors to hold counts 
    int N_cells = barcodes_by_cell_true_list.size();
    IntegerVector decoded_counts(N_barcodes, 0);
    IntegerVector successful_decoded_counts(N_barcodes, 0);
    
    for (int ci = 0; ci < N_cells; ci++) {
      // Get decoded counts and true barcode IDs for this cell
      IntegerVector true_barcodes_ids = barcodes_by_cell_true_list[ci];
      IntegerVector decoded_barcode = spots_bc_labels_list[ci];
      LogicalVector successfully_decoded = !is_na(decoded_barcode);
      int n_spots = true_barcodes_ids.size();
      for (int si = 0; si < n_spots; si++) {
        if (successfully_decoded[si]) {
          // For this spot, get its true barcode 
          int true_bc = true_barcodes_ids[si] - 1; // R is 1-indexed, convert to 0-indexed
          // Get the decoded barcode
          int this_read = decoded_barcode[si] - 1; // R is 1-indexed, convert to 0-indexed
          // Increment the counts
          decoded_counts[this_read]++;
          if (this_read == true_bc) {
            successful_decoded_counts[this_read]++;
          }
        }
      }
    }
    
    NumericVector PPV_true = as<NumericVector>(successful_decoded_counts) / as<NumericVector>(decoded_counts);
    return PPV_true;
    
  }

void run_MCMCSA_cpp(
  const VectorXi& codebook,
  const MatrixXd& summary_stats_genes,
  const MatrixXd& summary_stats_blanks,
  const MatrixXd& schedules,
  const VectorXd& initial_state
  ) {
    
  }

