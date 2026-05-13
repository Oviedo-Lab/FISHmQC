
// Rcpp
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <nlopt.hpp>
#include <boost/math/distributions/normal.hpp>
using namespace Rcpp;
using namespace Eigen;

struct Codebook {
  int N_bits;
  std::vector<uint64_t> barcodes;
  std::vector<string> species;
};

struct SpotSim {
  std::vector<uint64_t> barcodes_true;
  std::vector<uint64_t> barcodes_read; 
  std::vector<uint64_t> barcodes_corrected;
  std::vector<int> labels_true;
  std::vector<int> labels_read;
  std::vector<int> labels_corrected;
  std::vector<int> cell_ids;
};

/*
 * Next up: remake the set_true_spot_info and generate_simulated_spots functions; 
 * releated: set_noise_by_flip_rates and set_luminance_noise_correlation_matrix
 */

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
NumericMatrix generate_spots_ci(
    IntegerMatrix bc_counts,
    NumericMatrix codebook_gen_mat,
    NumericMatrix noise,
    NumericMatrix corr_mat,
    int ci, // "ci" = cell index (cell number), used to grab rows from bc_counts
    double threshold,
    bool decode_and_label
  ) {
    int N_barcodes = bc_counts.ncol();
    int N_bits = codebook_gen_mat.ncol();
    
    int total_spots = 0;
    for (int j = 0; j < N_barcodes; j++) {
      total_spots += bc_counts(ci, j);
    }
    
    NumericMatrix spots_ci(total_spots, N_bits);
    std::fill(spots_ci.begin(), spots_ci.end(), NA_REAL);
    
    int end_row = 0;
    
    Environment MASS("package:MASS");
    Function mvrnorm = MASS["mvrnorm"];
    
    for (int bci = 0; bci < N_barcodes; bci++) {
      int count = bc_counts(ci, bci);
      end_row += count;
      int start_row = end_row - count;
      
      if (count >= 1) {
        // Get this barcode's mean vector
        NumericVector this_bc(N_bits);
        for (int k = 0; k < N_bits; k++) {
          this_bc[k] = codebook_gen_mat(bci, k);
        }
        
        // Build noise diagonal matrix
        NumericMatrix noise_matrix(N_bits, N_bits);
        for (int i = 0; i < N_bits; i++) {
          for (int j = 0; j < N_bits; j++) {
            if (i == j) {
              noise_matrix(i, j) = noise(bci, i);
            } else {
              noise_matrix(i, j) = 0.0;
            }
          }
        }
        
        // Sigma_noised = noise_matrix %*% corr_mat %*% noise_matrix
        NumericMatrix temp_sigma(N_bits, N_bits);
        for (int i = 0; i < N_bits; i++) {
          for (int j = 0; j < N_bits; j++) {
            double sum = 0.0;
            for (int k = 0; k < N_bits; k++) {
              sum += noise_matrix(i, k) * corr_mat(k, j);
            }
            temp_sigma(i, j) = sum;
          }
        }
        
        NumericMatrix Sigma_noised(N_bits, N_bits);
        for (int i = 0; i < N_bits; i++) {
          for (int j = 0; j < N_bits; j++) {
            double sum = 0.0;
            for (int k = 0; k < N_bits; k++) {
              sum += temp_sigma(i, k) * noise_matrix(k, j);
            }
            Sigma_noised(i, j) = sum;
          }
        }
        
        // Sample from multivariate normal using MASS::mvrnorm
        NumericMatrix spot_seeds(count, N_bits);
        if (count > 1) {
          spot_seeds = mvrnorm(count, this_bc, Sigma_noised);
        } else {
          spot_seeds.row(0) = Rcpp::as<NumericVector>(mvrnorm(count, this_bc, Sigma_noised));
        }
        
        for (int row = 0; row < count; row++) {
          for (int col = 0; col < N_bits; col++) {
            double val = spot_seeds(row, col);
            if (decode_and_label) {
              spots_ci(start_row + row, col) = (val > threshold) ? 1.0 : 0.0;
            } else {
              spots_ci(start_row + row, col) = val;
            }
          }
        }
      }
    }
    
    return spots_ci;
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

