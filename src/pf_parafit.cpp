// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <algorithm>
#include <vector>
#include <numeric>
#include <mutex> // Include mutex for thread-safe logging
using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;

std::mutex mtx; // Define a global mutex for synchronized printing

// Worker for parallel computation
struct Parafit_Permute : public Worker {
  const arma::mat& assoA;
  const arma::mat& paraB;
  const arma::mat& hostC;
  const arma::ucube& perm_cube;
  arma::vec& trace_results;
  const bool verbose;
  const int print_n;

  Parafit_Permute(const arma::mat& assoA, const arma::mat& paraB, const arma::mat& hostC, arma::ucube& perm_cube, arma::vec& trace_results, bool verbose, int print_n)
    : assoA(assoA), paraB(paraB), hostC(hostC), perm_cube(perm_cube), trace_results(trace_results), verbose(verbose), print_n(print_n) {}

  // Parallel computation
  void operator()(std::size_t begin, std::size_t end) override {

    // Parallel iters
    for (std::size_t i = begin; i < end; ++i) {

      // Shuffle rows independently
      arma::umat current_indices = perm_cube.slice(i);
      int n_rows = current_indices.n_rows;
      int n_cols = current_indices.n_cols;
      arma::mat assoA_i(n_rows, n_cols);
      for (int row = 0; row < n_rows; ++row) {
        for (int col = 0; col < n_cols; ++col) {
          assoA_i(row, col) = assoA(row, current_indices(row, col));
        }
      }

      // Perform the matrix multiplication to create D
      arma::mat matD = hostC * assoA_i * paraB;

      // Compute the sum of squares of the elements in the result matrix
      trace_results[i] = arma::accu(arma::square(matD));

      if (verbose){
        if (i % print_n == 0){
          std::lock_guard<std::mutex> lock(mtx);
          Rcpp::Rcout << "." << "\n";
        }
      }

    }

  }

};

//' Create a permutation cube
//'
//' Permute the values of each row independently, in multiple slices (permutations).
//' First permutation is the original indices.
//'
//' @param input_matrix matrix for dimensions
//' @param permutations number of permutations
//' @param seed random seed value
//'
//' @return Array with matrix permutations
// [[Rcpp::export]]
arma::ucube permute_cube_indices(const arma::mat& input_matrix, int permutations, int seed) {

  int n_rows = input_matrix.n_rows;
  int n_cols = input_matrix.n_cols;
  arma::ucube output_cube(n_rows, n_cols, permutations);

  // Create a vector of column indices
  std::vector<unsigned int> original_indices(n_cols);
  std::iota(original_indices.begin(), original_indices.end(), 0);

  arma::umat first_slice(n_rows, n_cols);
  for (int i = 0; i < n_rows; ++i) {
    for (int j = 0; j < n_cols; ++j) {
      first_slice(i, j) = original_indices[j];
    }
  }

  output_cube.slice(0) = first_slice;

  // Create a random engine
  std::mt19937 gen(seed);

  for (int k = 1; k < permutations; ++k) {
    arma::umat slice(n_rows, n_cols);
    for (int i = 0; i < n_rows; ++i) {
      // Create a copy of the column indices for each row
      std::vector<unsigned int> current_indices = original_indices;
      // Shuffle the column indices
      std::shuffle(current_indices.begin(), current_indices.end(), gen);
      // Assign the permuted indices to the current row of the slice
      for (int j = 0; j < n_cols; ++j) {
        slice(i, j) = current_indices[j];
      }
    }
    output_cube.slice(k) = arma::conv_to<arma::umat>::from(slice);
  }

  return output_cube;
}

//' Run the parafit algorithm
//'
//' Run the parafit algorithm to obtain the global test statistics and obtain a test
//' of significance using permutation. Optionally, obtain tests of individual associations.
//'
//' @param assoA association matrix
//' @param paraB parasite principal coordinates
//' @param hostC transposed host principal coordinates
//' @param permutations number of permutations
//' @param seed random seed value
//' @param verbose whether to print interation number to console
//' @param print_n print interation number at every print_n
//'
//' @return list of results, with elements global and links (the latter if test_links = TRUE).
// [[Rcpp::export]]
arma::vec pf_parafit_cpp(const arma::mat& assoA, const arma::mat& paraB, const arma::mat& hostC, int permutations, bool verbose, int print_n, int seed) {

  //int num_cols = HP.n_cols;
  arma::vec trace_results(permutations);
  arma::ucube perm_cube = permute_cube_indices(assoA, permutations, seed);

  int n_rows = assoA.n_rows;
  int n_cols = assoA.n_cols;

  for(int i = 0; i < permutations; i++){

    // Shuffle rows independently
    arma::umat current_indices = perm_cube.slice(i);

    arma::mat assoA_i(n_rows, n_cols);
    for (int row = 0; row < n_rows; ++row) {
      for (int col = 0; col < n_cols; ++col) {
        assoA_i(row, col) = assoA(row, current_indices(row, col));
      }
    }

    // Perform the matrix multiplication to create D
    arma::mat matD = hostC * assoA_i * paraB;

    // Compute the sum of squares of the elements in the result matrix
    trace_results[i] = arma::accu(arma::square(matD));



    if (verbose){
      if (i % print_n == 0){
        Rcpp::Rcout << i << "\n";
      }
    }

  }

  return trace_results;
}

//' Run the parafit algorithm
//'
//' Run the parafit algorithm to obtain the global test statistics and obtain a test
//' of significance using permutation. Optionally, obtain tests of individual associations.
//'
//' @param assoA association matrix
//' @param paraB parasite principal coordinates
//' @param hostC transposed host principal coordinates
//' @param permutations number of permutations
//' @param seed random seed value
//' @param verbose whether to print iteration number to the console
//' @param print_n the level to print iteration numbers by
//'
//' @return list of results, with elements global and links (the latter if test_links = TRUE).
// [[Rcpp::export]]
arma::vec pf_parafit_parallel_cpp(const arma::mat& assoA, const arma::mat& paraB, const arma::mat& hostC, int permutations, int seed, const bool verbose, const int print_n) {

  arma::vec trace_results(permutations);
  arma::ucube perm_cube = permute_cube_indices(assoA, permutations, seed);

  Parafit_Permute parafit_permute(assoA, paraB, hostC, perm_cube, trace_results, verbose, print_n);
  RcppParallel::parallelFor(0, permutations, parafit_permute);

  return trace_results;

}




