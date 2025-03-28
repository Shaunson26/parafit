// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <mutex> // Include mutex for thread-safe logging
using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;

std::mutex mtx; // Define a global mutex for synchronized printing

// Worker for parallel computation
struct parafit_permute : public Worker {
  const arma::mat& assoA;
  const arma::mat& paraB;
  const arma::mat& hostC;
  const arma::umat& perm_matrix;
  arma::vec& trace_results;
  const bool verbose;
  const int print_n;

  parafit_permute(const arma::mat& assoA, const arma::mat& paraB, const arma::mat& hostC, arma::umat& perm_matrix, arma::vec& trace_results, bool verbose, int print_n)
    : assoA(assoA), paraB(paraB), hostC(hostC), perm_matrix(perm_matrix), trace_results(trace_results), verbose(verbose), print_n(print_n) {}

  // Parallel computation
  void operator()(std::size_t begin, std::size_t end) override {

    for (std::size_t i = begin; i < end; ++i) {

      // Perform the matrix multiplication to create D
      arma::uvec current_permutation = perm_matrix.col(i);
      arma::mat assoA_i = assoA.rows(current_permutation);
      arma::mat matD = hostC * assoA_i * paraB;

      // Compute the sum of squares of the elements in the result matrix
      trace_results[i] = arma::accu(arma::square(matD));

      if (verbose){
        if (i % print_n == 0){
          std::lock_guard<std::mutex> lock(mtx);
          Rcpp::Rcout << i << "\n";
        }
      }

    }

  }

};


//' Run the parafit algorithm
//'
//' Run the parafit algorithm to obtain the global test statistics and obtain a test
//' of significance using permutation. Optionally, obtain tests of individual associations.
//'
//' @param assoA association matrix
//' @param paraB parasite principal coordinates
//' @param hostC transposed host principal coordinates
//' @param permutations number of permutations
//' @param verbose whether to print iteration number to the console
//' @param print_n the level to print iteration numbers by
//'
//' @return list of results, with elements global and links (the latter if test_links = TRUE).
// [[Rcpp::export]]
arma::vec pf_parafit_parallel_cpp(const arma::mat& assoA, const arma::mat& paraB, const arma::mat& hostC, std::size_t permutations, const bool verbose, const int print_n) {

  arma::vec trace_results(permutations);

  int n_rows = assoA.n_rows;

  arma::umat perm_matrix(n_rows, permutations);

  for(std::size_t i = 0; i < permutations; ++i){

    arma::uvec perm_indices;

    if (i == 0){
      perm_indices = arma::linspace<arma::uvec>(0, n_rows - 1, n_rows);
    } else {
      perm_indices = arma::randperm(n_rows);
    }

    perm_matrix.col(i) = perm_indices;

  };

  parafit_permute do_parafit_permute(assoA, paraB, hostC, perm_matrix, trace_results, verbose, print_n);
  RcppParallel::parallelFor(0, permutations, do_parafit_permute);

  return trace_results;

}




