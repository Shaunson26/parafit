// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppParallel)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <algorithm>
#include <vector>
#include <numeric>
#include <random>
using namespace Rcpp;
using namespace arma;
using namespace RcppParallel;


// Worker for parallel computation
struct Parafit_Permute : public Worker {
  const arma::mat& assoA;
  const arma::mat& paraB;
  const arma::mat& hostC;
  arma::vec& trace_results;
  int seed;

  Parafit_Permute(const arma::mat& assoA, const arma::mat& paraB, const arma::mat& hostC, arma::vec& trace_results, int seed)
    : assoA(assoA), paraB(paraB), hostC(hostC), trace_results(trace_results), seed(seed) {}

  // Parallel computation
  void operator()(std::size_t begin, std::size_t end) override {

    std::mt19937 rng;
    int n_rows = assoA.n_rows;

    // Parallel iters
    for (std::size_t i = begin; i < end; ++i) {

      rng.seed(seed + i);  // Deterministic seed per iteration
      arma::mat assoA_i = assoA;

      if (i > 0) {

        for (int row = 0; row < n_rows; ++row){
          arma::rowvec row_data = assoA_i.row(row);  // copy row to a rowvec
          std::vector<double> row_vector(row_data.begin(), row_data.end()); // convert to std::vector
          std::shuffle(row_vector.begin(), row_vector.end(), rng);  // shuffle with thread-local RNG
          assoA_i.row(row) = arma::rowvec(row_vector);  // copy back shuffled data
        }

      }

      // Perform the matrix multiplication to create D
      arma::mat matD = hostC * assoA_i * paraB;

      // Compute the sum of squares of the elements in the result matrix
      trace_results[i] = arma::accu(arma::square(matD));

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
//' @param seed random seed value
//'
//' @return vector of trace statistics
// [[Rcpp::export]]
arma::vec pf_parafit_cpp(const arma::mat& assoA, const arma::mat& paraB, const arma::mat& hostC, int permutations, int seed) {

  arma::vec trace_results(permutations);
  Parafit_Permute parafit_permute(assoA, paraB, hostC, trace_results, seed);
  RcppParallel::parallelFor(0, permutations, parafit_permute);

  return trace_results;

}




