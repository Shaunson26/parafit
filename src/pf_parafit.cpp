#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Run the parafit algorithm
//'
//' Run the parafit algorithm to obtain the global test statistics and obtain a test
//' of significance using permutation. Optionally, obtain tests of individual associations.
//'
//' @param assoA association matrix
//' @param paraB parasite principal coordinates
//' @param hostC transposed host principal coordinates
//' @param permutations number of permutations
//' @param verbose whether to print interation number to console
//' @param print_n print interation number at every print_n
//'
//' @return list of results, with elements global and links (the latter if test_links = TRUE).
// [[Rcpp::export]]
arma::vec pf_parafit_cpp(const arma::mat& assoA, const arma::mat& paraB, const arma::mat& hostC, int permutations, bool verbose, int print_n) {

  //int num_cols = HP.n_cols;
  arma::vec trace_results(permutations);

  for(int i = 0; i < permutations; i++){

    arma::mat assoA_permuted = assoA;

    if (i > 0) {

      assoA_permuted = arma::shuffle(assoA_permuted, 0);

    }

    // Perform the matrix multiplication to create D
    arma::mat matD = hostC * assoA_permuted * paraB;

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

