#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' Obtain Principal Coordinate from a distance matrix
//'
//' Obtain principal coordinates from a distance matrix (triangle or square)
//'
//' @param x a distance matrix
//'
//' @return matrix of principal coordinates
// [[Rcpp::export]]
arma::mat pf_pcoa_cpp(const arma::mat& x) {
  int n = x.n_rows;

  // Compute squared distances
  arma::mat D2 = arma::pow(x, 2.0);

  // Compute the centering matrix
  arma::mat I = arma::eye<arma::mat>(n, n);
  arma::mat One = arma::ones<arma::mat>(n, n);
  arma::mat H = I - One / n;

  // Double-centering the squared distance matrix
  arma::mat B = -0.5 * H * D2 * H;

  // Eigen decomposition of B
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, B);

  // Select positive eigenvalues and corresponding eigenvectors
  arma::uvec positive_indices = arma::find(eigval > 0);
  arma::vec positive_eigval = eigval(positive_indices);
  arma::mat positive_eigvec = eigvec.cols(positive_indices);

  // Compute principal coordinates
  arma::mat principal_coords = positive_eigvec * arma::diagmat(arma::sqrt(positive_eigval));

  // Reverse the column order
  principal_coords = principal_coords.cols(arma::reverse(arma::regspace<arma::uvec>(0, positive_eigvec.n_cols - 1)));


  return principal_coords;
}
