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
Rcpp::List pf_pcoa_cpp(Rcpp::RObject x) {

  // Check if the input is a matrix
  if (!Rf_isMatrix(x)) {
    Rcpp::stop("pf_pcoa_cpp: Input must be a matrix, not a data frame or other object.");
  }

  arma::mat x_mat = Rcpp::as<arma::mat>(x);

  int n = x_mat.n_rows;

  // Compute squared distances
  arma::mat D2 = arma::pow(x_mat, 2.0);

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

  // Compute eigenvalues matrix
  //arma::vec relative_eigval = positive_eigval / arma::sum(positive_eigval);
  //arma::vec cumulative_eigval = arma::cumsum(relative_eigval);
  //arma::mat eigenvalues_matrix = arma::join_horiz(positive_eigval, relative_eigval, cumulative_eigval);

  // Set column names for eigenvalues matrix
  //Rcpp::NumericMatrix eigenvalues_R(eigenvalues_matrix.n_rows, eigenvalues_matrix.n_cols);
  //std::copy(eigenvalues_matrix.begin(), eigenvalues_matrix.end(), eigenvalues_R.begin());
  //Rcpp::CharacterVector colnames = {"eigenvalue", "relative_eigenvalue", "cumulative_eigenvalue"};
  //Rcpp::colnames(eigenvalues_R) = colnames;

  // Compute principal coordinates
  arma::mat principal_coords = positive_eigvec * arma::diagmat(arma::sqrt(positive_eigval));

  // Reverse the column order
  //principal_coords = principal_coords.cols(arma::reverse(arma::regspace<arma::uvec>(0, principal_coords.n_cols - 1)));

  //return principal_coords;
  return Rcpp::List::create(Rcpp::Named("values") = positive_eigval, Rcpp::Named("vectors") = principal_coords);
}
