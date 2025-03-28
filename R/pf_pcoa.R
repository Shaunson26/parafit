#' Obtain Principal Coordinate from a distance matrix
#'
#' Obtain principal coordinates from a distance matrix (triangle or square)
#'
#' @param x a distance matrix
#' @param correction description
#' @param use_cpp use ape pacakge or CPP
#'
#' @return matrix of principal coordinates
#'
#' @export
pf_pcoa <- function(x, correction = "none", use_cpp = FALSE) {

  if (use_cpp){

    x_pc <- pf_pcoa_cpp(as.matrix(x))

    return(
      list(
        pco = x_pc$values,
        vectors = x_pc$vectors
      )
    )

  }

  epsilon <- sqrt(.Machine$double.eps)

  x_pc <- ape::pcoa(x, correction = correction)

  if (x_pc$correction[2] == 1) {
    stopifnot(
      'Matrix has negative eigenvalues. Rerun with correction="lingoes" or correction="cailliez"' = min(x_pc$values[, 'Relative_eig']) > -epsilon
    )

    sum_eig_sq <- sum(x_pc$values[, 'Eigenvalues']^2)
    x_vectors <- x_pc$vectors

  }

  if (x_pc$correction[2] != 1) {
    sum_eig_sq <- sum(x_pc$values[, 'Corr_eig']^2)
    x_vectors <- x_pc$vectors.cor
  }

  list(
    pco = x_pc,
    vectors = x_vectors,
    sum_eig_sq = sum_eig_sq
  )

}
