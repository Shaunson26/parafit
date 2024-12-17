#' Obtain Principal Coordinate from a distance matrix
#'
#' Obtain principal coordinates from a distance matrix (triangle or square)
#'
#' @param x a distance matrix
#'
#' @return matrix of principal coordinates
#'
#' @export
pf_pcoa <- function(x){
  pf_pcoa_cpp(as.matrix(x))
}
