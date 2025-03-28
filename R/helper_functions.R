#' Wrangle PCOA
#'
#' Add additional objects to PCOA
#'
#' @param x POCA from ...
#' @return list
wrangle_pf_pcoa_cpp <- function(x){

  stopifnot(inherits(x, 'list'))

  x$values[,1] <- x$values[nrow(x$values):1,]
  eigen_relative <-  x$value[,1] / sum(x$value[,1])
  eigen_cumulative <- cumsum(eigen_relative)
  x$values <-
    data.frame(
      eigenvalue =  x$values[,1],
      relative_eigenvalue = eigen_relative,
      cumulative_eigenvalue = eigen_cumulative
    )

  x$vectors <- x$vectors[, ncol(x$vectors):1]

  x
}
