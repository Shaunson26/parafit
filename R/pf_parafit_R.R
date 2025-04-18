#' Parafit
#'
#' Run the parafit algorithm to obtain the global test statistics and obtain a test
#' of significance using permutation. Optionally, obtain tests of individual associations.
#'
#' @param assoA association matrix
#' @param paraB parasite principal coordinates
#' @param hostC transposed host principal coordinates
#' @param permutations number of permutations
#' @param seed random seed value
#'
#' @return vector of results, with elements global and links (the latter if test_links = TRUE)
#' @export
pf_parafit_r <- function(assoA, paraB, hostC, permutations, seed){

  set.seed(seed)

  trace_results = vector('numeric', permutations)

  for(i in 1:permutations){

    assoA_i = assoA

    if (i > 1){
      assoA_i <- t(apply(assoA_i, MARGIN = 1, sample))
    }

    matD <- hostC %*% assoA_i %*% paraB
    trace_results[i] <- sum(matD^2)

  }

  return(matrix(trace_results, ncol = 1))

}








