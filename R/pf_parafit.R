#' Run the parafit algorithm
#'
#' Run the parafit algorithm to obtain the global test statistics and obtain a test
#' of significance using permutation. Optionally, obtain tests of individual associations.
#'
#' @param host_pcoa host principal coordinates
#' @param parasite_pcoa parasite principal coordinates
#' @param associations association matrix
#' @param permutations number of permutations to conduct for significance testing
#' @param test_links test individual associations in associations matrix
#' @param seed seed for randomisation
#' @param parallel test links with parallel computing
#' @param cores number of cores if using parallel computing
#'
#' @return list of results, with elements global and links (the latter if test_links = TRUE).
#' @export
pf_parafit <- function(host_pcoa, parasite_pcoa, associations, permutations, test_links = FALSE, seed, parallel = F, cores = 1){

  stopifnot(
    'host_pcoa column number must match associations row number' = ncol(host_pcoa) == nrow(associations),
    'parasite_pcoa column number must match associations column number' = nrow(parasite_pcoa) == ncol(associations),
    'associations must be a matrix' = inherits(associations, 'matrix')
  )

  message('Running parafit ...')


  if (missing(seed)){
    seed = sample(1:1000, 1)
  }

  # Global
  set.seed(seed)

  global_trace = pf_parafit_cpp(assoA = associations, paraB = parasite_pcoa, hostC = host_pcoa, permutations = permutations)

  results <-
    list(
      global =
        list(
          stat = global_trace[1],
          p = sum(global_trace >= global_trace[1]) / length(global_trace),
          permutations = permutations,
          stat_perm = global_trace[,1]
        )
    )

  # Links
  if (test_links) {

    test_link = function(i){

      associations_0 <- associations

      associations_0[
        HP_link_inds[i, 'row'],
        HP_link_inds[i, 'col']
      ] <- 0

      set.seed(seed)

      pf_parafit_cpp(assoA = associations_0, paraB = parasite_pcoa, hostC = host_pcoa, permutations = permutations)

    }

    HP_link_inds <- which(associations > 0, arr.ind = TRUE)
    HP_link_inds <- HP_link_inds[ order(HP_link_inds[,'row']),]
    n_HP_links <- nrow(HP_link_inds)

    message('There are ', n_HP_links, ' links to be tested', sep = '', appendLF = FALSE)

    # link_traces <-
    #   sapply(seq(n_HP_links), function(i){
    #
    #     if (i == 1 || i %% 5 == 0){ message(i)}
    #
    #     associations_0 <- associations
    #
    #     associations_0[
    #       HP_link_inds[i, 'row'],
    #       HP_link_inds[i, 'col']
    #     ] <- 0
    #
    #     set.seed(seed)
    #
    #     pf_parafit_cpp(assoA = associations_0, paraB = parasite_pcoa, hostC = host_pcoa, permutations = permutations)
    #   })

    if (parallel){
      message(' ... in parallel')
      cl <- parallel::makeCluster(cores)
      parallel::clusterEvalQ(cl, devtools::load_all())
      #parallel::clusterEvalQ(cl, require(parafit))

      link_traces <- pbapply::pbsapply(seq(n_HP_links), FUN = test_link, cl = cl)
      parallel::stopCluster(cl)
    }

    if (!parallel){
      message('')
      link_traces <- sapply(seq(n_HP_links), test_link)
    }

    link_stats <-
      global_trace[,1] - link_traces

    link_pvalues <-
      apply(link_stats, MARGIN = 2, function(x) sum(x >= x[1]) / length(x))

    results$links <-
      data.frame(
        host = HP_link_inds[,'row'],
        parasite = HP_link_inds[,'col'],
        stat_1 = unname(link_stats[1,]),
        p_1 = link_pvalues
      )

  }

  message('Done')

  structure(
    .Data = results,
    class = 'pf_parafit'
  )
}

#' Print method for pf_parafit class
#'
#' @param x pf_parafit object
#' @export print.pf_parafit
#' @export
print.pf_parafit <- function(x){
  cat('Parafit\n\n')
  cat('ParaFitGlobal = ', x$global$stat, ', P value = ', round(x$global$p, 4), ' (', x$global$permutations, ')\n\n', sep = '')
  cat('Individual tests of links\n\n')
  if (is.null(x$links)) {
    cat('Not tested')
  } else {
  print(x$links)
  }
}
