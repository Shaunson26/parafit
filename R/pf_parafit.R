#' Run the parafit algorithm
#'
#' Run the parafit algorithm to obtain the global test statistics and obtain a test
#' of significance using permutation. Optionally, obtain tests of individual associations.
#'
#' @param host_pcoa host principal coordinates (samples in rows, axes in columns)
#' @param parasite_pcoa parasite principal coordinates (samples in rows, axes in columns)
#' @param associations association matrix (host in rows, parasite in columns)
#' @param permutations number of permutations to conduct for significance testing
#' @param test_links test individual associations in associations matrix
#' @param seed seed for randomisation
#' @param parallel test links with parallel computing
#' @param cores number of cores if using parallel computing
#' @param verbose whether to print iteration numbers to the console
#' @param .print_n what iteration values to print by
#' @param use_r use R version of parafit (slower)
#'
#' @return list of results, with elements global and links (the latter if test_links = TRUE).
#' @export
pf_parafit <- function(host_pcoa, parasite_pcoa, associations, permutations, test_links = FALSE, seed, parallel = F, cores, verbose = FALSE, .print_n = 100, use_r = FALSE){

  stopifnot(
    'host_pcoa must be a matrix' = inherits(host_pcoa, 'matrix'),
    'parasite_pcoa must be a matrix' = inherits(parasite_pcoa, 'matrix'),
    'associations must be a matrix' = inherits(associations, 'matrix'),
    'host_pcoa column number must match associations row number' = nrow(host_pcoa) == nrow(associations),
    'parasite_pcoa column number must match associations column number' = nrow(parasite_pcoa) == ncol(associations)
  )

  if (parallel){
    stopifnot('cores must be provided if parallel = TRUE' = !missing(cores))
  }

  if (missing(seed)){
    seed = sample(x = 1:1000, size =  1)
  }

  RcppParallel::setThreadOptions(numThreads = 1)
  parafit_fn <- pf_parafit_cpp

  if (use_r){
    parafit_fn <- pf_parafit_r
  }

  if (parallel){
    RcppParallel::setThreadOptions(numThreads = cores)
  }

  message('Running parafit ...')

  # Global ----
  set.seed(seed)

  global_trace <-
    parafit_fn(assoA = associations,
               paraB = parasite_pcoa,
               hostC = t(host_pcoa),
               permutations = permutations,
               seed = seed)

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

      print_i <- getOption('pf_link_print')
      print_i <- ifelse(is.null(print_i), 100, print_i)

      if (i %in% c(1, n_HP_links) | i %% print_i == 0) message(i)

      associations_0 <- associations

      associations_0[
        HP_link_inds[i, 'row'],
        HP_link_inds[i, 'col']
      ] <- 0

      set.seed(seed)

      parafit_fn(assoA = associations_0,
                 paraB = parasite_pcoa,
                 hostC = t(host_pcoa),
                 permutations = permutations,
                 seed = seed)

    }

    HP_link_inds <- which(associations > 0, arr.ind = TRUE)
    HP_link_inds <- HP_link_inds[ order(HP_link_inds[,'row']),]
    n_HP_links <- nrow(HP_link_inds)

    parallel_text <- ifelse(parallel, ' ... in parallel', '')

    message('There are ', n_HP_links, ' links to be tested', parallel_text, sep = '')

    message('Testing links ...')

    link_traces <- sapply(seq(n_HP_links), test_link)

    link_stats <-
      global_trace[,1] - link_traces

    link_pvalues <-
      apply(link_stats, MARGIN = 2, function(x) sum(x >= x[1]) / length(x))

    results$links <-
      data.frame(
        host = rownames(associations)[HP_link_inds[,'row']],
        parasite = colnames(associations)[HP_link_inds[,'col']],
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
#' @param ... passed on to print
#' @export print.pf_parafit
#' @export
print.pf_parafit <- function(x, ...){
  base::cat('Parafit\n\n')
  base::cat('ParaFitGlobal = ', x$global$stat, ', P value = ', base::round(x$global$p, 4), ' (', x$global$permutations, ')\n\n', sep = '')
  base::cat('Individual tests of links\n\n')
  if (is.null(x$links)) {
    cat('Not tested')
  } else {
    base::cat(base::nrow(x$links), 'links tested.\n\n')
    base::print(utils::head(x$links))
  }
}
