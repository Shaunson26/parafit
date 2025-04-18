#' Parafit
#'
#' Run the parafit algorithm to obtain the global test statistics and obtain a test
#' of significance using permutation. Optionally, obtain tests of individual associations.
#'
#' @param HP association matrix
#' @param para.D parasite principal coordinates
#' @param host.D transposed host principal coordinates
#' @param nperm number of permutations
#' @param seed random seed value
#' @param test.links =FALSE,
#' @param seed =NULL,
#' @param correction ="none",
#' @param silent =FALSE,
#' @param parallel = F,
#' @param cores = 1
#'
#' @return list of results, with elements global and links (the latter if test_links = TRUE)
#' @export
pf_parafit_ape <- function(host.D, para.D, HP, nperm=999, test.links=FALSE, seed=NULL, correction="none", silent=FALSE, parallel = F, cores = 1) {

  HP <- as.matrix(HP)
  host.D <- as.matrix(host.D)
  para.D <- as.matrix(para.D)

  HP_dim <- dim(HP)
  n.host = nrow(host.D)
  n.para = nrow(para.D)

  stopifnot('The number of host species in host.D does not match the row number in HP' = HP_dim[1] == n.host,
            'The number of parasite species in para.D does not match the column number in HP' = HP_dim[2] == n.para)

  message("n.hosts = ", n.host, ", n.parasites = ", n.para)

  epsilon <- sqrt(.Machine$double.eps)

  if(is.null(seed)) {
    runif(1)
    seed <- .Random.seed[trunc(runif(1,1,626))]
  }


  # Host distances ----
  message('Host distances')

  host.pc <- ape::pcoa(host.D, correction = correction)

  host.pc_correction_factor = host.pc$correction[2]

  if (host.pc_correction_factor == 1) {

    host.pc_relative_eig_min <- min(host.pc$values[,'Relative_eig'])

    if (host.pc_relative_eig_min < -epsilon) stop('Host D matrix has negative eigenvalues. Rerun with correction="lingoes" or correction="cailliez"')

    sum.host.values.sq <- sum(host.pc$values[,'Eigenvalues']^2)
    host.vectors <- host.pc$vectors

  }

  if (host.pc_correction_factor != 1) {

    sum.host.values.sq <- sum(host.pc$values[,'Corr_eig']^2)
    host.vectors <- host.pc$vectors.cor
  }


  # Parasite distances ----
  message('Parasite distances')

  para.pc <- ape::pcoa(para.D, correction = correction)
  para.pc_correction_factor = para.pc$correction[2]

  if(para.pc_correction_factor == 1) {

    para.pc_relative_eig_min <- min(para.pc$values[,'Relative_eig'])

    if (para.pc_relative_eig_min < -epsilon) stop('Parasite D matrix has negative eigenvalues. Rerun with correction="lingoes" or correction="cailliez"')
    sum.para.values.sq <- sum(para.pc$values[,'Eigenvalues']^2)
    para.vectors <- para.pc$vectors
  }

  if(para.pc_correction_factor != 1) {
    sum.para.values.sq <- sum(para.pc$values[,'Corr_eig']^2)
    para.vectors <- para.pc$vectors.cor
  }

  # Compute and test the global statistics ---
  message('Global statistic')

  p.per.h <- apply(HP, 1, sum)
  h.per.p <- apply(HP, 2, sum)

  mat.4 <- t(host.vectors) %*% HP %*% para.vectors
  global <- sum(mat.4^2)
  global.perm <- rep(NA, nperm)

  set.seed(seed)

  for(i in 1:nperm) {

    HP.perm <- apply(HP, 2, sample)
    mat.4.perm <- t(host.vectors) %*% HP.perm %*% para.vectors
    global.perm[i] <- sum(mat.4.perm^2)

  }

  p.global <- sum(1, global.perm >= global)/(nperm+1)

  out <-
    list(
      ParaFitGlobal = global,
      p.global = p.global,
      global.perm = global.perm,
      para.per.host = p.per.h,
      host.per.para = h.per.p,
      nperm = nperm
    )

  # Test individual H-P links
  if (test.links) {

    message('Testing links', appendLF = FALSE)

    tracemax <- max(sum.host.values.sq, sum.para.values.sq)

    # 1. Create the list of H-P pairs
    list.hp <- which(t(cbind(HP, rep(0, n.host))) > 0)
    host_inds <- (list.hp %/% (n.para + 1)) + 1
    para_inds <- list.hp %% (n.para+1)

    HP.list <-
      data.frame(
        Host = row.names(host.D)[host_inds],
        Parasite = row.names(para.D)[para_inds]
      )

    n.links <- length(list.hp)

    message(' (n = ', n.links, ')', appendLF = FALSE)

    test_link <- function(k){
      message('')
      message('  link: ', k)

      # 2. Compute reference values of link statistics
      HP.k <- HP
      HP.k[HP.list$Host[k], HP.list$Parasite[k]] <- 0

      mat.4.k <- t(host.vectors) %*% HP.k %*% para.vectors
      trace.k <- sum(mat.4.k^2)
      stat1.k <- global - trace.k
      stat2.k <- NA

      den <- tracemax - global
      if(den > epsilon) { stat2.k <- stat1.k/den }


      # 3. Test link statistics by permutations
      set.seed(seed)

      has_stat2 <- !is.na(stat2.k)
      stat1.perm = rep(NA, nperm)
      stat2.perm = rep(NA, nperm)

      for(i in 1:nperm) {

        HP.k.perm <- apply(HP.k, 2, sample)
        mat.4.k.perm <- t(host.vectors) %*% HP.k.perm %*% para.vectors
        trace.k.perm <- sum(mat.4.k.perm^2)
        stat1.perm[i] <- global.perm[i] - trace.k.perm

        if(has_stat2) {
          den <- tracemax - global.perm[i]
          if(den > epsilon) {
            stat2.perm[i] <- stat1.perm[i]/den
          }
        }

      }

      p.stat1.k <- sum(1, stat1.perm > stat1.k)/(nperm+1)
      p.stat2.k <- sum(1, stat2.perm > stat2.k)/(nperm+1)

      # data.frame(
      #   F1.stat = stat1.k,
      #   p.F1 = p.stat1.k,
      #   F2.stat = stat2.k,
      #   p.F2 = p.stat2.k
      # )
      #

      stat1.perm


    }


    if (parallel){
      message(' ... in parallel')
      cl <- parallel::makeCluster(cores)
      res <- pbapply::pblapply(1:n.links, FUN = test_link, cl = cl)
      parallel::stopCluster(cl)
    }

    if (!parallel){
      res <- lapply(1:n.links, test_link)
    }


    HP.list <-
      data.frame(
        HP.list,
        do.call(rbind.data.frame, res)
      )


    out$link.table <- HP.list

    #out$link.table <- do.call(cbind, res)

  }

  out
}
