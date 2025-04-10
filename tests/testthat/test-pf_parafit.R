test_that("parafit error checking works", {

  expect_error(
    pf_parafit(host_pcoa = pf_pcoa(gopher.D)),
    regexp = 'host_pcoa must be a matrix'
  )

  expect_error(
    pf_parafit(host_pcoa = pf_pcoa(gopher.D)$vectors,
               parasite_pcoa = pf_pcoa(lice.D)),
    regexp = 'parasite_pcoa must be a matrix'
  )

  expect_error(
    pf_parafit(host_pcoa = pf_pcoa(gopher.D)$vectors,
               parasite_pcoa = pf_pcoa(lice.D)$vectors,
               associations = as.data.frame(gopher.lice.links),
               permutations = 10,
               test_links = FALSE,
               parallel = F,
               seed = 1010),
    regexp = 'associations must be a matrix'
  )


  expect_error(
    pf_parafit(host_pcoa = pf_pcoa(gopher.D)$vectors[1:3,],
               parasite_pcoa = pf_pcoa(lice.D)$vectors,
               associations = gopher.lice.links),
    regexp = 'host_pcoa column number must match associations row number'
  )

  expect_error(
    pf_parafit(host_pcoa = pf_pcoa(gopher.D)$vectors,
               parasite_pcoa = pf_pcoa(lice.D)$vectors[1:3,],
               associations = gopher.lice.links),
    regexp = 'parasite_pcoa column number must match associations column number'
  )

  expect_error(
    pf_parafit(host_pcoa = pf_pcoa(gopher.D)$vectors,
               parasite_pcoa = pf_pcoa(lice.D)$vectors,
               associations = gopher.lice.links,
               permutations = 10,
               test_links = FALSE,
               parallel = T,
               seed = 1010),
    regexp = 'cores must be provided if parallel = TRUE'
  )



})

test_that("parafit single global works", {

  res <-
    suppressMessages(
      pf_parafit(host_pcoa = pf_pcoa(gopher.D)$vectors,
                 parasite_pcoa = pf_pcoa(lice.D)$vectors,
                 associations = gopher.lice.links,
                 permutations = 10,
                 test_links = FALSE,
                 parallel = F,
                 seed = 1010)
    )

  expect_true(all.equal(res$global$stat, 0.01389872, tolerance = 0.000001))
  expect_true(all.equal(res$global$p, 0.1, tolerance = 0.1))
  expect_true(res$global$permutations == 10)
  expect_true(
    all.equal(
      res$global$stat_perm[1:5],
      c(0.013898715, 0.007029701, 0.007882358, 0.008075136, 0.006659498),
      tolerance = 0.0000001)
  )
})

test_that("parafit single link works", {

  res <-
    suppressMessages(
      pf_parafit(host_pcoa = pf_pcoa(gopher.D)$vectors,
                 parasite_pcoa = pf_pcoa(lice.D)$vectors,
                 associations = gopher.lice.links,
                 permutations = 10,
                 test_links = TRUE,
                 parallel = F,
                 seed = 1010)
    )

  expect_true(all.equal(res$global$stat, 0.01389872, tolerance = 0.000001))
  expect_true(all.equal(res$global$p, 0.1, tolerance = 0.1))
  expect_true(res$global$permutations == 10)
  expect_true(
    all.equal(
      res$global$stat_perm[1:5],
      c(0.013898715, 0.007029701, 0.007882358, 0.008075136, 0.006659498),
      tolerance = 0.0000001)
  )

  expect_true(!is.null(res$links))
  expect_true(inherits(res$links, 'data.frame'))
  expect_true(nrow(res$links) == 17)
  expect_true(ncol(res$links) == 4)
  expect_equal(names(res$links), c("host","parasite","stat_1","p_1"))
  expect_true(
    all.equal(
      res$links$stat_1,
      c(0.0009312199, 0.0011611181, 0.0010501927, 0.0006711532, 0.0017178115,
        0.0010412160, 0.0016337474, 0.0019256509, 0.0015769140, 0.0012866890,
        0.0007419528, 0.0013933007, 0.0015419676, 0.0007215688, 0.0004458883,
        0.0006382121, 0.0008493063),
      tolerance = 0.0000001)
  )

  expect_true(
    all.equal(
      res$links$p_1,
      c(0.2, 0.1, 0.2, 0.5, 0.1, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.3, 0.4, 0.2, 0.2),
      tolerance = 0.1)
  )


})

test_that("parafit parallel global works", {

  res <-
    suppressMessages(
      pf_parafit(host_pcoa = pf_pcoa(gopher.D)$vectors,
                 parasite_pcoa = pf_pcoa(lice.D)$vectors,
                 associations = gopher.lice.links,
                 permutations = 10,
                 test_links = FALSE,
                 parallel = T,
                 cores = 4,
                 seed = 1010)
    )

  expect_true(all.equal(res$global$stat, 0.01389872, tolerance = 0.000001))
  expect_true(all.equal(res$global$p, 0.1, tolerance = 0.1))
  expect_true(res$global$permutations == 10)
  expect_true(
    all.equal(
      res$global$stat_perm,
      c(0.013898715, 0.007029701, 0.007882358, 0.008075136, 0.006659498,
        0.007376169, 0.008334728, 0.010105484, 0.007409536, 0.007762936),
      tolerance = 0.0000001)
  )
})

test_that("parafit parallel link works", {

  res <-
    suppressMessages(
      pf_parafit(host_pcoa = pf_pcoa(gopher.D)$vectors,
                 parasite_pcoa = pf_pcoa(lice.D)$vectors,
                 associations = gopher.lice.links,
                 permutations = 10,
                 test_links = TRUE,
                 parallel = T,
                 cores = 4,
                 seed = 1010)
    )

  expect_true(all.equal(res$global$stat, 0.01389872, tolerance = 0.000001))
  expect_true(all.equal(res$global$p, 0.1, tolerance = 0.1))
  expect_true(res$global$permutations == 10)
  expect_true(
    all.equal(
      res$global$stat_perm,
      c(0.013898715, 0.007029701, 0.007882358, 0.008075136, 0.006659498,
        0.007376169, 0.008334728, 0.010105484, 0.007409536, 0.007762936),
      tolerance = 0.0000001)
  )

  expect_true(!is.null(res$links))
  expect_true(inherits(res$links, 'data.frame'))
  expect_true(nrow(res$links) == 17)
  expect_true(ncol(res$links) == 4)
  expect_equal(names(res$links), c("host","parasite","stat_1","p_1"))
  expect_true(
    all.equal(
      res$links$stat_1,
      # as above
      c(0.0009312199, 0.0011611181, 0.0010501927, 0.0006711532, 0.0017178115,
        0.0010412160, 0.0016337474, 0.0019256509, 0.0015769140, 0.0012866890,
        0.0007419528, 0.0013933007, 0.0015419676, 0.0007215688, 0.0004458883,
        0.0006382121, 0.0008493063),
      tolerance = 0.0000001)
  )

  expect_true(
    all.equal(
      res$links$p_1,
      c(0.2, 0.1, 0.2, 0.5, 0.1, 0.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.3, 0.4, 0.2, 0.2),
      tolerance = 0.1)
  )


})

test_that("simplier results across single and parallel", {



  single <-
    pf_parafit_cpp(
      assoA = gopher.lice.links,
      pf_pcoa(lice.D)$vectors,
      t(pf_pcoa(gopher.D)$vectors),
      permutations = 999,
      seed = 1010,
      verbose = F,
      print_n = 0
    )

  parallel_1 <-
    pf_parafit_parallel_cpp(
      assoA = gopher.lice.links,
      pf_pcoa(lice.D)$vectors,
      t(pf_pcoa(gopher.D)$vectors),
      permutations = 999,
      seed = 1010,
      verbose = F,
      print_n = 0
    )

  parallel_2 <-
    pf_parafit_parallel_cpp(
      assoA = gopher.lice.links,
      pf_pcoa(lice.D)$vectors,
      t(pf_pcoa(gopher.D)$vectors),
      permutations = 999,
      seed = 1010,
      verbose = F,
      print_n = 0
    )

  res <-
    cbind(
      single, parallel_1, parallel_2
    )


  expect_true(all.equal(res[,1], res[,2]))
  expect_true(all.equal(res[,1], res[,3]))
  expect_true(all.equal(res[,2], res[,3]))

})
