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
                 permutations = 999,
                 test_links = FALSE,
                 parallel = F,
                 seed = 1010)
    )

  #saveRDS(res, testthat::test_path('expected-results-data', 'parafit-single-global.rds'))

  res_expected <-
    readRDS(testthat::test_path('expected-results-data', 'parafit-single-global.rds'))

  expect_true(
    all.equal(res$global$stat, res_expected$global$stat)
  )

  expect_true(
    all.equal(res$global$p, res_expected$global$p)
  )

  expect_true(
    all.equal(res$global$stat_perm, res$global$stat_perm)
  )
})

test_that("parafit single link works", {

  res <-
    suppressMessages(
      pf_parafit(host_pcoa = pf_pcoa(gopher.D)$vectors,
                 parasite_pcoa = pf_pcoa(lice.D)$vectors,
                 associations = gopher.lice.links,
                 permutations = 999,
                 test_links = TRUE,
                 parallel = F,
                 seed = 1010)
    )

  #saveRDS(res, testthat::test_path('expected-results-data', 'parafit-single-link.rds'))

  res_expected <-
    readRDS(testthat::test_path('expected-results-data', 'parafit-single-link.rds'))

  expect_true(
    all.equal(res$global$stat, res_expected$global$stat)
  )

  expect_true(
    all.equal(res$global$p, res_expected$global$p)
  )

  expect_true(
    all.equal(res$global$stat_perm, res$global$stat_perm)
  )

  expect_true(
    all.equal(res$links$stat_1, res_expected$links$stat_1)
  )

  expect_true(
    all.equal(res$links$p_1, res_expected$links$p_1)
  )

})

test_that("parafit parallel global works", {

  res <-
    suppressMessages(
      pf_parafit(host_pcoa = pf_pcoa(gopher.D)$vectors,
                 parasite_pcoa = pf_pcoa(lice.D)$vectors,
                 associations = gopher.lice.links,
                 permutations = 999,
                 test_links = FALSE,
                 parallel = T,
                 cores = 4,
                 seed = 1010)
    )

  #saveRDS(res, testthat::test_path('expected-results-data', 'parafit-parallel-global.rds'))

  res_expected <-
    readRDS(testthat::test_path('expected-results-data', 'parafit-parallel-global.rds'))

  expect_true(
    all.equal(res$global$stat, res_expected$global$stat)
  )

  expect_true(
    all.equal(res$global$p, res_expected$global$p)
  )

  expect_true(
    all.equal(res$global$stat_perm, res$global$stat_perm)
  )



})

test_that("parafit parallel link works", {

  res <-
    suppressMessages(
      pf_parafit(host_pcoa = pf_pcoa(gopher.D)$vectors,
                 parasite_pcoa = pf_pcoa(lice.D)$vectors,
                 associations = gopher.lice.links,
                 permutations = 999,
                 test_links = TRUE,
                 parallel = T,
                 cores = 4,
                 seed = 1010)
    )

  #saveRDS(res, testthat::test_path('expected-results-data', 'parafit-parallel-link.rds'))

  res_expected <-
    readRDS(testthat::test_path('expected-results-data', 'parafit-parallel-link.rds'))

  expect_true(
    all.equal(res$global$stat, res_expected$global$stat)
  )

  expect_true(
    all.equal(res$global$p, res_expected$global$p)
  )

  expect_true(
    all.equal(res$global$stat_perm, res$global$stat_perm)
  )

  expect_true(
    all.equal(res$links$stat_1, res_expected$links$stat_1)
  )

  expect_true(
    all.equal(res$links$p_1, res_expected$links$p_1)
  )

})

test_that("same results across single and parallel and r", {

  single <-
    suppressMessages(
      pf_parafit(host_pcoa = pf_pcoa(gopher.D)$vectors,
                 parasite_pcoa = pf_pcoa(lice.D)$vectors,
                 associations = gopher.lice.links,
                 permutations = 9999,
                 test_links = TRUE,
                 parallel = F,
                 cores = 1,
                 seed = 1010)
    )

  single_r <-
    suppressMessages(
      pf_parafit(host_pcoa = pf_pcoa(gopher.D)$vectors,
                 parasite_pcoa = pf_pcoa(lice.D)$vectors,
                 associations = gopher.lice.links,
                 permutations = 9999,
                 test_links = TRUE,
                 seed = 1010,
                 use_r = TRUE)
    )


  parallel_1 <-
    suppressMessages(
      pf_parafit(host_pcoa = pf_pcoa(gopher.D)$vectors,
                 parasite_pcoa = pf_pcoa(lice.D)$vectors,
                 associations = gopher.lice.links,
                 permutations = 9999,
                 test_links = TRUE,
                 parallel = T,
                 cores = 4,
                 seed = 1010)
    )

  parallel_2 <-
    suppressMessages(
      pf_parafit(host_pcoa = pf_pcoa(gopher.D)$vectors,
                 parasite_pcoa = pf_pcoa(lice.D)$vectors,
                 associations = gopher.lice.links,
                 permutations = 9999,
                 test_links = TRUE,
                 parallel = T,
                 cores = 4,
                 seed = 9999)
    )

  ## Single vs R
  expect_true(
    all.equal(single$global$stat, single_r$global$stat)
  )

  expect_true(
    # 0.0006 vs 0.0004
    all.equal(single$global$p, single_r$global$p) == "Mean relative difference: 0.3333333"
  )

  expect_true(
    all.equal(single$global$stat_perm, single_r$global$stat_perm) == "Mean relative difference: 0.1578976"
  )

  expect_true(
    all.equal(single$links$stat_1, single_r$links$stat_1)
  )

  # single_links <- single$links[ order(single$links$p_1), ]
  # single_r_links <- single_r$links[ order(single_r$links$p_1), ]
  # max(abs(single_links$p_1 - single_r_links$p_1)) # 0.3
  zz <- merge.data.frame(single_links, single_r_links, by = c('host', 'parasite'))
  zz[ order(zz$p_1.x), c('host', 'parasite', 'p_1.x', 'p_1.y')]

  ## Single vs parallel
  expect_true(
    all.equal(single$global$stat, parallel_1$global$stat)
  )

  expect_true(
    all.equal(single$global$p, parallel_1$global$p)
  )

  expect_true(
    all.equal(single$global$stat_perm, parallel_1$global$stat_perm)
  )

  expect_true(
    all.equal(single$links$stat_1, parallel_1$links$stat_1)
  )

  expect_true(
    all.equal(single$links$p_1, parallel_1$links$p_1)
  )


  # Parallel different seeds

  expect_true(
    all.equal(parallel_1$global$stat, parallel_2$global$stat)
  )

  expect_true(
    # 0.0006 vs 0.0004
    all.equal(parallel_1$global$p, parallel_2$global$p) == "Mean relative difference: 0.3333333"
  )

  expect_true(
    all.equal(parallel_1$global$stat_perm, parallel_2$global$stat_perm) == "Mean relative difference: 0.1580925"
  )

  expect_true(
    all.equal(parallel_1$links$stat_1, parallel_2$links$stat_1)
  )

  expect_true(
    all.equal(parallel_1$links$p_1, parallel_2$links$p_1) == "Mean relative difference: 0.01839303"
  )

  expect_true(
    all.equal(max(abs(parallel_1$links$p_1 - parallel_2$links$p_1)), 0.00650065, tolerance = 0.000001)
  )


})
