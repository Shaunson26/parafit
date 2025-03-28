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
      c(0.013898715, 0.008502097, 0.009143968, 0.008501108, 0.007330396),
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
      c(0.013898715, 0.008502097, 0.009143968, 0.008501108, 0.007330396),
      tolerance = 0.0000001)
  )

  expect_true(!is.null(res$links))
  expect_true(inherits(res$links, 'data.frame'))
  expect_true(nrow(res$links) == 17)
  expect_true(ncol(res$links) == 4)
  expect_equal(names(res$links), c("host","parasite","stat_1","p_1"))
  expect_true(
    all.equal(
      res$links$stat_1[1:5],
      c(0.0009312199, 0.0011611181, 0.0010501927, 0.0006711532, 0.0017178115),
      tolerance = 0.0000001)
  )

  expect_true(
    all.equal(
      res$links$p_1[1:5],
      c(0.1, 0.4, 0.1, 0.1, 0.1),
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
      res$global$stat_perm[1:5],
      #c(0.013898715, 0.008502097, 0.009143968, 0.008501108, 0.007330396),
      c(0.013898715, 0.007730290, 0.009128391, 0.006877452, 0.008071909),
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
      res$global$stat_perm[1:5],
      #c(0.013898715, 0.008502097, 0.009143968, 0.008501108, 0.007330396),
      c(0.013898715, 0.007730290, 0.009128391, 0.006877452, 0.008071909),
      tolerance = 0.0000001)
  )

  expect_true(!is.null(res$links))
  expect_true(inherits(res$links, 'data.frame'))
  expect_true(nrow(res$links) == 17)
  expect_true(ncol(res$links) == 4)
  expect_equal(names(res$links), c("host","parasite","stat_1","p_1"))
  expect_true(
    all.equal(
      res$links$stat_1[1:5],
      c(0.0009312199, 0.0011611181, 0.0010501927, 0.0006711532, 0.0017178115),
      tolerance = 0.0000001)
  )

  expect_true(
    all.equal(
      res$links$p_1[1:5],
      c(0.1, 0.4, 0.1, 0.1, 0.1),
      tolerance = 0.1)
  )


})
