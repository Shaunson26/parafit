test_that("parafit global works", {

  res <-
    suppressMessages(
      pf_parafit(host_pcoa = t(pf_pcoa(gopher.D)),
                 parasite_pcoa = pf_pcoa(lice.D),
                 associations = gopher.lice.links,
                 permutations = 10,
                 test_links = FALSE,
                 seed = 1010)
    )

  expect_true(all.equal(res$global$stat, 0.01389872, tolerance = 0.000001))
  expect_true(all.equal(res$global$p, 0.1, tolerance = 0.001))
  expect_true(res$global$permutations == 10)
  expect_true(
    all.equal(
      res$global$stat_perm[1:5],
      c(0.013898715, 0.008502097, 0.009143968, 0.008501108, 0.007330396),
      tolerance = 0.0000001)
  )
})

test_that("parafit link works", {

  res <-
    suppressMessages(
      pf_parafit(host_pcoa = t(pf_pcoa(gopher.D)),
                 parasite_pcoa = pf_pcoa(lice.D),
                 associations = gopher.lice.links,
                 permutations = 10,
                 test_links = TRUE,
                 seed = 1010)
    )

  expect_true(all.equal(res$global$stat, 0.01389872, tolerance = 0.000001))
  expect_true(all.equal(res$global$p, 0.1, tolerance = 0.001))
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

# test_that("parafit link parallel works", {
#
#   pbapply::pboptions(type = 'none')
#
#   res <-
#     suppressMessages(
#       pf_parafit(host_pcoa = t(pf_pcoa(gopher.D)),
#                  parasite_pcoa = pf_pcoa(lice.D),
#                  associations = gopher.lice.links,
#                  permutations = 99,
#                  test_links = TRUE,
#                  parallel = TRUE,
#                  cores = 2,
#                  seed = 1010)
#     )
#
#   pbapply::pboptions(type = 'text')
#
#   expect_true(all.equal(res$global$stat, 0.01389872, tolerance = 0.000001))
#   expect_true(all.equal(res$global$p, 0.01010101, tolerance = 0.00001))
#   expect_true(res$global$permutations == 99)
#   expect_true(
#     all.equal(
#       res$global$stat_perm[1:5],
#       c(0.013898715, 0.008502097, 0.009143968, 0.008501108, 0.007330396),
#       tolerance = 0.0000001)
#   )
#
#   expect_true(!is.null(res$links))
#   expect_true(inherits(res$links, 'data.frame'))
#   expect_true(nrow(res$links) == 17)
#   expect_true(ncol(res$links) == 4)
#   expect_equal(names(res$links), c("host","parasite","stat_1","p_1"))
#   expect_true(
#     all.equal(
#       res$links$stat_1[1:5],
#       c(0.0009312199, 0.0011611181, 0.0010501927, 0.0006711532, 0.0017178115),
#       tolerance = 0.0000001)
#   )
#
#   expect_true(
#     all.equal(
#       res$links$p_1[1:5],
#       c(0.03030303, 0.15151515, 0.03030303, 0.10101010, 0.01010101),
#       tolerance = 0.000001)
#   )
#
#
# })
