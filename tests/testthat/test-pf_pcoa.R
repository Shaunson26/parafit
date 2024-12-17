test_that("pf_pcoa works", {
  res <- pf_pcoa(gopher.D)
  expect_true(ncol(res) == 14)
  expect_true(nrow(res) == 15)

  expect_true(
    all.equal(
      res[1:5,1],
      c(0.06099355, 0.05965239, -0.10934852, -0.08247512, -0.11000241),
      tolerance = 0.0000001
    )
  )
})
