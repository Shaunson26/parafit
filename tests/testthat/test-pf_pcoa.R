test_that("pf_pcoa works", {

  res <- pf_pcoa(gopher.D)

  expect_true(ncol(res$vectors) == 14)
  expect_true(nrow(res$vectors) == 15)

  expect_true(
    all.equal(
      unname(res$vectors[1:5,1]),
      c(0.06099355, 0.05965239, -0.10934852, -0.08247512, -0.11000241),
      tolerance = 0.0000001
    )
  )

  expect_true(all.equal(res$sum_eig_sq, 0.009184943))

  # cpp
  res_cpp <- pf_pcoa_cpp(gopher.D)
  res_cpp <- wrangle_pf_pcoa_cpp(res_cpp)

  expect_true(
    all.equal(
      unname(res$pco$values[,c('Eigenvalues', 'Relative_eig','Cumul_eig')]),
      unname(res_cpp$values[,1:3])
    )
  )

  expect_true(all.equal(unname(abs(res$pco$vectors)), abs(res_cpp$vectors)))


  # Sponge
  expect_error(pf_pcoa(sponge_dist), regexp = "Matrix has negative eigenvalues")

  res <- pf_pcoa(sponge_dist, correction = 'lingoes')
  expect_true(ncol(res$vectors) == 164)
  expect_true(nrow(res$vectors) == 166)

  expect_true(
    all.equal(
      unname(res$vectors[1:5,1]),
      c(-0.01778891, -0.01917971, -0.01759496, -0.01298697, -0.01176871),
      tolerance = 0.000001
    )
  )

  expect_true(all.equal(res$sum_eig_sq, 3.992881, tolerance = 0.00001))

  # cpp
  expect_error(pf_pcoa_cpp(sponge_dist), regexp = "Input must be a matrix")
  res_cpp <- pf_pcoa_cpp(as.matrix(sponge_dist))
  res_cpp <- wrangle_pf_pcoa_cpp(res_cpp)

  expect_true(
    all.equal(
      unname(res$pco$values[1:82,c('Eigenvalues')]),
      unname(res_cpp$values[,1]),
      tolerance = 0.00000001
    )
  )

  expect_true(
    all.equal(
      unname(abs(res$pco$vectors[,])),
      abs(res_cpp$vectors[,1:73])
    )
  )

})
