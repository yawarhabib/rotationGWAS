ex <- rotationGWAS:::toy_trial(n = 50, n_markers = 200, seed = 4)

test_that("the rotation matrix whitens Omega and is the symmetric root", {
  r <- rotation(ex$means, ex$Omega, ex$relationship, ex$covariates)
  expect_lt(max(abs(r$W %*% ex$Omega %*% r$W - diag(r$n))), 1e-8)
  expect_lt(max(abs(r$W - t(r$W))), 1e-10)
  expect_lt(max(abs(r$W %*% r$W - solve(ex$Omega))), 1e-6)
  expect_equal(r$X[, 1], as.numeric(r$W %*% rep(1, r$n)))
  expect_equal(r$p, 3); expect_equal(colnames(r$X), c("intercept", "PC1", "PC2"))
})

test_that("inputs are matched by name: permuting any input leaves every result unchanged", {
  r <- rotation(ex$means, ex$Omega, ex$relationship, ex$covariates)
  set.seed(9); o1 <- sample(r$n); o2 <- sample(r$n); o3 <- sample(r$n)
  r2 <- rotation(ex$means, ex$Omega[o1, o1], ex$relationship[o2, o2], ex$covariates[o3, ])
  expect_equal(r2$ids, r$ids); expect_equal(r2$W, r$W); expect_equal(r2$y, r$y); expect_equal(r2$K, r$K)
  expect_equal(r2$fit$free$vg, r$fit$free$vg); expect_equal(r2$fit$one$vg, r$fit$one$vg)
  o4 <- sample(r$n)                                   # the order of the means sets the order of everything
  r3 <- rotation(ex$means[o4], ex$Omega, ex$relationship, ex$covariates)
  expect_equal(r3$ids, r$ids[o4]); expect_equal(r3$fit$free$deviance, r$fit$free$deviance)
})

test_that("rotated least squares equals generalised least squares with weight solve(Omega)", {
  r <- rotation(ex$means, ex$Omega, ex$relationship, ex$covariates)
  z <- ex$markers[r$ids, 17]
  Xr <- cbind(r$X, r$W %*% z); b_rot <- lm.fit(Xr, r$y)$coefficients
  Xu <- cbind(1, ex$covariates[r$ids, ], z); Oi <- solve(ex$Omega)
  b_gls <- solve(crossprod(Xu, Oi %*% Xu), crossprod(Xu, Oi %*% ex$means[r$ids]))
  expect_lt(max(abs(as.numeric(b_rot) - as.numeric(b_gls))), 1e-8)
})

test_that("the guards refuse what the rotation cannot take", {
  m <- ex$means; O <- ex$Omega; A <- ex$relationship; Q <- ex$covariates
  expect_error(rotation(unname(m), O, A, Q), "named")
  expect_error(rotation(m, O[-1, -1], A, Q), "lacks")
  expect_error(rotation(m, O, A[-(1:2), -(1:2)], Q), "lacks")
  expect_error(rotation(m, O, A, Q[-1, ]), "lack")
  expect_error(rotation(m, O, A, cbind(Q, one = 1)), "constant")
  expect_error(rotation(m, O, A, cbind(Q, PC1b = Q[, 1])), "collinear")
  On <- O; On[1, 2] <- On[1, 2] + 1; expect_error(rotation(m, On, A, Q), "symmetric")
  Oneg <- O; Oneg[1, 1] <- -1; expect_error(rotation(m, Oneg, A, Q), "positive definite")
  expect_error(rotation(m, O, A, Q, max_condition = 1), "ill-conditioned")
  Osing <- O; Osing[2, ] <- Osing[1, ]; Osing[, 2] <- Osing[, 1]; expect_error(rotation(m, Osing, A, Q))
  m2 <- m; m2[3] <- NA; expect_error(rotation(m2, O, A, Q), "NA")
  Aneg <- A; Aneg[1, 1] <- -5; expect_error(rotation(m, O, Aneg, Q), "semi-definite")
})

test_that("extra accessions in Omega or the relationship matrix are dropped with a message", {
  m <- ex$means[1:40]
  expect_message(r <- rotation(m, ex$Omega, ex$relationship, ex$covariates), "dropped")
  expect_equal(r$n, 40); expect_equal(r$ids, names(m))
})

test_that("the model without a marker is fitted both ways and the summary of Omega is right", {
  r <- rotation(ex$means, ex$Omega, ex$relationship, ex$covariates)
  expect_true(r$fit$free$ve > 0); expect_equal(r$fit$one$ve, 1); expect_true(is.finite(r$fit$one$deviance))
  expect_equal(r$omega[["asv"]], asv(ex$Omega)); expect_equal(r$omega[["diag_mean"]], mean(diag(ex$Omega)))
  expect_output(print(r), "Rotation of 50 adjusted means")
})

test_that("asv is the centred trace over n - 1 and is free of a constant added to every element", {
  M <- ex$Omega; n <- nrow(M); P <- diag(n) - 1 / n
  expect_equal(asv(M), sum(diag(P %*% M)) / (n - 1))
  expect_equal(asv(M + 7), asv(M))
})

test_that("a diagonal Omega is accepted as inverse-variance weighting, and the functions that need the full matrix refuse it", {
  D <- diag(diag(ex$Omega)); dimnames(D) <- dimnames(ex$Omega)
  expect_message(rd <- rotation(ex$means, D, ex$relationship, ex$covariates), "inverse-variance")
  expect_true(rd$diagonal)
  expect_output(print(rd), "inverse-variance weighting")
  r <- rotation(ex$means, ex$Omega, ex$relationship, ex$covariates)
  expect_false(r$diagonal)
  s_full <- genome_scan(r, ex$markers, scheme = "inverse_variance")       # the inverse-variance scheme of the full matrix ...
  s_diag <- genome_scan(rd, ex$markers, scheme = "rotated")               # ... is the rotation of its diagonal
  expect_equal(s_diag$t, s_full$t, tolerance = 1e-6)
  expect_error(efficiency(rd, ex$markers), "standard errors only")
  expect_error(estimation_share(rd), "standard errors only")
  expect_error(specification_test(rd), "standard errors only")
  expect_error(independent_effect(rd), "standard errors only")
})

test_that("attaching the package prints its one line, and it can be silenced", {
  expect_message(rotationGWAS:::.onAttach("", "rotationGWAS"), "the intercept is not a column of ones")
  expect_silent(suppressPackageStartupMessages(rotationGWAS:::.onAttach("", "rotationGWAS")))
})
