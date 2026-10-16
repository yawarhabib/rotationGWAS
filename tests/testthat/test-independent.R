# independent_effect(): the two-kernel engine against the one-kernel engine and the full-form deviance, the fit's
# optimality, the scan under the extended model, the boundary identity with the unweighted analysis, and a trial
# with a real independent effect.

ex <- rotationGWAS:::toy_trial(n = 120, n_markers = 80, seed = 2)
r <- rotation(ex$means, ex$Omega, ex$relationship, ex$covariates)
Kind <- r$W %*% r$W

full_form_dev <- function(y, X, V) {            # log|V| + log|X'V^-1 X| + y'Py, the paper's script 04
  Vi <- solve(V); XtViX <- crossprod(X, Vi %*% X)
  P <- Vi - Vi %*% X %*% solve(XtViX, crossprod(X, Vi))
  as.numeric(determinant(V, logarithm = TRUE)$modulus + determinant(XtViX, logarithm = TRUE)$modulus + crossprod(y, P %*% y))
}

test_that("the two-kernel deviance equals the one-kernel deviance when the second component is zero, and the full form", {
  s2 <- rotationGWAS:::reml_setup2(r$y, r$X, r$K, Kind)
  for (vv in list(c(3, 0.8), c(0.5, 2))) {
    expect_equal(rotationGWAS:::reml_deviance2(s2, vv[1], 0, vv[2], full = FALSE), rotationGWAS:::reml_deviance(r$setup, vv[1], vv[2], full = FALSE), tolerance = 1e-10)
    expect_equal(rotationGWAS:::reml_deviance2(s2, vv[1], 0.7, vv[2]), full_form_dev(r$y, r$X, vv[1] * r$K + 0.7 * Kind + diag(vv[2], r$n)), tolerance = 1e-8)
  }
})

test_that("the extended fit nests the additive one and is a local minimum of the full-form deviance", {
  ie <- independent_effect(r)
  expect_true(ie$fit$deviance <= r$fit$free$deviance + 1e-6)
  expect_equal(ie$lr, max(r$fit$free$deviance - ie$fit$deviance, 0))
  expect_equal(ie$delta_aic, 2 - ie$lr)
  f <- ie$fit
  d0 <- full_form_dev(r$y, r$X, f$vg * r$K + f$vi * Kind + diag(f$ve, r$n))
  expect_equal(d0, f$deviance, tolerance = 1e-8)
  for (k in 1:6) {                                   # no nearby point in the interior does better
    g <- c(f$vg, f$vi, f$ve) * exp(stats::rnorm(3, 0, 0.05)); g[c(f$vg, f$vi, f$ve) == 0] <- 0
    expect_true(full_form_dev(r$y, r$X, g[1] * r$K + g[2] * Kind + diag(g[3], r$n)) >= d0 - 1e-6)
  }
  expect_equal(ie$fit_one$ve, 1)
  expect_output(print(ie), "independent across accessions")
})

test_that("the scan under the extended model is the full-form P3D test at its fitted components", {
  ie <- independent_effect(r)
  s <- genome_scan(r, ex$markers, model = ie)
  f <- ie$fit; V <- f$vg * r$K + f$vi * Kind + diag(f$ve, r$n); Vi <- solve(V)
  P <- Vi - Vi %*% r$X %*% solve(crossprod(r$X, Vi %*% r$X), crossprod(r$X, Vi)); P <- (P + t(P)) / 2
  Zt <- r$W %*% ex$markers[r$ids, ]; df <- r$n - ncol(r$X) - 1
  zPz <- colSums(Zt * (P %*% Zt)); zPy <- as.numeric(crossprod(Zt, P %*% r$y)); yPy <- as.numeric(crossprod(r$y, P %*% r$y))
  b <- zPy / zPz; se <- sqrt((yPy - b * zPy) / df / zPz)
  expect_equal(s$beta, unname(b), tolerance = 1e-8); expect_equal(s$se, unname(se), tolerance = 1e-8)
  expect_equal(attr(s, "model"), "additive + independent effect")
  expect_equal(attr(s, "fits")$vi, f$vi)
  expect_error(genome_scan(r, ex$markers, scheme = "unweighted", model = ie), "rotation only")
  expect_error(genome_scan(r, ex$markers, model = r), "independent_effect")
})

test_that("at a residual variance of zero the extended model is the unweighted analysis, and the scan warns", {
  ie <- independent_effect(r)
  ie0 <- ie; ie0$fit$ve <- 0; ie0$fit$boundary <- "ve_zero"
  expect_warning(s0 <- genome_scan(r, ex$markers, model = ie0), "unweighted analysis")
  # the unweighted P3D test at the same ratio of the two components, in the full form on the unrotated data
  X0 <- cbind(1, r$covariates); V <- ie0$fit$vg * r$relationship + diag(ie0$fit$vi, r$n); Vi <- solve(V)
  P <- Vi - Vi %*% X0 %*% solve(crossprod(X0, Vi %*% X0), crossprod(X0, Vi))
  Z <- ex$markers[r$ids, ]; df <- r$n - ncol(X0) - 1
  zPz <- colSums(Z * (P %*% Z)); zPy <- as.numeric(crossprod(Z, P %*% r$means)); yPy <- as.numeric(crossprod(r$means, P %*% r$means))
  b <- zPy / zPz; se <- sqrt((yPy - b * zPy) / df / zPz)
  expect_equal(s0$t, unname(b / se), tolerance = 1e-6)
})

test_that("a real independent effect is found, and it brings the residual variance of the rotated model back towards one", {
  set.seed(7)
  ex2 <- rotationGWAS:::toy_trial(n = 150, n_markers = 60, seed = 3)
  means2 <- ex2$means + stats::rnorm(length(ex2$means), 0, sqrt(6))      # an effect independent across accessions, variance 6
  r2 <- rotation(means2, ex2$Omega, ex2$relationship, ex2$covariates)
  st <- specification_test(r2)
  ie <- independent_effect(r2)
  expect_true(st$residual > 1)                       # the additive model's residual absorbs the missing effect
  expect_true(ie$fit$vi > 0); expect_true(ie$lr > 0)
  expect_true(abs(ie$fit$ve - 1) < abs(st$residual - 1))   # the extended model's residual is nearer its target
  expect_true(ie$p_mixture < 0.5)
})

test_that("the scan under the extended model can hold the components of the fit with the residual variance at one", {
  ie <- independent_effect(r)
  s1 <- genome_scan(r, ex$markers, model = ie, fit = "residual_one")
  f <- ie$fit_one; V <- f$vg * r$K + f$vi * Kind + diag(1, r$n); Vi <- solve(V)
  P <- Vi - Vi %*% r$X %*% solve(crossprod(r$X, Vi %*% r$X), crossprod(r$X, Vi)); P <- (P + t(P)) / 2
  Zt <- r$W %*% ex$markers[r$ids, ]; df <- r$n - ncol(r$X) - 1
  zPz <- colSums(Zt * (P %*% Zt)); zPy <- as.numeric(crossprod(Zt, P %*% r$y)); yPy <- as.numeric(crossprod(r$y, P %*% r$y))
  b <- zPy / zPz; se <- sqrt((yPy - b * zPy) / df / zPz)
  expect_equal(s1$beta, unname(b), tolerance = 1e-8); expect_equal(s1$se, unname(se), tolerance = 1e-8)
  expect_equal(attr(s1, "fits")$ve, 1); expect_equal(attr(s1, "fits")$vi, f$vi); expect_match(attr(s1, "model"), "fixed at one")
  ie0 <- ie; ie0$fit$ve <- 0; ie0$fit$boundary <- "ve_zero"                 # the boundary of the free fit does not concern this scan
  expect_no_warning(s0 <- genome_scan(r, ex$markers, model = ie0, fit = "residual_one"))
  expect_equal(s0$t, s1$t)
  expect_error(genome_scan(r, ex$markers, fit = "residual_one"), "give `model`")
})
