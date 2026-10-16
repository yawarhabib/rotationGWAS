# The engine against numerical identities on synthetic data. Nothing here depends on any trial.

sim <- function(n = 60, p = 3, m = 400, vg = 4, ve = 1.5, seed = 1, rank_markers = m) {
  set.seed(seed)
  G <- matrix(rbinom(n * rank_markers, 2, 0.3), n, rank_markers)      # dosages 0/1/2
  C <- scale(G, center = TRUE, scale = FALSE)
  K <- tcrossprod(C) / ncol(C); K <- K / mean(diag(K))                  # centred: K 1 = 0, rank <= n - 1
  X <- cbind(1, matrix(rnorm(n * (p - 1)), n, p - 1))
  u <- as.numeric(t(chol(K + diag(1e-10, n))) %*% rnorm(n)) * sqrt(vg)
  y <- as.numeric(X %*% c(10, 1, -1)[seq_len(p)]) + u + rnorm(n, sd = sqrt(ve))
  list(y = y, X = X, K = K)
}

full_form <- function(y, X, K, vg, ve) {                                   # log|V| + log|X'V^-1 X| + y'Py
  V <- vg * K + diag(ve, length(y)); R <- chol(V)
  Vi <- chol2inv(R); XtViX <- crossprod(X, Vi %*% X); Rx <- chol(XtViX)
  P <- Vi - Vi %*% X %*% chol2inv(Rx) %*% crossprod(X, Vi)
  2 * sum(log(diag(R))) + 2 * sum(log(diag(Rx))) + as.numeric(crossprod(y, P %*% y))
}

test_that("the error-contrast deviance equals the full form plus log|X'X|", {
  s <- sim(); st <- reml_setup(s$y, s$X, s$K)
  for (vc in list(c(4, 1.5), c(0.3, 2), c(50, 0.01), c(1e-4, 1))) {
    expect_equal(reml_deviance(st, vc[1], vc[2], full = TRUE), full_form(s$y, s$X, s$K, vc[1], vc[2]), tolerance = 1e-8)
  }
})

test_that("the free fit is the REML optimum: y'Py = n - p and no lower deviance nearby", {
  s <- sim(); st <- reml_setup(s$y, s$X, s$K); f <- reml_fit_free(st)
  expect_equal(f$boundary, "none")
  expect_equal(sum(st$eta^2 / (f$vg * st$d + f$ve)), st$df, tolerance = 1e-8)
  for (fac in c(0.9, 0.99, 1.01, 1.1)) {
    expect_gte(reml_deviance(st, f$vg * fac, f$ve, full = FALSE), f$deviance_ec - 1e-9)
    expect_gte(reml_deviance(st, f$vg, f$ve * fac, full = FALSE), f$deviance_ec - 1e-9)
  }
})

test_that("the free fit agrees with rrBLUP::mixed.solve within its tolerance", {
  skip_if_not_installed("rrBLUP")
  s <- sim(n = 80, m = 600, seed = 3); st <- reml_setup(s$y, s$X, s$K); f <- reml_fit_free(st)
  ms <- rrBLUP::mixed.solve(s$y, X = s$X, Z = diag(length(s$y)), K = s$K)
  expect_equal(f$vg, ms$Vu, tolerance = 2e-3)
  expect_equal(f$ve, ms$Ve, tolerance = 2e-3)
  expect_gte(reml_deviance(st, ms$Vu, ms$Ve, full = FALSE), f$deviance_ec - 1e-8)   # ours is never above the minimum
})

test_that("the fit with the residual fixed at one has zero derivative at an interior optimum", {
  s <- sim(vg = 4, ve = 1, seed = 5); st <- reml_setup(s$y, s$X, s$K); f1 <- reml_fit_fixed(st, 1)
  expect_equal(f1$boundary, "none"); expect_equal(f1$ve, 1)
  h <- 1e-5 * f1$vg
  grad <- (reml_deviance(st, f1$vg + h, 1) - reml_deviance(st, f1$vg - h, 1)) / (2 * h)
  expect_lt(abs(grad) * f1$vg, 1e-5)
  expect_lt(f1$deviance_ec - reml_deviance(st, f1$vg * 1.05, 1, full = FALSE), 0)
})

test_that("a response in the column space of X puts the genetic variance on the boundary", {
  s <- sim(); y0 <- as.numeric(s$X %*% c(1, 2, 3)); st <- reml_setup(y0, s$X, s$K)
  expect_true(all(abs(st$eta) < 1e-10))
  f1 <- reml_fit_fixed(st, 1); expect_equal(f1$boundary, "vg_zero"); expect_equal(f1$vg, 0)
})

test_that("the trait scale leaves the ratio, the test statistic and the boundary flags unchanged", {
  s <- sim(seed = 7); st <- reml_setup(s$y, s$X, s$K); f <- reml_fit_free(st)
  for (c in c(0.01, 100)) {
    st2 <- reml_setup(c * s$y, s$X, s$K); f2 <- reml_fit_free(st2)
    expect_equal(f2$vg, c^2 * f$vg, tolerance = 1e-5); expect_equal(f2$ve, c^2 * f$ve, tolerance = 1e-5)
    expect_equal(f2$lambda, f$lambda, tolerance = 1e-5)
    expect_equal(f2$deviance_ec, f$deviance_ec + st$df * log(c^2), tolerance = 1e-8)
    g <- reml_fit_fixed(st, 1); g2 <- reml_fit_fixed(st2, c^2)
    expect_equal(g2$vg, c^2 * g$vg, tolerance = 1e-5)
    expect_equal(g2$deviance_ec - f2$deviance_ec, g$deviance_ec - f$deviance_ec, tolerance = 1e-8)
  }
})

test_that("a REML estimate of the residual at zero is reported as a boundary, not as a number", {
  s <- sim(seed = 17); st <- reml_setup(s$y, s$X, s$K); f <- reml_fit_free(st)   # this draw puts ve on the boundary
  expect_equal(f$boundary, "ve_zero"); expect_lt(f$ve, 1e-10 * f$vg)
  expect_true(is.finite(f$deviance))
})

test_that("a kernel of rank below n - p is handled without any floor or ones term", {
  s <- sim(n = 60, p = 3, rank_markers = 5, seed = 11); st <- reml_setup(s$y, s$X, s$K)
  expect_equal(sum(st$d > 1e-10 * max(st$d)), 5)
  for (vc in list(c(2, 1), c(0.1, 3))) expect_equal(reml_deviance(st, vc[1], vc[2]), full_form(s$y, s$X, s$K, vc[1], vc[2]), tolerance = 1e-8)
  f <- reml_fit_free(st); expect_true(is.finite(f$deviance)); expect_true(f$ve > 0)
})

test_that("the profile-likelihood interval has the right deviance at its limits and contains the estimate", {
  s <- sim(n = 80, m = 600, vg = 3, ve = 2, seed = 13); st <- reml_setup(s$y, s$X, s$K); f <- reml_fit_free(st)
  ci <- reml_profile_ve(st, f, level = 0.95)
  expect_true(ci$bounded); expect_lt(ci$lower, f$ve); expect_gt(ci$upper, f$ve)
  expect_equal(reml_fit_fixed(st, ci$upper)$deviance_ec - f$deviance_ec, ci$crit, tolerance = 1e-6)
  if (ci$lower > 0) expect_equal(reml_fit_fixed(st, ci$lower)$deviance_ec - f$deviance_ec, ci$crit, tolerance = 1e-6)
  ci90 <- reml_profile_ve(st, f, level = 0.90)
  expect_gte(ci90$lower, ci$lower); expect_lte(ci90$upper, ci$upper)
})

test_that("the scan projection reproduces z'Pz, z'Py and y'Py from the full matrices", {
  s <- sim(seed = 1); st <- reml_setup(s$y, s$X, s$K); f <- reml_fit_free(st)
  expect_equal(f$boundary, "none")                       # an interior ratio, so that the direct matrices are well conditioned
  V <- s$K + diag(f$lambda, length(s$y)); Vi <- solve(V)
  P <- Vi - Vi %*% s$X %*% solve(crossprod(s$X, Vi %*% s$X)) %*% crossprod(s$X, Vi)
  pr <- scan_projection(st, f$lambda)
  Z <- matrix(rbinom(length(s$y) * 7, 2, 0.4), length(s$y), 7); HZ <- pr$H %*% Z
  expect_equal(colSums(HZ^2), diag(crossprod(Z, P %*% Z)), tolerance = 1e-9)
  expect_equal(as.numeric(crossprod(HZ, pr$Hy)), as.numeric(crossprod(Z, P %*% s$y)), tolerance = 1e-9)
  expect_equal(pr$yPy, as.numeric(crossprod(s$y, P %*% s$y)), tolerance = 1e-9)
})

test_that("generalised least squares at given components reduces to ordinary least squares when vg = 0", {
  s <- sim(seed = 19); g <- gls_fixed(s$y, s$X, s$K, vg = 0, ve = 2.5)
  ols <- lm.fit(s$X, s$y)
  expect_equal(g$beta, as.numeric(ols$coefficients), tolerance = 1e-10)
  expect_equal(g$se, sqrt(2.5 * diag(solve(crossprod(s$X)))), tolerance = 1e-10)
})
