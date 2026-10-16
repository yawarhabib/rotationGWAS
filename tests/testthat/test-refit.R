# refit_markers(): the fit against its full-form definition, the selection by window, and the input rules.

ex <- rotationGWAS:::toy_trial(n = 120, n_markers = 80, seed = 2)
r <- rotation(ex$means, ex$Omega, ex$relationship, ex$covariates)
s <- genome_scan(r, ex$markers)

# the full-form refit: X = [X0, W z], vg by a one-dimensional REML search with ve = 1, then GLS at (vg, 1)
full_form_refit <- function(r, z) {
  n <- r$n; X <- cbind(r$X, as.numeric(r$W %*% z)); y <- r$y; K <- r$K
  S <- diag(n) - X %*% solve(crossprod(X), t(X)); ee <- eigen(S %*% K %*% S, symmetric = TRUE)
  keep <- seq_len(n - ncol(X)); d <- pmax(ee$values[keep], 0); eta <- as.numeric(crossprod(ee$vectors[, keep], y))
  dev <- function(lv) sum(log(exp(lv) * d + 1)) + sum(eta^2 / (exp(lv) * d + 1))
  vg <- exp(optimize(dev, interval = log(c(1e-4, 1e6)), tol = 1e-12)$minimum)
  Vi <- chol2inv(chol(vg * K + diag(n))); C <- solve(crossprod(X, Vi %*% X))
  b <- as.numeric(C %*% crossprod(X, Vi %*% y))[ncol(X)]; se <- sqrt(C[ncol(X), ncol(X)])
  c(vg = vg, beta = b, se = se, neglog10p = -log10(2 * pt(-abs(b / se), n - ncol(X))))
}

test_that("the refit equals its full-form definition for named markers", {
  mk <- c(ex$truth$causal, "m0003", "m0050")
  rf <- refit_markers(r, ex$markers, which = mk, scan = s)
  expect_equal(rf$marker, mk)
  for (i in seq_along(mk)) {
    ref <- full_form_refit(r, ex$markers[r$ids, mk[i]])
    # two optimisers on a deviance flat to ten decimals near its minimum agree on vg to about 1e-7 relative
    expect_equal(rf$genetic[i], ref[["vg"]], tolerance = 1e-5)
    expect_equal(rf$beta[i], ref[["beta"]], tolerance = 1e-6)
    expect_equal(rf$se[i], ref[["se"]], tolerance = 1e-6)
    expect_equal(rf$neglog10p[i], ref[["neglog10p"]], tolerance = 1e-6)
  }
  expect_equal(attr(rf, "df"), r$n - r$p - 1)
  expect_equal(rf$scan_neglog10p, s$neglog10p[match(mk, s$marker)])
  expect_equal(attr(rf, "threshold"), attr(s, "threshold"))
  expect_equal(rf$significant, rf$neglog10p > attr(s, "threshold"))
})

test_that("the window selects the markers of the rotated scan within it, and the refit decides close to the scan", {
  thr <- attr(s, "threshold")
  rf <- refit_markers(r, ex$markers, scan = s, window = 1.5)
  near <- s$marker[!is.na(s$neglog10p) & abs(s$neglog10p - thr) <= 1.5]
  expect_setequal(rf$marker, near)
  expect_equal(attr(rf, "window"), 1.5)
  expect_true(all(abs(rf$neglog10p - rf$scan_neglog10p) < 1))     # the refit and the scan differ by much less than a unit here
  expect_output(print(rf), "residual variance fixed at one")
  expect_message(refit_markers(r, ex$markers, scan = s, window = 0), "no marker within")
})

test_that("the input rules hold and a marker without variance gives NA", {
  expect_error(refit_markers(r, ex$markers), "name the markers")
  expect_error(refit_markers(r, ex$markers, which = "nosuch"), "not in markers")
  expect_error(refit_markers(r, ex$markers, scan = ex$means), "genome_scan")
  Z <- ex$markers; Z[, 2] <- 1
  rf <- refit_markers(r, Z, which = c("m0002", "m0003"), threshold = 3)
  expect_true(is.na(rf$neglog10p[1])); expect_false(rf$significant[1])
  expect_true(is.finite(rf$neglog10p[2]))
  rf2 <- refit_markers(r, ex$markers, which = 3)
  expect_equal(rf2$marker, colnames(ex$markers)[3]); expect_true(is.na(rf2$significant))
})
