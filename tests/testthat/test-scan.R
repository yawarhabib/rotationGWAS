# genome_scan(): the P3D test against its full-form definition, the identities between schemes, and the input rules.

ex <- rotationGWAS:::toy_trial(n = 40, n_markers = 60, seed = 2)
r <- rotation(ex$means, ex$Omega, ex$relationship, ex$covariates)

# the P3D statistic written out in full: V = vg K + ve I at the null fit, P = V^-1 - V^-1 X0 (X0' V^-1 X0)^-1 X0' V^-1,
# b = z'Py / z'Pz, s2 = (y'Py - b z'Py) / (n - p - 1), se = sqrt(s2 / z'Pz); the package uses the spectral form instead
full_form_scan <- function(y, X0, K, Zt, vg, ve) {
  n <- length(y); V <- vg * K + diag(ve, n); Vi <- solve(V)
  P <- Vi - Vi %*% X0 %*% solve(crossprod(X0, Vi %*% X0), crossprod(X0, Vi)); P <- (P + t(P)) / 2
  df <- n - ncol(X0) - 1
  zPz <- colSums(Zt * (P %*% Zt)); zPy <- as.numeric(crossprod(Zt, P %*% y)); yPy <- as.numeric(crossprod(y, P %*% y))
  b <- zPy / zPz; s2 <- (yPy - b * zPy) / df; se <- sqrt(s2 / zPz)
  data.frame(beta = b, se = se, neglog10p = -log10(2 * pt(-abs(b / se), df)))
}

test_that("the rotated scan equals the full-form P3D test at the fitted variance components", {
  s <- genome_scan(r, ex$markers)
  ref <- full_form_scan(r$y, r$X, r$K, r$W %*% ex$markers[r$ids, ], r$fit$free$vg, r$fit$free$ve)
  expect_equal(s$beta, ref$beta, tolerance = 1e-8)
  expect_equal(s$se, ref$se, tolerance = 1e-8)
  expect_equal(s$neglog10p, ref$neglog10p, tolerance = 1e-8)
  expect_equal(attr(s, "df"), r$n - r$p - 1)
  expect_equal(attr(s, "threshold"), -log10(0.05 / sum(!is.na(s$neglog10p))))
})

test_that("the approximate schemes equal their full-form tests at their own null fits", {
  s <- genome_scan(r, ex$markers, scheme = c("inverse_variance", "unweighted"))
  f <- attr(s, "fits")
  X0 <- cbind(1, r$covariates); Z <- ex$markers[r$ids, ]
  w <- 1 / sqrt(diag(r$Omega))
  ref_iv <- full_form_scan(r$means * w, X0 * w, r$relationship * tcrossprod(w), Z * w, f$vg[f$scheme == "inverse_variance"], f$ve[f$scheme == "inverse_variance"])
  ref_un <- full_form_scan(r$means, X0, r$relationship, Z, f$vg[f$scheme == "unweighted"], f$ve[f$scheme == "unweighted"])
  expect_equal(s$neglog10p[s$scheme == "inverse_variance"], ref_iv$neglog10p, tolerance = 1e-8)
  expect_equal(s$neglog10p[s$scheme == "unweighted"], ref_un$neglog10p, tolerance = 1e-8)
  expect_equal(s$beta[s$scheme == "unweighted"], ref_un$beta, tolerance = 1e-8)
})

test_that("with a diagonal Omega the rotation and inverse-variance weighting coincide, and with Omega = cI all three do", {
  Od <- diag(diag(ex$Omega)); dimnames(Od) <- dimnames(ex$Omega)
  rd <- rotation(ex$means, Od, ex$relationship, ex$covariates)
  sd <- genome_scan(rd, ex$markers, scheme = c("rotated", "inverse_variance"))
  expect_equal(sd$neglog10p[sd$scheme == "rotated"], sd$neglog10p[sd$scheme == "inverse_variance"], tolerance = 1e-8)
  Oc <- diag(mean(diag(ex$Omega)), nrow(ex$Omega)); dimnames(Oc) <- dimnames(ex$Omega)
  rc <- rotation(ex$means, Oc, ex$relationship, ex$covariates)
  sc <- genome_scan(rc, ex$markers, scheme = c("rotated", "inverse_variance", "unweighted"))
  # each scheme fits its own model without a marker, so the fitted ratios agree to the optimiser's tolerance, about 1e-7
  lp <- split(sc$neglog10p, sc$scheme)
  expect_equal(lp$rotated, lp$inverse_variance, tolerance = 1e-6)
  expect_equal(lp$rotated, lp$unweighted, tolerance = 1e-6)
  b <- split(sc$beta, sc$scheme)
  expect_equal(b$rotated, b$unweighted, tolerance = 1e-6)   # means and markers are scaled alike, so the effect per dose is unchanged
})

test_that("a marker without variance gets NA and does not count in the threshold; accession order does not matter", {
  Z <- ex$markers; Z[, 5] <- 1
  s <- genome_scan(r, Z)
  expect_true(is.na(s$neglog10p[5])); expect_false(s$significant[5])
  expect_equal(attr(s, "fits")$n_tested, ncol(Z) - 1)
  expect_equal(attr(s, "threshold"), -log10(0.05 / (ncol(Z) - 1)))
  perm <- sample(nrow(ex$markers))
  s2 <- genome_scan(r, ex$markers[perm, ])
  s1 <- genome_scan(r, ex$markers)
  expect_equal(s1$neglog10p, s2$neglog10p)
})

test_that("in a larger toy trial the causal marker is the strongest under the rotation and the summary reports the schemes", {
  ex2 <- rotationGWAS:::toy_trial(n = 200, n_markers = 100, seed = 2)
  r2 <- rotation(ex2$means, ex2$Omega, ex2$relationship, ex2$covariates)
  s <- genome_scan(r2, ex2$markers, scheme = c("rotated", "unweighted"))
  rot <- s[s$scheme == "rotated", ]
  expect_equal(rot$marker[which.max(rot$neglog10p)], ex2$truth$causal)
  expect_true(rot$significant[rot$marker == ex2$truth$causal])
  expect_equal(rot$beta[rot$marker == ex2$truth$causal], ex2$truth$effect, tolerance = 0.3)
  sm <- summary(s)
  expect_equal(sm$fits$scheme, c("rotated", "unweighted"))
  expect_true(all(sm$fits$n_tested == ncol(ex2$markers)))
  expect_true(is.finite(sm$fits$inflation[1]))
  expect_output(print(sm), "Genome scan of 200 accessions")
})

test_that("the input rules hold", {
  expect_error(genome_scan(r, ex$markers, scheme = "rotation"), "scheme must be")
  expect_error(genome_scan(r, ex$markers, alpha = 1), "alpha")
  Zna <- ex$markers; Zna[1, 1] <- NA
  expect_error(genome_scan(r, Zna), "NA")
  expect_error(genome_scan(r, ex$markers[-1, ]), "lack")
  expect_error(genome_scan(ex$means, ex$markers), "rotation")
})
