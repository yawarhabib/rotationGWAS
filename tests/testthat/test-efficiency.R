ex <- rotationGWAS:::toy_trial(n = 50, n_markers = 120, seed = 6)
r <- rotation(ex$means, ex$Omega, ex$relationship, ex$covariates)

sandwich <- function(z, X0, Sigma, Vwork) {          # the coded form of the paper's scripts, one marker at a time
  M <- solve(Vwork); Xm <- cbind(X0, z); bread <- solve(crossprod(Xm, M %*% Xm))
  (bread %*% crossprod(Xm, M %*% Sigma %*% M %*% Xm) %*% bread)[ncol(Xm), ncol(Xm)]
}

test_that("the vectorised efficiency equals the sandwich form, marker by marker", {
  e <- efficiency(r, ex$markers)
  sg <- attr(e, "sigma2_g"); X0 <- cbind(1, ex$covariates[r$ids, ]); A <- ex$relationship[r$ids, r$ids]; O <- ex$Omega[r$ids, r$ids]
  Sigma <- sg * A + O
  for (j in c(1, 17, 50, 120)) {
    z <- ex$markers[r$ids, j]
    v_rot <- sandwich(z, X0, Sigma, Sigma); v_inv <- sandwich(z, X0, Sigma, sg * A + diag(diag(O))); v_unw <- sandwich(z, X0, Sigma, sg * A + diag(mean(diag(O)), r$n))
    expect_equal(e$se_rotated[j], sqrt(v_rot), tolerance = 1e-10)
    expect_equal(e$inverse_variance[j], sqrt(v_inv / v_rot), tolerance = 1e-10)
    expect_equal(e$unweighted[j], sqrt(v_unw / v_rot), tolerance = 1e-10)
  }
  expect_true(all(e$inverse_variance >= 1 - 1e-12, na.rm = TRUE)); expect_true(all(e$unweighted >= 1 - 1e-12, na.rm = TRUE))
})

test_that("with Omega proportional to the identity the three schemes coincide in the scan, and the efficiency is refused as diagonal", {
  O <- diag(2, r$n); dimnames(O) <- list(r$ids, r$ids)
  expect_message(r0 <- rotation(ex$means, O, ex$relationship, ex$covariates), "inverse-variance")
  s <- genome_scan(r0, ex$markers, scheme = c("rotated", "inverse_variance", "unweighted"))
  t_rot <- s$t[s$scheme == "rotated"]; t_ivw <- s$t[s$scheme == "inverse_variance"]; t_unw <- s$t[s$scheme == "unweighted"]
  expect_equal(t_ivw, t_rot, tolerance = 1e-6); expect_equal(t_unw, t_rot, tolerance = 1e-6)
  expect_error(efficiency(r0, ex$markers), "standard errors only")
})

test_that("the marker coding (0/1/2 or -1/0/1), the row order and the chunk size do not change the factors", {
  e <- efficiency(r, ex$markers)
  e2 <- efficiency(r, ex$markers - 1); expect_equal(e2$inverse_variance, e$inverse_variance); expect_equal(e2$unweighted, e$unweighted)
  set.seed(3); o <- sample(r$n); e3 <- efficiency(r, ex$markers[o, ]); expect_equal(e3, e)
  e4 <- efficiency(r, ex$markers, chunk_size = 7); expect_equal(e4, e)
})

test_that("a constant marker is NA, a user genetic variance is used and recorded, and the share follows", {
  Z <- ex$markers; Z[, 5] <- 2
  e <- efficiency(r, Z); expect_true(is.na(e$inverse_variance[5])); expect_true(is.na(e$se_rotated[5]))
  e9 <- efficiency(r, ex$markers, sigma2_g = 9); expect_equal(attr(e9, "sigma2_g"), 9); expect_equal(attr(e9, "sigma2_g_source"), "supplied")
  expect_equal(attr(e, "sigma2_g"), r$fit$one$vg)
  s <- estimation_share(r, 9); expect_equal(s$share, asv(ex$Omega) / (9 * asv(ex$relationship) + asv(ex$Omega)))
  expect_equal(attr(e9, "share")$share, s$share)
  expect_output(print(summary(e)), "relative to the rotation"); expect_output(print(s), "Estimation share")
  expect_error(efficiency(r, ex$markers[-1, ]), "lack")
  Zna <- ex$markers; Zna[2, 3] <- NA; expect_error(efficiency(r, Zna), "NA")
})
