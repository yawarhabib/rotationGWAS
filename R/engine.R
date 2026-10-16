# The REML engine of rotationGWAS.
#
# One kernel, in the error-contrast spectral form. For the model
#   y = X b + u + e,   u ~ N(0, vg K),   e ~ N(0, ve I),
# let Qc be an orthonormal basis of the orthogonal complement of the column space
# of X (n x (n - p)), and Qc' K Qc = U diag(d) U'. With eta = U' Qc' y, minus
# twice the residual log-likelihood is, up to a constant,
#   dev(vg, ve) = sum(log(vg d + ve)) + sum(eta^2 / (vg d + ve)).
# This form needs no inversion of V, takes a kernel of any rank (a zero
# eigenvalue only contributes log(ve) + eta^2 / ve), and makes every evaluation
# O(n - p) after one decomposition. The full-form deviance of the paper's
# scripts, log|V| + log|X' V^-1 X| + y' P y, equals dev + log|X' X|.
#
# Nothing in this file carries a constant on the scale of the trait: brackets
# and starting values come from the eigenvalues d and the contrasts eta.

# ---------------------------------------------------------------------------
# Setup: checks, projection, decomposition
# ---------------------------------------------------------------------------

#' Set up the error-contrast spectral form of one model
#'
#' @param y numeric response of length n (the rotated means).
#' @param X numeric n x p matrix of fixed effects, full column rank (the rotated
#'   intercept, rotated covariates, and a rotated marker when one is tested).
#' @param K numeric n x n symmetric positive semi-definite kernel (the rotated
#'   relationship matrix); any rank.
#' @return a list of class `reml_setup` with the eigenvalues `d` (length
#'   n - p, non-negative), the contrasts `eta`, `n`, `p`, `logdet_XtX`, and the
#'   matrix `B` = U' Qc' (n - p by n) used by the genome scan.
#' @noRd
reml_setup <- function(y, X, K) {
  y <- as.numeric(y); X <- as.matrix(X); n <- length(y); p <- ncol(X)
  if (nrow(X) != n) stop("X must have one row per element of y")
  if (!is.matrix(K) || any(dim(K) != n)) stop("K must be an n x n matrix")
  if (anyNA(y) || anyNA(X) || anyNA(K)) stop("y, X and K must not contain NA")
  qx <- qr(X)
  if (qx$rank < p) stop("X is not of full column rank (rank ", qx$rank, " for ", p, " columns)")
  if (p >= n) stop("there must be more observations than fixed effects")
  Q  <- qr.Q(qx, complete = TRUE)
  Qc <- Q[, (p + 1):n, drop = FALSE]
  K  <- (K + t(K)) / 2
  M  <- crossprod(Qc, K %*% Qc)
  ee <- eigen((M + t(M)) / 2, symmetric = TRUE)
  d  <- ee$values
  neg <- -min(d) / max(abs(d), 1)
  if (neg > 1e-6) stop("K is not positive semi-definite (relative negative eigenvalue ", signif(neg, 3), ")")
  d  <- pmax(d, 0)
  B  <- crossprod(ee$vectors, t(Qc))          # (n - p) x n, rows orthonormal in R^n
  eta <- as.numeric(B %*% y)
  R <- qr.R(qx)
  structure(list(d = d, eta = eta, n = n, p = p, df = n - p,
                 logdet_XtX = 2 * sum(log(abs(diag(R)))), B = B, y = y, X = X, K = K),
            class = "reml_setup")
}

# ---------------------------------------------------------------------------
# Deviance
# ---------------------------------------------------------------------------

#' REML deviance of the error-contrast form
#'
#' @param s a `reml_setup`.
#' @param vg,ve variance components (vg >= 0, ve > 0 unless vg > 0).
#' @param full if TRUE, the full-form deviance log|V| + log|X'V^-1 X| + y'Py of
#'   the paper's scripts (error-contrast deviance plus log|X'X|).
#' @noRd
reml_deviance <- function(s, vg, ve, full = TRUE) {
  v <- vg * s$d + ve
  if (any(v <= 0)) return(Inf)
  dev <- sum(log(v)) + sum(s$eta^2 / v)
  if (full) dev + s$logdet_XtX else dev
}

# ---------------------------------------------------------------------------
# Free fit: both components estimated, profiled over the ratio
# ---------------------------------------------------------------------------

#' REML fit with both variance components estimated
#'
#' The deviance is profiled over the ratio lambda = ve / vg: for a given
#' lambda the REML estimate of vg is sum(eta^2 / (d + lambda)) / (n - p), so
#' the fit is a one-dimensional search in log(lambda). The bracket is set from
#' the eigenvalues, a coarse grid finds the basin, and `optimize` refines it.
#' The two ends of the bracket are the boundaries ve -> 0 and vg -> 0; the fit
#' reports which, if either, is the optimum.
#'
#' @return list: vg, ve, lambda, deviance (full form), deviance_ec, boundary
#'   (`"none"`, `"ve_zero"` or `"vg_zero"`), the profiled deviance function
#'   value at both boundaries, and `df`.
#' @noRd
reml_fit_free <- function(s) {
  d <- s$d; eta <- s$eta; m <- s$df
  if (all(d <= 0)) stop("the projected kernel is zero: the genetic variance is not identifiable")
  dpos <- d[d > 0]
  scale <- mean(dpos)
  prof <- function(loglam) {                   # profiled deviance in log(lambda), error-contrast form
    lam <- exp(loglam); w <- d + lam
    vg <- sum(eta^2 / w) / m
    m * log(vg) + sum(log(w)) + m
  }
  lo <- log(scale) - 30; hi <- log(scale) + 30  # lambda from 1e-13 to 1e13 times the eigenvalue scale
  grid <- seq(lo, hi, length.out = 121)
  pv <- vapply(grid, prof, numeric(1))
  i <- which.min(pv)
  a <- grid[max(1, i - 1)]; b <- grid[min(length(grid), i + 1)]
  o <- stats::optimize(prof, interval = c(a, b), tol = 1e-12)
  if (prof(o$minimum) > pv[i]) o <- list(minimum = grid[i], objective = pv[i])
  loglam <- o$minimum; lam <- exp(loglam)
  vg <- sum(eta^2 / (d + lam)) / m; ve <- lam * vg
  # boundaries: the profiled deviance at the ends of the bracket
  dev_lo <- prof(lo); dev_hi <- prof(hi)
  boundary <- "none"
  if (i <= 2 || dev_lo <= o$objective + 1e-8) boundary <- "ve_zero"
  if (i >= length(grid) - 1 || dev_hi <= o$objective + 1e-8) boundary <- "vg_zero"
  if (boundary == "ve_zero") { lam <- exp(lo); vg <- sum(eta^2 / (d + lam)) / m; ve <- lam * vg }
  if (boundary == "vg_zero") { ve <- sum(eta^2) / m; vg <- 0; lam <- Inf }
  dev_ec <- reml_deviance(s, vg, ve, full = FALSE)
  list(vg = vg, ve = ve, lambda = lam, deviance = dev_ec + s$logdet_XtX, deviance_ec = dev_ec,
       boundary = boundary, df = m, residual = "estimated")
}

# ---------------------------------------------------------------------------
# Fit with the residual variance fixed
# ---------------------------------------------------------------------------

#' REML fit with the residual variance fixed (at one for the rotated model)
#'
#' Only vg is estimated. The derivative of the deviance at vg = 0 decides
#' whether the estimate is on the boundary; otherwise the bracket is grown
#' from the moment estimate sum(eta^2) / sum(d) until the derivative changes
#' sign, and `optimize` in log(vg) finds the minimum.
#'
#' @param ve the fixed residual variance (1 for the rotated model).
#' @return list: vg, ve, deviance (full form), deviance_ec, boundary
#'   (`"none"` or `"vg_zero"`), df.
#' @noRd
reml_fit_fixed <- function(s, ve = 1) {
  d <- s$d; eta <- s$eta; m <- s$df
  if (ve <= 0) stop("ve must be positive")
  dd <- function(vg) { v <- vg * d + ve; sum(d / v) - sum(eta^2 * d / v^2) }   # d dev / d vg
  if (dd(0) >= 0) {
    vg <- 0; boundary <- "vg_zero"
  } else {
    v0 <- max(sum(eta^2) / sum(d), ve / max(d))      # moment start, never zero
    hi <- v0
    for (k in 1:200) { if (dd(hi) > 0) break; hi <- hi * 4 }
    if (dd(hi) <= 0) stop("no finite maximum of the likelihood in vg")
    lo <- hi / 4; for (k in 1:200) { if (dd(lo) < 0) break; lo <- lo / 4 }
    o <- stats::optimize(function(lv) reml_deviance(s, exp(lv), ve, full = FALSE), interval = log(c(lo, hi)), tol = 1e-12)
    vg <- exp(o$minimum); boundary <- "none"
  }
  dev_ec <- reml_deviance(s, vg, ve, full = FALSE)
  list(vg = vg, ve = ve, lambda = if (vg > 0) ve / vg else Inf, deviance = dev_ec + s$logdet_XtX, deviance_ec = dev_ec,
       boundary = boundary, df = m, residual = "fixed")
}

# ---------------------------------------------------------------------------
# Profile-likelihood interval for the residual variance
# ---------------------------------------------------------------------------

#' Profile-likelihood interval for the residual variance
#'
#' At each value of ve the genetic variance is re-optimised; the interval at
#' level `level` is the set of ve whose profile deviance lies within the
#' `level` quantile of chi-squared on one degree of freedom of the minimum.
#' The lower limit is reported as 0 when the deviance is still inside at
#' 1e-8 times the total variance; the upper limit is NA, with `bounded = FALSE`,
#' when it is still inside at 1e8 times the total variance.
#'
#' @param free the free fit of the same setup (from `reml_fit_free`).
#' @param crit the cut-off on the profile deviance; by default the `level`
#'   quantile of chi-squared on one degree of freedom, but the package passes
#'   the paper's rounded values, 2.71 for the 90% and 3.84 for the 95% interval,
#'   so that the tests and the intervals agree exactly.
#' @return list: lower, upper, level, crit, bounded, and `inside_one` (whether
#'   one lies in the interval), plus the profile deviance function.
#' @noRd
reml_profile_ve <- function(s, free = reml_fit_free(s), level = 0.95, crit = stats::qchisq(level, 1)) {
  dmin <- free$deviance_ec
  prof <- function(ve) reml_fit_fixed(s, ve)$deviance_ec - dmin - crit     # <= 0 inside the interval
  tot <- free$vg * mean(s$d) + free$ve                                     # the scale of the variances
  if (free$ve > 0 && prof(free$ve) > 1e-6) stop("the free fit is not the minimum of the profile")
  ve_hat <- if (free$ve > 0) free$ve else 1e-8 * tot
  # lower limit
  lo <- 1e-8 * tot
  lower <- if (prof(lo) <= 0) 0 else stats::uniroot(prof, c(lo, ve_hat), tol = 1e-10 * tot)$root
  # upper limit
  hi <- max(ve_hat, lo) * 2; bounded <- TRUE
  while (prof(hi) <= 0) { hi <- hi * 2; if (hi > 1e8 * tot) { bounded <- FALSE; break } }
  upper <- if (bounded) stats::uniroot(prof, c(max(ve_hat, lo), hi), tol = 1e-10 * tot)$root else NA_real_
  list(lower = lower, upper = upper, level = level, crit = crit, bounded = bounded,
       inside_one = (lower <= 1) && (is.na(upper) || upper >= 1), profile = prof)
}

# ---------------------------------------------------------------------------
# Generalised least squares at given variance components
# ---------------------------------------------------------------------------

#' Fixed effects and standard errors at given variance components
#'
#' V = vg K + ve I is inverted by Cholesky; the covariance of the estimates is
#' (X' V^-1 X)^-1 with no scale factor, as in the refit with the residual
#' fixed at one.
#' @noRd
gls_fixed <- function(y, X, K, vg, ve) {
  n <- length(y)
  V <- vg * K + diag(ve, n)
  R <- chol(V)
  Vi_X <- backsolve(R, forwardsolve(t(R), X))
  Vi_y <- backsolve(R, forwardsolve(t(R), y))
  XtViX <- crossprod(X, Vi_X)
  C <- solve(XtViX)
  beta <- as.numeric(C %*% crossprod(X, Vi_y))
  list(beta = beta, se = sqrt(diag(C)), cov = C)
}

#' The projection of the P3D scan at a given variance ratio
#'
#' With V proportional to K + lambda I, P = V^-1 - V^-1 X (X'V^-1 X)^-1 X'V^-1
#' equals Qc (Qc' V Qc)^-1 Qc' = B' diag(1 / (d + lambda)) B, so that for any
#' column z, z'Pz = sum((H z)^2) with H = diag(1 / sqrt(d + lambda)) B. The scan
#' uses H once per block of markers.
#' @return the (n - p) x n matrix H and the vector H y.
#' @noRd
scan_projection <- function(s, lambda) {
  w <- 1 / sqrt(s$d + lambda)
  H <- s$B * w                                 # rows scaled
  list(H = H, Hy = as.numeric(H %*% s$y), yPy = sum((s$eta * w)^2))
}
