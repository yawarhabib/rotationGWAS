# The REML engine for two kernels: the additive kernel and a second one (WW', the rotated image of an effect
# independent across accessions). The error-contrast form again, but without a joint diagonalisation: with Qc the
# orthonormal complement of X, M1 = Qc' K1 Qc, M2 = Qc' K2 Qc and eta = Qc' y, minus twice the residual
# log-likelihood is, up to a constant,
#   dev(vg, vi, ve) = log|vg M1 + vi M2 + ve I| + eta' (vg M1 + vi M2 + ve I)^-1 eta,
# evaluated by a Cholesky factorisation of an (n - p) x (n - p) matrix. The full-form deviance of the paper's script
# 04 is this plus log|X'X|, the same constant as for one kernel, so the two engines' deviances are comparable.

#' Set up the error-contrast form of a two-kernel model
#' @return a list of class `reml_setup2`: `M1`, `M2`, `eta`, `n`, `p`, `df`, `logdet_XtX`, `Qc`, `y`, `X`, `K1`, `K2`.
#' @noRd
reml_setup2 <- function(y, X, K1, K2) {
  y <- as.numeric(y); X <- as.matrix(X); n <- length(y); p <- ncol(X)
  if (nrow(X) != n) stop("X must have one row per element of y")
  for (K in list(K1, K2)) if (!is.matrix(K) || any(dim(K) != n)) stop("the kernels must be n x n matrices")
  if (anyNA(y) || anyNA(X) || anyNA(K1) || anyNA(K2)) stop("y, X and the kernels must not contain NA")
  qx <- qr(X)
  if (qx$rank < p) stop("X is not of full column rank")
  if (p >= n) stop("there must be more observations than fixed effects")
  Q <- qr.Q(qx, complete = TRUE); Qc <- Q[, (p + 1):n, drop = FALSE]
  sym <- function(M) (M + t(M)) / 2
  M1 <- sym(crossprod(Qc, K1 %*% Qc)); M2 <- sym(crossprod(Qc, K2 %*% Qc))
  R <- qr.R(qx)
  structure(list(M1 = M1, M2 = M2, eta = as.numeric(crossprod(Qc, y)), n = n, p = p, df = n - p,
                 logdet_XtX = 2 * sum(log(abs(diag(R)))), Qc = Qc, y = y, X = X, K1 = sym(K1), K2 = sym(K2)),
            class = "reml_setup2")
}

#' The error-contrast deviance of the two-kernel model (full form adds log|X'X|)
#' @noRd
reml_deviance2 <- function(s, vg, vi, ve, full = TRUE) {
  if (vg < 0 || vi < 0 || ve < 0) return(Inf)
  V <- vg * s$M1 + vi * s$M2 + diag(ve, s$df)
  R <- tryCatch(chol(V), error = function(e) NULL)
  if (is.null(R)) return(Inf)
  z <- forwardsolve(t(R), s$eta)
  dev <- 2 * sum(log(diag(R))) + sum(z^2)
  if (full) dev + s$logdet_XtX else dev
}

#' REML fit of the two-kernel model, free or with the residual fixed
#'
#' Nelder-Mead on the logarithms of the components, restarted from its own
#' solution until the deviance no longer decreases (script 04's rule), from
#' several starts: the one-kernel fit with a small second component, and the
#' second component carrying most of the variance. A component that ends
#' below `tol_boundary` times the total is reported on the boundary and set
#' to zero, the others re-optimised.
#'
#' @param s a `reml_setup2`.
#' @param ve `NULL` to estimate the residual variance, or its fixed value.
#' @param start optional starting values c(vg, vi, ve).
#' @return list: vg, vi, ve, deviance (full form), deviance_ec, boundary (a
#'   character vector, empty when interior: among "vg_zero", "vi_zero",
#'   "ve_zero"), restarts, df, residual ("estimated" or "fixed").
#' @noRd
reml_fit2 <- function(s, ve = NULL, start = NULL, tol_boundary = 1e-6) {
  fixed_ve <- !is.null(ve)
  tot <- sum(s$eta^2) / s$df                                    # the scale of the variances
  dev_of <- function(par, free) {                               # par: logs of the free components, in the order vg, vi, (ve)
    v <- exp(par)
    if (fixed_ve) reml_deviance2(s, v[1], v[2], ve, full = FALSE) else reml_deviance2(s, v[1], v[2], v[3], full = FALSE)
  }
  starts <- list()
  if (!is.null(start)) starts[[1]] <- log(pmax(start, 1e-8 * tot))
  m1 <- mean(diag(s$M1)); m2 <- mean(diag(s$M2))
  base <- if (fixed_ve) c(tot / m1 / 2, tot / m2 / 2) else c(tot / m1 / 3, tot / m2 / 3, tot / 3)
  starts[[length(starts) + 1]] <- log(base)
  starts[[length(starts) + 1]] <- log(base * c(1, 0.01, if (fixed_ve) NULL else 1))
  starts[[length(starts) + 1]] <- log(base * c(0.01, 1, if (fixed_ve) NULL else 1))
  best <- NULL; restarts_total <- 0
  for (st in starts) {
    cur <- st; val <- dev_of(cur); k <- 0
    repeat {                                                   # restart from the solution until no decrease
      o <- stats::optim(cur, dev_of, method = "Nelder-Mead", control = list(maxit = 5000, reltol = 1e-14))
      k <- k + 1
      if (o$value < val - 1e-9) { cur <- o$par; val <- o$value } else { cur <- o$par; val <- min(val, o$value); break }
      if (k > 50) break
    }
    restarts_total <- restarts_total + k
    if (is.null(best) || val < best$value) best <- list(par = cur, value = val)
  }
  v <- exp(best$par)
  comp <- if (fixed_ve) c(vg = v[1], vi = v[2], ve = ve) else c(vg = v[1], vi = v[2], ve = v[3])
  scale <- sum(comp[c("vg", "vi")] * c(m1, m2)) + comp[["ve"]]
  boundary <- names(comp)[comp < tol_boundary * scale & !(names(comp) == "ve" & fixed_ve)]
  if (length(boundary)) {                                       # set the boundary components to zero and re-optimise the rest
    comp[boundary] <- 0
    free_names <- setdiff(names(comp), c(boundary, if (fixed_ve) "ve"))
    if (length(free_names)) {
      f <- function(par) { cc <- comp; cc[free_names] <- exp(par); reml_deviance2(s, cc[["vg"]], cc[["vi"]], cc[["ve"]], full = FALSE) }
      o <- stats::optim(log(pmax(comp[free_names], 1e-8 * tot)), f, method = if (length(free_names) == 1) "Brent" else "Nelder-Mead",
                        lower = if (length(free_names) == 1) log(1e-12 * tot) else -Inf, upper = if (length(free_names) == 1) log(1e12 * tot) else Inf,
                        control = list(maxit = 5000, reltol = 1e-14))
      comp[free_names] <- exp(o$par)
    }
  }
  dev_ec <- reml_deviance2(s, comp[["vg"]], comp[["vi"]], comp[["ve"]], full = FALSE)
  list(vg = comp[["vg"]], vi = comp[["vi"]], ve = comp[["ve"]], deviance = dev_ec + s$logdet_XtX, deviance_ec = dev_ec,
       boundary = paste0(sub("^v", "v", boundary), if (length(boundary)) "_zero" else character(0)), restarts = restarts_total, df = s$df,
       residual = if (fixed_ve) "fixed" else "estimated")
}

#' Profile-likelihood interval for the residual variance of the two-kernel model
#' @noRd
reml_profile_ve2 <- function(s, free, level = 0.95, crit = stats::qchisq(level, 1)) {
  dmin <- free$deviance_ec
  prof <- function(ve) reml_fit2(s, ve = ve, start = c(max(free$vg, 1e-8), max(free$vi, 1e-8), ve))$deviance_ec - dmin - crit
  tot <- sum(s$eta^2) / s$df
  ve_hat <- if (free$ve > 0) free$ve else 1e-8 * tot
  lo <- 1e-8 * tot
  lower <- if (prof(lo) <= 0) 0 else stats::uniroot(prof, c(lo, ve_hat), tol = 1e-8 * tot)$root
  hi <- max(ve_hat, lo) * 2; bounded <- TRUE
  while (prof(hi) <= 0) { hi <- hi * 2; if (hi > 1e8 * tot) { bounded <- FALSE; break } }
  upper <- if (bounded) stats::uniroot(prof, c(max(ve_hat, lo), hi), tol = 1e-8 * tot)$root else NA_real_
  list(lower = lower, upper = upper, level = level, crit = crit, bounded = bounded,
       inside_one = (lower <= 1) && (is.na(upper) || upper >= 1))
}

#' The P3D projection of a general covariance V = vg K1 + vi K2 + ve I in the full form
#'
#' P = V^-1 - V^-1 X (X' V^-1 X)^-1 X' V^-1, returned as the factor H with
#' P = H' H (rows: the n - p non-null directions), so that the scan's block
#' products are the same as with one kernel.
#' @noRd
scan_projection2 <- function(s, vg, vi, ve) {
  V <- vg * s$K1 + vi * s$K2 + diag(ve, s$n)
  R <- chol((V + t(V)) / 2); Vi <- chol2inv(R)
  ViX <- Vi %*% s$X
  P <- Vi - ViX %*% solve(crossprod(s$X, ViX), t(ViX)); P <- (P + t(P)) / 2
  ee <- eigen(P, symmetric = TRUE)
  keep <- seq_len(s$df)
  lam <- pmax(ee$values[keep], 0)
  H <- t(ee$vectors[, keep, drop = FALSE]) * sqrt(lam)               # rows scaled: H = diag(sqrt(lam)) U'
  Hy <- as.numeric(H %*% s$y)
  list(H = H, Hy = Hy, yPy = sum(Hy^2))
}
