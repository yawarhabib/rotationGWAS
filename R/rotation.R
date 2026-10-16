# The rotation object: checked inputs, the rotation matrix, the rotated pieces, the model without a marker.

#' Rotate the adjusted means and the second-stage model
#'
#' Checks and aligns the second-stage inputs by accession name, builds the
#' rotation matrix \eqn{W = \Omega^{-1/2}} (the symmetric root), rotates the
#' means, the intercept, the covariates and the relationship matrix with it,
#' and fits the rotated model without a marker twice: with the residual
#' variance estimated, and with it fixed at one. After the rotation the
#' estimation errors of the means are independent with unit variance, so the
#' second stage has the form of an ordinary mixed-model GWAS (Piepho et al.
#' 2012). Every other function of the package takes the object returned here.
#'
#' @param means numeric vector of adjusted accession means (BLUEs, accessions
#'   fixed in the first stage), with the accession names as names.
#' @param Omega the variance-covariance matrix of the adjusted means from the
#'   first stage, an n x n matrix with the accession names as dimnames (the
#'   matrix itself, not its inverse and not a correlation matrix). A diagonal
#'   matrix, built from standard errors alone, is accepted: the rotation then
#'   reduces to inverse-variance weighting, the object says so, and the
#'   functions that need the full matrix (`efficiency()`,
#'   `estimation_share()`, `specification_test()`, `independent_effect()`)
#'   refuse it.
#' @param relationship the relationship matrix among accessions, n x n with
#'   dimnames, on its original scale; it may be singular (for example a centred
#'   VanRaden matrix).
#' @param covariates optional numeric matrix or data frame of covariates with
#'   the accession names as row names, for example principal components. The
#'   intercept is added by the package and must not be among the covariates.
#' @param max_condition the largest condition number of `Omega` accepted
#'   before the rotation is refused.
#' @return an object of class `rotation`, a list with the accession order
#'   (`ids`), `n`, `p`, `diagonal` (whether `Omega` is diagonal), the rotation matrix `W`, the rotated pieces
#'   (`y`, `X` with the rotated intercept as its first column, `K`), the aligned
#'   unrotated inputs (`means`, `Omega`, `relationship`, `covariates`), a
#'   numeric summary of `Omega` (`omega`), and the fits of the model without a
#'   marker (`fit`: `free` with both components estimated, `one` with the
#'   residual variance fixed at one).
#' @references Piepho, H.-P., J. Möhring, T. Schulz-Streeck and J. O. Ogutu,
#'   2012 A stage-wise approach for the analysis of multi-environment trials.
#'   Biometrical Journal 54: 844-860.
#' @examples
#' ex <- rotationGWAS:::toy_trial(n = 40, n_markers = 60, seed = 2)
#' r <- rotation(ex$means, ex$Omega, ex$relationship, ex$covariates)
#' r
#' @export
rotation <- function(means, Omega, relationship, covariates = NULL, max_condition = 1e10) {
  ids <- check_means(means); n <- length(ids)
  Omega <- check_square(Omega, ids, "Omega")
  relationship <- check_square(relationship, ids, "relationship", psd = TRUE)
  Q <- check_covariates(covariates, ids)
  diagonal <- all(abs(Omega[upper.tri(Omega)]) <= 1e-12 * max(abs(diag(Omega))))
  if (diagonal) message("Omega is diagonal (standard errors only): the rotation reduces to inverse-variance weighting; ",
                        "the efficiency, the estimation share, the specification test and the independent effect need the full matrix")
  # the rotation matrix from one symmetric eigen decomposition, which also serves as the guard
  ee <- eigen(Omega, symmetric = TRUE)
  lam <- ee$values
  if (min(lam) <= 0) stop("Omega is not positive definite (smallest eigenvalue ", signif(min(lam), 3), ")")
  if (max(lam) / min(lam) > max_condition)
    stop("Omega is too ill-conditioned to rotate (condition number ", signif(max(lam) / min(lam), 3), " > ", max_condition, ")")
  W <- ee$vectors %*% (t(ee$vectors) / sqrt(lam))
  W <- (W + t(W)) / 2
  y <- as.numeric(W %*% means[ids])
  X <- W %*% cbind(intercept = rep(1, n), Q)
  colnames(X) <- c("intercept", colnames(Q))
  K <- W %*% relationship %*% W; K <- (K + t(K)) / 2
  setup <- reml_setup(y, X, K)
  fit <- list(free = reml_fit_free(setup), one = reml_fit_fixed(setup, 1))
  structure(list(ids = ids, n = n, p = ncol(X), diagonal = diagonal, W = W, y = y, X = X, K = K,
                 means = means[ids], Omega = Omega, relationship = relationship, covariates = Q,
                 omega = omega_summary(Omega), fit = fit, setup = setup),
            class = "rotation")
}

#' @export
print.rotation <- function(x, digits = 4, ...) {
  cat("Rotation of", x$n, "adjusted means;", x$p, "fixed effects (the rotated intercept",
      if (x$p > 1) paste0("and ", x$p - 1, " covariate", if (x$p > 2) "s") else "", ")",
      if (isTRUE(x$diagonal)) "; Omega is diagonal, so this is inverse-variance weighting" else "", "\n")
  o <- x$omega
  cat(sprintf("Omega: diagonal from %s to %s (ratio %s); mean |correlation| %s, largest %s; condition number %s\n",
              signif(o[["diag_min"]], 3), signif(o[["diag_max"]], 3), signif(o[["diag_ratio"]], 3),
              signif(o[["mean_abs_cor"]], 3), signif(o[["max_abs_cor"]], 3), signif(o[["condition"]], 3)))
  f <- x$fit
  cat(sprintf("Model without a marker: genetic variance %s, residual variance %s (target one)%s\n",
              signif(f$free$vg, digits), signif(f$free$ve, digits),
              if (f$free$boundary != "none") paste0(" [boundary: ", f$free$boundary, "]") else ""))
  cat(sprintf("With the residual variance fixed at one: genetic variance %s%s\n",
              signif(f$one$vg, digits), if (f$one$boundary != "none") " [boundary]" else ""))
  invisible(x)
}

# ---------------------------------------------------------------------------
# input checks, all by name
# ---------------------------------------------------------------------------

check_means <- function(means) {
  if (!is.numeric(means) || is.matrix(means)) stop("means must be a numeric vector with accession names")
  ids <- names(means)
  if (is.null(ids) || any(ids == "")) stop("means must be named with the accession identifiers")
  if (anyDuplicated(ids)) stop("duplicated accession names in means: ", paste(unique(ids[duplicated(ids)]), collapse = ", "))
  if (anyNA(means)) stop("means contain NA for: ", paste(ids[is.na(means)], collapse = ", "))
  if (length(ids) < 3) stop("at least three accessions are needed")
  ids
}

check_square <- function(M, ids, what, psd = FALSE) {
  if (!is.matrix(M)) M <- as.matrix(M)
  if (nrow(M) != ncol(M)) stop(what, " must be a square matrix")
  rn <- rownames(M); cn <- colnames(M)
  if (is.null(rn) || is.null(cn)) stop(what, " must have the accession names as row and column names")
  if (!identical(rn, cn)) stop(what, " must have the same names on rows and columns")
  missing <- setdiff(ids, rn)
  if (length(missing)) stop(what, " lacks ", length(missing), " accession(s) of means: ", paste(utils::head(missing, 5), collapse = ", "), if (length(missing) > 5) ", ..." else "")
  if (length(rn) > length(ids)) message(what, ": ", length(rn) - length(ids), " accession(s) not in means dropped")
  M <- M[ids, ids, drop = FALSE]
  if (anyNA(M)) stop(what, " contains NA")
  asym <- max(abs(M - t(M))) / max(abs(M), 1e-300)
  if (asym > 1e-6) stop(what, " is not symmetric (relative asymmetry ", signif(asym, 3), ")")
  M <- (M + t(M)) / 2
  if (psd) {
    ev <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
    if (min(ev) < -1e-8 * max(abs(ev))) stop(what, " is not positive semi-definite (smallest eigenvalue ", signif(min(ev), 3), ")")
  }
  M
}

check_covariates <- function(covariates, ids) {
  n <- length(ids)
  if (is.null(covariates)) return(matrix(numeric(0), n, 0, dimnames = list(ids, NULL)))
  Q <- covariates
  if (is.data.frame(Q)) { rn <- rownames(Q); Q <- as.matrix(Q); rownames(Q) <- rn }
  if (!is.numeric(Q)) stop("covariates must be numeric")
  if (is.null(rownames(Q))) stop("covariates must have the accession names as row names")
  missing <- setdiff(ids, rownames(Q))
  if (length(missing)) stop("covariates lack ", length(missing), " accession(s) of means: ", paste(utils::head(missing, 5), collapse = ", "))
  Q <- Q[ids, , drop = FALSE]
  if (anyNA(Q)) stop("covariates contain NA")
  if (is.null(colnames(Q))) colnames(Q) <- paste0("covariate", seq_len(ncol(Q)))
  const <- apply(Q, 2, function(v) max(v) - min(v) == 0)
  if (any(const)) stop("covariate(s) ", paste(colnames(Q)[const], collapse = ", "), " are constant; the intercept is added by the package and must not be supplied")
  if (qr(cbind(1, Q))$rank < ncol(Q) + 1) stop("the covariates are collinear with each other or with the intercept")
  Q
}

# ---------------------------------------------------------------------------
# summaries of Omega
# ---------------------------------------------------------------------------

#' Average semivariance of a square matrix (Piepho 2019)
#'
#' Half the variance of a difference between two accessions, averaged over all
#' pairs; equal to the trace of the centred matrix divided by n - 1.
#' @param M a square numeric matrix.
#' @return a number.
#' @references Piepho, H.-P., 2019 A coefficient of determination (R2) for
#'   generalized linear mixed models. Biometrical Journal 61: 860-872.
#' @export
asv <- function(M) {
  M <- as.matrix(M); n <- nrow(M)
  if (n != ncol(M)) stop("M must be square")
  (sum(diag(M)) - sum(M) / n) / (n - 1)
}

omega_summary <- function(Omega) {
  d <- diag(Omega); R <- stats::cov2cor(Omega); off <- R[upper.tri(R)]
  ev <- eigen(Omega, symmetric = TRUE, only.values = TRUE)$values
  evr <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  c(diag_min = min(d), diag_max = max(d), diag_ratio = max(d) / min(d), diag_mean = mean(d), asv = asv(Omega),
    mean_abs_cor = mean(abs(off)), max_abs_cor = max(abs(off)), share_above_0.05 = mean(abs(off) > 0.05),
    condition = max(ev) / min(ev), condition_cor = max(evr) / min(evr))
}

# ---------------------------------------------------------------------------
# a synthetic example for the help pages and the tests (nothing from any trial)
# ---------------------------------------------------------------------------

#' A small synthetic second-stage dataset
#'
#' Simulates accessions with markers in two subpopulations, a VanRaden
#' relationship matrix, two principal components, and adjusted means whose
#' estimation errors have an unequal, correlated variance-covariance matrix
#' of the kind an unbalanced trial produces. For examples and tests only.
#' @param n accessions; `n_markers` markers; `seed` a seed.
#' @param vg,h2 genetic variance and the share of the genetic variance in the
#'   variance of a mean, which sets the scale of `Omega`.
#' @return list: `means`, `Omega`, `relationship`, `covariates`, `markers`
#'   (dosages 0/1/2, accessions in rows), `truth` (the causal marker and its
#'   effect).
#' @noRd
toy_trial <- function(n = 60, n_markers = 300, seed = 1, vg = 4, h2 = 0.8) {
  set.seed(seed)
  ids <- sprintf("acc%03d", seq_len(n))
  pop <- rep(c(0, 1), length.out = n)
  p1 <- stats::runif(n_markers, 0.1, 0.9); p2 <- pmin(pmax(p1 + stats::rnorm(n_markers, 0, 0.15), 0.05), 0.95)
  G <- t(sapply(seq_len(n), function(i) stats::rbinom(n_markers, 2, if (pop[i] == 0) p1 else p2)))
  dimnames(G) <- list(ids, sprintf("m%04d", seq_len(n_markers)))
  C <- scale(G, center = TRUE, scale = FALSE)
  A <- tcrossprod(C) / ncol(C); A <- A / mean(diag(A)); dimnames(A) <- list(ids, ids)
  pcs <- stats::prcomp(C)$x[, 1:2]; dimnames(pcs) <- list(ids, c("PC1", "PC2"))
  # an unbalanced-trial-like Omega: unequal precision and a shared component within groups of accessions
  grp <- rep(seq_len(ceiling(n / 6)), each = 6, length.out = n)
  Omega <- diag(stats::runif(n, 0.5, 2)) + 0.3 * outer(grp, grp, "==")
  Omega <- Omega * (vg * (1 - h2) / h2) / mean(diag(Omega)); dimnames(Omega) <- list(ids, ids)
  causal <- 17; beta <- 1.5
  u <- as.numeric(t(chol(A + diag(1e-8, n))) %*% stats::rnorm(n)) * sqrt(vg)
  e <- as.numeric(t(chol(Omega)) %*% stats::rnorm(n))
  means <- 20 + 0.5 * pcs[, 1] + beta * G[, causal] + u + e; names(means) <- ids
  list(means = means, Omega = Omega, relationship = A, covariates = pcs, markers = G,
       truth = list(causal = colnames(G)[causal], effect = beta, vg = vg))
}
