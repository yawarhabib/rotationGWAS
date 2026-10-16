# The answer before the scan: the marker-specific efficiency and the estimation share.

#' Marker-specific efficiency of approximate weighting
#'
#' For every marker, the standard error of the marker effect under
#' inverse-variance weighting and under the unweighted analysis relative to
#' its standard error under the rotation, predicted from the design without a
#' scan. The variance of the marker effect under a working covariance
#' \eqn{V_w} is evaluated under the true covariance
#' \eqn{\Sigma = \sigma_g^2 A + \Omega} by the sandwich expression, and the
#' rotation, whose working covariance is \eqn{\Sigma} itself, gives the
#' generalised least squares minimum. The working covariances are
#' \eqn{\sigma_g^2 A + \mathrm{diag}(\Omega)} for inverse-variance weighting
#' and \eqn{\sigma_g^2 A + cI} for the unweighted analysis, with \eqn{c} the
#' mean diagonal element of \eqn{\Omega}. Nothing is rotated or scanned here.
#'
#' @param x a `rotation` object.
#' @param markers numeric matrix of marker dosages (0, 1, 2), accessions in
#'   rows with the accession names as row names, markers in columns; complete
#'   (no NA).
#' @param sigma2_g the genetic variance on the scale of the relationship
#'   matrix. By default the estimate from the rotated model without a marker
#'   with the residual variance fixed at one (`x$fit$one$vg`).
#' @param chunk_size markers per block of matrix products.
#' @return a data frame of class `efficiency`, one row per marker:
#'   `marker`, `maf` (frequency of the less frequent allele), `se_rotated`
#'   (the predicted standard error under the rotation at the given variances),
#'   `inverse_variance` and `unweighted` (standard-deviation factors relative
#'   to the rotation, at or above one; their squares are the ratios on the
#'   variance scale). NA for a marker with no variance. Attributes
#'   `sigma2_g`, `sigma2_g_source`, `share` (the [estimation_share()]) and
#'   `n`, `p`.
#' @examples
#' ex <- rotationGWAS:::toy_trial(n = 40, n_markers = 60, seed = 2)
#' r <- rotation(ex$means, ex$Omega, ex$relationship, ex$covariates)
#' e <- efficiency(r, ex$markers)
#' summary(e)
#' @export
efficiency <- function(x, markers, sigma2_g = NULL, chunk_size = 5000) {
  stopifnot(inherits(x, "rotation")); needs_full_omega(x, "the efficiency of approximate weighting")
  Z <- check_markers(markers, x$ids)
  sg <- genetic_variance(x, sigma2_g)
  n <- x$n; X0 <- cbind(1, x$covariates)
  Sigma <- sg$value * x$relationship + x$Omega
  works <- list(rotated = Sigma, inverse_variance = sg$value * x$relationship + diag(diag(x$Omega)),
                unweighted = sg$value * x$relationship + diag(mean(diag(x$Omega)), n))
  proj <- lapply(works, function(V) {                  # P_M = M - M X0 (X0' M X0)^-1 X0' M with M = V^-1
    M <- chol2inv(chol(V)); MX <- M %*% X0
    P <- M - MX %*% solve(crossprod(X0, MX), t(MX)); P <- (P + t(P)) / 2
    list(P = P, B = P %*% Sigma %*% P)
  })
  m <- ncol(Z); v <- matrix(NA_real_, m, 3, dimnames = list(colnames(Z), names(works)))
  for (from in seq(1, m, by = chunk_size)) {
    j <- from:min(from + chunk_size - 1, m); Zj <- Z[, j, drop = FALSE]
    for (k in names(works)) {
      PZ <- proj[[k]]$P %*% Zj; den <- colSums(Zj * PZ)
      num <- if (k == "rotated") den else colSums(Zj * (proj[[k]]$B %*% Zj))
      v[j, k] <- ifelse(den > 1e-10 * max(den, 1e-300), num / den^2, NA_real_)
    }
  }
  tested <- is.finite(v[, "rotated"])
  out <- data.frame(marker = colnames(Z), maf = pmin(colMeans(Z) / 2, 1 - colMeans(Z) / 2),
                    se_rotated = sqrt(v[, "rotated"]),
                    inverse_variance = sqrt(v[, "inverse_variance"] / v[, "rotated"]),
                    unweighted = sqrt(v[, "unweighted"] / v[, "rotated"]), row.names = NULL)
  out[!tested, c("se_rotated", "inverse_variance", "unweighted")] <- NA_real_
  structure(out, class = c("efficiency", "data.frame"), sigma2_g = sg$value, sigma2_g_source = sg$source,
            share = estimation_share(x, sg$value), n = n, p = ncol(X0))
}

#' @export
summary.efficiency <- function(object, probs = c(0.5, 0.95, 1), ...) {
  q <- function(v) stats::quantile(v, probs, na.rm = TRUE)
  tab <- rbind(`inverse-variance` = q(object$inverse_variance), unweighted = q(object$unweighted))
  colnames(tab) <- paste0(c("median", "95%", "max")[seq_along(probs)], if (length(probs) > 3) "" else "")
  structure(list(table = tab, n_markers = sum(!is.na(object$inverse_variance)), n_untested = sum(is.na(object$inverse_variance)),
                 sigma2_g = attr(object, "sigma2_g"), sigma2_g_source = attr(object, "sigma2_g_source"), share = attr(object, "share")),
            class = "summary.efficiency")
}

#' @export
print.summary.efficiency <- function(x, digits = 4, ...) {
  cat("Standard error of the marker effect relative to the rotation (factor at or above one), over", x$n_markers, "markers",
      if (x$n_untested) paste0(" (", x$n_untested, " without variance, NA)") else "", "\n")
  print(round(x$table, digits))
  cat(sprintf("Genetic variance %s (%s); estimation share %s, on the average semivariance\n",
              signif(x$sigma2_g, digits), x$sigma2_g_source, signif(x$share$share, 3)))
  invisible(x)
}

#' @export
print.efficiency <- function(x, ...) { print(summary(x)); invisible(x) }

#' @export
`[.efficiency` <- function(x, ...) { out <- NextMethod(); class(out) <- "data.frame"; out }   # a subset is a plain data frame

#' The estimation share of the adjusted means
#'
#' The average semivariance (ASV) of the estimation error, \eqn{\Omega},
#' divided by the ASV of \eqn{\sigma_g^2 A + \Omega}: the part of the
#' variance of the adjusted means, on the genotype-difference basis, that
#' weighting can act on.
#'
#' @inheritParams efficiency
#' @return a list of class `estimation_share`: `asv_genetic`,
#'   `asv_estimation`, `share`, `sigma2_g`, `sigma2_g_source`.
#' @examples
#' ex <- rotationGWAS:::toy_trial(n = 40, n_markers = 60, seed = 2)
#' r <- rotation(ex$means, ex$Omega, ex$relationship, ex$covariates)
#' estimation_share(r)
#' @export
estimation_share <- function(x, sigma2_g = NULL) {
  stopifnot(inherits(x, "rotation")); needs_full_omega(x, "the estimation share")
  sg <- genetic_variance(x, sigma2_g)
  g <- sg$value * asv(x$relationship); e <- asv(x$Omega)
  structure(list(asv_genetic = g, asv_estimation = e, share = e / (g + e), sigma2_g = sg$value, sigma2_g_source = sg$source),
            class = "estimation_share")
}

#' @export
print.estimation_share <- function(x, digits = 4, ...) {
  cat(sprintf("Estimation share %s: ASV of the estimation error %s against %s for the genetic part (genetic variance %s, %s)\n",
              signif(x$share, 3), signif(x$asv_estimation, digits), signif(x$asv_genetic, digits), signif(x$sigma2_g, digits), x$sigma2_g_source))
  invisible(x)
}

# ---------------------------------------------------------------------------
# helpers shared with the scan
# ---------------------------------------------------------------------------

needs_full_omega <- function(x, what) {
  if (isTRUE(x$diagonal))
    stop("Omega is diagonal (standard errors only): ", what, " needs the full variance-covariance matrix of the adjusted means")
  invisible(TRUE)
}

genetic_variance <- function(x, sigma2_g) {
  if (is.null(sigma2_g)) {
    v <- x$fit$one$vg
    if (x$fit$one$boundary != "none") warning("the genetic variance of the model without a marker is on the boundary (zero)")
    list(value = v, source = "rotated model without a marker, residual variance fixed at one")
  } else {
    if (!is.numeric(sigma2_g) || length(sigma2_g) != 1 || !is.finite(sigma2_g) || sigma2_g < 0) stop("sigma2_g must be one non-negative number")
    list(value = sigma2_g, source = "supplied")
  }
}

check_markers <- function(markers, ids) {
  Z <- markers
  if (is.data.frame(Z)) { rn <- rownames(Z); Z <- as.matrix(Z); rownames(Z) <- rn }
  if (!is.matrix(Z) || !is.numeric(Z)) stop("markers must be a numeric matrix with accessions in rows")
  if (is.null(rownames(Z))) stop("markers must have the accession names as row names")
  missing <- setdiff(ids, rownames(Z))
  if (length(missing)) stop("markers lack ", length(missing), " accession(s) of means: ", paste(utils::head(missing, 5), collapse = ", "))
  Z <- Z[ids, , drop = FALSE]
  if (anyNA(Z)) stop("markers contain NA: impute missing calls before the analysis (a missing call cannot be handled after the rotation)")
  if (is.null(colnames(Z))) colnames(Z) <- paste0("marker", seq_len(ncol(Z)))
  storage.mode(Z) <- "double"
  Z
}
