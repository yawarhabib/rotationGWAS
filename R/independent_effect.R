# The independent effect: the rotated model extended by a genetic effect that is independent across accessions.

#' Fit an effect independent across accessions alongside the relationship matrix
#'
#' When the specification test finds the residual variance of the rotated
#' model above one, genetic covariance is missing from the model, if the
#' first stage is correctly specified. The
#' extension the paper examines adds an effect independent across accessions
#' (Oakey et al. 2007; Endelman 2023; Holland and Piepho 2024), whose kernel
#' is the identity before the rotation and \eqn{WW'} after it, the inverse of
#' the first-stage variance-covariance matrix. Without the rotation this
#' effect and the residual error have the same covariance matrix and only the
#' sum of their variances can be estimated; after it the two are separated.
#' The model without a marker is fitted with three variance components, by
#' REML with the residual variance estimated and with it fixed at one, and
#' compared with the additive model by the reduction in deviance. The scan can
#' then be run under this model with [genome_scan()] and its `model` argument.
#'
#' A residual variance at the boundary of zero means that all the estimation
#' error has been attributed to the independent effect; the model is then the
#' unweighted analysis of the adjusted means, and the scan under it is the
#' unweighted scan. A component is reported on the boundary when the fit puts
#' it below 1e-6 of the total variance; the paper's simulation summaries count
#' a residual variance below 0.001 as a boundary solution. The separation of
#' the independent effect from the residual is weak when the variance-
#' covariance matrix of the means is close to proportional to the identity.
#'
#' @param x a `rotation` object.
#' @return an object of class `independent_effect`: `fit` (the free fit:
#'   `vg`, `vi`, `ve`, `deviance`, `boundary`), `fit_one` (with the residual
#'   fixed at one), `additive` (the free fit of the additive model), `lr` (the
#'   reduction in deviance against the additive model), `p_mixture` (its
#'   p-value from the equal mixture of a point mass at zero and chi-squared on
#'   one degree of freedom, the variance being on the boundary under the null),
#'   `delta_aic` (AIC of the extended model minus the additive one), the
#'   profile-likelihood intervals for the residual variance under the extended
#'   model (`interval_90`, `interval_95`), `kernel` (\eqn{WW'}), and the setup
#'   the scan uses.
#' @references Oakey, H., A. P. Verbyla, B. R. Cullis, X. Wei and W. S.
#'   Pitchford, 2007 Joint modeling of additive and non-additive genetic line
#'   effects in single field trials. Theoretical and Applied Genetics 113:
#'   809-819. Endelman, J. B., 2023 Fully efficient, two-stage analysis of
#'   multi-environment trials with directional dominance and multi-trait
#'   genomic selection. Theoretical and Applied Genetics 136: 65. Holland,
#'   J. B., and H.-P. Piepho, 2024 Don't BLUP twice. G3 14: jkae250.
#' @examples
#' ex <- rotationGWAS:::toy_trial(n = 60, n_markers = 40, seed = 2)
#' r <- rotation(ex$means, ex$Omega, ex$relationship, ex$covariates)
#' ie <- independent_effect(r)
#' ie
#' @export
independent_effect <- function(x) {
  stopifnot(inherits(x, "rotation")); needs_full_omega(x, "the effect independent across accessions")
  Kind <- x$W %*% x$W; Kind <- (Kind + t(Kind)) / 2                     # WW' = Omega^-1
  s2 <- reml_setup2(x$y, x$X, x$K, Kind)
  add <- x$fit$free
  free <- reml_fit2(s2, start = c(add$vg, 0.1 * max(add$vg, 1e-8), max(add$ve, 1e-8)))
  one <- reml_fit2(s2, ve = 1, start = c(x$fit$one$vg, 0.1 * max(x$fit$one$vg, 1e-8), 1))
  lr <- max(add$deviance - free$deviance, 0)
  p_mix <- if (lr > 0) 0.5 * stats::pchisq(lr, 1, lower.tail = FALSE) else 1
  ci90 <- reml_profile_ve2(s2, free, 0.90, crit = 2.71); ci95 <- reml_profile_ve2(s2, free, 0.95, crit = 3.84)
  keep <- function(ci) list(lower = ci$lower, upper = ci$upper, bounded = ci$bounded, inside_one = ci$inside_one)
  structure(list(fit = free, fit_one = one, additive = add, lr = lr, p_mixture = p_mix, delta_aic = 2 - lr,
                 interval_90 = keep(ci90), interval_95 = keep(ci95), kernel = Kind, setup = s2, n = x$n, p = x$p),
            class = "independent_effect")
}

#' @export
print.independent_effect <- function(x, digits = 4, ...) {
  f <- x$fit; b <- f$boundary
  fmt <- function(ci) sprintf("[%s, %s]", signif(ci$lower, digits), if (ci$bounded) signif(ci$upper, digits) else "unbounded")
  cat(sprintf("Rotated model with the relationship matrix and an effect independent across accessions (kernel WW')\n"))
  cat(sprintf("Genetic variance %s, independent effect %s, residual variance %s (target one)%s\n",
              signif(f$vg, digits), signif(f$vi, digits), signif(f$ve, digits),
              if (length(b)) paste0(" [boundary: ", paste(b, collapse = ", "), "]") else ""))
  if ("ve_zero" %in% b) cat("  The residual variance is at the boundary of zero: with it, this model is the unweighted analysis of the adjusted means.\n")
  if ("vi_zero" %in% b) cat("  The independent effect is zero: the additive model stands.\n")
  cat(sprintf("Against the additive model (genetic variance %s, residual %s): deviance lower by %s (p = %s on the boundary mixture), AIC %s by %s\n",
              signif(x$additive$vg, digits), signif(x$additive$ve, digits), signif(x$lr, digits), signif(x$p_mixture, 3),
              if (x$delta_aic <= 0) "lower" else "higher", signif(abs(x$delta_aic), digits)))
  cat(sprintf("With the residual variance fixed at one: genetic variance %s, independent effect %s\n", signif(x$fit_one$vg, digits), signif(x$fit_one$vi, digits)))
  cat(sprintf("Profile-likelihood intervals for the residual variance under this model: 90%% %s, 95%% %s\n", fmt(x$interval_90), fmt(x$interval_95)))
  invisible(x)
}
