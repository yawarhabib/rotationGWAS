# The specification test: the residual variance of the rotated model against its target of one.

#' Test the residual variance of the rotated model against one
#'
#' After the rotation the estimation errors of the adjusted means have unit
#' variance, so the residual variance of the rotated model without a marker
#' has the target value one when the relationship matrix accounts for all
#' genetic covariance among the accessions. The test uses the increase in the
#' REML deviance when the residual variance is fixed at one, with the genetic
#' variance re-optimised. The one-sided test, the one the paper recommends,
#' takes the increase as zero when the estimate is at or below one and refers
#' it to 2.71, the 0.90 quantile of the chi-squared distribution on one degree
#' of freedom, which gives a one-sided test at the 5\% level; a residual
#' variance above one then indicates genetic covariance not represented by the
#' relationship matrix (see `independent_effect()`). The two-sided test, which
#' also detects a residual variance below one, refers the increase to 3.84.
#' Both can be read from the profile-likelihood intervals returned: the
#' one-sided test rejects when the 90\% interval lies above one, the two-sided
#' test when the 95\% interval excludes one; the intervals use the same
#' cut-offs as the tests, 2.71 and 3.84. The test presumes that the first
#' stage is correctly specified: a first stage that understates the estimation
#' error of the means also raises the residual variance.
#'
#' @param x a `rotation` object.
#' @return an object of class `specification_test`: the residual variance
#'   estimate (`residual`), the genetic variance with the residual estimated
#'   (`genetic`) and fixed at one (`genetic_at_one`), the increase in deviance
#'   (`increase`), the one-sided statistic (`statistic`, the increase or zero),
#'   the critical values (`crit_one_sided` 2.71, `crit_two_sided` 3.84), the
#'   decisions (`reject_one_sided`, `reject_two_sided`), p-values (`p_one_sided`,
#'   from the equal mixture of a point mass at zero and chi-squared on one
#'   degree of freedom; `p_two_sided`, chi-squared on one degree of freedom),
#'   the 90\% and 95\% profile-likelihood intervals (`interval_90`,
#'   `interval_95`, each with `lower`, `upper`, `bounded`), `profile` (a data
#'   frame of the increase in the REML deviance over a grid of values of the
#'   residual variance, drawn by `plot()`), `n`, `p` and the boundary status
#'   of the free fit.
#' @examples
#' ex <- rotationGWAS:::toy_trial(n = 120, n_markers = 80, seed = 2)
#' r <- rotation(ex$means, ex$Omega, ex$relationship, ex$covariates)
#' specification_test(r)
#' @export
specification_test <- function(x) {
  stopifnot(inherits(x, "rotation")); needs_full_omega(x, "the test of the residual variance")
  fr <- x$fit$free; f1 <- x$fit$one
  increase <- max(f1$deviance_ec - fr$deviance_ec, 0)
  above <- fr$ve > 1
  statistic <- if (above) increase else 0
  crit1 <- 2.71; crit2 <- 3.84
  p_one <- if (statistic > 0) 0.5 * stats::pchisq(statistic, 1, lower.tail = FALSE) else 1
  p_two <- stats::pchisq(increase, 1, lower.tail = FALSE)
  ci90 <- reml_profile_ve(x$setup, fr, 0.90, crit = crit1); ci95 <- reml_profile_ve(x$setup, fr, 0.95, crit = crit2)
  keep <- function(ci) list(lower = ci$lower, upper = ci$upper, bounded = ci$bounded, inside_one = ci$inside_one)
  # the profile of the deviance over the residual variance, for the plot: the genetic variance re-optimised at each point
  tot <- fr$vg * mean(x$setup$d) + fr$ve
  hi <- if (ci95$bounded) max(2, 1.3 * ci95$upper, 1.3 * fr$ve) else max(3, 3 * fr$ve)
  grid <- sort(unique(c(seq(0, hi, length.out = 81), 1))); grid[1] <- 1e-8 * tot      # one is always on the grid
  prof <- data.frame(ve = c(0, grid[-1]),
                     increase = pmax(vapply(grid, function(v) reml_fit_fixed(x$setup, v)$deviance_ec - fr$deviance_ec, numeric(1)), 0))
  structure(list(residual = fr$ve, genetic = fr$vg, genetic_at_one = f1$vg, increase = increase, statistic = statistic,
                 crit_one_sided = crit1, crit_two_sided = crit2,
                 reject_one_sided = statistic > crit1, reject_two_sided = increase > crit2,
                 p_one_sided = p_one, p_two_sided = p_two, interval_90 = keep(ci90), interval_95 = keep(ci95), profile = prof,
                 n = x$n, p = x$p, boundary = fr$boundary), class = "specification_test")
}

#' @export
print.specification_test <- function(x, digits = 4, ...) {
  fmt <- function(ci) sprintf("[%s, %s]", signif(ci$lower, digits), if (ci$bounded) signif(ci$upper, digits) else "unbounded")
  cat(sprintf("Residual variance of the rotated model without a marker: %s (target one; genetic variance %s, or %s with the residual fixed at one)%s\n",
              signif(x$residual, digits), signif(x$genetic, digits), signif(x$genetic_at_one, digits),
              if (x$boundary != "none") paste0(" [free fit on the boundary: ", x$boundary, "]") else ""))
  cat(sprintf("Increase in the REML deviance when the residual variance is fixed at one: %s\n", signif(x$increase, digits)))
  cat(sprintf("One-sided test (residual above one): statistic %s against %s, %s (p = %s)\n",
              signif(x$statistic, digits), x$crit_one_sided,
              if (x$reject_one_sided) "REJECTED: the residual variance is above one; genetic covariance is missing from the model, or the first stage understates the estimation error (see independent_effect())"
              else if (x$residual <= 1) "not rejected (the estimate is at or below one, statistic zero)" else "not rejected",
              signif(x$p_one_sided, 3)))
  cat(sprintf("Two-sided test: increase %s against %s, %s (p = %s)\n", signif(x$increase, digits), x$crit_two_sided,
              if (x$reject_two_sided) "REJECTED: the residual variance differs from one" else "not rejected", signif(x$p_two_sided, 3)))
  cat(sprintf("Profile-likelihood intervals for the residual variance: 90%% %s, 95%% %s\n", fmt(x$interval_90), fmt(x$interval_95)))
  invisible(x)
}
