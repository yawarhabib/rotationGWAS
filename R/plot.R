# Plot methods: the structure of Omega, the efficiency by allele frequency or along the genome, the scan as a Manhattan
# or a QQ plot, and the profile of the specification test. Base graphics only.

#' Plot methods for the package's objects
#'
#' @description
#' `plot()` on a `rotation` object shows the structure of the variance-covariance matrix of the adjusted means: the
#' absolute correlations as an image in the order of the accessions (`type = "correlation"`), or the variances along
#' the accessions (`type = "variance"`).
#'
#' On an `efficiency` object it shows the standard-error factor of inverse-variance weighting and of the unweighted
#' analysis for every marker, against the minor allele frequency (`type = "maf"`) or along the genome
#' (`type = "genome"`, which needs `map`).
#'
#' On a `genome_scan` object it draws a Manhattan plot per scheme (`type = "manhattan"`, positions from `map` or the
#' marker order without one) with the Bonferroni threshold, the significant markers in red and the `highlight`ed
#' markers circled; or a QQ plot of the observed against the expected \eqn{-\log_{10} p} of every scheme in the
#' object (`type = "qq"`), with the genomic inflation factors in the legend.
#'
#' On a `specification_test` object it draws the profile of the REML deviance of the residual variance: the increase
#' over its minimum as the residual variance moves, with the estimate, the target value one, the one-sided and
#' two-sided critical values and the 90\% and 95\% profile-likelihood intervals.
#'
#' @param x the object.
#' @param type which plot, see the description.
#' @param map optional data frame with the columns `marker`, `chromosome` and `position` (any unit), for the plots
#'   along the genome; the chromosomes are drawn in the order of their first appearance in `map`.
#' @param scheme for the Manhattan plot, the scheme or schemes to draw, one panel each; by default every scheme in
#'   the object.
#' @param highlight marker names to circle in the Manhattan plot.
#' @param main,xlab titles; `NULL` gives the defaults described above.
#' @param ... further graphical parameters for the main plotting call.
#' @return the data drawn, invisibly: the absolute correlation matrix or the variances; a data frame of the factors
#'   with their x coordinate; a data frame per scheme of the Manhattan coordinates, or a data frame per scheme of
#'   expected and observed \eqn{-\log_{10} p}; the profile.
#' @examples
#' ex <- rotationGWAS:::toy_trial(n = 60, n_markers = 200, seed = 2)
#' r <- rotation(ex$means, ex$Omega, ex$relationship, ex$covariates)
#' plot(r)
#' plot(efficiency(r, ex$markers))
#' s <- genome_scan(r, ex$markers, scheme = c("rotated", "unweighted"))
#' plot(s, type = "qq")
#' plot(s, scheme = "rotated", highlight = ex$truth$causal)
#' plot(specification_test(r))
#' @name plot-methods
NULL

#' @rdname plot-methods
#' @export
plot.rotation <- function(x, type = c("correlation", "variance"), main = NULL, ...) {
  type <- match.arg(type)
  n <- x$n
  if (type == "correlation") {
    R <- abs(stats::cov2cor(x$Omega))
    off <- R[upper.tri(R)]
    cols <- grDevices::hcl.colors(64, "YlOrRd", rev = TRUE)
    if (is.null(main)) main <- sprintf("|correlation| of the adjusted means: mean %.2f, largest %.2f", mean(off), max(off))
    graphics::image(seq_len(n), seq_len(n), R, col = cols, zlim = c(0, 1), xlab = "accession", ylab = "accession", main = main, ...)
    graphics::legend("topright", legend = c("1", "0.75", "0.5", "0.25", "0"), fill = cols[c(64, 48, 32, 16, 1)],
                     bty = "n", cex = 0.8, title = "|r|")
    invisible(R)
  } else {
    v <- diag(x$Omega)
    if (is.null(main)) main <- sprintf("variances of the adjusted means: %.3g to %.3g (ratio %.3g)", min(v), max(v), max(v) / min(v))
    graphics::plot(seq_len(n), v, type = "h", xlab = "accession", ylab = "variance of the adjusted mean", main = main,
                   ylim = c(0, max(v) * 1.05), ...)
    invisible(v)
  }
}

#' @rdname plot-methods
#' @export
plot.efficiency <- function(x, type = c("maf", "genome"), map = NULL, main = NULL, xlab = NULL, ...) {
  type <- match.arg(type)
  d <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  if (type == "maf") {
    xx <- d$maf; if (is.null(xlab)) xlab <- "minor allele frequency"; gp <- NULL
  } else {
    if (is.null(map)) stop("map must be given for type = \"genome\"")
    gp <- genome_positions(d$marker, map); xx <- gp$x; if (is.null(xlab)) xlab <- gp$label
  }
  ok <- !is.na(d$inverse_variance)
  top <- max(c(d$inverse_variance, d$unweighted), na.rm = TRUE)
  if (is.null(main)) main <- sprintf("standard error relative to the rotation (medians %.3f and %.3f)",
                                     stats::median(d$inverse_variance, na.rm = TRUE), stats::median(d$unweighted, na.rm = TRUE))
  graphics::plot(xx[ok], d$unweighted[ok], pch = 16, cex = 0.6, col = "firebrick", ylim = c(1, max(top, 1.05)),
                 xlab = xlab, ylab = "factor on the standard error", main = main, xaxt = if (is.null(gp)) "s" else "n", ...)
  if (!is.null(gp)) graphics::axis(1, at = gp$centres, labels = names(gp$centres))
  graphics::points(xx[ok], d$inverse_variance[ok], pch = 16, cex = 0.6, col = "royalblue")
  graphics::abline(h = 1, lty = 3)
  graphics::legend("topright", legend = c("unweighted", "inverse-variance weighting"), col = c("firebrick", "royalblue"), pch = 16, bty = "n", cex = 0.85)
  invisible(data.frame(marker = d$marker, x = xx, inverse_variance = d$inverse_variance, unweighted = d$unweighted, stringsAsFactors = FALSE))
}

#' @rdname plot-methods
#' @export
plot.genome_scan <- function(x, type = c("manhattan", "qq"), map = NULL, scheme = NULL, highlight = NULL, main = NULL, ...) {
  type <- match.arg(type)
  d <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  schemes <- unique(d$scheme)
  if (is.null(scheme)) scheme <- schemes
  if (!all(scheme %in% schemes)) stop("scheme not in the scan: ", paste(setdiff(scheme, schemes), collapse = ", "))
  thr <- attr(x, "threshold"); fits <- attr(x, "fits")
  if (type == "manhattan") {
    first <- d[d$scheme == scheme[1], ]
    gp <- genome_positions(first$marker, map)
    if (!is.null(highlight)) {
      h <- match(highlight, first$marker)
      if (anyNA(h)) stop("highlight marker(s) not in the scan: ", paste(highlight[is.na(h)], collapse = ", "))
    }
    if (length(scheme) > 1) { op <- graphics::par(mfrow = c(length(scheme), 1), mar = c(4, 4.2, 2.5, 1)); on.exit(graphics::par(op)) }
    out <- list()
    for (s in scheme) {
      ds <- d[d$scheme == s, ]; lp <- ds$neglog10p
      sig <- !is.na(lp) & lp > thr
      col <- c("grey25", "grey60")[(gp$chr_index %% 2) + 1]
      graphics::plot(gp$x, lp, pch = 16, cex = 0.6, col = col, xaxt = if (is.null(map)) "s" else "n", xlab = gp$label,
                     ylab = expression(-log[10] ~ p), main = if (is.null(main)) scheme_label(s) else main, ...)
      if (!is.null(map)) graphics::axis(1, at = gp$centres, labels = names(gp$centres))
      graphics::abline(h = thr, lty = 2)
      if (any(sig)) graphics::points(gp$x[sig], lp[sig], pch = 16, cex = 0.8, col = "firebrick")
      if (!is.null(highlight)) graphics::points(gp$x[h], lp[h], pch = 1, cex = 1.8, lwd = 2, col = "royalblue")
      out[[s]] <- data.frame(marker = ds$marker, x = gp$x, neglog10p = lp, significant = sig, stringsAsFactors = FALSE)
    }
    invisible(if (length(out) == 1) out[[1]] else out)
  } else {
    cols <- c(rotated = "black", inverse_variance = "royalblue", unweighted = "firebrick")
    out <- lapply(scheme, function(s) {
      obs <- d$neglog10p[d$scheme == s]; obs <- sort(obs[!is.na(obs)], decreasing = TRUE)
      data.frame(expected = -log10(stats::ppoints(length(obs))), observed = obs)
    })
    names(out) <- scheme
    lim <- c(0, max(unlist(lapply(out, function(q) max(q$expected, q$observed)))))
    if (is.null(main)) main <- "observed against expected -log10 p"
    graphics::plot(lim, lim, type = "n", xlab = expression(expected ~ -log[10] ~ p), ylab = expression(observed ~ -log[10] ~ p), main = main, ...)
    graphics::abline(0, 1, lty = 2)
    for (s in scheme) graphics::points(out[[s]]$expected, out[[s]]$observed, pch = 16, cex = 0.6, col = cols[[s]])
    infl <- fits$inflation[match(scheme, fits$scheme)]
    leg <- lapply(seq_along(scheme), function(i) bquote(.(scheme_label(scheme[i])) ~ (lambda == .(sprintf("%.3f", infl[i])))))
    graphics::legend("topleft", legend = as.expression(leg), col = cols[scheme], pch = 16, bty = "n", cex = 0.85)
    invisible(if (length(out) == 1) out[[1]] else out)
  }
}

#' @rdname plot-methods
#' @export
plot.specification_test <- function(x, main = NULL, ...) {
  p <- x$profile
  if (is.null(p)) stop("the object holds no profile")
  top <- min(max(p$increase, na.rm = TRUE), 3 * 3.84)
  if (is.null(main)) main <- sprintf("residual variance %.3g; one-sided test %s, two-sided %s", x$residual,
                                     if (x$reject_one_sided) "rejected" else "not rejected", if (x$reject_two_sided) "rejected" else "not rejected")
  graphics::plot(p$ve, p$increase, type = "l", lwd = 2, xlab = "residual variance of the rotated model", ylab = "increase in the REML deviance",
                 ylim = c(0, top), main = main, ...)
  graphics::abline(v = 1, col = "grey50")
  graphics::abline(h = x$crit_one_sided, lty = 3); graphics::abline(h = x$crit_two_sided, lty = 2)
  end <- function(ci) if (ci$bounded) ci$upper else max(p$ve)
  graphics::segments(x$interval_90$lower, x$crit_one_sided, end(x$interval_90), x$crit_one_sided, lwd = 3, col = "royalblue")
  graphics::segments(x$interval_95$lower, x$crit_two_sided, end(x$interval_95), x$crit_two_sided, lwd = 3, col = "firebrick")
  graphics::points(x$residual, 0, pch = 16)
  graphics::legend("top", legend = c("profile of the deviance", "estimate", "target value one", "2.71 and the 90% interval", "3.84 and the 95% interval"),
                   lty = c(1, NA, 1, 3, 2), pch = c(NA, 16, NA, NA, NA), lwd = c(2, NA, 1, 3, 3), col = c("black", "black", "grey50", "royalblue", "firebrick"),
                   bg = "white", box.col = "grey80", cex = 0.75)
  invisible(p)
}

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

scheme_label <- function(s) c(rotated = "rotated", inverse_variance = "inverse-variance weighting", unweighted = "unweighted")[[s]]

# cumulative genome coordinates of the markers from a map; without a map, the marker order
genome_positions <- function(markers, map) {
  if (is.null(map)) return(list(x = seq_along(markers), chr_index = rep(1L, length(markers)), centres = NULL, label = "marker"))
  if (!is.data.frame(map) || !all(c("marker", "chromosome", "position") %in% names(map)))
    stop("map must be a data frame with the columns marker, chromosome and position")
  i <- match(markers, map$marker)
  if (anyNA(i)) stop(sum(is.na(i)), " marker(s) not in map, for example ", paste(utils::head(markers[is.na(i)], 3), collapse = ", "))
  chr <- as.character(map$chromosome[i]); pos <- as.numeric(map$position[i])
  levels <- unique(as.character(map$chromosome))
  chr_index <- match(chr, levels)
  start <- tapply(pos, factor(chr, levels = levels), min); len <- tapply(pos, factor(chr, levels = levels), max) - start
  start[is.na(start)] <- 0; len[is.na(len)] <- 0
  gap <- 0.02 * sum(len)
  offset <- cumsum(c(0, utils::head(len + gap, -1)))
  x <- offset[chr_index] + (pos - start[chr_index])
  centres <- offset + len / 2; names(centres) <- levels
  list(x = as.numeric(x), chr_index = chr_index, centres = centres[!is.na(len) & len >= 0], label = "chromosome")
}
