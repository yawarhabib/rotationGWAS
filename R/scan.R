# The genome scan: the P3D single-marker test under the rotation, and under the two approximate schemes for comparison.

#' Genome scan of the rotated second stage
#'
#' Tests every marker one at a time with the variance components held at
#' their values from the model without a marker (the population parameters
#' previously determined approximation): the ratio of the residual to the
#' genetic variance is fixed, the residual scale is re-estimated at each
#' marker, and the marker effect is tested by a t statistic on
#' \eqn{n - p - 1} degrees of freedom. Under the rotation every marker column
#' is rotated with the same matrix as the means, so the scan is the ordinary
#' mixed-model GWAS of the rotated data. The two approximate schemes are
#' offered for comparison: inverse-variance weighting divides each mean by
#' its own standard error and treats the estimation errors as uncorrelated;
#' the unweighted analysis ignores the first-stage variance-covariance matrix
#' altogether. Each scheme fits its own model without a marker.
#'
#' The residual variance is estimated in the scan under every scheme, also
#' under the rotation, where its target is one; the test of the residual
#' variance is a separate step (`specification_test()`). Markers close to the significance threshold can
#' be refitted afterwards with the marker in the model and the residual
#' variance fixed at one (`refit_markers()`). When the test finds the residual
#' variance above one, the scan can be run under the model extended by an
#' effect independent across accessions, from [independent_effect()], through
#' `model`; that model is defined under the rotation only, and its variance
#' components are held at their fitted values in the same way, either from
#' its free fit or, with `fit = "residual_one"`, from its fit with the
#' residual variance fixed at one. The paper adds the effect to the model
#' when the test rejects but does not evaluate a scan under the extended
#' model; the option is offered here so that the scan can follow the model
#' the test led to.
#'
#' @param x a `rotation` object.
#' @param markers numeric matrix of marker dosages (0, 1, 2), accessions in
#'   rows with the accession names as row names, markers in columns; complete.
#' @param scheme the scheme or schemes to scan under: `"rotated"` (the
#'   default), `"inverse_variance"`, `"unweighted"`, or several of them.
#' @param alpha the genome-wide level of the Bonferroni threshold, applied to
#'   the number of markers with a test.
#' @param chunk_size markers per block of matrix products.
#' @param model optionally an `independent_effect` object: the scan then runs
#'   under the rotated model with the independent effect, at its fitted
#'   variance components (`scheme` must be `"rotated"`). A warning is given
#'   when that model's free fit has the residual variance at the boundary of
#'   zero, where the model is the unweighted analysis.
#' @param fit with `model`, which of its fits supplies the variance components
#'   held in the scan: the free fit (`"free"`, the default) or the fit with
#'   the residual variance fixed at one (`"residual_one"`).
#' @return a data frame of class `genome_scan` with one row per marker and
#'   scheme: `marker`, `scheme`, `maf`, `beta` (the marker effect per dose),
#'   `se`, `t`, `neglog10p` and `significant` (above the Bonferroni
#'   threshold). A marker without variance has `NA` in the test columns.
#'   Attributes: `fits`, a data frame with one row per scheme (the genetic,
#'   independent and residual variances of the model without a marker, the
#'   ratio of residual to genetic variance, the boundary status, the genomic
#'   inflation factor of the scan, the numbers of markers tested and
#'   significant), `threshold` (in \eqn{-\log_{10} p}), `alpha`, `n`, `p`,
#'   `df` and `model` (`"additive"`, `"additive + independent effect"` or
#'   `"additive + independent effect, residual variance fixed at one"`).
#' @examples
#' ex <- rotationGWAS:::toy_trial(n = 40, n_markers = 60, seed = 2)
#' r <- rotation(ex$means, ex$Omega, ex$relationship, ex$covariates)
#' s <- genome_scan(r, ex$markers, scheme = c("rotated", "unweighted"))
#' summary(s)
#' head(s[order(-s$neglog10p), ])
#' @export
genome_scan <- function(x, markers, scheme = "rotated", alpha = 0.05, chunk_size = 5000, model = NULL, fit = c("free", "residual_one")) {
  stopifnot(inherits(x, "rotation"))
  fit <- match.arg(fit)
  schemes <- c("rotated", "inverse_variance", "unweighted")
  if (!is.character(scheme) || !length(scheme) || !all(scheme %in% schemes))
    stop("scheme must be one or more of: ", paste(schemes, collapse = ", "))
  scheme <- unique(scheme)
  if (!is.numeric(alpha) || length(alpha) != 1 || !is.finite(alpha) || alpha <= 0 || alpha >= 1) stop("alpha must be one number in (0, 1)")
  if (!is.null(model)) {
    if (!inherits(model, "independent_effect")) stop("model must be an independent_effect object")
    if (!identical(scheme, "rotated")) stop("the model with the independent effect is defined under the rotation only; use scheme = \"rotated\"")
    if (model$n != x$n || !isTRUE(all.equal(model$setup$y, x$y))) stop("the model was fitted on another rotation object")
    if (fit == "free" && "ve_zero" %in% model$fit$boundary)
      warning("the residual variance of the model with the independent effect is at the boundary of zero; under it this scan is the unweighted analysis of the adjusted means")
    if (fit == "free" && "vi_zero" %in% model$fit$boundary) message("the independent effect is zero: this scan equals the scan under the additive model")
    if (fit == "residual_one" && "vi_zero" %in% model$fit_one$boundary)
      message("the independent effect is zero with the residual variance at one: this scan holds the additive model's components from its fit with the residual variance at one")
  } else if (fit != "free") stop("fit applies to the scan under the model with the independent effect; give `model`")
  Z <- check_markers(markers, x$ids)
  n <- x$n; m <- ncol(Z)
  maf <- pmin(colMeans(Z) / 2, 1 - colMeans(Z) / 2)
  X0 <- cbind(intercept = rep(1, n), x$covariates)
  res <- vector("list", length(scheme)); fits <- vector("list", length(scheme)); names(res) <- names(fits) <- scheme
  for (sch in scheme) {
    if (is.null(model)) {
      tr <- scheme_fit(x, sch, X0); st <- tr$setup; fr <- tr$fit
      pr <- scan_projection_safe(st, fr$lambda)
      HT <- if (is.null(tr$Tm)) pr$H else pr$H %*% tr$Tm        # the markers enter transformed like every other column
      df <- st$df - 1; vi <- NA_real_; boundary <- fr$boundary
    } else {
      fr <- if (fit == "free") model$fit else model$fit_one
      pr <- scan_projection2(model$setup, fr$vg, fr$vi, fr$ve)
      HT <- pr$H %*% x$W
      df <- model$setup$df - 1; vi <- fr$vi; boundary <- if (length(fr$boundary)) paste(fr$boundary, collapse = ",") else "none"
    }
    beta <- se <- tval <- lp <- rep(NA_real_, m)
    for (from in seq(1, m, by = chunk_size)) {
      j <- from:min(from + chunk_size - 1, m)
      Hz <- HT %*% Z[, j, drop = FALSE]
      zPz <- colSums(Hz^2); zPy <- as.numeric(crossprod(Hz, pr$Hy))
      ok <- zPz > 1e-8 * max(zPz, 1e-300)
      idx <- j[ok]
      b <- zPy[ok] / zPz[ok]
      s2 <- (pr$yPy - b * zPy[ok]) / df
      beta[idx] <- b; se[idx] <- sqrt(s2 / zPz[ok])
    }
    ok <- is.finite(beta)
    tval[ok] <- beta[ok] / se[ok]
    lp[ok] <- -log10(2 * stats::pt(-abs(tval[ok]), df))
    res[[sch]] <- data.frame(marker = colnames(Z), scheme = sch, maf = maf, beta = beta, se = se, t = tval, neglog10p = lp,
                             row.names = NULL, stringsAsFactors = FALSE)
    fits[[sch]] <- data.frame(scheme = sch, vg = fr$vg, vi = vi, ve = fr$ve, lambda = if (fr$vg > 0) fr$ve / fr$vg else Inf, boundary = boundary,
                              inflation = inflation_factor(lp), n_tested = sum(ok), n_significant = NA_integer_,
                              row.names = NULL, stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, res)
  n_tested <- fits[[1]]$n_tested
  if (n_tested == 0) stop("no marker has variance among the accessions")
  threshold <- -log10(alpha / n_tested)
  out$significant <- !is.na(out$neglog10p) & out$neglog10p > threshold
  fits <- do.call(rbind, fits)
  fits$n_significant <- as.integer(vapply(scheme, function(s) sum(out$significant[out$scheme == s]), numeric(1)))
  rownames(out) <- NULL; rownames(fits) <- NULL
  structure(out, class = c("genome_scan", "data.frame"), fits = fits, threshold = threshold, alpha = alpha,
            n = n, p = ncol(X0), df = x$setup$df - 1,
            model = if (is.null(model)) "additive" else if (fit == "free") "additive + independent effect" else "additive + independent effect, residual variance fixed at one")
}

#' @export
summary.genome_scan <- function(object, top = 5, ...) {
  fits <- attr(object, "fits")
  plain <- as.data.frame(unclass(object), stringsAsFactors = FALSE)          # a plain data frame, so that the pieces print as such
  best <- lapply(split(plain, plain$scheme), function(d) {
    d <- d[order(-d$neglog10p), , drop = FALSE]; utils::head(d[!is.na(d$neglog10p), c("marker", "beta", "se", "neglog10p")], top) })
  structure(list(fits = fits, threshold = attr(object, "threshold"), alpha = attr(object, "alpha"), model = attr(object, "model"),
                 n = attr(object, "n"), p = attr(object, "p"), df = attr(object, "df"), top = best), class = "summary.genome_scan")
}

#' @export
print.summary.genome_scan <- function(x, digits = 4, ...) {
  cat(sprintf("Genome scan of %d accessions, %d fixed effects, %d degrees of freedom per marker; Bonferroni threshold %s in -log10 p at alpha %s%s\n",
              x$n, x$p, x$df, signif(x$threshold, 4), x$alpha, if (identical(x$model, "additive")) "" else paste0("; model: ", x$model)))
  f <- x$fits
  tab <- data.frame(scheme = f$scheme, genetic = signif(f$vg, digits), independent = signif(f$vi, digits), residual = signif(f$ve, digits),
                    boundary = f$boundary, inflation = round(f$inflation, 4), tested = f$n_tested, significant = f$n_significant)
  if (all(is.na(tab$independent))) tab$independent <- NULL
  print(tab, row.names = FALSE)
  for (s in names(x$top)) {
    cat(sprintf("\nStrongest markers, %s:\n", s)); d <- x$top[[s]]
    d$beta <- signif(d$beta, digits); d$se <- signif(d$se, digits); d$neglog10p <- round(d$neglog10p, 3)
    print(d, row.names = FALSE)
  }
  invisible(x)
}

#' @export
print.genome_scan <- function(x, ...) { print(summary(x)); invisible(x) }

#' @export
`[.genome_scan` <- function(x, ...) { out <- NextMethod(); class(out) <- "data.frame"; out }   # a subset is a plain data frame

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

# the transformation of a scheme, its setup and its model without a marker; the rotation reuses the object's own
scheme_fit <- function(x, scheme, X0) {
  if (scheme == "rotated") return(list(Tm = x$W, setup = x$setup, fit = x$fit$free))
  if (scheme == "inverse_variance") {
    w <- 1 / sqrt(diag(x$Omega)); Tm <- diag(w, x$n)
    y <- as.numeric(x$means * w); X <- X0 * w; K <- x$relationship * tcrossprod(w)
  } else {
    Tm <- NULL; y <- as.numeric(x$means); X <- X0; K <- x$relationship
  }
  st <- reml_setup(y, X, K)
  list(Tm = Tm, setup = st, fit = reml_fit_free(st))
}

# the scan projection at the fitted ratio; at a genetic variance of zero the ratio is infinite and the projection is
# the plain residual projection, whose scale cancels in the t statistic
scan_projection_safe <- function(s, lambda) {
  if (is.infinite(lambda)) return(list(H = s$B, Hy = s$eta, yPy = sum(s$eta^2)))
  scan_projection(s, lambda)
}

# the genomic inflation factor of a scan from its -log10 p values, as the paper's script 05 computes it
inflation_factor <- function(lp) {
  pv <- 10^(-lp); pv <- pv[!is.na(pv) & pv > 0 & pv <= 1]
  if (!length(pv)) return(NA_real_)
  stats::qchisq(stats::median(pv), 1, lower.tail = FALSE) / stats::qchisq(0.5, 1)
}
