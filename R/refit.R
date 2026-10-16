# The refit: chosen markers fitted in the model with the residual variance fixed at one.

#' Refit chosen markers with the residual variance fixed at one
#'
#' For each named marker, the rotated model is fitted by REML with the marker
#' among the fixed effects and the residual variance fixed at its target of
#' one, the genetic variance re-estimated. The marker effect, its standard
#' error from \eqn{(X'V^{-1}X)^{-1}} without a scale factor, and a t test on
#' \eqn{n - p - 1} degrees of freedom are returned. With known variance
#' components this fit reproduces the single-stage analysis of the plot data,
#' which the genome scan, with its residual variance estimated and its
#' variance ratio held at the model without a marker, approximates. The
#' paper's rule is to refit the markers whose scan statistic lies within a
#' window of the significance threshold, on either side; given a
#' [genome_scan()] object those markers are selected here.
#'
#' @param x a `rotation` object.
#' @param markers numeric matrix of marker dosages (0, 1, 2), accessions in
#'   rows with the accession names as row names, markers in columns; complete.
#' @param which the markers to refit: a character vector of marker names (or
#'   column indices). If `NULL`, `scan` must be given and the markers within
#'   `window` of its threshold under the rotation are refitted.
#' @param scan optionally the `genome_scan` object of the same markers; its
#'   threshold is used for `significant`, its rotated statistic is carried
#'   along, and it selects the markers when `which` is `NULL`.
#' @param window the half-width, in \eqn{-\log_{10} p}, of the band around
#'   the threshold from which markers are selected when `which` is `NULL`.
#' @param threshold the significance threshold in \eqn{-\log_{10} p} when no
#'   `scan` is given; `NULL` leaves `significant` as `NA`.
#' @return a data frame of class `refit`, one row per marker: `marker`,
#'   `scan_neglog10p` (when a scan is given), `genetic` (the genetic variance
#'   of the refit), `beta`, `se`, `t`, `neglog10p`, `significant`, and
#'   `boundary` (whether the genetic variance went to zero). Attributes
#'   `threshold`, `window`, `df`, `n`, `p`.
#' @examples
#' ex <- rotationGWAS:::toy_trial(n = 120, n_markers = 80, seed = 2)
#' r <- rotation(ex$means, ex$Omega, ex$relationship, ex$covariates)
#' s <- genome_scan(r, ex$markers)
#' refit_markers(r, ex$markers, scan = s, window = 1.5)
#' refit_markers(r, ex$markers, which = ex$truth$causal, scan = s)
#' @export
refit_markers <- function(x, markers, which = NULL, scan = NULL, window = 1.5, threshold = NULL) {
  stopifnot(inherits(x, "rotation"))
  Z <- check_markers(markers, x$ids)
  if (!is.null(scan)) {
    if (!inherits(scan, "genome_scan")) stop("scan must be a genome_scan object")
    rot <- scan[scan$scheme == "rotated", ]
    if (!nrow(rot)) stop("the scan holds no rotated scheme; refit selection follows the rotated scan")
    if (!all(rot$marker %in% colnames(Z))) stop("the scan and markers do not hold the same markers")
    threshold <- attr(scan, "threshold")
  }
  if (is.null(which)) {
    if (is.null(scan)) stop("name the markers in `which`, or give a genome_scan to select them by `window`")
    if (!is.numeric(window) || length(window) != 1 || window < 0) stop("window must be one non-negative number")
    which <- rot$marker[!is.na(rot$neglog10p) & abs(rot$neglog10p - threshold) <= window]
    if (!length(which)) message("no marker within ", window, " of the threshold ", signif(threshold, 4), " under the rotation")
  }
  if (is.numeric(which)) which <- colnames(Z)[which]
  if (!is.character(which)) stop("which must name markers")
  missing <- setdiff(which, colnames(Z))
  if (length(missing)) stop("marker(s) not in markers: ", paste(utils::head(missing, 5), collapse = ", "))
  which <- unique(which)
  X0 <- x$X; y <- x$y; K <- x$K; W <- x$W; n <- x$n
  out <- lapply(which, function(mk) {
    z <- as.numeric(W %*% Z[, mk])
    X <- cbind(X0, marker = z)
    if (qr(X)$rank < ncol(X)) return(data.frame(marker = mk, genetic = NA_real_, beta = NA_real_, se = NA_real_, t = NA_real_,
                                                  neglog10p = NA_real_, boundary = NA, stringsAsFactors = FALSE))
    st <- reml_setup(y, X, K)
    f <- reml_fit_fixed(st, 1)
    g <- gls_fixed(y, X, K, f$vg, 1)
    k <- ncol(X); b <- g$beta[k]; s <- g$se[k]; tt <- b / s
    data.frame(marker = mk, genetic = f$vg, beta = b, se = s, t = tt, neglog10p = -log10(2 * stats::pt(-abs(tt), st$df)),
               boundary = f$boundary != "none", stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out)
  if (is.null(out)) out <- data.frame(marker = character(0), genetic = numeric(0), beta = numeric(0), se = numeric(0), t = numeric(0),
                                      neglog10p = numeric(0), boundary = logical(0))
  if (!is.null(scan)) out <- cbind(out[1], scan_neglog10p = rot$neglog10p[match(out$marker, rot$marker)], out[-1])
  out$significant <- if (is.null(threshold)) NA else !is.na(out$neglog10p) & out$neglog10p > threshold
  rownames(out) <- NULL
  structure(out, class = c("refit", "data.frame"), threshold = threshold, window = if (is.null(scan)) NA_real_ else window,
            df = n - x$p - 1, n = n, p = x$p)
}

#' @export
`[.refit` <- function(x, ...) { out <- NextMethod(); class(out) <- "data.frame"; out }   # a subset is a plain data frame

#' @export
print.refit <- function(x, digits = 4, ...) {
  thr <- attr(x, "threshold")
  cat(sprintf("Refit of %d marker%s with the residual variance fixed at one (t on %d degrees of freedom)%s\n",
              nrow(x), if (nrow(x) == 1) "" else "s", attr(x, "df"),
              if (!is.null(thr)) sprintf("; threshold %s in -log10 p, %d significant", signif(thr, 4), sum(x$significant, na.rm = TRUE)) else ""))
  d <- as.data.frame(unclass(x), stringsAsFactors = FALSE)
  for (cn in intersect(c("scan_neglog10p", "genetic", "beta", "se", "t", "neglog10p"), names(d))) d[[cn]] <- signif(d[[cn]], digits)
  print(d, row.names = FALSE)
  invisible(x)
}
