# The data shipped with the package: a simulated population and three rolling-panel trials, built by data-raw/make_panels.R.

#' A simulated population of 150 inbred lines with 1,998 markers
#'
#' The genotypes behind the three rolling-panel trials [panel_dense],
#' [panel_sparse] and [panel_misspecified]. Three subpopulations of eight
#' founders each; fifteen families of ten recombinant inbred lines, each
#' family a biparental cross of two founders of one subpopulation, every line
#' a mosaic of the two parental haplotypes with crossovers along each of ten
#' chromosomes, homozygous apart from two per cent residual heterozygosity.
#' The founder haplotypes carry subpopulation-specific allele frequencies and
#' linkage disequilibrium along the chromosome. Markers with a minor allele
#' frequency below 0.05 were removed. Everything is simulated; nothing comes
#' from any real trial.
#'
#' @format A list:
#' \describe{
#'   \item{markers}{integer matrix of dosages (0, 1, 2) of the counted allele,
#'     150 lines in rows, 1,998 markers in columns, with the line and marker
#'     names as dimnames.}
#'   \item{map}{data frame with the `marker`, its `chromosome` (1 to 10) and
#'     its `position` in centimorgans on a chromosome of 100 cM.}
#'   \item{relationship}{the relationship matrix of VanRaden (2008) computed
#'     from the markers and scaled to a mean diagonal of one, so that the
#'     genetic variance is the variance of a line's genetic value.}
#'   \item{covariates}{the first three principal components of the centred
#'     markers, with the line names as row names.}
#'   \item{family}{the family of each line, named by line.}
#'   \item{subpopulation}{the subpopulation of each line's founders, named by
#'     line.}
#' }
#' @source Simulated by `data-raw/make_panels.R` in the package sources.
#' @references VanRaden, P. M., 2008 Efficient methods to compute genomic
#'   predictions. Journal of Dairy Science 91: 4414-4423.
#' @examples
#' str(panel_genotypes$markers)
#' table(panel_genotypes$family, panel_genotypes$subpopulation)
"panel_genotypes"

#' Three simulated rolling-panel trials, analysed in their first stage
#'
#' Plot-level data of three multi-environment trials of the 150 lines of
#' [panel_genotypes], together with the adjusted means and their
#' variance-covariance matrix from the first stage, so that the second stage
#' of the package can be run on them directly and the first stage can be
#' repeated from the plots. The trials are rolling panels of the kind the
#' paper's scenario C uses: each line is grown in consecutive environments,
#' its window rolling evenly over the environments, with three replicates in
#' each environment.
#'
#' - `panel_dense`: 8 environments, three per line, 1,350 plots, 56 lines
#'   per environment on average. The adjusted means are estimated with nearly
#'   equal precision and are nearly uncorrelated, so the three weighting
#'   schemes agree.
#' - `panel_sparse`: 60 environments, three per line, 1,350 plots, 7.5 lines
#'   per environment on average, the same lines, replication and number of
#'   plots as the dense panel. The environment effects are estimated from few
#'   lines, so the adjusted means are imprecise and correlated and the
#'   weighting matters.
#' - `panel_misspecified`: the dense design with unequal replication, one to
#'   six environments per line (1,437 plots), and a genetic effect independent
#'   across lines added to the genetic value, with variance 6 against the
#'   genetic variance of 10 carried by the relationship matrix. The
#'   relationship matrix does not represent this effect, so the residual
#'   variance of the rotated model lies above one. The unequal replication
#'   helps to separate the effect from the residual after the rotation: its
#'   rotated kernel is the inverse of the variance-covariance matrix of the
#'   means, and the separation is weak when that matrix is close to
#'   proportional to the identity, as it is on the dense design with equal
#'   replication, where many lines share each environment.
#'
#' The plot values were simulated as the sum of an intercept of 50, an
#' environment effect (standard deviation 4), the line's genetic value, a
#' replicate-within-environment effect (variance 0.5), a line-by-environment
#' interaction (variance 3) and a plot error (variance 6). The genetic value
#' is a polygenic part with covariance 10 times the relationship matrix, the
#' effects of three causal markers on chromosomes 2, 5 and 8 (2.0, 1.2 and 0.6
#' per dose of the counted allele), a regression on the three principal
#' components, and, in the misspecified panel, the independent effect.
#'
#' The first stage fitted the paper's model to the plots with lme4, line and
#' environment as fixed effects and replicate within environment and
#' line-by-environment as random effects, and took the estimated marginal
#' means of the lines and their variance-covariance matrix from emmeans with
#' the Kenward-Roger adjustment (on these designs, which are balanced within
#' each line-by-environment cell, the adjustment leaves the matrix unchanged).
#'
#' @format Each panel is a list:
#' \describe{
#'   \item{plots}{data frame of the plot data: `accession`, `environment`,
#'     `replicate` (factors) and the response `y`.}
#'   \item{means}{the adjusted means of the lines, named by line.}
#'   \item{Omega}{their variance-covariance matrix from the first stage, with
#'     the line names as dimnames.}
#'   \item{design}{the number of environments, the environments per line, the
#'     replicates, the number of plots and the mean number of lines per
#'     environment.}
#'   \item{first_stage}{the model formula, its variance components and the
#'     Kenward-Roger inflation factor of the variances of the means.}
#'   \item{truth}{the generating values: `mu`, the `causal` markers and their
#'     `effect` per dose, the variances (`genetic_variance`,
#'     `independent_variance`, `interaction_variance`, `replicate_variance`,
#'     `error_variance`), `environment_sd`, `pc_coefficient`, the simulated
#'     genetic `value` of each line and the `environment_effect`s.}
#' }
#' @source Simulated by `data-raw/make_panels.R` in the package sources.
#' @examples
#' r <- rotation(panel_sparse$means, panel_sparse$Omega,
#'               panel_genotypes$relationship, panel_genotypes$covariates)
#' r
#' efficiency(r, panel_genotypes$markers)
#' @name panels
#' @aliases panel_dense panel_sparse panel_misspecified
NULL

#' @rdname panels
"panel_dense"

#' @rdname panels
"panel_sparse"

#' @rdname panels
"panel_misspecified"
