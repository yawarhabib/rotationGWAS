# rotationGWAS 0.9.0

First release, accompanying the manuscript "Exact weighting in two-stage genome-wide association analysis: when does it matter?" (Habib and Piepho, 2026).

- `rotation()`: the inputs aligned by accession name, the rotation matrix as the symmetric inverse square root of the variance-covariance matrix of the adjusted means, the means, the intercept, the covariates and the relationship matrix rotated with it, and the rotated model without a marker fitted with the residual variance estimated and fixed at one. A diagonal matrix (standard errors only) is accepted and labelled inverse-variance weighting.
- `efficiency()` and `estimation_share()`: the standard error of every marker's effect under inverse-variance weighting and under the unweighted analysis relative to the rotation, and the share of the variance of the means that is estimation error, both before any scan.
- `genome_scan()`: the P3D scan under the rotation and, for comparison, under the two approximate schemes; also under the model with the independent effect, with its components from the free fit or from the fit with the residual variance fixed at one.
- `specification_test()`: the residual variance of the rotated model against its target value of one, one-sided and two-sided, with profile-likelihood intervals.
- `independent_effect()`: the model extended by an effect independent across accessions, fitted by a two-kernel REML engine and compared with the additive model.
- `refit_markers()`: chosen markers, or those within a window of the threshold, refitted with the marker in the model and the residual variance fixed at one.
- Plot methods: the structure of the variance-covariance matrix of the means (`rotation`), the factors by allele frequency or along the genome (`efficiency`), Manhattan and QQ plots (`genome_scan`), and the profile of the deviance of the residual variance with the intervals (`specification_test`).
- Data: a simulated population of 150 inbred lines (`panel_genotypes`) and three rolling-panel trials analysed in their first stage (`panel_dense`, `panel_sparse`, `panel_misspecified`).
- Vignette: the steps of the paper's Conclusions on the shipped trials, starting with the first stage in lme4 and emmeans.
- One line on attaching the package, silenced by `suppressPackageStartupMessages()`.
