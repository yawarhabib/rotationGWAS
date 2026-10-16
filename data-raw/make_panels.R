# Builds the data shipped with rotationGWAS: one simulated population of inbred accessions with markers, and three
# rolling-panel field trials, each simulated at the plot level and analysed in its first stage. Everything here is
# simulated from the seed below; nothing comes from any real trial. Run from the package root; writes data/*.rda.
# The first stage needs lme4, emmeans and pbkrtest (all in Suggests).
#
# The dense and the sparse panel follow the synthetic designs of scenario C of the paper: every accession is grown
# in three consecutive environments with three replicates, the windows rolling evenly over the environments, so that
# the number of accessions, the environments per accession, the replication and the total number of plots are the
# same in both and only the number of accessions sharing an environment differs.
#   dense          8 environments, 56 accessions per environment on average
#   sparse        60 environments, 7.5 accessions per environment
#   misspecified   the dense design with unequal replication, one to six environments per accession, and an effect
#                  independent across accessions added to the genetic value, which the relationship matrix does not
#                  represent. The unequal replication helps to separate that effect from the residual after the
#                  rotation: its rotated kernel is the inverse of Omega, and the separation is weak when Omega is close
#                  to proportional to the identity, as on the dense design with equal replication, where many lines
#                  share each environment.
# The first-stage model is the paper's: accession and environment fixed, replicate within environment and
# accession-by-environment random; the adjusted means are the estimated marginal means of the accessions with the
# Kenward-Roger adjusted variance-covariance matrix.

set.seed(20261012)
stopifnot(file.exists("DESCRIPTION"))
suppressMessages({ library(lme4); library(emmeans) })
stopifnot(requireNamespace("pbkrtest", quietly = TRUE))

## ---------- 1. the population: 150 inbred lines in 15 families from 24 founders, 2,000 markers on 10 chromosomes ----------
# Three subpopulations of eight founders each. Every family is a biparental cross of two founders of one subpopulation,
# five crosses per subpopulation, ten recombinant inbred lines per cross: a line's haplotype is a mosaic of the two
# parental haplotypes with crossovers along each chromosome, and the line is homozygous for it apart from residual
# heterozygosity. The relationship matrix thus has family blocks within the subpopulation structure.
n <- 150; n_chr <- 10; per_chr <- 200; m <- n_chr * per_chr
n_founders <- 24; n_families <- 15; lines_per_family <- 10
ids <- sprintf("L%03d", seq_len(n))
map <- data.frame(marker = sprintf("chr%02d_%04d", rep(seq_len(n_chr), each = per_chr), rep(seq_len(per_chr), n_chr)),
                  chromosome = rep(seq_len(n_chr), each = per_chr),
                  position = as.numeric(unlist(lapply(seq_len(n_chr), function(cc) sort(round(runif(per_chr, 0.05, 99.95), 2))))),
                  stringsAsFactors = FALSE)                  # positions in cM, chromosomes of 100 cM

# founder haplotypes: allele frequencies with a base value per marker shifted in each subpopulation, and linkage
# disequilibrium from a marker copying the allele of its left neighbour with probability 0.7
founder_pop <- rep(1:3, each = n_founders / 3)
p0 <- runif(m, 0.15, 0.85)
pf <- pmin(pmax(sweep(matrix(rnorm(3 * m, 0, 0.15), 3, m), 2, p0, "+"), 0.03), 0.97)
founders <- matrix(0L, n_founders, m)
for (i in seq_len(n_founders)) {
  h <- rbinom(m, 1, pf[founder_pop[i], ])
  copy <- which(rbinom(m, 1, 0.7) == 1 & map$chromosome == c(0, map$chromosome[-m]))
  for (j in copy) h[j] <- h[j - 1]
  founders[i, ] <- h
}

# the crosses: two different founders of one subpopulation, five crosses per subpopulation
crosses <- do.call(rbind, lapply(1:3, function(k) {
  f <- which(founder_pop == k)
  t(replicate(n_families / 3, sort(sample(f, 2))))
}))
family <- rep(seq_len(n_families), each = lines_per_family)

# a recombinant inbred line: on each chromosome a mosaic of the two parental haplotypes with a Poisson number of
# crossovers (mean two, the inbreeding generations included) at uniform positions
ril <- function(h1, h2) {
  h <- integer(m)
  for (cc in seq_len(n_chr)) {
    j <- which(map$chromosome == cc); pos <- map$position[j]
    xo <- sort(runif(rpois(1, 2), 0, 100))
    seg <- findInterval(pos, xo)                          # which segment each marker lies in
    from1 <- (seg + rbinom(1, 1, 0.5)) %% 2 == 0         # alternate parents from a random start
    h[j] <- ifelse(from1, h1[j], h2[j])
  }
  h
}
hap <- t(sapply(seq_len(n), function(i) ril(founders[crosses[family[i], 1], ], founders[crosses[family[i], 2], ])))
markers <- 2L * hap
het <- matrix(rbinom(n * m, 1, 0.02) == 1, n, m)          # residual heterozygosity, two per cent of the calls
markers[het] <- 1L
dimnames(markers) <- list(ids, map$marker)
subpopulation <- factor(paste0("S", founder_pop[crosses[family, 1]]))

# polymorphic markers with a minor allele frequency of at least 0.05
freq <- colMeans(markers) / 2
keep <- pmin(freq, 1 - freq) >= 0.05
markers <- markers[, keep]; map <- map[keep, ]; rownames(map) <- NULL
m <- ncol(markers)

# the relationship matrix of VanRaden (2008), scaled to a mean diagonal of one, and the first three principal components
p <- colMeans(markers) / 2
C <- sweep(markers, 2, 2 * p)
relationship <- tcrossprod(C) / (2 * sum(p * (1 - p)))
relationship <- relationship / mean(diag(relationship))
dimnames(relationship) <- list(ids, ids)
covariates <- prcomp(C)$x[, 1:3]
dimnames(covariates) <- list(ids, c("PC1", "PC2", "PC3"))

panel_genotypes <- list(markers = markers, map = map, relationship = relationship, covariates = covariates,
                        family = setNames(sprintf("F%02d", family), ids), subpopulation = setNames(as.character(subpopulation), ids))

## ---------- 2. the generating model ----------
# three causal markers on different chromosomes, with a minor allele frequency between 0.2 and 0.4
cand <- which(pmin(p, 1 - p) >= 0.2 & pmin(p, 1 - p) <= 0.4)
causal <- vapply(c(2, 5, 8), function(cc) sample(cand[map$chromosome[cand] == cc], 1), integer(1))
truth <- list(mu = 50, causal = colnames(markers)[causal], effect = c(2.0, 1.2, 0.6),         # per dose of the counted allele
              genetic_variance = 10, independent_variance = 0, interaction_variance = 3,
              replicate_variance = 0.5, error_variance = 6, environment_sd = 4,
              pc_coefficient = c(0.4, -0.3, 0.2))
rootA <- with(eigen(relationship, symmetric = TRUE), vectors %*% diag(sqrt(pmax(values, 0))))

make_design <- function(n_env, env_per_acc = 3, n_rep = 3) {
  env_per_acc <- rep_len(env_per_acc, n)
  starts <- 1 + floor((seq_len(n) - 1) * (n_env - env_per_acc + 1) / n)
  rows <- do.call(rbind, lapply(seq_len(n), function(i)
    expand.grid(accession = ids[i], environment = starts[i]:(starts[i] + env_per_acc[i] - 1), replicate = seq_len(n_rep),
                stringsAsFactors = FALSE)))
  data.frame(accession = factor(rows$accession, levels = ids),
             environment = factor(sprintf("E%02d", rows$environment)),
             replicate = factor(rows$replicate))
}

simulate_plots <- function(design, tr) {
  n_env <- nlevels(design$environment)
  g <- as.numeric(rootA %*% rnorm(n)) * sqrt(tr$genetic_variance)
  ind <- if (tr$independent_variance > 0) rnorm(n, 0, sqrt(tr$independent_variance)) else rep(0, n)
  qtl <- as.numeric(markers[, tr$causal] %*% tr$effect)
  str <- as.numeric(covariates %*% tr$pc_coefficient)
  value <- g + ind + qtl + str                                   # the genetic value of each accession, structure included
  env_eff <- rnorm(n_env, 0, tr$environment_sd)
  er <- interaction(design$environment, design$replicate, drop = TRUE)
  ae <- interaction(design$accession, design$environment, drop = TRUE)
  rep_eff <- rnorm(nlevels(er), 0, sqrt(tr$replicate_variance))
  ge_eff <- rnorm(nlevels(ae), 0, sqrt(tr$interaction_variance))
  y <- tr$mu + value[as.integer(design$accession)] + env_eff[as.integer(design$environment)] +
    rep_eff[as.integer(er)] + ge_eff[as.integer(ae)] + rnorm(nrow(design), 0, sqrt(tr$error_variance))
  list(plots = data.frame(design, y = round(y, 4)), value = setNames(value, ids), environment_effect = env_eff)
}

## ---------- 3. the first stage ----------
first_stage <- function(plots) {
  m1 <- lmer(y ~ accession + environment + (1 | environment:replicate) + (1 | accession:environment), data = plots, REML = TRUE)
  stopifnot(length(m1@optinfo$conv$lme4$messages) == 0 || all(grepl("boundary", m1@optinfo$conv$lme4$messages)))
  emm <- emmeans(m1, ~ accession, mode = "kenward-roger", rg.limit = 2e4)
  s <- summary(emm)
  means <- setNames(as.numeric(s$emmean), as.character(s$accession))
  Omega <- as.matrix(vcov(emm)); dimnames(Omega) <- list(names(means), names(means))
  L <- emm@linfct
  naive <- as.matrix(L %*% as.matrix(vcov(m1)) %*% t(L))
  kr <- mean(diag(Omega) / diag(naive))
  stopifnot(kr >= 1 - 1e-8)                            # the adjustment was applied; on these designs, balanced within accession and environment, it leaves Omega unchanged
  vc <- as.data.frame(VarCorr(m1))
  list(means = means, Omega = Omega,
       variance_components = setNames(vc$vcov, vc$grp), kenward_roger_factor = kr,
       singular = isSingular(m1))
}

build_panel <- function(n_env, tr, label, env_per_acc = 3) {
  design <- make_design(n_env, env_per_acc)
  sim <- simulate_plots(design, tr)
  fs <- first_stage(sim$plots)
  per_env <- table(design$environment[design$replicate == "1"])
  cat(sprintf("%-13s %2d environments, %d plots, %.1f accessions per environment (%d to %d); residual %.3f, G x E %.3f, replicate %.3f; KR factor %.6f%s\n",
              label, n_env, nrow(sim$plots), mean(per_env), min(per_env), max(per_env),
              fs$variance_components[["Residual"]], fs$variance_components[["accession:environment"]],
              fs$variance_components[["environment:replicate"]], fs$kenward_roger_factor, if (fs$singular) " [singular]" else ""))
  list(plots = sim$plots, means = fs$means, Omega = fs$Omega,
       design = list(environments = n_env, environments_per_accession = as.integer(table(design$accession) / 3), replicates = 3L, plots = nrow(sim$plots),
                     accessions_per_environment = as.numeric(mean(per_env))),
       first_stage = list(model = "y ~ accession + environment + (1 | environment:replicate) + (1 | accession:environment)",
                          variance_components = fs$variance_components, kenward_roger_factor = fs$kenward_roger_factor),
       truth = c(tr, list(value = sim$value, environment_effect = sim$environment_effect)))
}

panel_dense <- build_panel(8, truth, "dense")
panel_sparse <- build_panel(60, truth, "sparse")
truth_mis <- truth; truth_mis$independent_variance <- 6
env_per_acc <- sample(1:6, n, replace = TRUE, prob = c(0.1, 0.2, 0.3, 0.2, 0.12, 0.08))     # unequal replication
panel_misspecified <- build_panel(8, truth_mis, "misspecified", env_per_acc)

## ---------- 4. save ----------
dir.create("data", showWarnings = FALSE)
save(panel_genotypes, file = "data/panel_genotypes.rda", compress = "xz", version = 2)
save(panel_dense, file = "data/panel_dense.rda", compress = "xz", version = 2)
save(panel_sparse, file = "data/panel_sparse.rda", compress = "xz", version = 2)
save(panel_misspecified, file = "data/panel_misspecified.rda", compress = "xz", version = 2)
cat(sprintf("\n%d accessions, %d markers kept of %d; causal markers %s with effects %s per dose\n", n, m, n_chr * per_chr,
            paste(truth$causal, collapse = ", "), paste(truth$effect, collapse = ", ")))
print(file.info(list.files("data", full.names = TRUE))[, "size", drop = FALSE])
