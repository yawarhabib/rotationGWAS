# The data shipped with the package: aligned by name, and showing what the help pages say they show.

G <- panel_genotypes$markers; A <- panel_genotypes$relationship; Q <- panel_genotypes$covariates

test_that("the population and the panels are aligned by name and well formed", {
  expect_equal(dim(G), c(150L, 1998L)); expect_true(is.integer(G)); expect_true(all(G %in% 0:2))
  expect_equal(rownames(G), rownames(A)); expect_equal(colnames(A), rownames(A)); expect_equal(rownames(Q), rownames(A))
  expect_equal(panel_genotypes$map$marker, colnames(G))
  expect_equal(abs(mean(diag(A)) - 1), 0, tolerance = 1e-12)
  expect_true(min(eigen(A, symmetric = TRUE, only.values = TRUE)$values) > -1e-8)
  expect_equal(as.integer(table(panel_genotypes$family)), rep(10L, 15))
  for (pn in list(panel_dense, panel_sparse, panel_misspecified)) {
    expect_equal(names(pn$means), rownames(A)); expect_equal(dimnames(pn$Omega), list(rownames(A), rownames(A)))
    expect_true(min(eigen(pn$Omega, symmetric = TRUE, only.values = TRUE)$values) > 0)
    expect_equal(levels(pn$plots$accession), rownames(A)); expect_equal(nrow(pn$plots), pn$design$plots)
    expect_equal(as.integer(table(pn$plots$accession) / 3), pn$design$environments_per_accession)
    expect_equal(nlevels(pn$plots$environment), pn$design$environments)
    expect_true(all(pn$truth$causal %in% colnames(G)))
  }
  expect_equal(nrow(panel_dense$plots), 1350L); expect_equal(nrow(panel_sparse$plots), 1350L)
  expect_equal(panel_misspecified$truth$independent_variance, 6); expect_equal(panel_dense$truth$independent_variance, 0)
})

test_that("the dense panel: the schemes agree, the residual variance is near one, the strongest causal marker is found", {
  r <- rotation(panel_dense$means, panel_dense$Omega, A, Q)
  e <- efficiency(r, G)
  expect_lt(max(e$unweighted, na.rm = TRUE), 1.01); expect_lt(max(e$inverse_variance, na.rm = TRUE), 1.01)
  st <- specification_test(r)
  expect_false(st$reject_one_sided); expect_false(st$reject_two_sided)
  s <- genome_scan(r, G, scheme = c("rotated", "unweighted"))
  rot <- s[s$scheme == "rotated", ]
  expect_s3_class(rot, "data.frame"); expect_false(inherits(rot, "genome_scan"))       # a subset is a plain data frame
  expect_true(rot$significant[rot$marker == panel_dense$truth$causal[1]])
  expect_equal(rot$marker[which.max(rot$neglog10p)], panel_dense$truth$causal[1])
})

test_that("the sparse panel: the same lines and plots, but imprecise correlated means and a visible loss under approximate weighting", {
  r <- rotation(panel_sparse$means, panel_sparse$Omega, A, Q)
  expect_gt(r$omega[["mean_abs_cor"]], 0.3); expect_gt(r$omega[["diag_ratio"]], 2)
  e <- efficiency(r, G)
  expect_gt(median(e$inverse_variance, na.rm = TRUE), 1.05); expect_gt(max(e$unweighted, na.rm = TRUE), 1.5)
  expect_gt(estimation_share(r)$share, 0.3)
  st <- specification_test(r)
  expect_false(st$reject_one_sided)
  s <- genome_scan(r, G)
  expect_true(s$significant[s$marker == panel_sparse$truth$causal[1]])
})

test_that("the misspecified panel: the residual variance is above one, and the independent effect is found", {
  r <- rotation(panel_misspecified$means, panel_misspecified$Omega, A, Q)
  st <- specification_test(r)
  expect_true(st$reject_one_sided); expect_gt(st$residual, 2)
  ie <- independent_effect(r)
  expect_gt(ie$fit$vi, 0); expect_lt(ie$p_mixture, 0.05); expect_lt(ie$delta_aic, 0)
  expect_lt(abs(ie$fit_one$vi - panel_misspecified$truth$independent_variance), 3)   # 6.4 against 6
  expect_true(ie$interval_95$lower <= 1 && ie$interval_95$upper >= 1)
})
