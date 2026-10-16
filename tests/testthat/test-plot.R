# The plot methods: they draw without error on a null device and return the data they drew.

ex <- rotationGWAS:::toy_trial(n = 60, n_markers = 120, seed = 2)
r <- rotation(ex$means, ex$Omega, ex$relationship, ex$covariates)
map <- data.frame(marker = colnames(ex$markers), chromosome = rep(1:3, each = 40), position = rep(seq(1, 40), 3), stringsAsFactors = FALSE)
grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)

test_that("the rotation plots return the correlations or the variances", {
  R <- plot(r)
  expect_equal(dim(R), c(r$n, r$n)); expect_equal(unname(diag(R)), rep(1, r$n)); expect_true(all(R >= 0 & R <= 1))
  v <- plot(r, type = "variance", main = "x")
  expect_equal(v, diag(r$Omega))
})

test_that("the efficiency plots return the factors with their coordinate, by frequency or along the genome", {
  e <- efficiency(r, ex$markers)
  d1 <- plot(e)
  expect_equal(nrow(d1), ncol(ex$markers)); expect_equal(d1$x, e$maf); expect_equal(d1$unweighted, e$unweighted)
  d2 <- plot(e, type = "genome", map = map)
  expect_equal(nrow(d2), ncol(ex$markers))
  expect_true(all(diff(d2$x[1:40]) > 0))                                   # positions increase within a chromosome
  expect_true(d2$x[41] > d2$x[40])                                         # and the next chromosome starts beyond the last
  expect_error(plot(e, type = "genome"), "map must be")
  expect_error(plot(e, type = "genome", map = map[-1, ]), "not in map")
})

test_that("the scan plots return the Manhattan coordinates per scheme and the QQ quantiles", {
  s <- genome_scan(r, ex$markers, scheme = c("rotated", "unweighted"))
  m1 <- plot(s, scheme = "rotated", map = map, highlight = ex$truth$causal)
  expect_equal(nrow(m1), ncol(ex$markers)); expect_equal(m1$neglog10p, s$neglog10p[s$scheme == "rotated"])
  expect_equal(m1$significant, s$significant[s$scheme == "rotated"])
  m2 <- plot(s)                                                            # both schemes, no map: the marker order
  expect_equal(names(m2), c("rotated", "unweighted")); expect_equal(m2$unweighted$x, seq_len(ncol(ex$markers)))
  expect_error(plot(s, scheme = "inverse_variance"), "not in the scan")
  expect_error(plot(s, highlight = "nothing"), "not in the scan")
  q <- plot(s, type = "qq")
  expect_equal(names(q), c("rotated", "unweighted"))
  expect_equal(q$rotated$observed, sort(s$neglog10p[s$scheme == "rotated"], decreasing = TRUE))
  expect_true(all(diff(q$rotated$expected) <= 0))
})

test_that("the specification test holds its profile, which the plot returns, with its minimum at the estimate", {
  st <- specification_test(r)
  p <- plot(st)
  expect_equal(p, st$profile); expect_equal(p$ve[1], 0); expect_true(all(p$increase >= 0))
  expect_true(max(p$ve) >= st$interval_95$upper)                            # the grid covers the 95% interval
  at <- which.min(abs(p$ve - st$residual))
  expect_lt(p$increase[at], 0.05)                                          # the profile is near its minimum at the estimate
  expect_lt(abs(p$increase[which.min(abs(p$ve - 1))] - st$increase), 0.05)   # and at one it matches the test's increase, to the grid
})
