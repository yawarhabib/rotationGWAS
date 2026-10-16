# specification_test(): the statistic from the two fits, the decisions against the intervals, and the two directions.

ex <- rotationGWAS:::toy_trial(n = 120, n_markers = 80, seed = 2)
r <- rotation(ex$means, ex$Omega, ex$relationship, ex$covariates)

test_that("the statistic is the increase in deviance at the residual fixed at one, zero when the estimate is at or below one", {
  st <- specification_test(r)
  expect_equal(st$increase, max(r$fit$one$deviance_ec - r$fit$free$deviance_ec, 0))
  expect_equal(st$statistic, if (r$fit$free$ve > 1) st$increase else 0)
  expect_equal(st$residual, r$fit$free$ve); expect_equal(st$genetic_at_one, r$fit$one$vg)
  expect_equal(st$crit_one_sided, 2.71); expect_equal(st$crit_two_sided, 3.84)
  expect_equal(st$p_two_sided, pchisq(st$increase, 1, lower.tail = FALSE))
  if (st$statistic > 0) expect_equal(st$p_one_sided, 0.5 * pchisq(st$statistic, 1, lower.tail = FALSE)) else expect_equal(st$p_one_sided, 1)
})

test_that("the decisions agree with the profile intervals (one-sided with the 90% interval above one, two-sided with the 95% interval excluding one)", {
  for (scale in c(0.25, 1, 4)) {
    Om <- ex$Omega * scale
    rs <- rotation(ex$means, Om, ex$relationship, ex$covariates)
    st <- specification_test(rs)
    # the increase at one is within 2.71 of the 90% interval's criterion only up to the interval's root tolerance
    expect_equal(st$reject_two_sided, !st$interval_95$inside_one)
    expect_equal(st$reject_one_sided, st$interval_90$lower > 1)
    expect_true(st$interval_90$lower >= st$interval_95$lower)
    expect_true(is.na(st$interval_95$upper) || st$interval_90$upper <= st$interval_95$upper)
  }
})

test_that("an Omega too small by a factor four is caught one-sided; one too large by a factor four is not, but is two-sided", {
  small <- rotation(ex$means, ex$Omega / 4, ex$relationship, ex$covariates)   # the rotated residual is then about four
  st_small <- specification_test(small)
  expect_true(st_small$residual > 1); expect_true(st_small$reject_one_sided); expect_true(st_small$reject_two_sided)
  large <- rotation(ex$means, ex$Omega * 4, ex$relationship, ex$covariates)   # the rotated residual is then about a quarter
  st_large <- specification_test(large)
  expect_true(st_large$residual < 1); expect_equal(st_large$statistic, 0); expect_false(st_large$reject_one_sided)
  expect_equal(st_large$p_one_sided, 1)
  expect_output(print(st_small), "REJECTED")
  expect_output(print(st_large), "statistic zero")
})
