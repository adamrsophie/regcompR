test_that("scoring rules match closed forms and scoringRules", {
  set.seed(1); d <- data.frame(x = rnorm(80)); d$y <- rnorm(80, 1 + d$x, 1.3); d$c <- rpois(80, exp(.5 + .3 * d$x))
  f <- regcomp(y ~ x, family = "normal", data = d); pm <- predict(f, type = "parameter")
  expect_equal(unname(scoring_rule(f, rule = "log")), unname(-dnorm(d$y, pm[, "mu"], pm[, "sigma"], log = TRUE)))
  expect_equal(scoring_rule(f, d[1:5, ], rule = "log"), scoring_rule(f, rule = "log")[1:5])
  skip_if_not_installed("scoringRules")
  expect_equal(unname(scoring_rule(f, rule = "crps")), unname(scoringRules::crps_norm(d$y, pm[, "mu"], pm[, "sigma"])), tolerance = 1e-7)
  fp <- regcomp(c ~ x, family = "poisson", data = d)
  expect_equal(unname(scoring_rule(fp, rule = "crps")), unname(scoringRules::crps_pois(d$c, predict(fp, type = "mean"))), tolerance = 1e-8)
})

test_that("cross validation, comparison and custom families run", {
  set.seed(1); d <- data.frame(x = rnorm(100)); d$c <- rpois(100, exp(.5 + .3 * d$x))
  fp <- regcomp(c ~ x, family = "poisson", data = d)
  cv <- cv_regcomp(fp, K = 4)
  expect_true(all(is.finite(cv$scores)))
  cm <- compare_models(a = fp, b = regcomp(c ~ 1, family = "poisson", data = d))
  expect_equal(nrow(cm), 2)
  fe <- custom_family("my_exp", list(par_spec("rate", "positive", role = "rate")),
    lpdf = function(y, p) dexp(y, p$rate, log = TRUE), support = function(p) list(lower = 0, upper = Inf), check = list(rate = 2))
  de <- data.frame(y = rexp(200, 2))
  expect_equal(unname(exp(coef(regcomp(y ~ 1, family = fe, data = de)))), 1 / mean(de$y), tolerance = 1e-6)
  expect_error(custom_family("bad", list(par_spec("rate", "positive")), lpdf = function(y, p) dexp(y, p$rate, log = TRUE) + 1,
    support = function(p) list(lower = 0, upper = Inf), check = list(rate = 1)), "integrates")
})
