test_that("HC0 to HC3 equal sandwich::vcovHC for GLMs", {
  skip_if_not_installed("sandwich")
  set.seed(1); d <- data.frame(x = rnorm(100), z = rnorm(100)); d$y <- rpois(100, exp(0.5 + 0.3 * d$x))
  d$b <- rbinom(100, 1, plogis(d$x))
  f <- regcomp(y ~ x + z, family = "poisson", data = d); g <- glm(y ~ x + z, family = poisson, data = d)
  fb <- regcomp(b ~ x + z, family = "bernoulli", data = d); gb <- glm(b ~ x + z, family = binomial, data = d)
  for (ty in c("HC0", "HC1", "HC2", "HC3")) {
    expect_equal(unname(vcov(f, type = ty)), unname(sandwich::vcovHC(g, type = ty)), tolerance = 1e-6, info = ty)
    expect_equal(unname(vcov(fb, type = ty)), unname(sandwich::vcovHC(gb, type = ty)), tolerance = 1e-6, info = ty)
  }
})

test_that("the Vuong statistic uses the pointwise log likelihood ratios", {
  set.seed(2); d <- data.frame(x = rnorm(300)); d$y <- rnbinom(300, mu = exp(1 + 0.3 * d$x), size = 2)
  a <- regcomp(y ~ x, family = "negative_binomial_nb2", data = d); b <- regcomp(y ~ x, family = "poisson", data = d)
  m <- log_lik(a) - log_lik(b)
  v <- vuong_test(a, b)
  expect_equal(v["Raw", "statistic"], sqrt(300) * mean(m) / sd(m))
  expect_equal(v["AIC-corrected", "statistic"], (sum(m) - 1) / (sqrt(300) * sd(m)))
  expect_equal(v$preferred[1], "model 1")
})

test_that("bootstrap, power simulation and best subset run and are reproducible", {
  set.seed(3); d <- data.frame(x1 = rnorm(80), x2 = rnorm(80), x3 = rnorm(80)); d$y <- rpois(80, exp(0.5 + 0.4 * d$x1))
  f <- regcomp(y ~ x1, family = "poisson", data = d)
  b1 <- boot_regcomp(f, R = 30); b2 <- boot_regcomp(f, R = 30)
  expect_equal(b1$t, b2$t)
  expect_equal(b1$R, 30)
  expect_lt(abs(log(b1$table$std_error[2] / sqrt(vcov(f)[2, 2]))), 0.5)
  expect_equal(nrow(boot_regcomp(f, R = 10, type = "case")$t), 10)
  p <- power_sim(f, "lambda:x1", nsim = 20)
  expect_true(p$power >= 0 && p$power <= 1 && p$nsim == 20)
  expect_error(power_sim(f, "lambda:nope"), "Unknown coefficient")
  bs <- best_subset(regcomp(y ~ x1 + x2 + x3, family = "poisson", data = d), criterion = "BIC")
  expect_equal(nrow(bs), 5)
  ff <- regcomp(y ~ x1 + x2 + x3, family = "poisson", data = d)
  expect_equal(regcompR:::info_criterion(ff, "AICc"), AIC(ff) + 2 * 4 * 5 / (80 - 4 - 1))
  expect_equal(regcompR:::info_criterion(ff, "HQ"), -2 * as.numeric(logLik(ff)) + 2 * 4 * log(log(80)))
})
