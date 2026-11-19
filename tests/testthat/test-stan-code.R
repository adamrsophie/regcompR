test_that("Stan code covers ordinal, multinomial, structured and multiple membership models", {
  set.seed(1); n <- 60; x <- rnorm(n)
  d <- data.frame(x, y = cut(x + rlogis(n), c(-Inf, -1, .3, Inf), ordered_result = TRUE),
                  yn = factor(sample(letters[1:3], n, TRUE)), g = factor(rep(1:10, each = 6)),
                  times = factor(rep(1:6, 10)), tnum = rep(1:6, 10), z = rnorm(n),
                  s1 = factor(sample(5, n, TRUE)), s2 = factor(sample(5, n, TRUE)), w = runif(n))
  expect_match(stan_code(y ~ x, family = "cumulative_logit", data = d)$code, "ordered\\[Kc_cuts\\] cuts")
  expect_match(stan_code(yn ~ x, family = "multinomial_logit", data = d)$code, "categorical_logit_lpmf")
  expect_match(stan_code(z ~ x + ar1(0 + times | g), family = "normal", data = d)$code, "rho_1\\^abs")
  expect_match(stan_code(z ~ x + ou(tnum | g), family = "normal", data = d)$code, "range_1")
  sc <- stan_code(z ~ x + (1 | mm(s1, s2, weights = cbind(w, 1 - w))), family = "normal", data = d)
  expect_match(sc$code, "W_1 \\*")
  expect_equal(dim(sc$data$W_1), c(n, 5L))
})

test_that("JAGS runner fits a Poisson GLMM", {
  skip_on_cran(); skip_if_not_installed("R2jags")
  skip_if(Sys.which("jags") == "" && !file.exists("/opt/homebrew/bin/jags"), "JAGS not installed")
  set.seed(2); g <- rep(1:15, each = 8); x <- rnorm(120)
  d <- data.frame(y = rpois(120, exp(.5 + .4 * x + rnorm(15, 0, .5)[g])), x, g = factor(g))
  fj <- suppressWarnings(regcomp(y ~ x + (1 | g), family = "poisson", data = d, engine = "bayes",
    control = regcomp_control(runner = "jags", chains = 2, iter_warmup = 500, iter_sampling = 500)))
  fm <- regcomp(y ~ x + (1 | g), family = "poisson", data = d)
  expect_equal(unname(coef(fj)), unname(coef(fm)), tolerance = 0.1)
})
