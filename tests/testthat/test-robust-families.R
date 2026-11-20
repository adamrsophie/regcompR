test_that("Huber and stochastic frontier laws are normalized and match references", {
  for (cs in list(list("huber", list(mu = 0.3, sigma = 1.2, k = 1.345)),
                  list("stochastic_frontier", list(mu = 1, sigma_v = 0.5, sigma_u = 0.8)))) {
    f <- get_family(cs[[1]]); d <- do.call(distribution, c(cs[1], cs[[2]])); y <- c(-2, 0.1, 0.9, 3)
    if (cs[[1]] == "stochastic_frontier") skip_if_not_installed("sn")
    expect_equal(log_prob(d, y), f$reference$lpdf(y, lapply(cs[[2]], rep_len, 4)), tolerance = 1e-7)
    expect_equal(integrate(function(x) prob(d, x), -Inf, Inf)$value, 1, tolerance = 1e-6)
    expect_equal(integrate(function(x) x * prob(d, x), -Inf, Inf)$value, mean(d), tolerance = 1e-6)
    expect_equal(integrate(function(x) (x - mean(d))^2 * prob(d, x), -Inf, Inf)$value, variance(d), tolerance = 1e-6)
    expect_equal(quantile(d, cdf(d, 0.7)), 0.7, tolerance = 1e-6)
  }
})

test_that("Huber regression with the rlm scale gives the Huber M estimate", {
  skip_if_not_installed("MASS")
  set.seed(2); x <- rnorm(200); y <- 1 + 2 * x + rt(200, 2)
  r <- MASS::rlm(y ~ x, k = 1.345)
  f <- regcomp(y ~ x, family = "huber", data = data.frame(x, y), fixed = list(sigma = r$s))
  expect_equal(unname(coef(f)), unname(coef(r)), tolerance = 1e-4)
})

test_that("stochastic frontier regression maximizes the skew normal likelihood", {
  skip_if_not_installed("sn")
  set.seed(3); n <- 300; x <- runif(n); y <- 1 + 0.5 * x + rnorm(n, 0, 0.2) - abs(rnorm(n, 0, 0.4))
  f <- regcomp(y ~ x, family = "stochastic_frontier", data = data.frame(x, y))
  nll <- function(p) -sum(sn::dsn(y, p[1] + p[2] * x, sqrt(exp(2 * p[3]) + exp(2 * p[4])), -exp(p[4] - p[3]), log = TRUE))
  o <- optim(c(1, .5, log(.2), log(.4)), nll, method = "BFGS", control = list(reltol = 1e-12))
  expect_equal(unname(coef(f)), o$par, tolerance = 1e-4)
  expect_equal(as.numeric(logLik(f)), -o$value, tolerance = 1e-7)
})

test_that("partial proportional odds equals clm with nominal effects", {
  skip_if_not_installed("ordinal")
  set.seed(1)
  d <- data.frame(x = rnorm(400), z = rnorm(400)); eta <- 0.8 * d$x; u <- runif(400)
  d$y <- factor(1 + (u > plogis(-1 - eta)) + (u > plogis(1 + 0.5 * d$z - eta)))
  f <- regcomp(y ~ x, family = partial_odds(~ z), data = d)
  m <- ordinal::clm(y ~ x, nominal = ~ z, data = d)
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(m)), tolerance = 1e-7)
  expect_equal(unname(coef(f)[c("cuts:1|2", "cuts:2|3", "d1:z", "d2:z", "mu:x")]), unname(coef(m)), tolerance = 1e-4)
})

test_that("independent multinomial probit probabilities equal the normal orthant integrals", {
  skip_if_not_installed("mvtnorm")
  f <- get_family("multinomial_probit"); e <- c(0.3, -0.5, 1)
  expect_equal(f$lpdf(1:4, list(eta = e)), f$reference$lpdf(1:4, list(eta = e)), tolerance = 1e-7)
  expect_equal(sum(exp(f$lpdf(1:4, list(eta = e)))), 1, tolerance = 1e-10)
  set.seed(4); n <- 300; x <- rnorm(n)
  y <- max.col(cbind(0, 0.5 + x, -0.3 - x) + matrix(rnorm(3 * n), n))
  fit <- regcomp(y ~ x, family = "multinomial_probit", data = data.frame(x, y))
  expect_length(coef(fit), 4)
})
