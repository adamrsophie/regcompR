test_that("summary parameterizations equal direct reparameterized gamma fits", {
  set.seed(1); d <- data.frame(x = runif(150)); d$y <- rgamma(150, shape = 3, scale = exp(1 + d$x) / 2)
  X <- cbind(1, d$x)
  maps <- list(mode = function(m, s) m / (1 - s^2), geometric_mean = function(m, s) m / (s^2 * exp(digamma(1 / s^2))),
               harmonic_mean = function(m, s) m / (1 - s^2))
  for (ty in names(maps)) {
    f <- regcomp(y ~ x, family = summary_param("gamma", type = ty), data = d)
    nll <- function(p) { s <- exp(p[3]); -sum(dgamma(d$y, 1 / s^2, scale = maps[[ty]](exp(X %*% p[1:2]), s) * s^2, log = TRUE)) }
    o <- suppressWarnings(optim(c(1, 1, log(0.5)), nll, method = "BFGS", control = list(reltol = 1e-12)))
    expect_equal(unname(coef(f)), o$par, tolerance = 1e-5, info = ty)
    expect_equal(as.numeric(logLik(f)), -o$value, tolerance = 1e-8, info = ty)
  }
  expect_error(summary_param("normal", type = "geometric_mean"), class = "regcomp_unsupported")
})

test_that("uplift and Qini equal glm predictions and tools4uplift", {
  set.seed(1); n <- 600
  d <- data.frame(t = rbinom(n, 1, 0.5), x = rnorm(n))
  d$y <- rbinom(n, 1, plogis(-0.5 + 0.3 * d$x + d$t * (0.2 + 0.8 * d$x)))
  f <- regcomp(y ~ t * x, family = "bernoulli", data = d)
  g <- glm(y ~ t * x, family = binomial, data = d)
  u <- uplift(f, "t")
  expect_equal(u, unname(predict(g, transform(d, t = 1), type = "response") - predict(g, transform(d, t = 0), type = "response")), tolerance = 1e-5)
  q <- qini(u, d$t, d$y)
  expect_equal(nrow(q), 10)
  skip_if_not_installed("tools4uplift")
  d$u <- u; pu <- tools4uplift::PerformanceUplift(d, "t", "y", "u", nb.group = 10)
  expect_equal(as.matrix(as.data.frame(unclass(q))), as.matrix(as.data.frame(unclass(pu))), ignore_attr = TRUE)
  expect_equal(attr(q, "qini"), tools4uplift::QiniArea(pu))
})

test_that("Skovgaard's r* is close to the exact t test in a normal linear model", {
  set.seed(2); d <- data.frame(x = rnorm(15)); d$z <- 1 + 0.3 * d$x + rnorm(15)
  s <- skovgaard_test(regcomp(z ~ x, family = "normal", data = d), "mu:x", R = 2000)
  pt <- summary(lm(z ~ x, d))$coefficients[2, 4]
  expect_lt(abs(s["r*", "p_value"] - pt), abs(s["r", "p_value"] - pt))
  expect_lt(abs(s["r*", "p_value"] - pt), 0.01)
})
