test_that("a nonlinear formula equals nls", {
  set.seed(1)
  d <- data.frame(x = runif(80, 0, 5))
  d$y <- 10 * exp(-0.6 * d$x) + rnorm(80, 0, 0.4)
  f <- regcomp(y ~ nl(a * exp(-b * x)), a ~ 1, b ~ 1, family = "normal", data = d,
               control = regcomp_control(start = list(a = 5, b = 1)))
  n <- nls(y ~ a * exp(-b * x), data = d, start = list(a = 5, b = 1))
  expect_equal(unname(coef(f)[1:2]), unname(coef(n)), tolerance = 1e-5)
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(n)), tolerance = 1e-6)
  nd <- data.frame(x = c(0, 1, 2))
  expect_equal(unname(predict(f, newdata = nd)), unname(predict(n, newdata = nd)), tolerance = 1e-5)
  # Standard errors from the observed information of the same likelihood.
  nll <- function(p) -sum(dnorm(d$y, p[1] * exp(-p[2] * d$x), exp(p[3]), log = TRUE))
  H <- stats::optimHess(coef(f), nll)
  expect_equal(unname(sqrt(diag(vcov(f)))), unname(sqrt(diag(solve(H)))), tolerance = 1e-4)
  expect_error(icc(f), class = "regcomp_unsupported")
})

test_that("nonlinear parameters take covariates and a non identity link", {
  set.seed(2)
  d <- data.frame(x = runif(200, 0, 4), z = rbinom(200, 1, 0.5))
  d$y <- rpois(200, (2 + 3 * d$z) * (1 - exp(-0.8 * d$x)))
  f <- regcomp(y ~ nl(log(a) + log1p(-exp(-k * x))), a ~ z, k ~ 1, family = "poisson", data = d,
               control = regcomp_control(start = list(a = 2, k = 1)))
  nll <- function(p) { m <- (p[1] + p[2] * d$z) * (1 - exp(-p[3] * d$x))
    if (any(m <= 0)) Inf else -sum(dpois(d$y, m, log = TRUE)) }
  o <- optim(c(2, 1, 1), nll, method = "BFGS", control = list(reltol = 1e-12))
  expect_equal(unname(coef(f)), o$par, tolerance = 1e-4)
  expect_equal(as.numeric(logLik(f)), -o$value, tolerance = 1e-6)
})

test_that("a nonlinear mixed model equals nlmer", {
  skip_if_not_installed("lme4")
  O <- as.data.frame(Orange); O$Tree <- factor(as.character(O$Tree))
  f <- regcomp(circumference ~ nl(Asym / (1 + exp((xmid - age) / scal))), Asym ~ 1 + (1 | Tree), xmid ~ 1, scal ~ 1,
               family = "normal", data = O, control = regcomp_control(start = list(Asym = 200, xmid = 725, scal = 350)))
  m <- lme4::nlmer(circumference ~ SSlogis(age, Asym, xmid, scal) ~ Asym | Tree, data = Orange,
                   start = c(Asym = 200, xmid = 725, scal = 350))
  expect_equal(unname(coef(f)[1:3]), unname(lme4::fixef(m)), tolerance = 1e-4)
  expect_equal(as.numeric(logLik(f)), as.numeric(logLik(m)), tolerance = 1e-5)
  vc <- as.data.frame(lme4::VarCorr(m))$sdcor
  expect_equal(exp(unname(c(coef(f)[4], f$theta))), vc[2:1], tolerance = 1e-3)
})

test_that("nonlinear formulas check their input", {
  d <- data.frame(x = 1:10, y = exp(-(1:10) / 3))
  ctl <- regcomp_control(start = list(a = 1, b = 1))
  expect_error(regcomp(y ~ nl(a * exp(-b * x)), a ~ 1, b ~ 1, family = "normal", data = d), "starting values")
  expect_error(regcomp(y ~ nl(sigma * exp(-b * x)), sigma ~ 1, b ~ 1, family = "normal", data = d, control = ctl), "rename sigma")
  expect_error(regcomp(y ~ nl(a * exp(-b * x)), a ~ 1, b ~ 1, family = "normal", data = d, engine = "tmb", control = ctl),
               class = "regcomp_unsupported")
  expect_error(regcomp(y ~ nl(x * exp(-b * x)), x ~ 1, b ~ 1, family = "normal", data = d, control = ctl), "data columns")
})
