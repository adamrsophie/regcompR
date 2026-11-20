# Models from other packages, checked against those packages.

test_that("Dirichlet regression matches DirichletReg", {
  skip_if_not_installed("DirichletReg")
  set.seed(4); n <- 200; x <- rnorm(n)
  Y <- cbind(rgamma(n, exp(0.5 + 0.3 * x)), rgamma(n, exp(0.2 - 0.2 * x)), rgamma(n, 1.5)); Y <- Y / rowSums(Y)
  dd <- data.frame(x = x); dd$Y <- Y
  f <- regcomp(Y ~ x, family = "dirichlet", data = dd)
  r <- DirichletReg::DirichReg(Yd ~ x, data.frame(Yd = I(DirichletReg::DR_data(Y)), x = x))
  expect_equal(as.numeric(logLik(f)), r$logLik, tolerance = 1e-6)
  expect_equal(unname(coef(f)), unname(unlist(r$coefficients)), tolerance = 1e-4)
})

test_that("quasi covariance equals glm quasibinomial and quasipoisson", {
  set.seed(1); d <- data.frame(x = rnorm(200), m = 20)
  d$s <- VGAM::rbetabinom(200, 20, plogis(d$x), 0.2)
  f <- regcomp(s ~ x, family = "binomial", size = m, data = d)
  expect_equal(unname(vcov(f, type = "quasi")), unname(vcov(glm(cbind(s, m - s) ~ x, quasibinomial, d))), tolerance = 1e-4)
  d$y <- rnbinom(200, mu = exp(1 + d$x), size = 2)
  fp <- regcomp(y ~ x, family = "poisson", data = d)
  expect_equal(summary(fp, vcov = "quasi")$table$std_error, unname(sqrt(diag(vcov(glm(y ~ x, quasipoisson, d))))), tolerance = 1e-4)
})

test_that("shifted lognormal regression is the maximum likelihood fit", {
  set.seed(2); x <- rnorm(400); y <- 0.3 + rlnorm(400, -0.5 + 0.3 * x, 0.4)
  f <- regcomp(y ~ x, family = "shifted_lognormal", data = data.frame(x, y))
  o <- optim(c(-0.5, 0.3, log(0.4), 0.2), function(b) -sum(ifelse(y > b[4], dlnorm(y - b[4], b[1] + b[2] * x, exp(b[3]), log = TRUE), -1e10)),
             control = list(maxit = 5000, reltol = 1e-12))
  expect_equal(as.numeric(logLik(f)), -o$value, tolerance = 1e-5)
})

test_that("Poisson regression on event time risk sets reproduces coxph", {
  skip_if_not_installed("survival")
  d <- na.omit(survival::lung[, c("time", "status", "age", "sex")]); d$status <- d$status - 1
  pp <- person_period(d, "time", "status", cuts = "events")
  f <- regcomp(event ~ 0 + factor(period) + age + sex, family = "poisson", data = pp)
  cx <- survival::coxph(survival::Surv(time, status) ~ age + sex, d, ties = "breslow")
  expect_equal(unname(coef(f)[c("lambda:age", "lambda:sex")]), unname(coef(cx)), tolerance = 1e-5)
  expect_equal(unname(sqrt(diag(vcov(f)))[c("lambda:age", "lambda:sex")]), unname(sqrt(diag(vcov(cx)))), tolerance = 1e-5)
  expect_error(person_period(d, "time", "status", cuts = c(5, 10)), "event time")
})
