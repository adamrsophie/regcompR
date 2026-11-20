# Regression tests for defects found in external review.

test_that("censoring intervals are intersected with the truncation region", {
  d <- data.frame(lo = 0, up = 3)
  f <- regcomp(cens(lo, up, trunc_upper = 1) ~ 1, family = "exponential", data = d)
  expect_equal(as.numeric(logLik(f)), 0, tolerance = 1e-8)
  z <- obs_loglik(get_family("exponential"), 2, list(lambda = 1), unclass(cens(2, Inf, trunc_upper = 1)))
  expect_equal(z, -Inf)
  d <- data.frame(lo = c(-2, -1, 0, 1), up = c(1, 2, 3, 4))
  m <- regcomp(cens(lo, up, trunc_lower = 0, trunc_upper = 2) ~ 1, data = d, family = "normal", fixed = list(sigma = 1))
  mu <- fitted_distribution(m)$pars$mu
  set.seed(9); u <- runif(4)
  fl <- pnorm(pmax(d$lo, 0), mu); fu <- pnorm(pmin(d$up, 2), mu)
  ref <- (fl + u * (fu - fl) - pnorm(0, mu)) / (pnorm(2, mu) - pnorm(0, mu))
  expect_equal(unname(residuals(m, type = "pit", seed = 9)), unname(ref), tolerance = 1e-8)
})

test_that("the Laplace objective rejects non maxima and is stable across warm starts", {
  f <- custom_family("bimodal_location", pars = list(par_spec("mu", role = "location")),
    lpdf = function(y, p) log((dnorm(y, p$mu - 4) + dnorm(y, p$mu + 4)) / 2))
  s <- build_spec(y ~ (1 | g), list(), f, data.frame(y = 0, g = "a"))
  st <- new.env(); st$u <- 0
  expect_equal(laplace_objective(s, 0, 0, st, regcomp_control()), -Inf)
  d <- data.frame(y = c(0, 0, 5, 12, 20, 35), g = gl(3, 2))
  s <- build_spec(y ~ (1 | g), list(), get_family("poisson"), d)
  st <- new.env(); st$u <- rep(0, s$Q)
  v <- replicate(3, laplace_objective(s, log(4), 0, st, regcomp_control()))
  expect_equal(v, rep(v[1], 3), tolerance = 1e-10)
})

test_that("prediction replays fitted bases and contrasts", {
  d <- data.frame(x = 1:10, y = c(4, 6, 7, 6, 6, 7, 10, 11, 13, 17))
  nd <- data.frame(x = c(3, 4, 8))
  f <- regcomp(y ~ poly(x, 2), data = d, family = "normal")
  expect_equal(unname(predict(f, nd)), unname(predict(f)[nd$x]), tolerance = 1e-8)
  expect_equal(unname(predict(f, nd[1, , drop = FALSE])), unname(predict(f)[3]), tolerance = 1e-8)
  s <- regcomp(y ~ scale(x), data = d, family = "normal")
  expect_equal(unname(predict(s, nd[1, , drop = FALSE])), unname(predict(s)[3]), tolerance = 1e-8)
  d <- data.frame(y = c(1, 2, 2, 4, 4, 5), f = factor(rep(letters[1:3], 2)))
  m <- regcomp(y ~ f, family = "normal", data = d)
  old <- options(contrasts = c("contr.sum", "contr.poly")); on.exit(options(old))
  expect_equal(unname(predict(m, d)), unname(predict(m)), tolerance = 1e-8)
})

test_that("fixed values must satisfy the family constraint", {
  d <- data.frame(y = c(0, 0, 1, 4, 9, 12))
  expect_error(regcomp(y ~ 1, data = d, family = "generalized_poisson", fixed = list(lambda = 1.2)),
               "fixed parameter")
})

test_that("Stan data rejects observations outside the truncation region", {
  dat <- data.frame(y = c(.2, 1, 3))
  expect_error(stan_code(cens(y, trunc_upper = 2) ~ 1, family = "normal", data = dat), "truncation region")
  s <- stan_code(y ~ 1, family = truncate_dist("normal", lower = 0), data = data.frame(y = c(.2, 1, 3)))
  expect_match(s$code, "negative_infinity", fixed = TRUE)
})

test_that("marginal predictions summarize the mixture and new groups get their own effects", {
  skip_if_not_installed("lme4")
  data("sleepstudy", package = "lme4")
  m <- regcomp(Reaction ~ Days + (1 | Subject), data = sleepstudy, family = "normal")
  sig <- fitted_distribution(m)$pars$sigma[1]; tau <- unname(exp(m$theta))
  nd <- data.frame(Days = 0, Subject = factor(c("308", "new")))
  v <- predict(m, nd, type = "variance", re = "marginal", nmarginal = 1000, seed = 1)
  expect_equal(unname(v), unname(rep(sig^2 + tau^2, 2)), tolerance = 0.1)
  expect_equal(unname(predict(m, nd, type = "mean", re = "conditional"))[2], unname(coef(m)[1]), tolerance = 1e-8)
  set.seed(55)
  s <- simulate(m, nsim = 1000, newdata = nd, re = "new")
  expect_equal(unname(apply(s, 1, var)), unname(rep(sig^2 + tau^2, 2)), tolerance = 0.15)
})

test_that("ordinal masses and truncated laws keep their mass far in the tails", {
  d <- distribution("cumulative_logit", mu = -50, cuts = c(-1, 1))
  expect_equal(log_prob(d, 2), log_diff_exp(plogis(49, lower.tail = FALSE, log.p = TRUE),
                                            plogis(51, lower.tail = FALSE, log.p = TRUE)), tolerance = 1e-10)
  b <- cens(c(10, -41), c(11, -40))
  got <- obs_loglik(get_family("normal"), b[, "lower"], list(mu = c(0, 0), sigma = c(1, 1)), unclass(b))
  expect_equal(got[1], log_diff_exp(pnorm(10, lower.tail = FALSE, log.p = TRUE), pnorm(11, lower.tail = FALSE, log.p = TRUE)))
  expect_equal(got[2], log_diff_exp(pnorm(-40, log.p = TRUE), pnorm(-41, log.p = TRUE)))
  n <- distribution(truncate_dist("normal", lower = 10, upper = 11), mu = 0, sigma = 1)
  expect_equal(log_prob(n, 10.5), dnorm(10.5, log = TRUE) - got[1], tolerance = 1e-10)
  expect_true(quantile(n, 0.5) > 10 && quantile(n, 0.5) < 11)
  expect_equal(log_prob(distribution("truncated_normal", mu = 0, sigma = 1), 0.3), dnorm(0.3, log = TRUE))
  zt <- distribution(truncate_dist("poisson", lower = 0.5), lambda = 2)
  expect_equal(sum(prob(zt, 0:60)), 1, tolerance = 1e-10)
  expect_equal(prob(zt, 0), 0)
})
