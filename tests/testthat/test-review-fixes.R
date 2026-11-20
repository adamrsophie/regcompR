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

test_that("generated Stan programs pass stanc for operator and special families", {
  skip_if_not_installed("cmdstanr")
  stanc <- tryCatch(file.path(cmdstanr::cmdstan_path(), "bin", "stanc"), error = function(e) "")
  skip_if_not(file.exists(stanc))
  ok <- function(s) { f <- tempfile(fileext = ".stan"); writeLines(s$code, f)
    z <- suppressWarnings(system2(stanc, c("--name=probe", f, paste0("--o=", tempfile())), stdout = TRUE, stderr = TRUE))
    is.null(attr(z, "status")) }
  d <- data.frame(y = c(1, 2, 2, 3), n = rep(3L, 4), o = ordered(c(1, 2, 3, 2)))
  expect_true(ok(stan_code(o ~ 1, family = "cumulative_probit", data = d)))
  expect_true(ok(stan_code(y ~ 0, family = "poisson", data = d)))
  expect_true(ok(stan_code(y ~ 1, family = mixture("normal", "normal"), data = d)))
  expect_true(ok(stan_code(y ~ 1, family = "expectile", data = d)))
  expect_true(ok(stan_code(y ~ 1, family = "nig", data = d)))
  expect_true(ok(stan_code(y ~ 1, family = "sichel", data = d)))
  expect_true(ok(stan_code(y ~ 1, N = n, family = "zipf_finite", data = d)))
  expect_true(ok(stan_code(y ~ 1, m1 = n, m2 = n, n = n, family = "fisher_noncentral_hypergeometric", data = d)))
  expect_true(ok(stan_code(y ~ 1, size = n, family = hurdle("binomial"), data = d)))
  expect_true(ok(stan_code(y ~ 1, size = n, family = truncate_dist("binomial", 1), data = d)))
  expect_error(stan_code(y ~ 1, family = "poisson", data = data.frame(y = c(0.2, 1))), "integer")
  c2 <- stan_code(cens(y, y + 2) ~ 1, family = "poisson", data = data.frame(y = c(0.2, 1.2)))
  expect_equal(c2$data$y, c(1L, 2L)); expect_equal(c2$data$y_upper, c(2L, 3L))
})

test_that("group lasso reaches the weighted group optimum", {
  set.seed(2)
  d <- data.frame(x = rnorm(80), z = 12 * rnorm(80)); d$y <- .3 + 2 * d$x + .2 * d$z + rnorm(80)
  m <- regcomp(y ~ x + z, data = d, family = "normal", fixed = list(sigma = 1),
               penalty = penalty("group_lasso", lambda = 4, groups = c("one", "one")))
  X <- model.matrix(~ x + z, d); s <- apply(X[, -1], 2, sd)
  value <- function(b) sum((d$y - X %*% b)^2) / 2 + 4 * sqrt(2) * sqrt(sum((s * b[-1])^2))
  ref <- optim(coef(m), value, method = "BFGS", control = list(reltol = 1e-14, maxit = 5000))
  expect_lt(value(coef(m)), ref$value + 1e-5)
})

test_that("weights, missing rows, fixed data and dot formulas are handled", {
  d <- data.frame(y = c(0, 1, 3, -1), w = c(1, 1, 1, 0))
  expect_equal(unname(exp(coef(regcomp(y ~ 1, family = "poisson", data = d, weights = w)))), 4 / 3, tolerance = 1e-6)
  d$w[4] <- -1; d$y[4] <- 2
  expect_error(regcomp(y ~ 1, family = "poisson", data = d, weights = w), "nonnegative")
  d <- data.frame(y = c(1, 2, NA, 4, 12, 8), x = c(0, 1, 2, 3, NA, 5))
  expect_equal(unname(coef(regcomp(y ~ x, family = "poisson", data = d))),
               unname(coef(glm(y ~ x, data = d, family = poisson()))), tolerance = 1e-6)
  d <- data.frame(y = c(0, 1, 2, 3, 5, 6), x = 1:6)
  f <- regcomp(y ~ x, family = "binomial", fixed = list(size = 6), data = d)
  expect_equal(unname(predict(f, d[1:2, ])), unname(predict(f)[1:2]))
  d <- data.frame(y = c(0, 1, 1, 2, 4, 5), x = 1:6, trials = c(2, 4, 5, 5, 6, 7), w = c(1, 2, 1, 1, 3, 1))
  f <- regcomp(y ~ ., family = "binomial", size = trials, data = d, weights = w)
  expect_equal(names(coef(f)), c("prob:(Intercept)", "prob:x"))
  expect_false(anyNA(predict(f, transform(d, extra = x^2))))
})

test_that("distribution functions respect support, atoms, tails and structured parameters", {
  expect_equal(cdf(distribution("chi", k = 3), c(-Inf, -1)), c(0, 0))
  expect_equal(cdf(distribution("reciprocal", a = 1, b = 2), 3), 1)
  d <- distribution("poisson_binomial", probs = c(.3, .6, .8))
  expect_equal(cdf(d, 0:3), cumsum(prob(d, 0:3))); expect_equal(quantile(d, .5), 2)
  expect_equal(quantile(distribution("normal", mu = 0, sigma = 1), -750, log_p = TRUE, lower_tail = FALSE),
               qnorm(-750, log.p = TRUE, lower.tail = FALSE), tolerance = 1e-8)
  expect_equal(cdf(distribution("phase_type", alpha = c(1, 0), S = matrix(c(-2, 0, 2, -2), 2)), 1), pgamma(1, 2, 2))
  g <- distribution("generalized_beta1", scale = 2, a = -1, p = 2, q = 3)
  expect_equal(prob(g, 3), dbeta(2 / 3, 2, 3) * 2 / 9)
  expect_equal(log_prob(distribution("gamma_count", lambda = 60, alpha = 1), 0:2), dpois(0:2, 60, log = TRUE))
  expect_equal(quantile(distribution("gauss_kuzmin"), c(.5, .9)), c(2, 13))
  expect_equal(mean(distribution("beta_negative_binomial", mu = 5, sigma = 2, nu = .7)), 5)
  expect_equal(prob(distribution("continuous_categorical", eta = rep(1000, 8)), rep(1 / 8, 8)), gamma(8), tolerance = 1e-8)
  expect_equal(prob(distribution("dirichlet", alpha = c(2, 2)), c(.4, .4)), 0)
  expect_equal(prob(distribution("wishart", nu = 4, Sigma = diag(2)), diag(c(1, -1))), 0)
  z <- regcomp(y ~ 1, data = data.frame(y = c(0, 0, .2, .5, .9)), family = "rectified_gaussian")
  r <- residuals(z, type = "pit", seed = 9)
  expect_false(isTRUE(all.equal(r[1], r[2])))
})
