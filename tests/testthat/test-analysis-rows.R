# One evaluated set of analysis rows shared by fitting, external vectors,
# cross validation and clustered covariances.

test_that("dot formulas and transformed covariates drop the same rows as glm", {
  dn <- data.frame(y = c(1, 3, 2, 5, 5, 8, 6, 9), x = 1:8); dn$x[3] <- NA
  ref <- unname(coef(glm(y ~ ., family = poisson(), data = dn)))
  expect_equal(unname(coef(regcomp(y ~ x, family = "poisson", data = dn))), ref, tolerance = 1e-6)
  expect_equal(unname(coef(regcomp(y ~ ., family = "poisson", data = dn))), ref, tolerance = 1e-6)
  dx <- data.frame(y = c(1, 3, 2, 5, 5, 8, 6, 9), x = c(1, 2, -1, 4, 5, 6, 7, 8))
  f <- regcomp(y ~ log(x), family = "poisson", data = dx)
  expect_equal(nobs(f), 7)
  expect_equal(unname(coef(f)), unname(suppressWarnings(coef(glm(y ~ log(x), family = poisson(), data = dx)))), tolerance = 1e-6)
  expect_equal(f$rows, c(1:2, 4:8))
})

test_that("external weights and data arguments are aligned before subsetting", {
  dw <- data.frame(y = c(1, 3, 2, 5, 5, 8, 6, 9), x = 1:8); wv <- 1:8
  expect_equal(unname(coef(regcomp(y ~ x, family = "poisson", data = dw, weights = wv, subset = 1:6))),
               unname(coef(glm(y ~ x, family = poisson(), data = dw, weights = wv, subset = 1:6))), tolerance = 1e-6)
  expect_error(regcomp(y ~ 1, family = "binomial", size = c(2, 3), data = data.frame(y = c(0, 1, 0, 1, 2, 1))), "one value or 6")
  d <- data.frame(y = c(.5, 1, 2, 3, 2, 1), x = 1:6)
  sc <- stan_code(y ~ x, family = "normal", data = d, weights = 2)
  expect_equal(sc$data$w, rep(2, 6))
  expect_equal(unname(coef(regcomp(y ~ x, family = "normal", data = d, weights = 2))),
               unname(coef(regcomp(y ~ x, family = "normal", data = d, weights = rep(2, 6)))), tolerance = 1e-8)
})

test_that("cross validation splits the stored analysis rows and external vectors", {
  dat <- data.frame(x = seq(-1, 1, length.out = 24), sz = rep(c(100, 500, 1000, 200), 6))
  dat$y <- c(1,2,0,1,3,2,1,2,3,2,2,1,2,1,0,3,1,2,3,2,1,0,2,1)
  sz_external <- dat$sz; folds <- rep(1:3, 8)
  a <- cv_regcomp(regcomp(y ~ x, size = sz_external, family = "binomial", data = dat), folds = folds)
  b <- cv_regcomp(regcomp(y ~ x, size = sz, family = "binomial", data = dat), folds = folds)
  expect_equal(a$scores, b$scores, tolerance = 1e-8)
  wext <- rep(c(1, 2), 12); dat$w <- wext
  a <- cv_regcomp(regcomp(y ~ x, size = sz, family = "binomial", data = dat, weights = wext), folds = folds)
  b <- cv_regcomp(regcomp(y ~ x, size = sz, family = "binomial", data = dat, weights = w), folds = folds)
  expect_equal(a$scores, b$scores, tolerance = 1e-8)
  dm <- data.frame(x = seq(-1, 1, length.out = 24), y = rep(c(1, 2, 0, 3, 1, 2), 4)); dm$y[5] <- NA
  fm <- regcomp(y ~ x, family = "poisson", data = dm)
  cv <- cv_regcomp(fm, K = 3)
  expect_true(is.finite(cv$estimate)); expect_length(cv$scores, nobs(fm))
  # Folds over the original rows are restricted to the analysis rows.
  expect_equal(cv_regcomp(fm, folds = rep(1:3, 8))$scores, cv_regcomp(fm, folds = rep(1:3, 8)[-5])$scores)
})

test_that("cluster formulas use the fitted analysis rows", {
  set.seed(2003)
  dat <- data.frame(x = seq(-1, 1, length.out = 24), g = rep(1:8, each = 3))
  dat$y <- rpois(24, exp(1 + .3 * dat$x))
  fit <- regcomp(y ~ x, family = "poisson", data = dat, subset = 1:18)
  expect_equal(vcov(fit, type = "sandwich", cluster = ~g), vcov(fit, type = "sandwich", cluster = dat$g[1:18]), tolerance = 1e-9)
  dat$y[2] <- NA
  fit2 <- regcomp(y ~ x, family = "poisson", data = dat)
  expect_equal(vcov(fit2, type = "sandwich", cluster = ~g), vcov(fit2, type = "sandwich", cluster = dat$g[-2]), tolerance = 1e-9)
})

test_that("Bayesian scores use the posterior predictive law; CRPS divergence is infinite", {
  base <- regcomp(y ~ 1, family = "normal", data = data.frame(y = -2:2), fixed = list(sigma = 1))
  b <- base; b$inference <- "bayesian"; b$u_draws <- NULL
  b$draws <- matrix(c(-3, 3), 2, dimnames = list(NULL, names(coef(base))))
  ea <- function(mu, sd) 2 * sd * dnorm(mu / sd) + mu * (2 * pnorm(mu / sd) - 1)
  crps <- mean(ea(c(-3, 3), 1)) - mean(outer(c(-3, 3), c(-3, 3), function(a, b) ea(a - b, sqrt(2)))) / 2
  expect_equal(unname(scoring_rule(b, data.frame(y = 0), rule = "crps")), crps, tolerance = 1e-6)
  expect_equal(unname(scoring_rule(b, data.frame(y = 0), rule = "dss")), log(10), tolerance = 1e-8)
  oc <- regcomp(y ~ 1, family = "cumulative_logit", data = data.frame(y = rep(1:2, each = 5)))
  ob <- oc; ob$inference <- "bayesian"; ob$u_draws <- NULL
  ob$draws <- matrix(c(-1, 3), 2, dimnames = list(NULL, names(coef(oc))))
  p1 <- mean(plogis(c(-1, 3)))
  expect_equal(unname(scoring_rule(ob, data.frame(y = 1), rule = "brier")), 2 * (1 - p1)^2, tolerance = 1e-8)
  expect_equal(unname(scoring_rule(ob, data.frame(y = 1), rule = "rps")), (1 - p1)^2, tolerance = 1e-8)
  fit <- regcomp(y ~ 1, xm = 1, family = "pareto", data = data.frame(y = exp(2:6)))
  expect_equal(unname(scoring_rule(fit, data.frame(y = 1), rule = "crps")), Inf)
})

test_that("binary Brier, comparison guards, new groups and RNG-safe plots", {
  d <- data.frame(y = c(0, 1, 1, 0, 1, 1), x = 1:6)
  fb <- regcomp(y ~ x, family = "bernoulli", data = d)
  expect_equal(unname(scoring_rule(fb, rule = "brier")), unname((d$y - predict(fb))^2), tolerance = 1e-10)
  db <- data.frame(y = c(1, 3, 2, 4), n = 5, x = 1:4)
  fbin <- regcomp(y ~ x, family = "binomial", size = n, data = db)
  expect_equal(unname(scoring_rule(fbin, rule = "brier")), unname((db$y / 5 - predict(fbin) / 5)^2), tolerance = 1e-10)
  dp <- data.frame(y = c(1, 3, 2, 5, 4), x = 1:5)
  expect_error(compare_models(regcomp(y ~ x, family = "poisson", data = dp),
                              regcomp(y ~ x, family = "poisson", data = dp, subset = 1:3)), "same observations")
  expect_error(compare_models(regcomp(y ~ x, family = "normal", data = dp),
                              regcomp(log(y) ~ x, family = "normal", data = dp)), "same observations")
  expect_true(supports("cumulative_logit", "rng"))
  set.seed(2001)
  dg <- data.frame(g = factor(rep(letters[1:8], each = 8)))
  dg$y <- rep(seq(-3, 3, length.out = 8), each = 8) + rnorm(64, sd = .4)
  fit <- regcomp(y ~ 1 + (1 | g), family = "normal", data = dg, fixed = list(sigma = .4))
  ref <- unname(.4^2 + exp(fit$theta)^2)
  got <- predict(fit, data.frame(dummy = 1), type = "variance", re = "marginal", nmarginal = 3000, seed = 44)
  expect_lt(abs(unname(got) - ref) / ref, .1)
  set.seed(5); before <- .Random.seed
  pdf(NULL); plot(distribution("multivariate_normal", mu = c(0, 0), Sigma = diag(2))); dev.off()
  expect_identical(.Random.seed, before)
})

test_that("no export masks a function of the base R packages", {
  ex <- getNamespaceExports("regcompR")
  for (p in c("base", "stats", "utils", "graphics", "grDevices", "methods"))
    expect_identical(intersect(ex, getNamespaceExports(p)), character(0), info = p)
})
