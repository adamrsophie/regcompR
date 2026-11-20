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
