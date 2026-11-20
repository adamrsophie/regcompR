# Model specification: factor by smooths, grouped FLAC, penalized
# identification and penalty metadata.

test_that("factor by smooths keep every level and match gamm", {
  skip_if_not_installed("mgcv"); skip_if_not_installed("nlme")
  set.seed(98)
  for (L in 2:3) {
    lv <- letters[1:L]
    d <- data.frame(x = rep(seq(0, 1, length.out = 50), L), g = factor(rep(lv, each = 50)))
    d$y <- (as.integer(d$g) - 1.5) * sin(2 * pi * d$x) + rnorm(nrow(d), 0, .08)
    f <- regcomp(y ~ g + s(x, by = g, k = 6), family = "normal", data = d)
    ref <- mgcv::gamm(y ~ g + s(x, by = g, k = 6), data = d, method = "ML")
    expect_equal(as.numeric(logLik(f)), as.numeric(logLik(ref$lme)), tolerance = 1e-5)
    nd <- data.frame(x = c(.1, .3, .8), g = factor(c(lv[L], "a", lv[2]), levels = lv))
    expect_equal(unname(predict(f, nd)), unname(as.numeric(predict(ref$gam, nd))), tolerance = 1e-4)
  }
  d$g2 <- factor(d$g, levels = rev(levels(d$g)))
  f2 <- regcomp(y ~ g2 + s(x, by = g2, k = 6), family = "normal", data = d)
  expect_equal(as.numeric(logLik(f2)), as.numeric(logLik(f)), tolerance = 1e-6)
})

test_that("FLAC gives the same fit for grouped and expanded binomial data", {
  chk <- function(d) {
    ex <- do.call(rbind, lapply(seq_len(nrow(d)), function(i) data.frame(x = d$x[i], w = d$w[i],
      y = c(rep(1, d$y[i]), rep(0, d$size[i] - d$y[i])))))
    for (type in c("firth", "flic", "flac")) {
      a <- regcomp(y ~ x, family = "binomial", data = d, size = size, weights = w, penalty = penalty(type))
      b <- regcomp(y ~ x, family = "bernoulli", data = ex, weights = w, penalty = penalty(type))
      expect_equal(unname(coef(a)), unname(coef(b)), tolerance = 1e-5)
    }
  }
  chk(data.frame(x = c(-2, -1, 0, 1, 2), y = c(0, 1, 4, 8, 9), size = 10, w = 1))
  chk(data.frame(x = c(-2, -1, 0, 1, 2), y = c(0, 1, 2, 6, 3), size = c(3, 7, 5, 10, 4), w = c(1, 2, 1, 3, 1)))
})

test_that("shrinkage penalties identify p > n designs; unpenalized fits still refuse them", {
  set.seed(45)
  hi <- data.frame(y = rnorm(5), matrix(rnorm(40), 5, 8))
  f <- regcomp(y ~ ., family = "normal", fixed = list(sigma = 1), data = hi, penalty = penalty("ridge", lambda = 2, standardize = FALSE))
  X <- model.matrix(~ ., hi[-1])
  expect_equal(unname(coef(f)), unname(drop(solve(crossprod(X) + diag(c(0, rep(2, 8))), crossprod(X, hi$y)))), tolerance = 1e-6)
  fl <- regcomp(y ~ ., family = "normal", fixed = list(sigma = 1), data = hi, penalty = penalty("lasso", lambda = 2))
  expect_lt(fl$diagnostics$kkt_residual, 1e-5)
  expect_error(regcomp(y ~ ., family = "normal", fixed = list(sigma = 1), data = hi), "rank deficient")
  expect_error(regcomp(y ~ ., family = "normal", fixed = list(sigma = 1), data = hi, penalty = penalty("firth")), "rank deficient")
})

test_that("stored penalties, objectives and proximal diagnostics agree", {
  set.seed(121)
  dp <- data.frame(x = rnorm(40), z = rnorm(40)); dp$y <- .3 + 2 * dp$x - dp$z + rnorm(40)
  for (pen in list(penalty("group_lasso", lambda = 3, groups = c("one", "one")), penalty("lasso", lambda = 3),
                   penalty("elastic_net", lambda = 3, alpha = .5))) {
    a <- regcomp(y ~ x + z, family = "normal", fixed = list(sigma = 1), data = dp, penalty = pen)
    b <- unname(coef(a))
    expect_equal(a$optimizer$objective, -as.numeric(logLik(a)) + a$penalty$value(b), tolerance = 1e-8)
    expect_lt(a$diagnostics$kkt_residual, 1e-5)
  }
  bs <- b <- NULL
  a <- regcomp(y ~ x + z, family = "normal", fixed = list(sigma = 1), data = dp, penalty = penalty("group_lasso", lambda = 3, groups = c("one", "one")))
  bs <- unname(coef(a))[a$penalty$idx] * a$penalty$scale
  expect_equal(a$penalty$value(unname(coef(a))), 3 * sqrt(2) * sqrt(sum(bs^2)))
  expect_warning(a1 <- regcomp(y ~ x + z, family = "normal", fixed = list(sigma = 1), data = dp,
                               penalty = penalty("lasso", lambda = 3), control = regcomp_control(prox_maxit = 3)), "converge")
  expect_equal(a1$optimizer$iterations, 3)
})
