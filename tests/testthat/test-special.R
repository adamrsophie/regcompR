test_that("exact Poisson matches poisson.test, glm and enumeration", {
  y <- c(3, 11); e <- c(1000, 1200)
  ep <- exact_poisson(y, c(0, 1), exposure = e); pt <- poisson.test(y[2:1], e[2:1])
  expect_equal(exp(ep$estimate), unname(pt$estimate), tolerance = 1e-5)
  expect_equal(ep$p_value, pt$p.value, tolerance = 1e-8)
  expect_equal(exp(ep$conf_int), as.numeric(pt$conf.int), tolerance = 1e-5)
  set.seed(2); st <- rep(1:4, each = 6); x <- rep(0:2, 8); ex <- runif(24, 1, 3)
  yy <- rpois(24, ex * exp(.3 * x + st / 4))
  g <- glm(yy ~ x + factor(st) + offset(log(ex)), poisson)
  expect_equal(exact_poisson(yy, x, strata = st, exposure = ex)$estimate, unname(coef(g)["x"]), tolerance = 1e-5)
})

test_that("person period expansion and multiple membership", {
  pp <- person_period(data.frame(t = c(2, 3, 1), d = c(1, 0, 1)), "t", "d")
  expect_equal(pp$period, c(1, 2, 1, 2, 3, 1)); expect_equal(pp$event, c(0, 1, 0, 0, 0, 1))
  set.seed(5); n <- 200; L <- 15; s1 <- sample(L, n, TRUE); w1 <- runif(n, .3, 1)
  d <- data.frame(y = rnorm(n) + rnorm(L)[s1], s1 = factor(s1), w1 = w1, w2 = 1 - w1)
  f1 <- regcomp(y ~ 1 + (1 | mm(s1, s1, weights = cbind(w1, w2))), family = "normal", data = d)
  f2 <- regcomp(y ~ 1 + (1 | s1), family = "normal", data = d)
  expect_equal(f1$logLik, f2$logLik, tolerance = 1e-6)
})

test_that("Heckman selection likelihood matches the bivariate normal integral", {
  skip_if_not_installed("mvtnorm")
  f <- get_family("heckman_selection")
  Y <- cbind(c(0, 1, 1), c(0, .4, -1.2))
  p <- list(mu = c(.2, .3, -.1), gamma = c(.5, -.2, 1), sigma = c(1.2, .8, 1.5), rho = c(.4, -.6, .3))
  expect_equal(f$lpdf(Y, p), f$reference$lpdf(Y, p), tolerance = 1e-8)
})
