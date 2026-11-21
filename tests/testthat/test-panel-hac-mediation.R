test_that("the within instrumental variable estimator equals plm", {
  skip_if_not_installed("plm")
  set.seed(1); d <- data.frame(id = rep(1:40, each = 5), t = rep(1:5, 40)); a <- rnorm(40)
  d$z1 <- rnorm(200); d$z2 <- rnorm(200); u <- rnorm(200); d$w <- rnorm(200)
  d$x <- d$z1 + 0.5 * d$z2 + a[d$id] + u + rnorm(200); d$y <- 1 + 0.5 * d$x + 0.3 * d$w + a[d$id] + u
  pd <- plm::pdata.frame(d, index = c("id", "t"))
  for (ef in c("individual", "twoways")) {
    f <- panel_within(y ~ x + w | z1 + z2 + w, data = d, id = "id", time = "t", effect = ef)
    p <- plm::plm(y ~ x + w | z1 + z2 + w, data = pd, model = "within", effect = ef)
    expect_equal(coef(f), coef(p), tolerance = 1e-10, info = ef)
    expect_equal(unname(vcov(f)), unname(vcov(p)), tolerance = 1e-10, info = ef)
  }
})

test_that("the HAC covariance equals sandwich::NeweyWest and boundary p values halve", {
  skip_if_not_installed("sandwich")
  set.seed(2); d <- data.frame(x = rnorm(150), w = rnorm(150)); d$y <- rbinom(150, 1, plogis(d$x))
  g <- glm(y ~ x + w, family = binomial, data = d, control = glm.control(epsilon = 1e-12)); f <- regcomp(y ~ x + w, family = "bernoulli", data = d)
  expect_equal(unname(vcov(f, type = "HAC", lag = 3)), unname(sandwich::NeweyWest(g, lag = 3, prewhite = FALSE, adjust = FALSE)), tolerance = 1e-4)
  gg <- data.frame(g = factor(rep(1:30, each = 6)), x = rnorm(180)); gg$y <- rpois(180, exp(0.2 + 0.3 * gg$x + rnorm(30, 0, .4)[gg$g]))
  f0 <- regcomp(y ~ x, family = "poisson", data = gg); f1 <- regcomp(y ~ x + (1 | g), family = "poisson", data = gg)
  expect_equal(anova(f0, f1, boundary = TRUE)[2, "Pr(>Chisq)"], anova(f0, f1)[2, "Pr(>Chisq)"] / 2)
})

test_that("shrinkage uses the heuristic factor and re-estimates the intercept", {
  set.seed(1); d <- data.frame(x1 = rnorm(60), x2 = rnorm(60), x3 = rnorm(60))
  d$y <- rbeta(60, plogis(0.3 * d$x1) * 10, (1 - plogis(0.3 * d$x1)) * 10)
  f <- regcomp(y ~ x1 + x2 + x3, family = "beta", data = d); f0 <- regcomp(y ~ 1, family = "beta", data = d)
  lr <- 2 * (f$logLik - f0$logLik)
  s <- shrink_regcomp(f)
  expect_equal(attr(s, "shrinkage"), (lr - 3) / lr, tolerance = 1e-6)
  expect_equal(unname(coef(s)[2:4]), unname(coef(f)[2:4]) * attr(s, "shrinkage"))
  b <- as.numeric(d$y > 0.5); fb <- regcomp(b ~ x1 + x2 + x3, family = "bernoulli", data = d)
  g <- glm(b ~ x1 + x2 + x3, binomial, d); off <- 0.8 * drop(model.matrix(g)[, -1] %*% coef(g)[-1])
  expect_equal(unname(coef(shrink_regcomp(fb, factor = 0.8))[1]), unname(coef(glm(b ~ offset(off), binomial, d))), tolerance = 1e-5)
})

test_that("mediation effects in linear models equal the product of coefficients", {
  set.seed(3); n <- 400
  d <- data.frame(t = rbinom(n, 1, 0.5), x = rnorm(n)); d$m <- rnorm(n, 0.5 * d$t + 0.3 * d$x); d$y <- rnorm(n, 0.4 * d$m + 0.2 * d$t)
  mm <- regcomp(m ~ t + x, family = "normal", data = d); my <- regcomp(y ~ t + m + x, family = "normal", data = d)
  r <- mediate_regcomp(mm, my, "t", "m", sims = 1000)
  ab <- coef(mm)[["mu:t"]] * coef(my)[["mu:m"]]
  expect_lt(abs(r["ACME (average)", "estimate"] - ab), 0.01)
  expect_lt(abs(r["ADE (average)", "estimate"] - coef(my)[["mu:t"]]), 0.01)
})
