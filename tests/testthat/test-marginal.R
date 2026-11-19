test_that("gee matches geepack and quasi likelihood matches glm", {
  set.seed(1); G <- 60; m <- 5
  d <- data.frame(g = rep(1:G, each = m), t = rep(1:m, G), x = rnorm(G * m))
  d$y <- rpois(G * m, exp(.2 + .4 * d$x + rnorm(G, 0, .5)[d$g]))
  q <- gee(y ~ x, data = d, id = seq_len(nrow(d)), family = quasipoisson())
  gq <- glm(y ~ x, quasipoisson, d)
  expect_equal(q$scale, summary(gq)$dispersion, tolerance = 1e-4)
  expect_equal(unname(sqrt(diag(q$vcov_naive))), unname(sqrt(diag(vcov(gq)))), tolerance = 1e-5)
  skip_if_not_installed("geepack")
  for (cs in c("exchangeable", "ar1", "unstructured")) {
    a <- gee(y ~ x, data = d, id = g, family = poisson(), corstr = cs, waves = t)
    b <- geepack::geeglm(y ~ x, data = d, id = g, family = poisson(), corstr = cs, waves = t)
    expect_equal(unname(coef(a)), unname(coef(b)), tolerance = 1e-4)
    expect_equal(unname(sqrt(diag(a$vcov))), unname(summary(b)$coefficients[, 2]), tolerance = 1e-4)
  }
})

test_that("cluster sandwich matches sandwich::vcovCL", {
  skip_if_not_installed("sandwich")
  set.seed(1); d <- data.frame(g = rep(1:30, each = 5), x = rnorm(150))
  d$y <- rpois(150, exp(0.3 + 0.5 * d$x + rnorm(30, 0, 0.5)[d$g]))
  fit <- regcomp(y ~ x, family = "poisson", data = d)
  gm <- glm(y ~ x, poisson, d)
  expect_equal(unname(vcov(fit, type = "sandwich", cluster = ~ g)),
               unname(sandwich::vcovCL(gm, cluster = ~ g, type = "HC0")), tolerance = 1e-5)
})
