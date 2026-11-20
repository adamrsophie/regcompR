test_that("far tail CDFs of map_airy and normal_laplace stay finite and proper", {
  for (id in c("map_airy", "normal_laplace")) {
    f <- get_family(id); p <- lapply(f$pars, function(q) 1)
    v <- family_cdf(f, c(-1e4, -50, 0, 50, 1e4), p)
    expect_true(all(is.finite(v)) && all(diff(v) >= 0) && v[1] < 1e-6 && v[5] > 1 - 1e-6)
  }
  f <- get_family("normal_laplace"); p <- list(mu = 0, sigma = 1, alpha = 1, beta = 1)
  dens <- function(y) exp(f$lpdf(y, lapply(p, rep, length(y))))
  expect_equal(family_cdf(f, 0.7, p), integrate(dens, -Inf, 0.7)$value, tolerance = 1e-7)
  expect_true(is.finite(scoring_rule(regcomp(y ~ 1, family = "map_airy", data = data.frame(y = c(-1, 0.3, 1.2, 2))),
                                     rule = "crps")[1]))
  expect_warning(m <- mean(distribution("mittag_leffler", alpha = 0.7, scale = 1)), "does not exist")
  expect_true(is.na(m))
})

test_that("multivariate hypergeometric takes the number of draws", {
  d <- distribution("multivariate_hypergeometric", K = c(3, 5, 2), size = 4)
  set.seed(1); r <- rng(d, 50)
  expect_true(all(rowSums(r) == 4))
  expect_equal(log_prob(d, rbind(c(1, 1, 1))), -Inf)
  expect_equal(mean(d), c(1.2, 2, 0.8))
  expect_error(distribution("multivariate_hypergeometric", K = c(3, 5, 2), size = 11), "exceeds")
})

test_that("Hessians near a joint constraint boundary match the closed form", {
  f <- get_family("generalized_log_series")
  # d2/dtheta2 of log P(y): -y/theta^2 - (beta - 1) y/(1 - theta)^2 + (log(1 - theta) + 1)/g^2,
  # g = (1 - theta) log(1 - theta); d2/dtheta dbeta = -y/(1 - theta).
  for (th in c(0.3, 0.4999999)) {
    g <- (1 - th) * log(1 - th)
    H <- hessian(distribution(f, theta = th, beta = 2), 2)[1, , ]
    expect_equal(H["theta", "theta"], -2 / th^2 - 2 / (1 - th)^2 + (log(1 - th) + 1) / g^2, tolerance = 1e-3)
    expect_equal(H["theta", "beta"], -2 / (1 - th), tolerance = 1e-3)
  }
})

test_that("marginal scores mix over random effects and grouped CV clusters its standard error", {
  set.seed(1); G <- 20
  d <- data.frame(g = factor(rep(1:G, each = 6)), x = rnorm(G * 6))
  d$y <- rpois(G * 6, exp(0.5 + 0.4 * d$x + rnorm(G, 0, 0.8)[d$g]))
  f <- regcomp(y ~ x + (1 | g), family = "poisson", data = d)
  m <- scoring_rule(f, rule = "log", re = "marginal", seed = 1)
  expect_identical(m, scoring_rule(f, rule = "log", re = "marginal", seed = 1))
  expect_lt(mean(m), mean(scoring_rule(f, rule = "log", re = "zero")))
  f0 <- regcomp(y ~ x, family = "poisson", data = d)
  expect_equal(scoring_rule(f0, rule = "dss", re = "marginal"), scoring_rule(f0, rule = "dss"))
  cv <- cv_regcomp(f, K = 4, group = "g", rule = "log")
  expect_identical(cv$re, "marginal")
  tot <- tapply(cv$scores - cv$estimate, d$g, sum)
  expect_equal(cv$se, sqrt(G / (G - 1) * sum(tot^2)) / nrow(d))
  expect_silent(cv_regcomp(f, folds = rep(c("a", "b", "c"), length.out = nrow(d)), rule = "log"))
})

test_that("capability labels separate closed forms from family specific algorithms", {
  expect_identical(family_capabilities("normal")$cdf, "analytic")
  expect_identical(family_capabilities("conway_maxwell_poisson")$log_prob, "algorithm")
  expect_identical(family_capabilities("poisson_lognormal")$log_prob, "algorithm")
  expect_identical(family_capabilities("generalized_poisson")$cdf, "numeric")
})
