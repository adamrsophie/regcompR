test_that("registry loads many families", {
  expect_gt(length(families()), 200)
})

test_that("closed form families match their reference implementations", {
  cases <- list(
    gumbel = list(mu = 0.2, beta = 1.3),
    gev = list(mu = 0, sigma = 1, xi = 0.2),
    hyperbolic_secant = list(mu = 0.1, sigma = 0.9),
    von_mises = list(mu = 0.4, kappa = 2.3),
    wrapped_cauchy = list(mu = 0.4, rho = 0.6),
    expectile = list(mu = 0, sigma = 1, tau = 0.3)
  )
  y <- c(0.3, 0.9, 1.7)
  for (id in names(cases)) {
    f <- get_family(id)
    if (!is.null(f$reference$package) && !f$reference$package %in% c("definition") &&
        !requireNamespace(f$reference$package, quietly = TRUE)) next
    d <- do.call(distribution, c(list(id), cases[[id]]))
    ref <- f$reference$lpdf(y, lapply(cases[[id]], rep_len, length(y)))
    expect_equal(log_prob(d, y), ref, tolerance = 1e-8, info = id)
  }
})

test_that("densities integrate to one", {
  d <- distribution("gamma_gompertz", b = 0.5, s = 2, beta = 1.5)
  expect_equal(integrate(function(x) prob(d, x), 0, Inf)$value, 1, tolerance = 1e-6)
  w <- distribution(wrap("exponential"), lambda = 0.7)
  expect_equal(integrate(function(x) prob(w, x), -pi, pi)$value, 1, tolerance = 1e-6)
})

test_that("generator identities hold", {
  y <- c(0.4, 1.1, 2.5)
  a <- distribution(exponentiated("exponential"), lambda = 0.8, a = 1.7)
  b <- distribution("generalized_exponential", alpha = 1.7, lambda = 0.8)
  expect_equal(log_prob(a, y), log_prob(b, y), tolerance = 1e-10)
})

test_that("multivariate laws are normalized", {
  d <- distribution("negative_multinomial", x0 = 2.5, probs = c(0.2, 0.3))
  g <- as.matrix(expand.grid(0:150, 0:150))
  expect_equal(sum(exp(log_prob(d, g))), 1, tolerance = 1e-8)
  S <- matrix(c(2, 0.5, 0.5, 1), 2)
  m <- distribution("multivariate_normal", mu = c(0, 1), Sigma = S)
  y <- rbind(c(0.2, 0.4))
  expect_equal(log_prob(m, y), -log(2 * pi) - 0.5 * log(det(S)) - 0.5 * drop(t(y[1, ] - c(0, 1)) %*% solve(S) %*% (y[1, ] - c(0, 1))))
  l <- distribution("lkj", eta = 2.5)
  expect_equal(log_prob(l, list(matrix(c(1, 0.4, 0.4, 1), 2))), dbeta(0.7, 2.5, 2.5, log = TRUE) - log(2))
})
