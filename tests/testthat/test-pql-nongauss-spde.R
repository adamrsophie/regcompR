test_that("PQL equals MASS::glmmPQL", {
  skip_if_not_installed("MASS"); skip_if_not_installed("nlme")
  set.seed(1); d <- data.frame(g = factor(rep(1:20, each = 8)), x = rnorm(160))
  d$y <- rpois(160, exp(0.3 + 0.4 * d$x + rnorm(20, 0, 0.5)[d$g]))
  f <- pql_regcomp(y ~ x + (1 | g), family = poisson(), data = d)
  m <- MASS::glmmPQL(y ~ x, random = ~ 1 | g, family = poisson(), data = d, verbose = FALSE)
  expect_equal(unname(coef(f)[1:2]), unname(nlme::fixef(m)), tolerance = 1e-5)
  expect_equal(exp(unname(f$theta)), as.numeric(nlme::VarCorr(m)[1, 2]), tolerance = 1e-4)
  expect_error(logLik(f), class = "regcomp_unsupported")
})

test_that("non-normal random intercepts give the exact marginal likelihood", {
  set.seed(1); d <- data.frame(g = factor(rep(1:25, each = 4)), x = rnorm(100))
  d$y <- rpois(100, exp(0.3 + 0.4 * d$x) * rgamma(25, 2, 2)[d$g])
  f <- re_nongauss(y ~ x, ~ g, family = "poisson", data = d, re = "gamma")
  X <- cbind(1, d$x)
  # Gamma frailty: the negative multinomial likelihood in closed form.
  nll <- function(p) { mu <- exp(X %*% p[1:2]); k <- exp(-p[3]); Y <- tapply(d$y, d$g, sum); M <- tapply(mu, d$g, sum)
    -(sum(dpois(d$y, mu, log = TRUE)) + sum(lgamma(k + Y) - lgamma(k) + k * log(k) - (k + Y) * log(k + M) + M)) }
  o <- optim(c(0.3, 0.4, log(0.5)), nll, method = "BFGS", control = list(reltol = 1e-12))
  expect_equal(unname(coef(f)), o$par, tolerance = 1e-4)
  expect_equal(f$logLik, -o$value, tolerance = 1e-6)
  ft <- re_nongauss(y ~ x, ~ g, family = "poisson", data = d, re = "laplace", nq = 30); p <- coef(ft)
  lt <- sum(sapply(levels(d$g), function(l) { i <- d$g == l
    log(integrate(function(b) sapply(b, function(bb) prod(dpois(d$y[i], exp(p[1] + p[2] * d$x[i] + bb)))) *
      exp(-abs(b) / exp(p[3])) / (2 * exp(p[3])), -Inf, Inf, rel.tol = 1e-10)$value) }))
  # The Laplace density has a kink at zero, so quadrature is less exact.
  expect_equal(ft$logLik, lt, tolerance = 1e-5)
})

test_that("SPDE fields give the exact normal marginal likelihood", {
  skip_if_not_installed("fmesher"); skip_if_not_installed("mvtnorm")
  set.seed(1); n <- 60; d <- data.frame(sx = runif(n), sy = runif(n), x = rnorm(n)); d$y <- 1 + 0.5 * d$x + sin(3 * d$sx) + rnorm(n, 0, 0.3)
  mesh <- fmesher::fm_mesh_2d(loc = as.matrix(d[, 1:2]), max.edge = c(0.25, 0.5), cutoff = 0.05)
  f <- regcomp(y ~ x + spde(sx + sy | 1, mesh), family = "normal", data = d)
  A <- as.matrix(fmesher::fm_basis(mesh, as.matrix(d[, 1:2]))); X <- cbind(1, d$x); b <- f$spec$blocks[[1]]
  p <- c(coef(f), f$theta)
  ll <- mvtnorm::dmvnorm(d$y, X %*% p[1:2], exp(2 * p[3]) * diag(n) + A %*% regcompR:::block_cov(b, p[4:5]) %*% t(A), log = TRUE)
  expect_equal(as.numeric(logLik(f)), ll, tolerance = 1e-6)
  expect_error(regcomp(y ~ x + spde(sx + sy | 1, mesh), family = "normal", data = d, engine = "tmb"), class = "regcomp_unsupported")
})

test_that("normalized Stan code keeps constants and renormalizes bounded priors", {
  d <- data.frame(g = factor(rep(1:5, each = 4)), t = rep(1:4, 5), y = rpois(20, 2))
  code <- regcompR:::stan_normalize(stan_code(y ~ ar1(0 + factor(t) | g), family = "poisson", data = d)$code)
  expect_false(grepl(" ~ ", code, fixed = TRUE))
  expect_match(code, "student_t_lccdf(0 | 3, 0, 2.5)", fixed = TRUE)
  expect_match(code, "log((1) - (-1))", fixed = TRUE)
})
