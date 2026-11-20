# Measure, atom and lattice contract: point masses of mixed laws, left
# limits, non unit lattices, stable tail masses and constraint propagation.

test_that("mixed laws combine atoms and densities under one measure", {
  dm <- distribution(mixture(inflate("normal"), "normal"), mu.1 = 0, sigma.1 = 1, zi.1 = .2,
                     mu.2 = 0, sigma.2 = 1, w2 = 0)
  expect_equal(prob(dm, 0), 0.1)
  expect_equal(cdf(dm, 0) - cdf(dm, -1e-10), 0.1, tolerance = 1e-6)
  expect_equal(prob(dm, 1), 0.5 * 0.8 * dnorm(1) + 0.5 * dnorm(1))
  dn <- distribution(inflate(inflate("beta")), mu = .5, phi = 2, zi = .2, zi_ = .3)
  expect_equal(prob(dn, 0), .3 + .7 * .2)
  expect_no_error(inflate(inflate("normal")))
  dh <- distribution(hurdle(inflate("normal")), mu = 0, sigma = 1, zi = .2, hu = .3)
  tot <- prob(dh, 0) + integrate(function(x) prob(dh, x), -Inf, 0)$value + integrate(function(x) prob(dh, x), 0, Inf)$value
  expect_equal(tot, 1, tolerance = 1e-8)
  dc <- distribution(convolve_dist(inflate("normal"), "normal"), mu.1 = 0, sigma.1 = 1, zi.1 = .2, mu.2 = 0, sigma.2 = 1)
  expect_equal(prob(dc, 0), .2 * dnorm(0) + .8 * dnorm(0, sd = sqrt(2)), tolerance = 1e-8)
  expect_error(convolve_dist(inflate("normal"), inflate("normal")), class = "regcomp_unsupported")
  dd <- distribution(discretize(inflate("exponential")), lambda = 1, zi = .2)
  expect_equal(prob(dd, 0), .2 + .8 * pexp(1))
  expect_equal(sum(prob(dd, 0:100)), 1)
  set.seed(3)
  expect_equal(mean(rng(dd, 20000) == 0), .2 + .8 * pexp(1), tolerance = .02)
  # Intrinsic mixed families report their atoms.
  tw <- distribution("tweedie", mu = 1, phi = 1, power = .5)
  expect_equal(regcompR:::family_atom(tw$family, 0, tw$pars), exp(-1 / 0.5))
  expect_equal(regcompR:::family_atom(get_family("rectified_gaussian"), 0, list(mu = 0, sigma = 1)), 0.5)
})

test_that("derived families keep base constraints and lattices", {
  expect_error(distribution("generalized_poisson", theta = 2, lambda = 1.2), class = "regcomp_domain_error")
  expect_error(distribution(exponentiated("generalized_poisson"), theta = 2, lambda = 1.2, a = 1), class = "regcomp_domain_error")
  expect_error(distribution(truncate_dist("generalized_poisson", lower = 1), theta = 2, lambda = 1.2), class = "regcomp_domain_error")
  expect_error(distribution(mixture("generalized_poisson", "poisson"), theta.1 = 2, lambda.1 = 1.2, lambda.2 = 1, w2 = 0),
               class = "regcomp_domain_error")
  dr <- distribution(round_dist("normal", digits = 1), mu = 0, sigma = 1)
  drg <- distribution(exponentiated(round_dist("normal", digits = 1)), mu = 0, sigma = 1, a = 1)
  expect_equal(prob(drg, c(-.3, .1, .2)), prob(dr, c(-.3, .1, .2)))
  expect_error(hurdle("poisson", at = .5), "lattice")
  expect_error(inflate("poisson", at = .5), "lattice")
  expect_no_error(inflate(round_dist("normal", digits = 1), at = .1))
})

test_that("interval masses keep representable tail probabilities", {
  d1 <- distribution(discretize("normal"), mu = 0, sigma = 1)
  expect_equal(log_prob(d1, 9), log(pnorm(9, lower.tail = FALSE) - pnorm(10, lower.tail = FALSE)), tolerance = 1e-8)
  dr <- distribution(round_dist("normal", digits = 1), mu = 0, sigma = 1)
  expect_equal(log_prob(dr, 9), log(pnorm(8.95, lower.tail = FALSE) - pnorm(9.05, lower.tail = FALSE)), tolerance = 1e-8)
  expect_equal(log_prob(dr, -9), log_prob(dr, 9), tolerance = 1e-8)
})

test_that("tiny rate hurdles keep their moments and sampler", {
  for (lam in c(1e-6, 1e-12, 1e-18)) {
    d <- distribution(hurdle("poisson"), lambda = lam, hu = .3)
    expect_equal(mean(d), .7 * (1 + lam / 2), tolerance = 1e-6)
    expect_equal(variance(d), .21, tolerance = 1e-5)
  }
  set.seed(1)
  x <- rng(distribution(hurdle("poisson"), lambda = 1e-18, hu = .3), 2000)
  expect_true(all(x %in% 0:1))
  expect_equal(mean(x == 0), .3, tolerance = .04)
})

test_that("lattice laws fit, score and give PIT values on their lattice", {
  fam <- round_dist("normal", digits = 1)
  expect_no_error(regcomp(y ~ 1, family = fam, data = data.frame(y = rep(c(-.1, 0, .1), 10)), fixed = list(sigma = .1)))
  dat <- data.frame(y = rep(0, 30))
  fit <- regcomp(y ~ 1, family = fam, data = dat, fixed = list(sigma = .1))
  fd <- fitted_distribution(fit, data.frame(y = 0))
  grid <- (-20:20) / 10
  expect_equal(unname(scoring_rule(fit, data.frame(y = 0), rule = "crps")), .1 * sum((cdf(fd, grid) - (0 <= grid))^2), tolerance = 1e-8)
  fd_all <- fitted_distribution(fit)
  flo <- cdf(fd_all, dat$y - .1); fup <- cdf(fd_all, dat$y)
  set.seed(123); expected <- flo + runif(nrow(dat)) * (fup - flo)
  expect_equal(unname(residuals(fit, type = "pit", seed = 123)), unname(expected), tolerance = 1e-8)
  # Unit lattice CRPS against its definition.
  pf <- regcomp(y ~ 1, family = "poisson", data = data.frame(y = c(1, 2, 4)))
  lam <- exp(coef(pf)[[1]]); k <- 0:200
  expect_equal(unname(scoring_rule(pf, data.frame(y = 2), rule = "crps")), sum((ppois(k, lam) - (2 <= k))^2), tolerance = 1e-8)
})

test_that("continuous PIT is the CDF at every scale and uses no random numbers", {
  dat <- data.frame(y = (-2:2) * 1e-8)
  fit <- regcomp(y ~ 1, family = "normal", data = dat, fixed = list(sigma = 1e-8))
  pp <- fitted_distribution(fit)$pars
  expect_equal(unname(residuals(fit, type = "pit", seed = 123)), pnorm(dat$y, pp$mu, pp$sigma), tolerance = 1e-12)
  set.seed(762); before <- .Random.seed
  invisible(residuals(fit, type = "pit"))
  expect_identical(.Random.seed, before)
  # A real atom is still randomized over exactly its mass.
  tob <- regcomp(y ~ 1, family = censor("normal", lower = 0), data = data.frame(y = c(0, 0, .5, 1, 2)))
  u <- residuals(tob, type = "pit", seed = 2)[1:2]
  F0 <- cdf(fitted_distribution(tob, data.frame(y = 0)), 0)
  expect_true(all(u > 0 & u <= F0))
})

test_that("structured parameters and responses respect their domains", {
  expect_error(distribution("matrix_normal", M = matrix(0, 1, 1), U = matrix(-1, 1, 1), V = matrix(1, 1, 1)), class = "regcomp_domain_error")
  expect_error(distribution("matrix_normal", M = matrix(0, 2, 1), U = diag(3), V = matrix(1, 1, 1)), class = "regcomp_domain_error")
  expect_error(distribution("wishart", nu = .5, Sigma = diag(2)), class = "regcomp_domain_error")
  expect_no_error(distribution("wishart", nu = 1.5, Sigma = diag(2)))
  expect_error(distribution("multivariate_normal", mu = c(0, 0), Sigma = matrix(c(1, 2, 2, 1), 2)), class = "regcomp_domain_error")
  ds <- distribution("stiefel_uniform", n = 3, k = 2)
  expect_false(validate_response(ds, matrix(0, 3, 2)))
  expect_equal(prob(ds, matrix(0, 3, 2)), 0)
  expect_true(validate_response(ds, diag(3)[, 1:2]))
  dmo <- distribution("bivariate_marshall_olkin_exponential", lambda1 = 1, lambda2 = 1, lambda12 = 1)
  expect_false(validate_response(dmo, c(-1, -2)))
  expect_equal(prob(dmo, c(-1, -2)), 0)
})

test_that("cure laws keep F + S = 1 above a finite uncured support", {
  dcu <- distribution(cure("uniform"), a = 0, b = 1, cure = .2)
  t <- c(.3, 1, 2)
  expect_equal(cdf(dcu, t) + survival(dcu, t), rep(1, 3))
  expect_equal(cdf(dcu, 1), .8)
  expect_equal(quantile(dcu, .9), Inf)
})

test_that("slicing keeps shared vector and matrix parameters", {
  ord <- distribution("cumulative_logit", mu = c(-1, 0, 1), cuts = c(-1, 1))
  expect_equal(unname(prob(ord[1], 3)), unname(prob(ord, rep(3, 3))[1]), tolerance = 1e-12)
  expect_identical(ord[2:3]$pars$cuts, ord$pars$cuts)
  mv <- distribution("multivariate_normal", mu = c(0, 1), Sigma = diag(2))
  expect_identical(mv[1]$pars, mv$pars)
  d <- distribution("normal", mu = 1:3, sigma = 2)
  expect_equal(d[c(3, 1)]$pars, list(mu = c(3L, 1L), sigma = c(2, 2)))
})

test_that("multinomial and Ewens laws take an explicit sample size", {
  set.seed(1)
  d <- distribution("multinomial", probs = c(.4, .6), size = 2)
  expect_true(all(rowSums(rng(d, 5)) == 2))
  expect_equal(prob(d, c(1, 1)), .48)
  expect_equal(prob(d, c(2, 1)), 0)
  expect_equal(prob(distribution("multinomial", probs = c(.4, .6)), c(3, 1)), dmultinom(c(3, 1), prob = c(.4, .6)))
  expect_error(rng(distribution("multinomial", probs = c(.4, .6)), 2), "size")
  expect_error(rng(d, 4, size = 3), "Unused")
  z <- rng(distribution("ewens", theta = 1, size = 3), 4)
  expect_true(all(drop(z %*% seq_len(ncol(z))) == 3))
})
