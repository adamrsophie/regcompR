# Stan generation: log CDF templates, joint constraints, unique helpers and
# prior matching.

stanc_path <- function() {
  p <- Sys.getenv("STANC", unset = "")
  if (!nzchar(p)) { h <- Sys.glob(file.path(path.expand("~"), ".cmdstan", "cmdstan-*", "bin", "stanc")); p <- if (length(h)) h[length(h)] else "" }
  p
}
expect_stanc <- function(code) {
  sc <- stanc_path(); if (!nzchar(sc) || !file.exists(sc)) return(invisible())
  f <- tempfile(fileext = ".stan"); writeLines(code, f)
  out <- suppressWarnings(system2(sc, c("--o=/dev/null", shQuote(f)), stdout = TRUE, stderr = TRUE))
  expect_null(attr(out, "status"), info = paste(out, collapse = "\n"))
}

test_that("log logistic truncation uses a log CDF", {
  a <- stan_code(y ~ 1, family = truncate_dist("log_logistic", upper = 2), data = data.frame(y = c(.5, .7, 1, 1.5)))
  expect_match(a$code, "log(loglogistic_cdf(2 | scale, shape))", fixed = TRUE)
  expect_stanc(a$code)
})

test_that("joint parameter constraints are enforced in Stan", {
  a <- stan_code(y ~ 1, family = "generalized_log_series", data = data.frame(y = c(1, 2, 3, 2)))
  expect_match(a$code, "if (!(theta[i] * beta[i] < 1)) target += negative_infinity();", fixed = TRUE)
  expect_stanc(a$code)
  m <- stan_code(y ~ 1, family = mixture("generalized_poisson", "poisson"), data = data.frame(y = c(0, 1, 2, 3)))
  expect_match(m$code, "negative_infinity", fixed = TRUE)
  expect_stanc(m$code)
})

test_that("composed operators give unique, valid Stan helpers", {
  d <- data.frame(y = c(0, 1, 2, 3))
  for (f in list(hurdle(mixture("poisson", "poisson")), inflate(mixture("poisson", "poisson")),
                 mixture(hurdle("poisson", at = 0), hurdle("poisson", at = 1)))) {
    a <- stan_code(y ~ 1, family = f, data = d)
    expect_false(grepl("real [a-z]+\\.[0-9]", a$code))
    expect_stanc(a$code)
  }
  a <- stan_code(y ~ 1, family = mixture(truncate_dist("normal", 0, 1), truncate_dist("normal", 2, 3)),
                 data = data.frame(y = c(.3, .7, 2.3, 2.7)))
  fns <- regmatches(a$code, gregexpr("real (rc_[A-Za-z0-9_]+)\\(", a$code))[[1]]
  expect_equal(anyDuplicated(fns), 0)
  expect_stanc(a$code)
})

test_that("priors match dotted parameters, scalar correlations and nothing else", {
  d <- data.frame(y = c(-1, .2, 3, 4.1, -.5, 3.3))
  a <- stan_code(y ~ 1, family = mixture("normal", "normal"), data = d,
                 priors = prior("normal(123, 0.2)", par = "mu.1", class = "Intercept"))
  expect_match(a$code, "b_mu_1[1] ~ normal(123, 0.2);", fixed = TRUE)
  expect_error(stan_code(y ~ 1, family = "normal", data = d, priors = prior("normal(0, 1)", par = "mux", class = "Intercept")),
               "match nothing")
  set.seed(1); dl <- data.frame(y = rnorm(24), t = factor(rep(1:3, 8)), g = factor(rep(1:8, each = 3)))
  for (st in c("ar1", "cs")) {
    f <- stats::as.formula(sprintf("y ~ 1 + %s(0 + t | g)", st))
    a <- stan_code(f, family = "normal", data = dl, priors = prior("uniform(0.2, 0.9)", class = "cor"))
    expect_match(a$code, "rho_1 ~ uniform(0.2, 0.9);", fixed = TRUE)
    expect_error(stan_code(f, family = "normal", data = dl, priors = prior("lkj(2)", class = "cor")), "lkj")
  }
  a <- stan_code(y ~ 1 + (0 + t | g), family = "normal", data = dl, priors = prior("lkj(4)", class = "cor"))
  expect_match(a$code, "lkj_corr_cholesky(4)", fixed = TRUE)
  expect_error(stan_code(y ~ 1 + (1 | g), family = "normal", data = dl, priors = prior("lkj(2)", class = "cor")), "match nothing")
})
