test_that("Stan code covers ordinal, multinomial, structured and multiple membership models", {
  set.seed(1); n <- 60; x <- rnorm(n)
  d <- data.frame(x, y = cut(x + rlogis(n), c(-Inf, -1, .3, Inf), ordered_result = TRUE),
                  yn = factor(sample(letters[1:3], n, TRUE)), g = factor(rep(1:10, each = 6)),
                  times = factor(rep(1:6, 10)), tnum = rep(1:6, 10), z = rnorm(n),
                  s1 = factor(sample(5, n, TRUE)), s2 = factor(sample(5, n, TRUE)), w = runif(n))
  expect_match(stan_code(y ~ x, family = "cumulative_logit", data = d)$code, "ordered\\[Kc_cuts\\] cuts")
  expect_match(stan_code(yn ~ x, family = "multinomial_logit", data = d)$code, "categorical_logit_lpmf")
  expect_match(stan_code(z ~ x + ar1(0 + times | g), family = "normal", data = d)$code, "rho_1\\^abs")
  expect_match(stan_code(z ~ x + ou(tnum | g), family = "normal", data = d)$code, "range_1")
  sc <- stan_code(z ~ x + (1 | mm(s1, s2, weights = cbind(w, 1 - w))), family = "normal", data = d)
  expect_match(sc$code, "W_1 \\*")
  expect_equal(dim(sc$data$W_1), c(n, 5L))
})
