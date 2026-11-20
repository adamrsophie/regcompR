grid_W <- function(k) {
  g <- expand.grid(i = 1:k, j = 1:k); A <- nrow(g)
  W <- outer(1:A, 1:A, function(a, b) as.numeric(abs(g$i[a] - g$i[b]) + abs(g$j[a] - g$j[b]) == 1))
  dimnames(W) <- list(paste0("r", 1:A), paste0("r", 1:A)); W
}

test_that("relatedness, CAR, Leroux and intrinsic CAR effects give the exact normal marginal likelihood", {
  skip_if_not_installed("mvtnorm")
  set.seed(1); W <- grid_W(5); A <- nrow(W)
  d <- data.frame(region = factor(rep(rownames(W), each = 3), levels = rownames(W)), x = rnorm(3 * A))
  phi <- as.vector(t(chol(solve(diag(rowSums(W)) - 0.9 * W))) %*% rnorm(A))
  d$y <- 1 + 0.5 * d$x + phi[as.integer(d$region)] + rnorm(3 * A, 0, 0.5)
  K <- solve(diag(rowSums(W)) - 0.5 * W); dimnames(K) <- dimnames(W)
  Z <- outer(as.integer(d$region), 1:A, `==`) * 1; X <- cbind(1, d$x)
  for (st in c("relmat", "car", "leroux", "icar")) {
    M <- if (st == "relmat") K else W
    f <- regcomp(stats::as.formula(sprintf("y ~ x + %s(region | 1, M)", st)), family = "normal", data = d)
    b <- f$spec$blocks[[1]]; np <- length(f$theta)
    nll <- function(p) -mvtnorm::dmvnorm(d$y, X %*% p[1:2], exp(2 * p[3]) * diag(nrow(d)) +
      Z %*% regcompR:::block_cov(b, p[-(1:3)]) %*% t(Z), log = TRUE)
    o <- optim(c(1, 0.5, log(0.5), rep(0, np)), nll, method = "BFGS", control = list(reltol = 1e-12, maxit = 1000))
    expect_equal(as.numeric(logLik(f)), -o$value, tolerance = 1e-6, info = st)
    expect_equal(unname(c(coef(f), f$theta)), o$par, tolerance = 1e-3, info = st)
  }
  expect_error(regcomp(y ~ x + car(region | 1, W[-1, -1]), family = "normal", data = d), "not rows")
  expect_error(regcomp(y ~ x + car(region | 1, W), family = "normal", data = d, engine = "tmb"), class = "regcomp_unsupported")
})

test_that("a separable space time effect gives the exact normal marginal likelihood", {
  skip_if_not_installed("mvtnorm")
  set.seed(1); W <- grid_W(3); A <- nrow(W)
  d <- expand.grid(time = 1:4, region = rownames(W)); d <- rbind(d, d); d$x <- rnorm(nrow(d))
  u <- rnorm(A * 4); d$y <- 1 + 0.5 * d$x + u[(match(d$region, rownames(W)) - 1) * 4 + d$time] + rnorm(nrow(d), 0, 0.5)
  f <- regcomp(y ~ x + spacetime(region + time | 1, W), family = "normal", data = d)
  b <- f$spec$blocks[[1]]; Z <- b$Xr; X <- cbind(1, d$x)
  nll <- function(p) -mvtnorm::dmvnorm(d$y, X %*% p[1:2], exp(2 * p[3]) * diag(nrow(d)) +
    Z %*% regcompR:::block_cov(b, p[4:6]) %*% t(Z), log = TRUE)
  o <- optim(c(1, 0.5, log(0.5), 0, 0, 0), nll, method = "BFGS", control = list(reltol = 1e-12, maxit = 2000))
  expect_equal(as.numeric(logLik(f)), -o$value, tolerance = 1e-6)
})

test_that("Stan code for areal effects uses the known matrices", {
  W <- grid_W(3); d <- data.frame(region = factor(rep(rownames(W), 2)), y = rpois(18, 2))
  sc <- stan_code(y ~ leroux(region | 1, W), family = "poisson", data = d)
  expect_match(sc$code, "inverse_spd", fixed = TRUE)
  expect_equal(sc$data$Wa_1, unname(W))
  expect_match(stan_code(y ~ icar(region | 1, W), family = "poisson", data = d)$code, "Kch_1", fixed = TRUE)
})
