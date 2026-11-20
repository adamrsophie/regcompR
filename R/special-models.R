# Special regression methods from the count, GLM and multilevel literature
# that are not plain likelihood fits of one family: discrete time hazard
# data expansion, exact conditional Poisson regression and quantile
# regression for counts.

#' Person period data for discrete time hazard models
#'
#' Expands one row per subject into one row per subject and period at risk
#' (Allison 1982; Singer and Willett 2003; Goldstein 2011, ch. 9). The new
#' column `event` is 1 in the period of the event and 0 otherwise, and
#' `period` counts the periods. A binary regression of `event` on period
#' effects and covariates is the discrete time hazard model: the `logit`
#' link gives the proportional odds model and the `cloglog` link the
#' grouped proportional hazards model; random effects such as `(1 | id)`
#' give frailty and multilevel survival models.
#'
#' @param data Data frame with one row per subject.
#' @param time Name of the column with the last observed period (a
#'   positive integer).
#' @param status Name of the column with 1 for an event in period `time`
#'   and 0 for censoring after it.
#' @param start First period (default 1).
#' @return The expanded data frame with columns `period` and `event`.
#' @export
#' @examples
#' d <- data.frame(id = 1:3, t = c(2, 3, 1), d = c(1, 0, 1), x = c(0.5, -1, 2))
#' pp <- person_period(d, "t", "d")
#' pp
#' fit <- regcomp(event ~ factor(period) + x, family = "bernoulli",
#'                links = list(prob = "cloglog"), data = person_period(
#'                  data.frame(t = rpois(200, 3) + 1, d = rbinom(200, 1, 0.7),
#'                             x = rnorm(200)), "t", "d"))
person_period <- function(data, time, status, start = 1) {
  tt <- data[[time]]
  if (any(tt < start | tt != round(tt), na.rm = TRUE)) stop_input("time must contain integers not below start.")
  reps <- tt - start + 1
  out <- data[rep(seq_len(nrow(data)), reps), , drop = FALSE]
  out$period <- unlist(lapply(reps, function(r) start - 1 + seq_len(r)))
  last <- cumsum(reps)
  out$event <- 0
  out$event[last] <- data[[status]]
  rownames(out) <- NULL
  out
}

#' Exact conditional Poisson regression
#'
#' Conditional inference for one coefficient of a Poisson log linear model
#' (Cox 1970; Hirji 2006), as in Stata's `expoisson` discussed by Hilbe
#' (2014, sec. 9.4). The counts are conditioned on their totals within the
#' strata, which removes the intercept and the stratum effects; the
#' conditional distribution of `T = sum(x * y)` is computed exactly by
#' polynomial convolution. Returns the conditional maximum likelihood
#' estimate (which equals the unconditional Poisson estimate with stratum
#' intercepts when it is finite), the median unbiased estimate when `T` is
#' at an end of its support, the exact two sided p value of `beta = 0`
#' (probability method), and the exact confidence interval.
#'
#' @param y Counts.
#' @param x Integer covariate whose coefficient is tested.
#' @param strata Optional factor of strata (nuisance effects).
#' @param exposure Optional exposure (offset on the rate scale).
#' @param level Confidence level.
#' @return A list with `estimate`, `type` ("CMLE" or "MUE"), `p_value`,
#'   `conf_int`, `statistic` (observed `T`) and the conditional
#'   distribution `support` and `counts` (log weights).
#' @export
#' @examples
#' exact_poisson(y = c(0, 1, 3, 4), x = c(0, 0, 1, 1), exposure = c(10, 12, 9, 11))
exact_poisson <- function(y, x, strata = NULL, exposure = NULL, level = 0.95) {
  n <- length(y)
  if (any(x != round(x))) stop_input("x must be integer valued; rescale it first.")
  if (is.null(strata)) strata <- rep(1, n)
  if (is.null(exposure)) exposure <- rep(1, n)
  xmin <- min(x); xs <- x - xmin
  # Log weights of T over its support from the multinomial generating
  # function (sum_i p_i z^x_i)^N within each stratum.
  logconv <- function(a, b) {
    out <- rep(-Inf, length(a) + length(b) - 1)
    for (i in which(is.finite(a))) {
      idx <- i + seq_along(b) - 1
      out[idx] <- log_sum_exp_pair(out[idx], a[i] + b)
    }
    out
  }
  total <- 0
  for (s in unique(strata)) {
    i <- strata == s; N <- sum(y[i]); if (N == 0) next
    pz <- rep(-Inf, max(xs[i]) + 1)
    for (j in which(i)) pz[xs[j] + 1] <- log_sum_exp_pair(pz[xs[j] + 1], log(exposure[j]))
    pz <- pz - log_sum_exp(pz)
    acc <- 0; base <- pz; e <- N
    while (e > 0) { if (e %% 2 == 1) acc <- logconv(acc, base); e <- e %/% 2; if (e > 0) base <- logconv(base, base) }
    total <- logconv(total, acc)
  }
  support <- seq_along(total) - 1
  tobs <- sum(xs * y)
  keep <- is.finite(total); support <- support[keep]; lw <- total[keep]
  cond <- function(b) { l <- lw + b * support; l - log_sum_exp(l) }
  ti <- which(support == tobs)
  pval <- { lp <- cond(0); min(1, sum(exp(lp[lp <= lp[ti] + 1e-7]))) }
  a <- (1 - level) / 2
  if (length(support) == 1) return(list(estimate = NA_real_, type = "none", p_value = 1, conf_int = c(-Inf, Inf),
    level = level, statistic = tobs + xmin * sum(y), support = support + xmin * sum(y), log_weights = lw))
  upper_tail <- function(b) sum(exp(cond(b)[support >= tobs]))
  lower_tail <- function(b) sum(exp(cond(b)[support <= tobs]))
  rng <- c(-50, 50)
  lo <- if (tobs == min(support)) -Inf else stats::uniroot(function(b) upper_tail(b) - a, rng, extendInt = "upX", tol = 1e-12)$root
  hi <- if (tobs == max(support)) Inf else stats::uniroot(function(b) lower_tail(b) - a, rng, extendInt = "downX", tol = 1e-12)$root
  if (tobs > min(support) && tobs < max(support)) {
    est <- stats::uniroot(function(b) sum(exp(cond(b)) * support) - tobs, rng, extendInt = "upX", tol = 1e-12)$root
    type <- "CMLE"
  } else {
    # Endpoint median unbiased estimate: the nontrivial inclusive tail at 1/2.
    est <- if (tobs == min(support)) stats::uniroot(function(b) lower_tail(b) - 0.5, rng, extendInt = "downX", tol = 1e-12)$root
      else stats::uniroot(function(b) upper_tail(b) - 0.5, rng, extendInt = "upX", tol = 1e-12)$root
    type <- "MUE"
  }
  list(estimate = est, type = type, p_value = pval, conf_int = c(lo, hi), level = level,
       statistic = tobs + xmin * sum(y), support = support + xmin * sum(y), log_weights = lw)
}

log_sum_exp_pair <- function(a, b) {
  m <- pmax(a, b)
  ifelse(is.infinite(m) & m < 0, -Inf, m + log(exp(a - m) + exp(b - m)))
}

#' Quantile regression for counts
#'
#' The jittering method of Machado and Santos Silva (2005), discussed by
#' Hilbe (2014, sec. 9.5): `z = y + U` with uniform noise, the transform
#' `log(z - tau)` (with a small floor below `tau`) is fitted by linear
#' quantile regression, and estimates are averaged over `m` jittered
#' samples. The conditional quantile of the count is
#' `ceiling(tau + exp(x' beta)) - 1`.
#'
#' @param formula Formula `count ~ predictors`.
#' @param data Data frame.
#' @param tau Quantile level in (0, 1).
#' @param m Number of jittered samples.
#' @param seed Seed for the jittering.
#' @return A list with `coefficients`, their `vcov` (between plus mean
#'   within jitter variance), `tau` and a `predict(newdata)` function giving
#'   the count quantiles.
#' @export
#' @examples
#' if (requireNamespace("quantreg", quietly = TRUE)) {
#'   set.seed(1)
#'   d <- data.frame(x = rnorm(300))
#'   d$y <- rpois(300, exp(1 + 0.4 * d$x))
#'   q <- quantile_count(y ~ x, d, tau = 0.5, m = 20)
#'   q$coefficients
#' }
quantile_count <- function(formula, data, tau = 0.5, m = 100, seed = 1) {
  if (!requireNamespace("quantreg", quietly = TRUE))
    stop_unsupported("Quantile regression for counts needs the 'quantreg' package; it")
  mf <- stats::model.frame(formula, data)
  y <- stats::model.response(mf); X <- stats::model.matrix(stats::terms(mf), mf)
  old <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, globalenv()))
  set.seed(seed)
  B <- matrix(0, m, ncol(X)); V <- 0
  for (j in seq_len(m)) {
    z <- y + stats::runif(length(y))
    tz <- log(pmax(z - tau, 1e-5))
    f <- quantreg::rq(tz ~ X - 1, tau = tau)
    B[j, ] <- f$coefficients
    V <- V + summary(f, se = "ker", covariance = TRUE)$cov
  }
  beta <- colMeans(B); names(beta) <- colnames(X)
  vc <- V / m + stats::cov(B); dimnames(vc) <- list(colnames(X), colnames(X))
  tt <- stats::delete.response(stats::terms(mf))
  list(coefficients = beta, vcov = vc, tau = tau, m = m,
       predict = function(newdata) { Xn <- stats::model.matrix(tt, stats::model.frame(tt, newdata))
         ceiling(tau + exp(drop(Xn %*% beta))) - 1 })
}
