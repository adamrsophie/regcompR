# Predictions, fitted distributions, simulation and pointwise likelihood.

new_spec <- function(fit, newdata, response = FALSE) {
  if (is.null(newdata)) return(fit$spec)
  spec <- fit$spec
  main <- spec$formulas[[primary_par(spec$fam, spec$rpars)]]
  ns <- build_spec(main, new_spec_dots(spec), spec$fam, as.data.frame(newdata),
                   fixed = c(spec$fixed), env = environment(main), prev = spec, response = response)
  ns$links <- spec$links
  ns
}

# Linear predictors for a coefficient vector and random effects.
pred_eta <- function(spec, beta, u, re) {
  uu <- if (re == "conditional" && length(u)) u else NULL
  compute_eta(spec, beta, uu)
}

pars_from <- function(spec, beta, u, re) eta_to_pars(spec, pred_eta(spec, beta, u, re))

#' Fitted distributions
#'
#' Returns the conditional distribution of the response for each row of
#' `newdata` (or the data used for fitting) as a vectorized distribution
#' object, so every method of [distribution()] (quantiles, probabilities,
#' moments, plots) applies to a fitted model. For Bayesian fits the
#' parameters are posterior means of each distribution parameter; use
#' [predict.regcomp_fit()] for full posterior summaries.
#'
#' @param fit A fitted model.
#' @param newdata Optional data frame.
#' @param re Random effects: `"conditional"` uses the estimated random
#'   effects (new levels get zero), `"zero"` sets all random effects to
#'   zero.
#' @return A `regcomp_dist` object.
#' @export
#' @examples
#' set.seed(1)
#' d <- data.frame(x = runif(100))
#' d$y <- rgamma(100, shape = 4, rate = 4 / exp(1 + d$x))
#' fit <- regcomp(y ~ x, family = "gamma", data = d)
#' fd <- fitted_distribution(fit, data.frame(x = c(0, 1)))
#' quantile(fd, 0.9)
fitted_distribution <- function(fit, newdata = NULL, re = c("conditional", "zero")) {
  re <- match.arg(re)
  spec <- new_spec(fit, newdata)
  if (fit$inference == "bayesian") {
    S <- nrow(fit$draws)
    acc <- NULL
    for (s in seq_len(S)) {
      p <- pars_from(spec, fit$draws[s, spec$beta_names], if (length(fit$u_draws)) fit$u_draws[s, ] else NULL, re)
      acc <- if (is.null(acc)) p else Map(`+`, acc, p)
    }
    pars <- lapply(acc, `/`, S)
  } else pars <- pars_from(spec, fit$coefficients, fit$u, re)
  structure(list(family = spec$fam, pars = pars), class = "regcomp_dist")
}

quantity_from_pars <- function(fam, pars, type, probs, at) {
  d <- structure(list(family = fam, pars = pars), class = "regcomp_dist")
  switch(type,
    mean = suppressWarnings(mean(d)),
    variance = suppressWarnings(variance(d)),
    sd = sqrt(suppressWarnings(variance(d))),
    median = family_quantile(fam, rep(0.5, length(d)), recycle_pars(pars)),
    quantile = sapply(probs, function(p) family_quantile(fam, rep(p, length(d)), recycle_pars(pars))),
    probability = sapply(at, function(a) exp(eval_lpdf(fam, rep(a, length(d)), pars))),
    density = sapply(at, function(a) exp(eval_lpdf(fam, rep(a, length(d)), pars))),
    cdf = sapply(at, function(a) family_cdf(fam, rep(a, length(d)), recycle_pars(pars))),
    survival = sapply(at, function(a) 1 - family_cdf(fam, rep(a, length(d)), recycle_pars(pars))),
    stop_input(sprintf("Unknown prediction type '%s'.", type)))
}

#' Predictions from a fitted model
#'
#' @param object A fitted model.
#' @param newdata Optional data frame.
#' @param type What to predict: `"parameter"` (all distribution parameters),
#'   `"link"` (linear predictors), `"mean"`, `"variance"`, `"sd"`,
#'   `"median"`, `"quantile"` (at `probs`), `"probability"` (mass or
#'   density at `at`), `"density"`, `"cdf"` or `"survival"` (at `at`).
#'   For ordinal and nominal families `"mean"` and `"probability"` give
#'   the matrix of category probabilities (all categories by default).
#' @param probs Probabilities for `type = "quantile"`.
#' @param at Response values for `"probability"`, `"density"`, `"cdf"` and
#'   `"survival"`.
#' @param interval `"none"`, or `"confidence"` for uncertainty of the
#'   predicted quantity. Frequentist fits use parametric simulation from the
#'   asymptotic normal distribution of the estimates (`ndraws` draws);
#'   Bayesian fits use posterior draws and return credible intervals.
#' @param level Interval level.
#' @param re `"conditional"`, `"zero"` or `"marginal"`. `"marginal"`
#'   averages the quantity over the random effect distribution by Monte
#'   Carlo integration (`nmarginal` draws), giving population averaged
#'   predictions.
#' @param ndraws Number of simulation draws for confidence intervals.
#' @param nmarginal Number of random effect draws for `re = "marginal"`.
#' @param seed Optional seed for the simulation.
#' @param ... Unused.
#' @return A numeric vector, a matrix (several `probs` or `at` values, or
#'   `type = "parameter"`), or a data frame with `estimate`, `lower` and
#'   `upper` when an interval is requested.
#' @export
predict.regcomp_fit <- function(object, newdata = NULL,
                                type = c("mean", "parameter", "link", "variance", "sd",
                                         "median", "quantile", "probability", "density",
                                         "cdf", "survival"),
                                probs = c(0.1, 0.5, 0.9), at = 0,
                                interval = c("none", "confidence"), level = 0.95,
                                re = c("conditional", "zero", "marginal"),
                                ndraws = 1000, nmarginal = 200, seed = NULL, ...) {
  type <- match.arg(type); interval <- match.arg(interval); re <- match.arg(re)
  if (!is.null(seed)) set.seed(seed)
  spec <- new_spec(object, newdata)
  fam <- spec$fam
  if (fam$variate_type == "categorical") {
    if (type == "mean") type <- "probability"
    if (missing(at)) at <- seq_along(spec$ylev %||% seq_len(max(object$spec$y, na.rm = TRUE)))
    if (!is.null(object$spec$ylev)) attr(at, "labels") <- object$spec$ylev[at]
  }
  one <- function(beta, u, theta) {
    if (type == "link") return(do.call(cbind, pred_eta(spec, beta, u, if (re == "marginal") "zero" else re)))
    if (type == "parameter") return(do.call(cbind, pars_from(spec, beta, u, if (re == "marginal") "zero" else re)))
    if (re == "marginal" && spec$Q > 0) {
      acc <- 0
      for (m in seq_len(nmarginal)) {
        ur <- draw_re(spec, theta)
        acc <- acc + quantity_from_pars(fam, pars_from(spec, beta, ur, "conditional"), type, probs, at)
      }
      return(acc / nmarginal)
    }
    quantity_from_pars(fam, pars_from(spec, beta, u, re), type, probs, at)
  }
  bayes <- object$inference == "bayesian"
  if (!bayes && interval == "none") {
    out <- one(object$coefficients, object$u, object$theta)
    return(label_pred(out, type, probs, at))
  }
  if (bayes) {
    S <- nrow(object$draws)
    vals <- lapply(seq_len(S), function(s) one(object$draws[s, spec$beta_names],
      if (length(object$u_draws)) object$u_draws[s, ] else NULL,
      object$draws[s, names(object$theta)]))
  } else {
    V <- object$vcov
    if (is.null(V) || anyNA(V)) stop_unsupported("Confidence intervals without a valid covariance matrix")
    est <- c(object$coefficients, object$theta)
    R <- tryCatch(chol(V[seq_along(est), seq_along(est)]), error = function(e) NULL)
    if (is.null(R)) stop_numerical("The covariance matrix is not positive definite.")
    np <- length(object$coefficients)
    vals <- lapply(seq_len(ndraws), function(s) {
      z <- est + drop(stats::rnorm(length(est)) %*% R)
      one(z[seq_len(np)], object$u, z[-seq_len(np)])
    })
  }
  arr <- simplify2array(vals)
  a <- (1 - level) / 2
  dims <- length(dim(arr))
  mar <- seq_len(dims - 1)
  est <- if (bayes) apply(arr, mar, mean) else one(object$coefficients, object$u, object$theta)
  if (interval == "none") return(label_pred(est, type, probs, at))
  lo <- apply(arr, mar, stats::quantile, a, na.rm = TRUE)
  hi <- apply(arr, mar, stats::quantile, 1 - a, na.rm = TRUE)
  if (is.null(dim(est))) {
    return(data.frame(estimate = est, lower = lo, upper = hi))
  }
  list(estimate = label_pred(est, type, probs, at), lower = label_pred(lo, type, probs, at),
       upper = label_pred(hi, type, probs, at))
}

label_pred <- function(x, type, probs, at) {
  if (is.matrix(x)) {
    if (type == "quantile") colnames(x) <- paste0("q", format(probs))
    if (type %in% c("probability", "density", "cdf", "survival"))
      colnames(x) <- paste0(type, "(", if (is.null(attr(at, "labels"))) format(at) else attr(at, "labels"), ")")
  }
  x
}

draw_re <- function(spec, theta) {
  u <- numeric(spec$Q)
  for (j in seq_along(spec$blocks)) {
    b <- spec$blocks[[j]]
    S <- block_cov(b, theta[spec$theta_idx[[j]]])
    R <- chol(S)
    u[spec$u_idx[[j]]] <- as.vector(t(matrix(stats::rnorm(b$q * b$L), b$L, b$q) %*% R))
  }
  u
}

#' Simulate responses from a fitted model
#'
#' Draws new responses from the fitted conditional distributions. Bayesian
#' fits draw from the posterior predictive distribution (one posterior draw
#' per simulated data set).
#'
#' @param object A fitted model.
#' @param nsim Number of simulated data sets.
#' @param seed Optional seed.
#' @param newdata Optional data frame.
#' @param re `"conditional"` keeps the estimated random effects; `"new"`
#'   draws new random effects from their estimated distribution.
#' @param ... Unused.
#' @return A matrix with one column per simulation.
#' @export
simulate.regcomp_fit <- function(object, nsim = 1, seed = NULL, newdata = NULL,
                                 re = c("conditional", "new"), ...) {
  re <- match.arg(re)
  if (!is.null(seed)) set.seed(seed)
  spec <- new_spec(object, newdata)
  sims <- sapply(seq_len(nsim), function(s) {
    if (object$inference == "bayesian") {
      i <- sample.int(nrow(object$draws), 1)
      beta <- object$draws[i, spec$beta_names]
      theta <- object$draws[i, names(object$theta)]
      u <- if (length(object$u_draws)) object$u_draws[i, ] else NULL
    } else { beta <- object$coefficients; theta <- object$theta; u <- object$u }
    if (re == "new" && spec$Q > 0) u <- draw_re(spec, theta)
    p <- pars_from(spec, beta, u, "conditional")
    family_rng(spec$fam, spec$n, p)
  })
  if (is.null(dim(sims))) sims <- matrix(sims, ncol = nsim)
  sims
}

#' Pointwise log likelihood
#'
#' Returns the log likelihood contribution of each observation (including
#' censoring and truncation). For Bayesian fits returns a draws by
#' observations matrix suitable for 'loo'.
#'
#' @param object A fitted model.
#' @param newdata Optional data frame containing the response.
#' @param ... Unused.
#' @return Numeric vector or matrix.
#' @export
log_lik <- function(object, newdata = NULL, ...) {
  spec <- if (is.null(newdata)) object$spec else {
    new_spec(object, newdata, response = TRUE)
  }
  f <- function(beta, u) obs_loglik(spec$fam, spec$y, pars_from(spec, beta, u, "conditional"), spec$cmat)
  if (object$inference == "bayesian") {
    t(sapply(seq_len(nrow(object$draws)), function(s)
      f(object$draws[s, spec$beta_names], if (length(object$u_draws)) object$u_draws[s, ] else NULL)))
  } else f(object$coefficients, object$u)
}

#' @export
fitted.regcomp_fit <- function(object, ...) predict(object, type = "mean")
