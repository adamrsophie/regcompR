#' Penalized quasi-likelihood for generalized linear mixed models
#'
#' Fits a GLMM by penalized quasi-likelihood (Breslow and Clayton 1993):
#' repeated linear mixed model fits, by maximum likelihood, to the working
#' response `z = eta + (y - mu) deta/dmu` with variance weights
#' `(dmu/deta)^2 / V(mu)`, as `MASS::glmmPQL()`. PQL is not a likelihood
#' method (no log likelihood, AIC or likelihood ratio tests); it is biased
#' for binary data with few observations per cluster, where the Laplace
#' fit of [regcomp()] is preferable. It handles quasi families.
#'
#' @param formula A formula with lme4 style random effects, `y ~ x + (1 | g)`.
#' @param family A [stats::family()] object (such as `binomial()`,
#'   `poisson()`, `quasipoisson()`).
#' @param data A data frame.
#' @param niter Maximum number of iterations.
#' @param tol Convergence tolerance on the linear predictor.
#' @return The final working normal mixed model fit (a `regcomp_fit`),
#'   whose `mu` coefficients are the fixed effects on the link scale and
#'   whose `theta` are the random effect parameters, with class
#'   `regcomp_pql` added and the number of iterations in `iterations`.
#' @export
#' @examples
#' set.seed(1)
#' d <- data.frame(g = factor(rep(1:20, each = 8)), x = rnorm(160))
#' d$y <- rpois(160, exp(0.3 + 0.4 * d$x + rnorm(20, 0, 0.5)[d$g]))
#' pql_regcomp(y ~ x + (1 | g), family = poisson(), data = d)
pql_regcomp <- function(formula, family, data, niter = 10, tol = 1e-6) {
  if (is.character(family)) family <- get(family, mode = "function")
  if (is.function(family)) family <- family()
  if (is.null(family$family)) stop_input("family must be a stats family object, such as binomial().")
  fixed <- formula; fixed[[3]] <- drop_bars(formula[[3]])
  if (is.null(fixed[[3]])) fixed[[3]] <- 1
  g0 <- stats::glm(fixed, family = family, data = data)
  d <- data[as.integer(rownames(stats::model.frame(g0))), , drop = FALSE]
  y <- g0$y; eta <- g0$linear.predictors; off <- if (is.null(g0$offset)) 0 else g0$offset
  wf <- formula; wf[[2]] <- as.name(".pql_z")
  ctl <- regcomp_control()
  for (it in seq_len(niter)) {
    mu <- family$linkinv(eta); me <- family$mu.eta(eta)
    d$.pql_z <- eta - off + (y - mu) / me
    d$.pql_lo <- -0.5 * log(me^2 / family$variance(mu))
    fit <- regcomp(wf, sigma ~ 1 + offset(.pql_lo), family = "normal", data = d, control = ctl)
    ctl$start <- NULL
    old <- eta
    eta <- as.numeric(stats::predict(fit, type = "mean", re = "conditional")) + off
    if (sum((eta - old)^2) < tol * sum(eta^2)) break
  }
  fit$iterations <- it; fit$pql_family <- family
  class(fit) <- c("regcomp_pql", class(fit))
  fit
}

#' @export
print.regcomp_pql <- function(x, digits = 4, ...) {
  cat(sprintf("<regcomp_pql> %s (%s link) GLMM by penalized quasi-likelihood, %d iterations, n = %d\n",
              x$pql_family$family, x$pql_family$link, x$iterations, x$nobs))
  V <- x$vcov[x$spec$beta_idx$mu, x$spec$beta_idx$mu, drop = FALSE]
  b <- x$coefficients[x$spec$beta_idx$mu]; se <- sqrt(diag(V))
  print(round(data.frame(estimate = b, std_error = se, z = b / se, row.names = sub("^mu:", "", names(b))), digits))
  cat("Random effect parameters:\n"); print(round(x$theta, digits))
  invisible(x)
}

#' @export
logLik.regcomp_pql <- function(object, ...) stop_unsupported("A log likelihood for a penalized quasi-likelihood fit")
