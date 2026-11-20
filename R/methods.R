# Methods for fitted models.

#' @export
print.regcomp_fit <- function(x, ...) {
  cat(sprintf("<regcomp_fit> %s regression, engine '%s' (%s)\n",
              x$family$display_name, x$engine, x$inference))
  if (x$approximation != "none") cat("  marginal likelihood approximation:", x$approximation, "\n")
  for (k in x$spec$rpars) cat(sprintf("  %s (%s link): %s\n", k, x$spec$links[[k]]$name,
                                      paste(deparse(rhs_of(x$spec$formulas[[k]])), collapse = "")))
  if (x$inference == "frequentist") {
    cat(sprintf("  log likelihood %.3f on %d parameters, n = %d\n", x$logLik, x$df, x$nobs))
    if (!isTRUE(x$diagnostics$converged)) cat("  WARNING: the optimizer did not report convergence\n")
  } else {
    cat(sprintf("  %d posterior draws, runner '%s'\n", nrow(x$draws), x$runner))
  }
  invisible(x)
}

#' @export
coef.regcomp_fit <- function(object, ...) {
  if (object$inference == "bayesian") colMeans(object$draws[, object$spec$beta_names, drop = FALSE])
  else object$coefficients
}

#' @export
logLik.regcomp_fit <- function(object, ...) {
  if (object$inference == "bayesian")
    stop_unsupported("A maximized log likelihood for a Bayesian fit; use log_lik() for pointwise values")
  structure(object$logLik, df = object$df, nobs = object$nobs, class = "logLik")
}

#' @export
nobs.regcomp_fit <- function(object, ...) object$nobs

#' Random effect estimates
#'
#' Returns the conditional modes (frequentist fits) or posterior means
#' (Bayesian fits) of the random effects, one data frame per term, and the
#' estimated standard deviations and correlations.
#'
#' @param object A fitted model.
#' @param ... Unused.
#' @return A list with `effects` (list of data frames) and `varcorr`.
#' @export
#' @examples
#' set.seed(3)
#' g <- gl(20, 10)
#' d <- data.frame(g, x = rnorm(200))
#' d$y <- rpois(200, exp(0.3 + 0.4 * d$x + rnorm(20, 0, 0.5)[g]))
#' fit <- regcomp(y ~ x + (1 | g), family = "poisson", data = d)
#' ranef_regcomp(fit)$varcorr
ranef_regcomp <- function(object, ...) {
  spec <- object$spec
  u <- if (object$inference == "bayesian") colMeans(object$u_draws) else object$u
  th <- if (object$inference == "bayesian") colMeans(object$draws[, names(object$theta), drop = FALSE]) else object$theta
  eff <- list(); vc <- list()
  for (j in seq_along(spec$blocks)) {
    b <- spec$blocks[[j]]
    uj <- u[spec$u_idx[[j]]]
    m <- matrix(uj, ncol = b$q, byrow = TRUE, dimnames = list(b$levels, b$terms))
    eff[[b$label]] <- as.data.frame(m)
    S <- block_cov(b, th[spec$theta_idx[[j]]])
    vc[[b$label]] <- list(sd = sqrt(diag(S)), cor = stats::cov2cor(S))
  }
  list(effects = eff, varcorr = vc)
}

#' Summarize a fitted model
#'
#' Frequentist fits report estimates, Wald standard errors, z statistics and
#' p values (labeled as Wald; penalized fits show no p values). Bayesian fits
#' report posterior means, standard deviations, quantiles, R-hat and
#' effective sample sizes. The two are different objects and are never
#' presented as interchangeable.
#'
#' @param object A fitted model.
#' @param level Interval level.
#' @param ... Unused.
#' @return An object of class `summary.regcomp_fit`.
#' @export
summary.regcomp_fit <- function(object, level = 0.95, ...) {
  if (object$inference == "bayesian") {
    dr <- object$draws
    a <- (1 - level) / 2
    tab <- data.frame(
      mean = colMeans(dr), sd = apply(dr, 2, stats::sd),
      lower = apply(dr, 2, stats::quantile, a), upper = apply(dr, 2, stats::quantile, 1 - a),
      rhat = object$rhat[colnames(dr)], ess_bulk = object$ess[colnames(dr)])
    out <- list(table = tab, fit = object, level = level, kind = "posterior")
  } else {
    est <- c(object$coefficients, object$theta)
    se <- rep(NA_real_, length(est))
    if (!is.null(object$vcov)) { dv <- diag(object$vcov); dv[dv < 0] <- NA; se <- sqrt(dv)[seq_along(est)] }
    z <- est / se
    q <- stats::qnorm(1 - (1 - level) / 2)
    tab <- data.frame(estimate = est, std_error = se, lower = est - q * se, upper = est + q * se,
                      z = z, p_value = 2 * stats::pnorm(-abs(z)))
    if (!is.null(object$penalty) && !object$penalty$type %in% pen_valid_cov) tab$p_value <- NA_real_
    out <- list(table = tab, fit = object, level = level, kind = "wald")
  }
  class(out) <- "summary.regcomp_fit"
  out
}

#' @export
print.summary.regcomp_fit <- function(x, digits = 4, ...) {
  print(x$fit)
  cat(if (x$kind == "wald") sprintf("\nWald estimates with %g%% confidence intervals:\n", 100 * x$level)
      else sprintf("\nPosterior summary with %g%% credible intervals:\n", 100 * x$level))
  print(round(as.matrix(x$table), digits))
  if (x$kind == "wald") {
    d <- x$fit$diagnostics
    msg <- c(if (!is.null(d$max_abs_gradient)) sprintf("max |gradient| %.2e", d$max_abs_gradient),
             if (!is.null(d$hessian_positive_definite))
               sprintf("Hessian %s", if (d$hessian_positive_definite) "positive definite" else "NOT positive definite"))
    if (length(msg)) cat(paste(msg, collapse = "; "), "\n")
    if (length(d$boundary_variance)) cat("Variance components near zero:", d$boundary_variance, "\n")
  }
  invisible(x)
}

#' @export
confint.regcomp_fit <- function(object, parm, level = 0.95, ...) {
  s <- summary(object, level = level)$table
  out <- as.matrix(s[, c("lower", "upper")])
  if (!missing(parm)) out <- out[parm, , drop = FALSE]
  out
}

#' @export
AIC.regcomp_fit <- function(object, ..., k = 2) ic_table(list(object, ...), k, "AIC", substitute(list(object, ...)))

#' @export
BIC.regcomp_fit <- function(object, ...) ic_table(list(object, ...), NULL, "BIC", substitute(list(object, ...)))

# One value for one model; a data frame of df and criterion for several, as in stats.
ic_table <- function(fits, k, name, calls) {
  if (any(vapply(fits, function(f) identical(f$inference, "bayesian"), TRUE)))
    stop_unsupported(sprintf("%s for a Bayesian fit; use log_lik() with an information criterion such as LOO", name))
  val <- vapply(fits, function(f) -2 * f$logLik + (if (is.null(k)) log(f$nobs) else k) * f$df, 0)
  if (length(fits) == 1) return(val)
  out <- data.frame(df = vapply(fits, function(f) f$df, 0), val)
  names(out)[2] <- name
  rownames(out) <- vapply(as.list(calls)[-1], deparse1, "")
  out
}

#' @export
formula.regcomp_fit <- function(x, ...) x$spec$formulas
