#' Priors for Bayesian fits
#'
#' Priors live on the coefficients of each linear predictor (on the link
#' scale), on random effect standard deviations and on random effect
#' correlations. A prior is written as a distribution call such as
#' `"normal(0, 2.5)"`; the same specification is translated to 'Stan' and to
#' the BUGS language used by 'JAGS' and 'OpenBUGS'.
#'
#' Supported prior distributions: `normal(mean, sd)`,
#' `student_t(df, location, scale)`, `cauchy(location, scale)`,
#' `logistic(location, scale)`, `exponential(rate)`, `gamma(shape, rate)`,
#' `uniform(lower, upper)` and, for correlations, `lkj(eta)`. Standard
#' deviation priors are truncated to positive values.
#'
#' Defaults: `normal(0, 5)` for intercepts and coefficients,
#' `student_t(3, 0, 2.5)` for standard deviations and `lkj(2)` for
#' correlations ('Stan'; BUGS engines use an inverse Wishart with identity
#' scale for correlated random effects, which is noted in the fit).
#'
#' @param dist Prior distribution as a string.
#' @param par Distribution parameter the prior applies to (`NULL` for all).
#' @param class `"b"` (coefficients other than the intercept), `"Intercept"`,
#'   `"sd"` or `"cor"`.
#' @param coef Optional coefficient name (for `class = "b"`) or random effect
#'   term label (for `"sd"` and `"cor"`).
#' @return An object of class `regcomp_prior`; combine several with `c()`.
#' @export
#' @examples
#' p <- c(prior("normal(0, 1)", par = "mu", class = "b"),
#'        prior("student_t(3, 0, 1)", class = "sd"))
#' p
prior <- function(dist, par = NULL, class = c("b", "Intercept", "sd", "cor"),
                  coef = NULL) {
  class <- match.arg(class)
  pd <- parse_prior(dist)
  structure(list(list(dist = pd$name, args = pd$args, par = par, class = class,
                      coef = coef, text = dist)), class = "regcomp_prior")
}

#' @export
c.regcomp_prior <- function(...) {
  structure(unlist(lapply(list(...), unclass), recursive = FALSE), class = "regcomp_prior")
}

#' @export
print.regcomp_prior <- function(x, ...) {
  for (p in x) cat(sprintf("%-28s class %-9s par %-6s coef %s\n", p$text, p$class,
                           if (is.null(p$par)) "(all)" else p$par,
                           if (is.null(p$coef)) "(all)" else p$coef))
  invisible(x)
}

prior_dists <- c(normal = 2, student_t = 3, cauchy = 2, logistic = 2, exponential = 1,
                 gamma = 2, uniform = 2, lkj = 1)

parse_prior <- function(s) {
  e <- tryCatch(str2lang(s), error = function(e) NULL)
  if (!is.call(e)) stop_input(sprintf("Cannot parse the prior '%s'.", s))
  nm <- as.character(e[[1]])
  if (!nm %in% names(prior_dists)) stop_input(sprintf("Unsupported prior distribution '%s'.", nm))
  args <- vapply(as.list(e)[-1], function(a) as.numeric(eval(a, baseenv())), 0)
  if (length(args) != prior_dists[[nm]])
    stop_input(sprintf("The prior %s needs %d arguments.", nm, prior_dists[[nm]]))
  bad <- !all(is.finite(args)) || switch(nm,
    normal = , cauchy = , logistic = args[2] <= 0, student_t = args[1] <= 0 || args[3] <= 0,
    exponential = , lkj = args[1] <= 0, gamma = any(args <= 0), uniform = args[1] >= args[2], FALSE)
  if (bad) stop_input(sprintf("The prior '%s' has invalid parameters.", s))
  list(name = nm, args = args)
}

find_prior <- function(priors, class, par, coef = NULL) {
  # The most specific matching entry wins (coefficient, then parameter,
  # then class); among equally specific entries the last one wins.
  hit <- NULL; best <- -1
  for (p in priors) {
    if (p$class != class) next
    if (!is.null(p$par) && !identical(p$par, par)) next
    if (!is.null(p$coef) && !identical(p$coef, coef)) next
    sp <- 2 * (!is.null(p$coef)) + (!is.null(p$par))
    if (sp >= best) { hit <- p; best <- sp }
  }
  if (is.null(hit)) {
    hit <- switch(class,
      b = list(dist = "normal", args = c(0, 5)),
      Intercept = list(dist = "normal", args = c(0, 5)),
      sd = list(dist = "student_t", args = c(3, 0, 2.5)),
      cor = list(dist = "lkj", args = 2),
      cuts = list(dist = "normal", args = c(0, 5)))
  }
  hit
}

prior_stan <- function(p, target) {
  a <- p$args
  switch(p$dist,
    normal = sprintf("%s ~ normal(%g, %g);", target, a[1], a[2]),
    student_t = sprintf("%s ~ student_t(%g, %g, %g);", target, a[1], a[2], a[3]),
    cauchy = sprintf("%s ~ cauchy(%g, %g);", target, a[1], a[2]),
    logistic = sprintf("%s ~ logistic(%g, %g);", target, a[1], a[2]),
    exponential = sprintf("%s ~ exponential(%g);", target, a[1]),
    gamma = sprintf("%s ~ gamma(%g, %g);", target, a[1], a[2]),
    uniform = sprintf("%s ~ uniform(%g, %g);", target, a[1], a[2]),
    lkj = sprintf("%s ~ lkj_corr_cholesky(%g);", target, a[1]))
}

prior_bugs <- function(p, target, positive = FALSE, dialect = "jags") {
  a <- p$args
  trunc <- if (!positive) "" else if (dialect == "jags") " T(0,)" else " I(0,)"
  switch(p$dist,
    normal = sprintf("%s ~ dnorm(%g, %g)%s", target, a[1], 1 / a[2]^2, trunc),
    student_t = sprintf("%s ~ dt(%g, %g, %g)%s", target, a[2], 1 / a[3]^2, a[1], trunc),
    cauchy = sprintf("%s ~ dt(%g, %g, 1)%s", target, a[1], 1 / a[2]^2, trunc),
    logistic = sprintf("%s ~ dlogis(%g, %g)%s", target, a[1], 1 / a[2], trunc),
    exponential = sprintf("%s ~ dexp(%g)", target, a[1]),
    gamma = sprintf("%s ~ dgamma(%g, %g)", target, a[1], a[2]),
    uniform = sprintf("%s ~ dunif(%g, %g)", target, a[1], a[2]),
    stop_unsupported(sprintf("The prior '%s' in the BUGS language", p$dist)))
}
