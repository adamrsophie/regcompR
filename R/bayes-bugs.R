# BUGS language program generation for JAGS and OpenBUGS.

bugs_inv <- function(link, x) {
  switch(link,
    identity = x, log = sprintf("exp(%s)", x), softplus = sprintf("log(1 + exp(%s))", x),
    sqrt = sprintf("pow(%s, 2)", x), inverse = sprintf("1 / (%s)", x),
    logit = sprintf("1 / (1 + exp(-(%s)))", x), probit = sprintf("phi(%s)", x),
    cloglog = sprintf("1 - exp(-exp(%s))", x), cauchit = sprintf("(arctan(%s) / 3.141592653589793 + 0.5)", x),
    loglog = sprintf("exp(-exp(-(%s)))", x), tan_half = sprintf("(2 * arctan(%s))", x),
    logshift1 = sprintf("(1 + exp(%s))", x), logshift2 = sprintf("(2 + exp(%s))", x),
    stop_unsupported(sprintf("The link '%s' in the BUGS language", link)))
}

#' Generate the BUGS program of a model
#'
#' Returns the BUGS language program run by the `"bayes"` engine with the
#' `"jags"` or `"openbugs"` runners, and its data list. Families with a
#' native BUGS distribution of the same parameterization use it directly;
#' other families with a closed form log density use the zeros trick
#' (a Poisson observation of zero with mean `C - log f(y)`), which gives the
#' same posterior. Families whose density needs series, special functions
#' outside the BUGS language or numerical integration are not available in
#' BUGS and raise `regcomp_unsupported`; use a 'Stan' runner for them.
#'
#' @inheritParams stan_code
#' @param dialect `"jags"` or `"openbugs"` (they differ in truncation syntax).
#' @return A list with `code` and `data`.
#' @export
#' @examples
#' d <- data.frame(y = rpois(20, 3), x = rnorm(20))
#' cat(bugs_code(y ~ x, family = "poisson", data = d)$code)
bugs_code <- function(formula, ..., family, data, weights = NULL, links = NULL,
                      fixed = NULL, priors = NULL, parameterization = NULL,
                      dialect = c("jags", "openbugs")) {
  dialect <- match.arg(dialect)
  fam <- get_family(family, parameterization)
  dots <- eval(substitute(alist(...)))
  spec <- build_spec(formula, dots, fam, as.data.frame(data), weights = substitute(weights),
                     links = links, fixed = fixed, env = parent.frame())
  list(code = bugs_program(spec, priors, dialect), data = bugs_data(spec))
}

bugs_program <- function(spec, priors = NULL, dialect = "jags") {
  fam <- spec$fam
  bg <- fam$bugs
  if (is.null(bg)) stop_unsupported("A BUGS language implementation", fam$id)
  if (!is.null(spec$cmat)) stop_unsupported("Censored or truncated responses in the BUGS engines", fam$id)
  useweights <- any(spec$weights != 1)
  if (useweights && is.null(bg$loglik))
    stop_unsupported("Weights with this family in the BUGS engines", fam$id)
  L <- c("model {", "  for (i in 1:N) {")
  for (k in spec$rpars) {
    e <- sprintf("inprod(X_%s[i, ], b_%s[]) + off_%s[i]", k, k, k)
    for (j in which(vapply(spec$blocks, `[[`, "", "par") == k)) {
      b <- spec$blocks[[j]]
      e <- if (b$type == "smooth") paste0(e, sprintf(" + inprod(Zs_%d[i, ], s_%d[])", j, j))
           else paste0(e, sprintf(" + inprod(Xr_%d[i, ], r_%d[g_%d[i], ])", j, j, j))
    }
    L <- c(L, sprintf("    eta_%s[i] <- %s", k, e),
           sprintf("    %s[i] <- %s", k, bugs_inv(spec$links[[k]]$name, sprintf("eta_%s[i]", k))))
  }
  for (k in names(spec$fixed)) L <- c(L, sprintf("    %s[i] <- fix_%s", k, k))
  pm <- stats::setNames(lapply(names(fam$pars), function(k) sprintf("%s[i]", k)), names(fam$pars))
  if (!useweights && !is.null(bg$native)) {
    L <- c(L, sprintf("    y[i] ~ %s", stan_fill(bg$native, "y[i]", pm)))
  } else {
    L <- c(L, sprintf("    ll[i] <- %s", stan_fill(bg$loglik, "y[i]", pm)),
           "    zeros[i] ~ dpois(10000 - w[i] * ll[i])")
  }
  L <- c(L, "  }")
  for (k in spec$rpars) {
    cn <- colnames(spec$preds[[k]]$X)
    for (j in seq_along(cn)) {
      cls <- if (cn[j] == "(Intercept)") "Intercept" else "b"
      L <- c(L, "  " %+% prior_bugs(find_prior(priors, cls, k, cn[j]), sprintf("b_%s[%d]", k, j), dialect = dialect))
    }
  }
  for (j in seq_along(spec$blocks)) {
    b <- spec$blocks[[j]]
    if (b$type == "smooth") {
      L <- c(L, sprintf("  for (l in 1:L_%d) { s_%d[l] ~ dnorm(0, 1 / pow(sd_%d, 2)) }", j, j, j),
             "  " %+% prior_bugs(find_prior(priors, "sd", b$par, b$label), sprintf("sd_%d", j), TRUE, dialect))
    } else if (b$corr) {
      L <- c(L, sprintf("  for (l in 1:L_%d) { r_%d[l, 1:q_%d] ~ dmnorm(zero_%d[], Omega_%d[, ]) }", j, j, j, j, j),
             sprintf("  Omega_%d[1:q_%d, 1:q_%d] ~ dwish(I_%d[, ], q_%d + 1)", j, j, j, j, j),
             sprintf("  Sigma_%d[1:q_%d, 1:q_%d] <- inverse(Omega_%d[, ])", j, j, j, j))
    } else {
      L <- c(L, sprintf("  for (e in 1:q_%d) {", j),
             sprintf("    for (l in 1:L_%d) { r_%d[l, e] ~ dnorm(0, 1 / pow(sd_%d[e], 2)) }", j, j, j),
             "    " %+% prior_bugs(find_prior(priors, "sd", b$par, b$label), sprintf("sd_%d[e]", j), TRUE, dialect),
             "  }")
    }
  }
  paste(c(L, "}"), collapse = "\n")
}

bugs_data <- function(spec) {
  out <- stan_data(spec)
  keep <- c("N", "y", "w", grep("^(X_|off_|fix_|Xr_|g_|L_|q_|Zs_)", names(out), value = TRUE),
            names(spec$data_pars))
  out <- out[intersect(names(out), keep)]
  if (!is.null(spec$fam$bugs$loglik) && (is.null(spec$fam$bugs$native) || any(spec$weights != 1)))
    out$zeros <- rep(0, spec$n)
  else out$w <- NULL
  for (j in seq_along(spec$blocks)) {
    b <- spec$blocks[[j]]
    if (b$type != "smooth" && b$corr) {
      out[[paste0("zero_", j)]] <- rep(0, b$q)
      out[[paste0("I_", j)]] <- diag(b$q)
    }
  }
  out
}
