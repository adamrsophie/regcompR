# Bootstrap, power by simulation, the Vuong test and best subset selection.

# Refit a model to new data (same call otherwise).
refit_data <- function(object, data) {
  cl <- object$call; env <- environment(object$spec$formula) %||% parent.frame()
  cl$data <- data; cl$subset <- NULL
  if (!is.null(object$model_args$weights)) cl$weights <- object$model_args$weights
  ctl <- if (!is.null(cl$control)) eval(cl$control, env) else regcomp_control()
  ctl$hessian <- FALSE; ctl$restarts <- 0
  if (!is.null(object$spec$nl)) ctl$start <- lapply(stats::setNames(nm = object$spec$nl$pars),
    function(k) unname(object$coefficients[object$spec$beta_idx[[k]]]))
  cl$control <- ctl
  tryCatch(suppressWarnings(eval(cl, env)), error = function(e) NULL)
}

# The analysis data with the response replaced by simulated values.
with_response <- function(object, data, y) {
  lhs <- object$spec$formula[[2]]
  if (!is.name(lhs)) stop_unsupported("Simulation based refits when the response is an expression (store it as a column)")
  nm <- as.character(lhs)
  if (is.factor(data[[nm]])) y <- factor(levels(data[[nm]])[y], levels = levels(data[[nm]]))
  data[[nm]] <- y
  data
}

#' Bootstrap a fitted model
#'
#' The parametric bootstrap simulates responses from the fitted model (with
#' new random effects, as `lme4::bootMer()` by default) and refits; the
#' case bootstrap resamples rows, or whole clusters with `cluster`. Returns
#' the replicates of `statistic` with bootstrap standard errors and
#' percentile intervals.
#'
#' @param object A fitted model (frequentist).
#' @param R Number of bootstrap replicates.
#' @param type `"parametric"` or `"case"`.
#' @param statistic Function of a fit returning a numeric vector; default
#'   the coefficients.
#' @param cluster For the case bootstrap, a one sided formula or vector
#'   naming the clusters to resample.
#' @param level Interval level.
#' @param seed Random seed.
#' @return An object of class `regcomp_boot` with `t0`, the replicate
#'   matrix `t` and a `table` of estimates, bootstrap standard errors and
#'   percentile intervals. Failed refits are dropped and counted.
#' @export
#' @examples
#' d <- data.frame(x = rnorm(60)); d$y <- rpois(60, exp(0.5 + 0.3 * d$x))
#' f <- regcomp(y ~ x, family = "poisson", data = d)
#' boot_regcomp(f, R = 20)
boot_regcomp <- function(object, R = 199, type = c("parametric", "case"), statistic = stats::coef,
                         cluster = NULL, level = 0.95, seed = 1) {
  type <- match.arg(type)
  if (object$inference != "frequentist") stop_unsupported("Bootstrapping a Bayesian fit")
  local_seed(seed)
  d <- object$model_data
  t0 <- statistic(object)
  if (type == "case" && !is.null(cluster)) {
    cl <- if (inherits(cluster, "formula")) eval(cluster[[2]], d, environment(cluster)) else cluster
    if (length(cl) != nrow(d)) stop_input("cluster must have one value per analysis row.")
    groups <- split(seq_len(nrow(d)), as.character(cl))
  }
  sims <- if (type == "parametric") simulate(object, nsim = R, re = "new")
  reps <- lapply(seq_len(R), function(r) {
    nd <- if (type == "parametric") with_response(object, d, sims[, r]) else if (is.null(cluster)) d[sample.int(nrow(d), replace = TRUE), , drop = FALSE] else {
      pick <- sample(names(groups), length(groups), replace = TRUE)
      do.call(rbind, lapply(seq_along(pick), function(i) { x <- d[groups[[pick[i]]], , drop = FALSE]
        if (!is.null(cluster) && inherits(cluster, "formula") && is.name(cluster[[2]])) x[[as.character(cluster[[2]])]] <- paste0(pick[i], ".", i); x }))
    }
    f <- refit_data(object, nd)
    if (is.null(f)) NULL else tryCatch(statistic(f), error = function(e) NULL)
  })
  ok <- !vapply(reps, is.null, TRUE)
  tt <- do.call(rbind, reps[ok]); colnames(tt) <- names(t0)
  a <- (1 - level) / 2
  tab <- data.frame(estimate = t0, std_error = apply(tt, 2, stats::sd),
                    lower = apply(tt, 2, stats::quantile, a), upper = apply(tt, 2, stats::quantile, 1 - a), row.names = names(t0))
  structure(list(t0 = t0, t = tt, table = tab, R = sum(ok), failed = sum(!ok), type = type), class = "regcomp_boot")
}

#' @export
print.regcomp_boot <- function(x, digits = 4, ...) {
  cat(sprintf("<regcomp_boot> %s bootstrap, %d replicates%s\n", x$type, x$R,
              if (x$failed) sprintf(" (%d failed refits dropped)", x$failed) else ""))
  print(round(x$table, digits))
  invisible(x)
}

#' Power by simulation
#'
#' Simulates responses from the fitted model, optionally with chosen
#' coefficient values (the effect size to detect) and a new design
#' (`newdata`, for another sample size), refits, and reports the share of
#' simulations whose Wald test of `term` rejects at `alpha`, with a
#' binomial confidence interval. This is the approach of the simr package.
#'
#' @param object A fitted model (frequentist).
#' @param term Coefficient name as in `coef(object)` (for example
#'   `"mu:x"`).
#' @param nsim Number of simulations.
#' @param alpha Test level.
#' @param coefficients Optional named vector of coefficient values to
#'   simulate from (other coefficients keep their estimates).
#' @param newdata Optional design to simulate on (must contain every
#'   covariate and the response column, whose values are replaced).
#' @param seed Random seed.
#' @return A one row data frame with power, its 95% interval, the number of
#'   successful refits and of failures.
#' @export
#' @examples
#' d <- data.frame(x = rnorm(50)); d$y <- rpois(50, exp(0.5 + 0.3 * d$x))
#' f <- regcomp(y ~ x, family = "poisson", data = d)
#' power_sim(f, "lambda:x", nsim = 20, coefficients = c("lambda:x" = 0.2))
power_sim <- function(object, term, nsim = 100, alpha = 0.05, coefficients = NULL, newdata = NULL, seed = 1) {
  if (object$inference != "frequentist") stop_unsupported("Power simulation from a Bayesian fit")
  if (!term %in% names(object$coefficients)) stop_input(sprintf("Unknown coefficient '%s'.", term))
  local_seed(seed)
  sim_obj <- object
  if (length(coefficients)) {
    bad <- setdiff(names(coefficients), names(object$coefficients))
    if (length(bad)) stop_input(sprintf("Unknown coefficient%s: %s.", if (length(bad) > 1) "s" else "", paste(bad, collapse = ", ")))
    sim_obj$coefficients[names(coefficients)] <- coefficients
  }
  d <- if (is.null(newdata)) object$model_data else as.data.frame(newdata)
  sims <- simulate(sim_obj, nsim = nsim, newdata = newdata, re = "new")
  p <- vapply(seq_len(nsim), function(s) {
    f <- refit_data(object, with_response(object, d, sims[, s]))
    if (is.null(f)) return(NA_real_)
    v <- tryCatch(stats::vcov(object = update_cov(f))[term, term], error = function(e) NA_real_)
    if (!is.finite(v) || v <= 0) return(NA_real_)
    2 * stats::pnorm(-abs(f$coefficients[[term]] / sqrt(v)))
  }, 0)
  ok <- is.finite(p); k <- sum(p[ok] < alpha); n <- sum(ok)
  ci <- if (n) stats::binom.test(k, n)$conf.int else c(NA, NA)
  data.frame(term = term, power = if (n) k / n else NA_real_, lower = ci[1], upper = ci[2], nsim = n, failed = sum(!ok))
}

# The covariance of a refit made with hessian = FALSE.
update_cov <- function(f) {
  if (!is.null(f$vcov)) return(f)
  spec <- f$spec
  if (spec$Q > 0) {
    np <- length(f$coefficients); st <- new.env(); st$u <- f$u
    ctl <- regcomp_control()
    ob <- function(p) laplace_objective(spec, p[seq_len(np)], p[-seq_len(np)], st, ctl)
    H <- num_hessian(ob, c(unname(f$coefficients), unname(f$theta)), s = c(coef_scale(spec, f$coefficients), rep(1, length(f$theta))))
    V <- safe_inverse(-H, c(coef_scale(spec, f$coefficients), rep(1, length(f$theta))))
  } else {
    bs <- coef_scale(spec, f$coefficients)
    H <- num_jacobian(function(b) ml_gradient(spec, b), unname(f$coefficients), bs)
    V <- safe_inverse(-(H + t(H)) / 2, bs)
  }
  nn <- c(spec$beta_names, names(f$theta))[seq_len(nrow(V))]
  dimnames(V) <- list(nn, nn); f$vcov <- V
  f
}

#' Vuong test of non nested models
#'
#' Compares two models fitted to the same observations through their
#' pointwise log likelihood ratios `m_i`: `z = sqrt(n) mean(m) / sd(m)`,
#' with the AIC and BIC corrected versions that subtract `(k1 - k2)` and
#' `(k1 - k2) log(n) / 2` from the summed ratio (Vuong 1989). Positive
#' values favor the first model. The raw statistic equals `pscl::vuong()`;
#' the corrections count every estimated parameter (pscl counts only the
#' regression coefficients).
#'
#' @param object1,object2 Fitted models on the same rows.
#' @return A data frame with the three statistics and one sided p values.
#' @export
#' @examples
#' d <- data.frame(x = rnorm(200)); d$y <- rnbinom(200, mu = exp(1 + 0.3 * d$x), size = 2)
#' vuong_test(regcomp(y ~ x, family = "negative_binomial_nb2", data = d),
#'            regcomp(y ~ x, family = "poisson", data = d))
vuong_test <- function(object1, object2) {
  if (object1$inference != "frequentist" || object2$inference != "frequentist") stop_unsupported("The Vuong test for a Bayesian fit")
  same_response(list(object1, object2), "The Vuong test")
  if (object1$spec$Q > 0 || object2$spec$Q > 0)
    stop_unsupported("The Vuong test for models with random effects or smooths (pointwise log likelihoods are conditional)")
  m <- log_lik(object1) - log_lik(object2)
  n <- length(m)
  dk <- object1$df - object2$df
  adj <- c(Raw = 0, `AIC-corrected` = dk, `BIC-corrected` = dk * log(n) / 2)
  z <- (sum(m) - adj) / (sqrt(n) * stats::sd(m))
  data.frame(statistic = z, p_value = ifelse(z > 0, stats::pnorm(z, lower.tail = FALSE), stats::pnorm(z)),
             preferred = ifelse(z > 0, "model 1", "model 2"), row.names = names(adj))
}

#' Best subset selection
#'
#' Fits every subset of the fixed effect terms of one parameter (keeping
#' random effects) and ranks them by AIC, AICc, BIC or the Hannan-Quinn
#' criterion. Practical for up to about 15 candidate terms.
#'
#' @param object A fitted model.
#' @param scope Optional one sided formula of candidate terms (default the
#'   terms of the fit).
#' @param criterion `"AIC"`, `"AICc"`, `"BIC"` or `"HQ"`.
#' @param par Parameter whose formula is searched (default the first).
#' @param top Number of best models to return.
#' @return A data frame of the best models (terms, criterion, df), with
#'   the best fit as attribute `"fit"`.
#' @export
#' @examples
#' d <- data.frame(x1 = rnorm(80), x2 = rnorm(80), x3 = rnorm(80))
#' d$y <- rpois(80, exp(0.5 + 0.4 * d$x1))
#' best_subset(regcomp(y ~ x1 + x2 + x3, family = "poisson", data = d), criterion = "BIC")
best_subset <- function(object, scope = NULL, criterion = c("AIC", "AICc", "BIC", "HQ"), par = NULL, top = 5) {
  no_nl(object, "Best subset selection")
  criterion <- match.arg(criterion)
  if (object$inference != "frequentist") stop_unsupported("Best subset selection of a Bayesian fit")
  par <- par %||% object$spec$rpars[1]
  universe <- if (is.null(scope)) fixed_terms_of(object, par) else attr(stats::terms(scope), "term.labels")
  if (length(universe) > 15) stop_input("Best subset selection is limited to 15 candidate terms; use step_regcomp().")
  subsets <- unlist(lapply(0:length(universe), function(k) utils::combn(universe, k, simplify = FALSE)), recursive = FALSE)
  fits <- lapply(subsets, function(tt) refit_terms(object, tt, par))
  v <- vapply(fits, function(f) if (is.null(f)) Inf else info_criterion(f, criterion), 0)
  o <- order(v)[seq_len(min(top, length(v)))]
  out <- data.frame(terms = vapply(subsets[o], function(tt) if (length(tt)) paste(tt, collapse = " + ") else "(none)", ""),
                    criterion = v[o], df = vapply(fits[o], function(f) if (is.null(f)) NA_real_ else f$df, 0))
  names(out)[2] <- criterion
  attr(out, "fit") <- fits[[o[1]]]
  out
}

# Information criteria; AICc = AIC + 2k(k + 1)/(n - k - 1) and HQ = -2 logLik + 2 k log(log(n)).
info_criterion <- function(f, criterion) {
  k <- f$df; n <- f$nobs; ll <- f$logLik
  switch(criterion, AIC = -2 * ll + 2 * k, BIC = -2 * ll + log(n) * k,
         AICc = -2 * ll + 2 * k + 2 * k * (k + 1) / max(n - k - 1, 1e-8), HQ = -2 * ll + 2 * k * log(log(n)))
}

#' Skovgaard's modified directed likelihood ratio
#'
#' Higher order inference on one coefficient: the signed root of the
#' likelihood ratio statistic `r` and Skovgaard's (1996) modification
#' `r* = r + log(u / r) / r + log(C) / r`, whose normal approximation has
#' relative error of third order. The sample space derivatives are
#' approximated by covariances of the score and of the log likelihood
#' ratio over `R` responses simulated from the fit (Skovgaard 2001), as in
#' `likelihoodAsy::rstar()`. For models without random effects.
#'
#' @param object A fitted model (frequentist, no random effects).
#' @param term Coefficient name as in `coef(object)`.
#' @param value Value of the coefficient under the null hypothesis.
#' @param R Number of simulated samples.
#' @param seed Random seed.
#' @return A data frame with the Wald statistic, `r` and `r*` and their two
#'   sided normal p values.
#' @references Skovgaard (1996), An explicit large-deviation approximation
#'   to one-parameter tests. Bernoulli, 2, 145-165. Skovgaard (2001),
#'   Likelihood asymptotics. Scand. J. Statist., 28, 3-32.
#' @export
#' @examples
#' d <- data.frame(x = rnorm(15)); d$y <- rpois(15, exp(0.5 + 0.4 * d$x))
#' f <- regcomp(y ~ x, family = "poisson", data = d)
#' skovgaard_test(f, "lambda:x", R = 500)
skovgaard_test <- function(object, term, value = 0, R = 1000, seed = 1) {
  if (object$inference != "frequentist") stop_unsupported("Skovgaard's statistic for a Bayesian fit")
  spec <- object$spec
  if (spec$Q > 0 || length(spec$aux)) stop_unsupported("Skovgaard's statistic for models with random effects, smooths or thresholds")
  k <- match(term, names(object$coefficients))
  if (is.na(k)) stop_input(sprintf("Unknown coefficient '%s'.", term))
  bh <- unname(object$coefficients); p <- length(bh); bs <- coef_scale(spec, object$coefficients)
  info <- function(sp, b) { H <- num_jacobian(function(x) ml_gradient(sp, x), b, bs); -(H + t(H)) / 2 }
  # Constrained fit with the coefficient held at value.
  b0 <- bh; b0[k] <- value
  full <- function(z) { b <- b0; b[-k] <- z; b }
  o <- stats::nlminb(bh[-k], function(z) -ml_objective(spec, full(z)), function(z) -ml_gradient(spec, full(z))[-k],
                     scale = 1 / bs[-k], control = list(eval.max = 1000, iter.max = 1000, rel.tol = 1e-12))
  bt <- full(o$par)
  r <- sign(bh[k] - value) * sqrt(max(2 * (ml_objective(spec, bh) - ml_objective(spec, bt)), 0))
  jh <- info(spec, bh); jt <- info(spec, bt)
  local_seed(seed)
  sims <- simulate(object, nsim = R)
  U <- t(vapply(seq_len(R), function(s) {
    sp <- spec; sp$y <- if (is.matrix(spec$y)) matrix(sims[, s], nrow(spec$y)) else sims[, s]
    c(ml_gradient(sp, bh), ml_gradient(sp, bt), ml_objective(sp, bh) - ml_objective(sp, bt))
  }, numeric(2 * p + 1)))
  C <- stats::cov(U)
  S <- C[1:p, p + 1:p, drop = FALSE]; ih <- C[1:p, 1:p, drop = FALSE]; q <- C[2 * p + 1, 1:p]
  iinv <- qr.solve(ih, tol = 1e-20)
  SS2 <- t(S) %*% iinv %*% jh; SS1 <- drop(q %*% iinv %*% jh)
  num <- SS1[k] - if (p > 1) drop(SS2[-k, k] %*% solve(t(SS2[-k, -k, drop = FALSE]), SS1[-k])) else 0
  u <- num * sqrt(solve(jh)[k, k])
  Cp <- if (p > 1) det(SS2[-k, -k, drop = FALSE]) / sqrt(det(jt[-k, -k, drop = FALSE]) * det(jh[-k, -k, drop = FALSE])) else 1
  rs <- r + log(Cp) / r + log(u / r) / r
  w <- (bh[k] - value) * sqrt(1 / solve(jh)[k, k])
  st <- c(Wald = w, r = r, `r*` = rs)
  data.frame(statistic = st, p_value = 2 * stats::pnorm(-abs(st)), row.names = names(st))
}
