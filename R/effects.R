# Average marginal effects, average predictions, R squared and ICC.

# A function of the coefficients (and random effects) giving one predicted
# quantity per row of nd.
pred_fn <- function(object, nd, type, par, re) {
  spec <- new_spec(object, nd)
  function(beta, u) {
    if (type == "link") return(pred_eta(spec, beta, u, re)[[par]])
    p <- pars_from(spec, beta, u, re)
    if (type == "parameter") return(p[[par]])
    quantity_from_pars(spec$fam, p, "mean", NULL, 0)
  }
}

# Estimates of q(beta, u) (a vector) with delta method (frequentist) or
# posterior (Bayesian) standard errors and intervals.
effect_inference <- function(object, q, level, vcov = "model", ...) {
  a <- (1 - level) / 2
  if (object$inference == "bayesian") {
    bn <- object$spec$beta_names
    D <- t(vapply(seq_len(nrow(object$draws)), function(s)
      q(object$draws[s, bn], if (length(object$u_draws)) object$u_draws[s, ] else NULL), q(object$draws[1, bn], NULL)))
    D <- matrix(D, nrow = nrow(object$draws))
    return(data.frame(estimate = colMeans(D), std_error = apply(D, 2, stats::sd), statistic = NA_real_, p_value = NA_real_,
                      lower = apply(D, 2, stats::quantile, a), upper = apply(D, 2, stats::quantile, 1 - a)))
  }
  b <- object$coefficients; u <- object$u
  est <- q(b, u)
  J <- vapply(seq_along(b), function(j) { h <- 1e-5 * max(1, abs(b[j])); bp <- b; bm <- b
    bp[j] <- b[j] + h; bm[j] <- b[j] - h; (q(bp, u) - q(bm, u)) / (2 * h) }, est)
  J <- matrix(J, nrow = length(est))
  se <- sqrt(pmax(rowSums((J %*% stats::vcov(object, type = vcov, ...)) * J), 0))
  z <- stats::qnorm(1 - a)
  data.frame(estimate = est, std_error = se, statistic = est / se, p_value = 2 * stats::pnorm(-abs(est / se)),
             lower = est - z * se, upper = est + z * se)
}

effect_data <- function(object, newdata) {
  d <- newdata %||% object$model_data
  if (is.null(d)) stop_input("Pass newdata: the fit stores no analysis rows.")
  as.data.frame(d)
}

#' Average marginal effects and contrasts
#'
#' For a numeric covariate, the average over the rows of the derivative of
#' the predicted quantity with respect to the covariate (the average
#' marginal effect, or average slope); for a factor, logical or character
#' covariate, the average over the rows of the difference between the
#' predictions with every row set to each level and with every row set to
#' the reference level (g-computation; the average treatment effect under
#' the usual assumptions). These are the `avg_slopes()` and
#' `avg_comparisons()` of marginaleffects and the average marginal effects
#' of Stata's `margins, dydx()`.
#'
#' Frequentist fits use the delta method with the chosen covariance
#' (random effects are held at their estimates); Bayesian fits summarize the
#' effect over the posterior draws.
#'
#' @param object A fitted model.
#' @param variables Covariates; by default every covariate of the model.
#' @param type `"mean"` (the predicted mean), `"parameter"` (a distribution
#'   parameter) or `"link"` (its linear predictor).
#' @param par Distribution parameter for `type = "parameter"` or `"link"`.
#' @param newdata Rows to average over; the analysis rows by default.
#' @param re `"conditional"` (with the estimated random effects) or
#'   `"zero"` (random effects set to zero).
#' @param level Interval level.
#' @param vcov,... Covariance choice for frequentist fits (see
#'   [vcov.regcomp_fit()]).
#' @return A data frame with one row per covariate (or contrast).
#' @export
#' @examples
#' d <- data.frame(x = rnorm(300), g = factor(sample(c("a", "b", "c"), 300, TRUE)))
#' d$y <- rbinom(300, 1, plogis(-0.5 + 0.8 * d$x + (d$g == "b")))
#' f <- regcomp(y ~ x + g, family = "bernoulli", data = d)
#' avg_effects(f)
avg_effects <- function(object, variables = NULL, type = c("mean", "parameter", "link"), par = NULL,
                        newdata = NULL, re = c("conditional", "zero"), level = 0.95, vcov = "model", ...) {
  type <- match.arg(type); re <- match.arg(re)
  par <- par %||% object$spec$rpars[1]
  d <- effect_data(object, newdata)
  if (is.null(variables)) {
    variables <- unique(unlist(lapply(object$spec$rpars, function(k) all.vars(rhs_of(object$spec$formulas[[k]])))))
    variables <- setdiff(intersect(variables, names(d)), all.vars(object$spec$formulas[[1]][[2]]))
    grp <- unlist(lapply(object$spec$blocks, function(b) all.vars(str2lang(b$group))))
    variables <- setdiff(variables, grp)
  }
  out <- lapply(variables, function(v) {
    x <- d[[v]]
    if (is.numeric(x) && length(unique(x)) > 2 || is.numeric(x) && !all(x %in% 0:1)) {
      h <- 1e-4 * max(stats::sd(x), 1e-8)
      up <- d; up[[v]] <- x + h; dn <- d; dn[[v]] <- x - h
      fu <- pred_fn(object, up, type, par, re); fd <- pred_fn(object, dn, type, par, re)
      r <- effect_inference(object, function(b, u) mean((fu(b, u) - fd(b, u)) / (2 * h)), level, vcov, ...)
      return(cbind(term = v, contrast = "dY/dX", r))
    }
    lev <- if (is.factor(x)) levels(x) else sort(unique(x))
    set <- function(l) { nd <- d; nd[[v]] <- if (is.factor(x)) factor(l, levels = lev) else l; nd }
    f0 <- pred_fn(object, set(lev[1]), type, par, re)
    do.call(rbind, lapply(lev[-1], function(l) { f1 <- pred_fn(object, set(l), type, par, re)
      cbind(term = v, contrast = paste(l, "-", lev[1]),
            effect_inference(object, function(b, u) mean(f1(b, u) - f0(b, u)), level, vcov, ...)) }))
  })
  out <- do.call(rbind, out); rownames(out) <- NULL
  out
}

#' Marginal means (average predictions)
#'
#' The prediction averaged over rows, overall or with every row set to each
#' value of a covariate (counterfactual or g-computation marginal means, as
#' `avg_predictions(variables = )` of marginaleffects). Unlike the
#' estimated marginal means of emmeans, which average over a balanced
#' reference grid, these average over the observed rows.
#'
#' @inheritParams avg_effects
#' @param variable Covariate to set, or `NULL` for the overall average.
#' @param at Values of the covariate; its levels (or mean and the mean plus
#'   and minus one standard deviation for a numeric covariate) by default.
#' @return A data frame with one row per value.
#' @export
#' @examples
#' d <- data.frame(x = rnorm(300), g = factor(sample(c("a", "b"), 300, TRUE)))
#' d$y <- rpois(300, exp(0.2 + 0.5 * d$x + 0.4 * (d$g == "b")))
#' marginal_means(regcomp(y ~ x + g, family = "poisson", data = d), "g")
marginal_means <- function(object, variable = NULL, at = NULL, type = c("mean", "parameter", "link"), par = NULL,
                            newdata = NULL, re = c("conditional", "zero"), level = 0.95, vcov = "model", ...) {
  type <- match.arg(type); re <- match.arg(re)
  par <- par %||% object$spec$rpars[1]
  d <- effect_data(object, newdata)
  if (is.null(variable)) {
    f <- pred_fn(object, d, type, par, re)
    return(effect_inference(object, function(b, u) mean(f(b, u)), level, vcov, ...))
  }
  x <- d[[variable]]
  at <- at %||% if (is.numeric(x)) mean(x) + c(-1, 0, 1) * stats::sd(x) else if (is.factor(x)) levels(x) else sort(unique(x))
  out <- do.call(rbind, lapply(at, function(a) { nd <- d; nd[[variable]] <- if (is.factor(x)) factor(a, levels(x)) else a
    f <- pred_fn(object, nd, type, par, re)
    effect_inference(object, function(b, u) mean(f(b, u)), level, vcov, ...) }))
  out <- cbind(stats::setNames(data.frame(at), variable), out); rownames(out) <- NULL
  out
}

# Refit with the main formula reduced to an intercept (keeping its random
# effect terms when keep_re) and every other parameter intercept only.
null_fit <- function(object, keep_re = FALSE) {
  fo <- object$spec$formulas[[object$spec$rpars[1]]]
  tl <- attr(stats::terms(fo), "term.labels"); re <- if (keep_re) tl[grepl("\\|", tl)]
  nf <- stats::reformulate(c("1", if (length(re)) paste0("(", re, ")")), response = fo[[2]], env = environment(fo))
  cl <- as.list(object$call); env <- environment(fo)
  isf <- vapply(seq_along(cl), function(i) i > 1 && (is.null(names(cl)) || names(cl)[i] %in% c("", "formula")) &&
    inherits(tryCatch(eval(cl[[i]], env), error = function(e) NULL), "formula"), TRUE)
  cl <- cl[!isf | seq_along(cl) == which(isf)[1]]; cl[[which(isf)[1]]] <- nf
  cl$penalty <- NULL; cl$fixed <- NULL
  eval(as.call(cl), env)
}

# Distribution specific observation level variance on the link scale
# (Nakagawa, Johnson and Schielzeth 2017).
obs_level_variance <- function(object) {
  fam <- object$family$id; k <- object$spec$rpars[1]; lk <- object$spec$links[[k]]$name
  if (fam %in% c("normal", "gaussian") && lk == "identity") {
    s <- unique(fitted_distribution(object)$pars$sigma); if (length(s) > 1) s <- mean(s^2) else s <- s^2
    return(s)
  }
  if (fam %in% c("bernoulli", "binomial"))
    return(switch(lk, logit = pi^2 / 3, probit = 1, cloglog = pi^2 / 6,
                  stop_unsupported(sprintf("Nakagawa R2 with the %s link", lk))))
  if (fam %in% c("poisson", "negative_binomial_nb2") && lk == "log") {
    # Expected count from the intercept only model with the same random
    # effects, as in MuMIn::r.squaredGLMM and performance.
    n0 <- null_fit(object, keep_re = TRUE)
    lam <- exp(unname(n0$coefficients[1]) + 0.5 * random_variance(n0))
    extra <- if (fam == "poisson") 0 else 1 / mean(fitted_distribution(object)$pars$phi)
    return(log1p(1 / lam + extra))
  }
  stop_unsupported(sprintf("Nakagawa R2 for the '%s' family", fam))
}

# Fixed and random effect variances of the main parameter's linear predictor.
nakagawa_parts <- function(object) {
  if (object$inference != "frequentist") stop_unsupported("R2 and ICC for a Bayesian fit")
  spec <- object$spec; k <- spec$rpars[1]
  var_f <- stats::var(drop(spec$preds[[k]]$X %*% object$coefficients[spec$beta_idx[[k]]]))
  list(fixed = var_f, random = random_variance(object), residual = obs_level_variance(object))
}

# Mean over rows of the random effect variance of the main parameter's
# linear predictor, diag(Z Lambda Z') (Johnson 2014 for random slopes).
random_variance <- function(object) {
  spec <- object$spec; k <- spec$rpars[1]
  if (spec$Q == 0 || is.null(spec$Z[[k]])) return(0)
  covs <- lapply(seq_along(spec$blocks), function(j) { b <- spec$blocks[[j]]
    S <- block_cov(b, object$theta[spec$theta_idx[[j]]])
    if (b$type != "re" || b$par != k) S <- S * 0
    Matrix::kronecker(Matrix::Diagonal(b$L), Matrix::Matrix(S)) })
  Z <- spec$Z[[k]]
  mean(Matrix::rowSums((Z %*% Matrix::bdiag(covs)) * Z))
}

#' R squared and intraclass correlation
#'
#' `r2()` gives likelihood based pseudo R squared measures against the
#' intercept only model (McFadden; Cox and Snell; Nagelkerke) and, for the
#' families where it is defined, the marginal and conditional R squared of
#' Nakagawa, Johnson and Schielzeth (2017): the fixed effect variance, and
#' the fixed plus random effect variance, of the main parameter's linear
#' predictor over their sum with the distribution specific observation level
#' variance (normal: sigma^2; binary and binomial: pi^2 / 3 for logit, 1 for
#' probit, pi^2 / 6 for cloglog; Poisson and NB2 with log link: the
#' lognormal approximation from the intercept only model with the same
#' random effects). `icc()` gives the adjusted and conditional ICC from the
#' same variances.
#'
#' @param object A frequentist fit.
#' @return `r2()`: a named vector; `icc()`: a named vector with the
#'   adjusted and conditional ICC.
#' @export
#' @examples
#' d <- data.frame(x = rnorm(300), g = factor(rep(1:30, each = 10)))
#' d$y <- rpois(300, exp(0.5 + 0.4 * d$x + rnorm(30, 0, 0.5)[d$g]))
#' f <- regcomp(y ~ x + (1 | g), family = "poisson", data = d)
#' r2(f)
#' icc(f)
r2 <- function(object) {
  if (object$inference != "frequentist") stop_unsupported("R2 for a Bayesian fit (use cross validated scores)")
  n0 <- null_fit(object)
  ll <- object$logLik; ll0 <- n0$logLik; n <- object$nobs
  cs <- 1 - exp(2 * (ll0 - ll) / n)
  out <- c(mcfadden = 1 - ll / ll0, cox_snell = cs, nagelkerke = cs / (1 - exp(2 * ll0 / n)))
  nk <- tryCatch(nakagawa_parts(object), regcomp_unsupported = function(e) NULL)
  if (!is.null(nk)) { tot <- nk$fixed + nk$random + nk$residual
    out <- c(out, marginal = nk$fixed / tot, conditional = (nk$fixed + nk$random) / tot) }
  out
}

#' @rdname r2
#' @export
icc <- function(object) {
  nk <- nakagawa_parts(object)
  if (nk$random == 0) stop_input("The main parameter has no random effects.")
  c(adjusted = nk$random / (nk$random + nk$residual), conditional = nk$random / (nk$fixed + nk$random + nk$residual))
}
