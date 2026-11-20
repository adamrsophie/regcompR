# Likelihood ratio and Wald tests, simple slopes and variance inflation.

same_response <- function(fits, what) {
  if (any(vapply(fits, function(f) identical(f$inference, "bayesian"), TRUE)))
    stop_unsupported(sprintf("%s for Bayesian fits", what))
  y1 <- fits[[1]]$spec$y
  for (f in fits[-1]) if (NROW(f$spec$y) != NROW(y1) || !isTRUE(all.equal(unclass(f$spec$y), unclass(y1), check.attributes = FALSE)))
    stop_input(sprintf("%s needs fits to the same observations of the same response.", what))
}

#' Likelihood ratio tests of nested models
#'
#' Compares fits of the same response on the same rows, ordered by their
#' number of parameters, each against the previous one. The models must be
#' nested; that is not checked. When the smaller model sets a variance
#' component to zero the null value is on the boundary and the chi square
#' p value is conservative (about twice the correct one for one variance;
#' Self and Liang 1987).
#'
#' @param object,... Fitted models.
#' @return A data frame of class `anova` with degrees of freedom, log
#'   likelihood, AIC, the likelihood ratio statistic, its degrees of
#'   freedom and p value.
#' @export
#' @examples
#' d <- data.frame(x = rnorm(200), z = rnorm(200)); d$y <- rpois(200, exp(0.5 + 0.3 * d$x))
#' f1 <- regcomp(y ~ x, family = "poisson", data = d)
#' f2 <- regcomp(y ~ x + z, family = "poisson", data = d)
#' anova(f1, f2)
anova.regcomp_fit <- function(object, ...) {
  fits <- c(list(object), list(...))
  if (length(fits) < 2) stop_input("anova() of regcompR fits compares two or more nested models.")
  same_response(fits, "A likelihood ratio test")
  nm <- vapply(as.list(substitute(list(object, ...)))[-1], deparse1, "")
  o <- order(vapply(fits, `[[`, 0, "df")); fits <- fits[o]
  ll <- vapply(fits, `[[`, 0, "logLik"); df <- vapply(fits, `[[`, 0, "df")
  lr <- c(NA, 2 * diff(ll)); ddf <- c(NA, diff(df))
  out <- data.frame(df = df, logLik = ll, AIC = -2 * ll + 2 * df, Chisq = lr, `Chi Df` = ddf,
                    `Pr(>Chisq)` = stats::pchisq(lr, ddf, lower.tail = FALSE), check.names = FALSE, row.names = nm[o])
  structure(out, heading = "Likelihood ratio tests\n", class = c("anova", "data.frame"))
}

#' Wald tests of linear hypotheses
#'
#' Tests `L beta = rhs` with the statistic
#' `(L b - rhs)' (L V L')^-1 (L b - rhs)`, chi square with `nrow(L)`
#' degrees of freedom, where `V` is the chosen covariance.
#'
#' @param object A frequentist fit.
#' @param L The hypothesis: a character vector of coefficient names (all
#'   zero), a named numeric vector (one linear combination, for example
#'   `c("mu:x" = 1, "mu:z" = -1)` for equal coefficients), or a matrix with
#'   one row per combination and coefficient names as column names.
#' @param rhs Right hand side, recycled to the number of rows of `L`.
#' @param vcov `"model"`, `"sandwich"` or `"quasi"` (see
#'   [vcov.regcomp_fit()]).
#' @param ... Passed to [vcov.regcomp_fit()] (such as `cluster`).
#' @return A data frame with the estimates of the combinations, the
#'   statistic, its degrees of freedom and the p value.
#' @export
#' @examples
#' d <- data.frame(x = rnorm(200), z = rnorm(200)); d$y <- rpois(200, exp(0.5 + 0.3 * d$x + 0.3 * d$z))
#' f <- regcomp(y ~ x + z, family = "poisson", data = d)
#' wald_test(f, c("lambda:x" = 1, "lambda:z" = -1))
#' wald_test(f, c("lambda:x", "lambda:z"))
wald_test <- function(object, L, rhs = 0, vcov = "model", ...) {
  if (object$inference != "frequentist") stop_unsupported("A Wald test for a Bayesian fit (use the posterior draws)")
  b <- object$coefficients; bn <- names(b)
  if (is.character(L)) L <- matrix(diag(1, length(L)), length(L), dimnames = list(NULL, L))
  else if (is.null(dim(L))) L <- matrix(L, 1, dimnames = list(NULL, names(L)))
  if (is.null(colnames(L)) || !all(colnames(L) %in% bn))
    stop_input(sprintf("L must name coefficients of the fit: %s.", paste(setdiff(colnames(L), bn), collapse = ", ")))
  M <- matrix(0, nrow(L), length(bn), dimnames = list(NULL, bn)); M[, colnames(L)] <- L
  V <- stats::vcov(object, type = vcov, ...)
  est <- drop(M %*% b) - rep_len(rhs, nrow(M))
  W <- drop(crossprod(est, solve(M %*% V %*% t(M), est)))
  data.frame(estimate = paste(signif(est, 4), collapse = ", "), chisq = W, df = nrow(M),
             p_value = stats::pchisq(W, nrow(M), lower.tail = FALSE))
}

#' Simple slopes under moderation
#'
#' The slope of a numeric covariate `x` on the linear predictor of one
#' distribution parameter at chosen values of a moderator, averaged over
#' the analysis rows with the moderator set to each value (the average
#' marginal effect on the link scale), with delta method standard errors
#' and Wald tests. With a linear interaction `x + z + x:z` this is
#' `b_x + b_xz z`. Random effects do not enter the slope.
#'
#' @param object A frequentist fit.
#' @param x Name of the focal numeric covariate.
#' @param moderator Name of the moderator.
#' @param at Values of the moderator; by default its mean and the mean plus
#'   and minus one standard deviation (its levels for a factor).
#' @param par Distribution parameter; the first regression parameter by
#'   default.
#' @param vcov,... Covariance choice passed to [vcov.regcomp_fit()].
#' @return A data frame with one row per moderator value.
#' @export
#' @examples
#' d <- data.frame(x = rnorm(300), z = rnorm(300))
#' d$y <- rpois(300, exp(0.5 + 0.3 * d$x + 0.2 * d$x * d$z))
#' simple_slopes(regcomp(y ~ x * z, family = "poisson", data = d), "x", "z")
simple_slopes <- function(object, x, moderator, at = NULL, par = NULL, vcov = "model", ...) {
  no_nl(object, "Simple slopes")
  if (object$inference != "frequentist") stop_unsupported("Simple slopes for a Bayesian fit")
  d <- object$model_data
  if (is.null(d) || !is.numeric(d[[x]])) stop_input("x must be a numeric covariate of the analysis data.")
  par <- par %||% object$spec$rpars[1]
  z <- d[[moderator]]
  at <- at %||% if (is.numeric(z)) mean(z) + c(-1, 0, 1) * stats::sd(z) else levels(factor(z))
  V <- stats::vcov(object, type = vcov, ...); idx <- object$spec$beta_idx[[par]]
  h <- 1e-4 * max(1, stats::sd(d[[x]]))
  X <- function(nd) colMeans(new_spec(object, nd)$preds[[par]]$X)
  rows <- lapply(at, function(a) {
    nd <- d; nd[[moderator]] <- if (is.factor(z)) factor(a, levels(z)) else a
    up <- nd; up[[x]] <- up[[x]] + h
    cvec <- (X(up) - X(nd)) / h
    est <- sum(cvec * object$coefficients[idx]); se <- sqrt(drop(t(cvec) %*% V[idx, idx] %*% cvec))
    data.frame(moderator = a, slope = est, std_error = se, statistic = est / se, p_value = 2 * stats::pnorm(-abs(est / se)))
  })
  out <- do.call(rbind, rows); names(out)[1] <- moderator
  out
}

#' Variance inflation factors
#'
#' Generalized variance inflation factors (Fox and Monette 1992) of the
#' terms of each distribution parameter's regression, from the correlations
#' of the coefficient estimates other than the intercept (as in
#' `car::vif()`). Values above
#' about 5 or 10 indicate strong collinearity. `gvif_adj` is
#' `GVIF^(1 / (2 df))`, comparable across terms with several columns.
#'
#' @param object A fitted model.
#' @return A data frame with one row per term and parameter.
#' @export
#' @examples
#' d <- data.frame(x = rnorm(100)); d$z <- d$x + rnorm(100, 0, 0.3); d$y <- rnorm(100)
#' vif(regcomp(y ~ x + z, family = "normal", data = d))
vif <- function(object) {
  no_nl(object, "Variance inflation factors")
  V <- stats::vcov(object)
  out <- lapply(object$spec$rpars, function(k) {
    X <- object$spec$preds[[k]]$X; asg <- attr(X, "assign")
    keep <- colnames(X) != "(Intercept)"
    if (sum(keep) < 2) return(NULL)
    idx <- object$spec$beta_idx[[k]][keep]
    R <- stats::cov2cor(V[idx, idx, drop = FALSE])
    tl <- attr(stats::terms(object$spec$formulas[[k]]), "term.labels"); tl <- tl[!grepl("\\|", tl)]
    a <- if (is.null(asg) || max(asg) > length(tl)) seq_len(sum(keep)) else asg[keep]
    lab <- if (identical(a, seq_len(sum(keep)))) colnames(X)[keep] else tl
    do.call(rbind, lapply(unique(a), function(t) { j <- a == t
      g <- det(R[j, j, drop = FALSE]) * det(R[!j, !j, drop = FALSE]) / det(R)
      data.frame(parameter = k, term = lab[t], df = sum(j), gvif = g, gvif_adj = g^(1 / (2 * sum(j)))) }))
  })
  do.call(rbind, out)
}

#' Johnson-Neyman regions of significance
#'
#' The values of a numeric moderator at which the simple slope of `x` (see
#' [simple_slopes()]) is significantly different from zero (Johnson and
#' Neyman 1936; Bauer and Curran 2005). Slopes are computed on a grid over
#' the observed range of the moderator and the bounds are found by linear
#' interpolation of the z statistic. The level is per moderator value, not
#' adjusted for looking at many values.
#'
#' @inheritParams simple_slopes
#' @param level Significance level.
#' @param n Grid size.
#' @param plot Draw the slope with its confidence band against the
#'   moderator, shading the region of significance.
#' @return A list with the `bounds` (moderator values where significance
#'   changes), the `range` of the moderator and the `grid` of slopes,
#'   invisibly when plotted.
#' @export
#' @examples
#' d <- data.frame(x = rnorm(300), z = rnorm(300))
#' d$y <- rpois(300, exp(0.5 + 0.3 * d$x + 0.2 * d$x * d$z))
#' johnson_neyman(regcomp(y ~ x * z, family = "poisson", data = d), "x", "z")
johnson_neyman <- function(object, x, moderator, par = NULL, level = 0.05, n = 101, vcov = "model", plot = FALSE, ...) {
  no_nl(object, "Johnson-Neyman intervals")
  z <- object$model_data[[moderator]]
  if (!is.numeric(z)) stop_input("The moderator must be numeric.")
  g <- simple_slopes(object, x, moderator, at = seq(min(z), max(z), length.out = n), par = par, vcov = vcov, ...)
  crit <- stats::qnorm(1 - level / 2); s <- abs(g$statistic) - crit; zz <- g[[moderator]]
  i <- which(diff(sign(s)) != 0)
  bounds <- zz[i] - s[i] * (zz[i + 1] - zz[i]) / (s[i + 1] - s[i])
  g$significant <- s > 0
  out <- list(bounds = bounds, range = range(z), grid = g)
  if (!plot) return(out)
  lo <- g$slope - crit * g$std_error; hi <- g$slope + crit * g$std_error
  graphics::plot(zz, g$slope, type = "n", ylim = range(lo, hi), xlab = moderator, ylab = paste("slope of", x),
                 main = "Johnson-Neyman plot")
  graphics::polygon(c(zz, rev(zz)), c(lo, rev(hi)), col = grDevices::adjustcolor("#1b6ca8", 0.2), border = NA)
  graphics::lines(zz, g$slope, lwd = 2, col = "#1b6ca8"); graphics::abline(h = 0, lty = 2)
  graphics::rug(zz[g$significant], col = "#1b6ca8", lwd = 2)
  if (length(bounds)) graphics::abline(v = bounds, lty = 3)
  invisible(out)
}

#' Stepwise selection of terms
#'
#' Forward, backward or bidirectional selection of the terms of one
#' distribution parameter's formula by AIC or BIC, as in `stats::step()`
#' and StepReg. Selection makes the reported standard errors, intervals and
#' p values of the final model too optimistic (Harrell 2015, sec. 4.3);
#' prefer penalized estimation (see [penalty()]) or cross validation for
#' prediction, and treat stepwise results as exploratory.
#'
#' @param object A frequentist fit: the starting model.
#' @param scope A one sided formula with the largest set of terms considered
#'   for `par`; by default the terms of the starting model (backward only).
#' @param direction `"both"`, `"backward"` or `"forward"`.
#' @param criterion `"AIC"` or `"BIC"`.
#' @param par Distribution parameter whose formula is searched; the main
#'   one by default. Random effect terms are kept.
#' @param trace Print each step.
#' @return The selected fit, with the path of steps in `attr(, "steps")`.
#' @export
#' @examples
#' d <- data.frame(x1 = rnorm(200), x2 = rnorm(200), x3 = rnorm(200))
#' d$y <- rpois(200, exp(0.5 + 0.4 * d$x1))
#' f0 <- regcomp(y ~ 1, family = "poisson", data = d)
#' fs <- step_regcomp(f0, scope = ~ x1 + x2 + x3, direction = "forward", criterion = "BIC")
#' formula(fs)
#' attr(fs, "steps")
step_regcomp <- function(object, scope = NULL, direction = c("both", "backward", "forward"),
                         criterion = c("AIC", "BIC"), par = NULL, trace = FALSE) {
  no_nl(object, "Stepwise selection")
  direction <- match.arg(direction); criterion <- match.arg(criterion)
  if (object$inference != "frequentist") stop_unsupported("Stepwise selection of a Bayesian fit")
  par <- par %||% object$spec$rpars[1]
  ic <- function(f) -2 * f$logLik + (if (criterion == "AIC") 2 else log(f$nobs)) * f$df
  fixed_terms <- function(f) { tl <- attr(stats::terms(f$spec$formulas[[par]]), "term.labels"); tl[!grepl("\\|", tl)] }
  universe <- if (is.null(scope)) fixed_terms(object) else attr(stats::terms(scope), "term.labels")
  refit <- function(f, terms) {
    fo <- f$spec$formulas[[par]]; tl <- attr(stats::terms(fo), "term.labels"); re <- tl[grepl("\\|", tl)]
    rhs <- c(if (attr(stats::terms(fo), "intercept") == 0) "0", terms, if (length(re)) paste0("(", re, ")"))
    nf <- stats::reformulate(if (length(rhs)) rhs else "1", response = if (length(fo) == 3) fo[[2]] else NULL, env = environment(fo))
    cl <- f$call; i <- which(vapply(as.list(cl), function(e) tryCatch(identical(eval(e, environment(fo)), fo) ||
      (inherits(e, "formula") && identical(deparse(e), deparse(fo))), error = function(e) FALSE), TRUE))
    if (!length(i)) stop_unsupported(sprintf("Stepwise selection when the formula of '%s' is not written in the call", par))
    cl[[i[1]]] <- nf
    tryCatch(eval(cl, environment(fo)), error = function(e) NULL)
  }
  cur <- object; cur_ic <- ic(cur); steps <- data.frame(step = "start", criterion = cur_ic)
  repeat {
    have <- fixed_terms(cur); cand <- list()
    if (direction != "forward") for (t in have) cand[[paste("-", t)]] <- setdiff(have, t)
    if (direction != "backward") for (t in setdiff(universe, have)) cand[[paste("+", t)]] <- c(have, t)
    if (!length(cand)) break
    fits <- lapply(cand, function(tt) refit(cur, tt))
    v <- vapply(fits, function(f) if (is.null(f)) Inf else ic(f), 0)
    if (min(v) >= cur_ic - 1e-8) break
    b <- which.min(v); cur <- fits[[b]]; cur_ic <- v[b]
    steps <- rbind(steps, data.frame(step = names(cand)[b], criterion = cur_ic))
    if (trace) cat(names(cand)[b], criterion, "=", format(cur_ic), "\n")
  }
  names(steps)[2] <- criterion
  attr(cur, "steps") <- steps
  cur
}
