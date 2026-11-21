# Uplift modeling: individual treatment effects and Qini evaluation.

#' Individual uplift predictions
#'
#' The predicted mean response with the treatment set to its treated value
#' minus the prediction with it set to the control value, for each row: the
#' individual uplift (conditional average treatment effect) of a model with
#' treatment by covariate interactions.
#'
#' @param object A fitted model.
#' @param treat Name of the treatment variable (0/1, logical or a two level
#'   factor).
#' @param newdata Optional data (default the model data).
#' @param treated,control Values of `treat` for treatment and control
#'   (default 1 and 0, or the second and first factor levels).
#' @return A numeric vector of uplifts.
#' @export
#' @examples
#' set.seed(1); n <- 400
#' d <- data.frame(t = rbinom(n, 1, 0.5), x = rnorm(n))
#' d$y <- rbinom(n, 1, plogis(-0.5 + 0.3 * d$x + d$t * (0.2 + 0.8 * d$x)))
#' f <- regcomp(y ~ t * x, family = "bernoulli", data = d)
#' u <- uplift(f, "t")
#' qini(u, d$t, d$y)
uplift <- function(object, treat, newdata = NULL, treated = NULL, control = NULL) {
  d <- if (is.null(newdata)) object$model_data else as.data.frame(newdata)
  v <- d[[treat]]
  if (is.null(v)) stop_input(sprintf("'%s' is not a column of the data.", treat))
  if (is.factor(v)) {
    treated <- treated %||% levels(v)[2]; control <- control %||% levels(v)[1]
    d1 <- d; d1[[treat]] <- factor(treated, levels = levels(v)); d0 <- d; d0[[treat]] <- factor(control, levels = levels(v))
  } else {
    treated <- treated %||% (if (is.logical(v)) TRUE else 1); control <- control %||% (if (is.logical(v)) FALSE else 0)
    d1 <- d; d1[[treat]] <- treated; d0 <- d; d0[[treat]] <- control
  }
  as.numeric(stats::predict(object, d1, type = "mean")) - as.numeric(stats::predict(object, d0, type = "mean"))
}

#' Qini curve and coefficient
#'
#' Ranks observations by predicted uplift and, for the top `k / groups`
#' share, computes the incremental number of responders
#' `R_T - R_C N_T / N_C` and the observed uplift in each group (Radcliffe
#' 2007). The Qini coefficient is the area between the Qini curve
#' (incremental responders as a share of all treated, by share of treated)
#' and the line of random targeting. Equal to
#' `tools4uplift::PerformanceUplift()` and `QiniArea()`.
#'
#' @param uplift Predicted uplifts (for example from [uplift()]).
#' @param treat Treatment indicator (1 treated, 0 control).
#' @param outcome Binary outcome (1 response).
#' @param groups Number of groups.
#' @return An object of class `regcomp_qini`: a data frame with columns
#'   `cum_per`, `T_Y1`, `T_n`, `C_Y1`, `C_n`, `incremental_Y1`,
#'   `inc_uplift` and `uplift` (percent), with the Qini coefficient as
#'   attribute `"qini"`.
#' @export
qini <- function(uplift, treat, outcome, groups = 10) {
  n <- length(uplift)
  if (length(treat) != n || length(outcome) != n) stop_input("uplift, treat and outcome must have the same length.")
  if (groups < 2) stop_input("groups must be at least 2.")
  treat <- as.numeric(treat); outcome <- as.numeric(outcome)
  # Ranks rounded to two decimals, as tools4uplift does.
  r <- round(rank(-uplift, ties.method = "max") / n, 2)
  g <- integer(n)
  for (i in seq_len(groups)) g[r > (i - 1) / groups & r <= i / groups] <- i
  cum <- t(vapply(seq_len(groups), function(i) { s <- g <= i
    c(sum(treat[s] == 1 & outcome[s] == 1), sum(treat[s] == 1), sum(treat[s] == 0 & outcome[s] == 1), sum(treat[s] == 0)) }, numeric(4)))
  inc <- cum[, 1] - cum[, 3] * cum[, 2] / cum[, 4]
  up <- vapply(seq_len(groups), function(i) { s <- g == i
    mean(outcome[s & treat == 1]) - mean(outcome[s & treat == 0]) }, 0)
  out <- data.frame(cum_per = seq_len(groups) / groups, T_Y1 = cum[, 1], T_n = cum[, 2], C_Y1 = cum[, 3], C_n = cum[, 4],
                    incremental_Y1 = inc, inc_uplift = inc / cum[groups, 2] * 100, uplift = up * 100)
  x <- out$T_n / out$T_n[groups]
  area <- sum((out$inc_uplift + c(0, out$inc_uplift[-groups])) / 2 * diff(c(0, x)))
  structure(out, qini = area - out$inc_uplift[groups] / 2, class = c("regcomp_qini", "data.frame"))
}

#' @export
print.regcomp_qini <- function(x, digits = 4, ...) {
  print(round(as.data.frame(unclass(x))[, names(x)], digits), row.names = FALSE)
  cat(sprintf("Qini coefficient: %s\n", format(attr(x, "qini"), digits = digits)))
  invisible(x)
}

#' @export
plot.regcomp_qini <- function(x, ...) {
  graphics::plot(c(0, x$cum_per), c(0, x$inc_uplift), type = "b", xlab = "Share targeted",
                 ylab = "Incremental responders (% of treated)", main = "Qini curve", ...)
  graphics::lines(c(0, 1), c(0, x$inc_uplift[nrow(x)]), lty = 2)
  invisible(x)
}

#' Causal mediation analysis
#'
#' Average causal mediation effects (ACME), average direct effects (ADE),
#' the total effect and the proportion mediated from a mediator model and
#' an outcome model, by the quasi-Bayesian Monte Carlo method of Imai,
#' Keele and Tingley (2010), as `mediation::mediate()`: coefficients are
#' drawn from their asymptotic normal distribution, mediators are simulated
#' from the mediator model under treatment and control, and outcomes are
#' predicted under each combination. Any families (and distributional
#' models) without random effects.
#'
#' @param med Fitted mediator model (with the treatment as covariate).
#' @param out Fitted outcome model (with treatment and mediator).
#' @param treat Name of the treatment variable.
#' @param mediator Name of the mediator variable (the response of `med`).
#' @param treated,control Treatment values compared.
#' @param sims Number of Monte Carlo draws.
#' @param level Interval level.
#' @param seed Random seed.
#' @return A data frame with estimates, intervals and two sided Monte Carlo
#'   p values for ACME and ADE under control and treatment, their
#'   averages, the total effect and the proportion mediated (the median of
#'   its draws).
#' @references Imai, Keele and Tingley (2010), A general approach to causal
#'   mediation analysis. Psychological Methods, 15, 309-334.
#' @export
#' @examples
#' set.seed(1); n <- 300
#' d <- data.frame(t = rbinom(n, 1, 0.5), x = rnorm(n))
#' d$m <- rnorm(n, 0.5 * d$t + 0.3 * d$x); d$y <- rpois(n, exp(0.2 + 0.3 * d$m + 0.2 * d$t))
#' mm <- regcomp(m ~ t + x, family = "normal", data = d)
#' my <- regcomp(y ~ t + m + x, family = "poisson", data = d)
#' mediate_regcomp(mm, my, "t", "m", sims = 200)
mediate_regcomp <- function(med, out, treat, mediator, treated = 1, control = 0, sims = 1000, level = 0.95, seed = 1) {
  for (f in list(med, out)) {
    if (f$inference != "frequentist") stop_unsupported("Mediation analysis of a Bayesian fit")
    if (f$spec$Q > 0) stop_unsupported("Mediation analysis with random effects")
  }
  local_seed(seed)
  d <- out$model_data
  draw <- function(f) { V <- stats::vcov(f); b <- f$coefficients
    sweep(matrix(stats::rnorm(sims * length(b)), sims) %*% chol(V), 2, b, `+`) }
  Bm <- draw(med); Bo <- draw(out)
  at <- function(t) { x <- d; x[[treat]] <- if (is.factor(d[[treat]])) factor(t, levels = levels(d[[treat]])) else t; x }
  d1 <- at(treated); d0 <- at(control)
  setm <- function(x, m) { x[[mediator]] <- if (is.factor(d[[mediator]])) factor(levels(d[[mediator]])[m], levels = levels(d[[mediator]])) else m; x }
  eff <- t(vapply(seq_len(sims), function(s) {
    fm <- med; fm$coefficients[] <- Bm[s, ]; fo <- out; fo$coefficients[] <- Bo[s, ]
    m1 <- simulate(fm, newdata = d1)[, 1]; m0 <- simulate(fm, newdata = d0)[, 1]
    y <- function(x, m) as.numeric(stats::predict(fo, setm(x, m), type = "mean"))
    y11 <- y(d1, m1); y10 <- y(d1, m0); y01 <- y(d0, m1); y00 <- y(d0, m0)
    c(acme0 = mean(y01 - y00), acme1 = mean(y11 - y10), ade0 = mean(y10 - y00), ade1 = mean(y11 - y01), total = mean(y11 - y00))
  }, numeric(5)))
  eff <- cbind(eff, acme = (eff[, 1] + eff[, 2]) / 2, ade = (eff[, 3] + eff[, 4]) / 2)
  eff <- cbind(eff, prop = eff[, "acme"] / eff[, "total"])
  a <- (1 - level) / 2
  # Means of the draws; the median for the proportion mediated, a ratio
  # (as mediation does).
  est <- colMeans(eff); est["prop"] <- stats::median(eff[, "prop"])
  data.frame(estimate = est, lower = apply(eff, 2, stats::quantile, a), upper = apply(eff, 2, stats::quantile, 1 - a),
             p_value = apply(eff, 2, function(v) min(1, 2 * min(mean(v <= 0), mean(v >= 0)))),
             row.names = c("ACME (control)", "ACME (treated)", "ADE (control)", "ADE (treated)", "Total effect",
                           "ACME (average)", "ADE (average)", "Prop. mediated"))
}
