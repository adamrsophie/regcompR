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
