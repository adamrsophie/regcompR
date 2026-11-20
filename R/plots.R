# Base graphics plots for distributions and fitted models.

palette_rc <- function(n) {
  base <- c("#1b6ca8", "#d1495b", "#2e8b57", "#edae49", "#7b4b94", "#00798c",
            "#8c564b", "#4d4d4d")
  rep_len(base, n)
}

plot_range <- function(d, xlim) {
  if (!is.null(xlim)) return(xlim)
  fam <- d$family
  q <- tryCatch(c(min(quantile(d, 0.001)), max(quantile(d, 0.999))),
                error = function(e) c(-5, 5))
  s <- fam$support(d$pars)
  lo <- max(min(q), min(s$lower)); hi <- min(max(q), max(s$upper))
  if (!is.finite(lo)) lo <- -5
  if (!is.finite(hi)) hi <- lo + 10
  if (lo == hi) hi <- lo + 1
  c(lo, hi)
}

#' Plot a distribution
#'
#' Draws the density (continuous) or probability mass (discrete), the
#' distribution function, the survival function, the hazard or the quantile
#' function of a distribution object. Vectorized distribution objects are
#' overlaid (up to `max_curves`).
#'
#' @param x A distribution object.
#' @param type `"density"`, `"cdf"`, `"survival"`, `"hazard"` or
#'   `"quantile"`.
#' @param xlim Optional range of the response.
#' @param n Number of evaluation points for continuous laws.
#' @param max_curves Maximum number of overlaid distributions.
#' @param col Colors.
#' @param add Add to an existing plot.
#' @param legend Show a legend when several curves are drawn.
#' @param main,xlab,ylab Titles.
#' @param ... Passed to the plotting functions.
#' @return The input, invisibly.
#' @export
#' @examples
#' plot(distribution("gamma", mu = c(1, 2, 4), sigma = 0.5))
#' plot(distribution("negative_binomial_nb2", mu = 5, phi = 2))
#' plot(distribution("weibull", shape = 1.5, scale = 2), type = "hazard")
plot.regcomp_dist <- function(x, type = c("density", "cdf", "survival", "hazard", "quantile"),
                              xlim = NULL, n = 401, max_curves = 8, col = NULL, add = FALSE,
                              legend = TRUE, main = NULL, xlab = NULL, ylab = NULL, ...) {
  type <- match.arg(type)
  fam <- x$family
  if (fam$variate_type != "scalar") return(plot_multivariate(x, ...))
  m <- min(length(x), max_curves)
  if (is.null(col)) col <- palette_rc(m)
  sub <- function(i) structure(list(family = fam, pars = lapply(x$pars, function(v) rep_len(v, length(x))[i])),
                               class = "regcomp_dist")
  if (type == "quantile") {
    p <- seq(0.001, 0.999, length.out = n)
    vals <- sapply(seq_len(m), function(i) quantile(sub(i), p))
    if (!add) graphics::matplot(p, vals, type = if (fam$discrete) "s" else "l", lty = 1, col = col,
                                main = if (is.null(main)) fam$display_name else main,
                                xlab = if (is.null(xlab)) "probability" else xlab,
                                ylab = if (is.null(ylab)) "quantile" else ylab, ...)
    else graphics::matlines(p, vals, lty = 1, col = col, ...)
    return(invisible(x))
  }
  rng <- plot_range(sub(seq_len(m)), xlim)
  xs <- if (fam$discrete) seq(ceiling(rng[1]), floor(rng[2])) else seq(rng[1], rng[2], length.out = n)
  if (fam$discrete && length(xs) > 2000) xs <- unique(round(seq(rng[1], rng[2], length.out = 2000)))
  f <- function(i) {
    d <- sub(i)
    switch(type,
      density = prob(d, xs),
      cdf = cdf(d, xs),
      survival = survival(d, xs),
      hazard = prob(d, xs) / pmax(if (fam$discrete) survival(d, xs - 1) else survival(d, xs), 1e-300))
  }
  vals <- sapply(seq_len(m), f)
  vals <- matrix(vals, ncol = m)
  vals[!is.finite(vals)] <- NA
  ylab0 <- switch(type, density = if (fam$discrete) "probability" else "density",
                  cdf = "F(y)", survival = "S(y)", hazard = "hazard")
  main0 <- if (is.null(main)) fam$display_name else main
  if (fam$discrete && type %in% c("density", "hazard")) {
    off <- if (m > 1) (seq_len(m) - (m + 1) / 2) * 0.6 / m else 0
    if (!add) graphics::plot(range(xs), c(0, max(vals, na.rm = TRUE)), type = "n", main = main0,
                             xlab = if (is.null(xlab)) "y" else xlab,
                             ylab = if (is.null(ylab)) ylab0 else ylab, ...)
    for (i in seq_len(m)) graphics::segments(xs + off[i], 0, xs + off[i], vals[, i], col = col[i], lwd = 2)
  } else {
    if (!add) graphics::matplot(xs, vals, type = if (fam$discrete) "s" else "l", lty = 1, lwd = 2,
                                col = col, main = main0, xlab = if (is.null(xlab)) "y" else xlab,
                                ylab = if (is.null(ylab)) ylab0 else ylab, ...)
    else graphics::matlines(xs, vals, type = if (fam$discrete) "s" else "l", lty = 1, lwd = 2, col = col)
  }
  if (legend && m > 1 && !add) {
    labs <- vapply(seq_len(m), function(i) paste(names(x$pars), "=",
      vapply(x$pars, function(v) format(rep_len(v, length(x))[i], digits = 3), ""), collapse = ", "), "")
    graphics::legend("topright", legend = labs, col = col, lwd = 2, bty = "n", cex = 0.8)
  }
  invisible(x)
}

plot_multivariate <- function(x, ...) {
  fam <- x$family
  if (!is.null(fam$extra$plot)) return(invisible(fam$extra$plot(x, ...)))
  set.seed(1)
  s <- family_rng(fam, 500, x$pars)
  if (is.matrix(s) && ncol(s) >= 2) {
    graphics::pairs(s[, seq_len(min(4, ncol(s))), drop = FALSE], pch = 16, cex = 0.4,
                    col = grDevices::adjustcolor("#1b6ca8", 0.4),
                    main = paste(fam$display_name, "(random draws)"))
  } else graphics::hist(s, main = fam$display_name)
  invisible(x)
}

#' Plot a fitted model
#'
#' Diagnostic and interpretation plots:
#' * `"qq"`: normal QQ plot of randomized quantile residuals with a
#'   simulation envelope.
#' * `"pit"`: histogram of PIT values (uniform under a correct model).
#' * `"residuals"`: quantile residuals against the fitted mean.
#' * `"rootogram"`: hanging rootogram of observed and expected frequencies
#'   (discrete responses).
#' * `"ppc"`: observed response distribution against simulated replicates
#'   (posterior predictive check for Bayesian fits).
#' * `"effects"`: a distribution parameter or the mean against one
#'   covariate, other covariates at their mean or most frequent level, with
#'   a confidence or credible band.
#' * `"fitted"`: fitted distributions of a few observations.
#' * `"trace"`: trace plots of posterior draws (Bayesian fits).
#' * `"ranef"`: caterpillar plot of the random effects of each grouping
#'   term with approximate 95% intervals (conditional standard deviations
#'   from the Laplace Hessian, or posterior quantiles).
#'
#' For ordinal and nominal families `"effects"` draws the category
#' probabilities against the covariate.
#'
#' @param x A fitted model.
#' @param which Plot type(s).
#' @param variable Covariate for `"effects"`.
#' @param par Distribution parameter for `"effects"`, or `"mean"`.
#' @param nsim Number of simulated replicates for `"qq"` and `"ppc"`.
#' @param rows Observations for `"fitted"`.
#' @param ask Ask before each new page when several plots are drawn.
#' @param ... Passed to plotting functions.
#' @return The input, invisibly.
#' @export
#' @examples
#' set.seed(1)
#' d <- data.frame(x = runif(200))
#' d$y <- rpois(200, exp(1 + d$x))
#' fit <- regcomp(y ~ x, family = "poisson", data = d)
#' plot(fit, which = "rootogram")
#' plot(fit, which = "effects", variable = "x")
plot.regcomp_fit <- function(x, which = c("qq", "pit", "residuals", "rootogram"),
                             variable = NULL, par = "mean", nsim = 50, rows = 1:4,
                             ask = FALSE, ...) {
  which <- match.arg(which, c("qq", "pit", "residuals", "rootogram", "ppc", "effects",
                              "fitted", "trace", "ranef"), several.ok = TRUE)
  if (!x$family$discrete) which <- setdiff(which, "rootogram")
  if (length(which) > 1 && ask) { op <- grDevices::devAskNewPage(TRUE); on.exit(grDevices::devAskNewPage(op)) }
  for (w in which) switch(w,
    qq = plot_qq(x, nsim, ...), pit = plot_pit(x, ...), residuals = plot_resid(x, ...),
    rootogram = plot_rootogram(x, ...), ppc = plot_ppc(x, nsim, ...),
    effects = plot_effects(x, variable, par, ...), fitted = plot(fitted_distribution(x)[rows], ...),
    ranef = plot_ranef(x, ...),
    trace = plot_trace(x, ...))
  invisible(x)
}

#' @export
`[.regcomp_dist` <- function(x, i) {
  # Only the observation axis is sliced: parameters that vary by observation
  # (regression parameters and data arguments; every parameter of a plain
  # scalar law). Shared thresholds, mean vectors and matrices stay whole.
  fam <- x$family
  obs <- vapply(fam$pars, function(p) fam$variate_type == "scalar" || isTRUE(p$regression) || isTRUE(p$data), TRUE)
  obs <- obs & vapply(x$pars[names(fam$pars)], function(v) is.null(dim(v)), TRUE)
  n <- max(c(1L, lengths(x$pars[names(fam$pars)][obs])))
  pars <- x$pars
  for (k in names(fam$pars)[obs]) pars[[k]] <- rep_len(pars[[k]], n)[i]
  structure(list(family = fam, pars = pars), class = "regcomp_dist")
}

plot_qq <- function(fit, nsim = 50, ...) {
  r <- residuals(fit, "quantile")
  r <- r[is.finite(r)]
  n <- length(r)
  th <- stats::qnorm(stats::ppoints(n))
  env <- apply(matrix(stats::rnorm(n * nsim), n), 2, sort)
  band <- apply(env, 1, stats::quantile, c(0.025, 0.975))
  graphics::plot(th, sort(r), pch = 16, cex = 0.6, col = "#1b6ca8", xlab = "standard normal quantiles",
                 ylab = "quantile residuals", main = "Quantile residuals", ...)
  graphics::polygon(c(th, rev(th)), c(band[1, ], rev(band[2, ])), col = grDevices::adjustcolor("grey", 0.4), border = NA)
  graphics::points(th, sort(r), pch = 16, cex = 0.6, col = "#1b6ca8")
  graphics::abline(0, 1, col = "#d1495b", lwd = 2)
}

plot_pit <- function(fit, ...) {
  u <- residuals(fit, "pit")
  graphics::hist(u, breaks = seq(0, 1, by = 0.1), col = "#9ecae1", border = "white",
                 main = "PIT histogram", xlab = "PIT", freq = FALSE, ...)
  graphics::abline(h = 1, col = "#d1495b", lwd = 2, lty = 2)
}

plot_resid <- function(fit, ...) {
  r <- residuals(fit, "quantile")
  m <- suppressWarnings(predict(fit, type = "mean"))
  if (is.matrix(m)) m <- drop(m %*% seq_len(ncol(m)))   # expected category score
  graphics::plot(m, r, pch = 16, cex = 0.6, col = grDevices::adjustcolor("#1b6ca8", 0.6),
                 xlab = "fitted mean", ylab = "quantile residual", main = "Residuals against fitted", ...)
  graphics::abline(h = c(-2, 0, 2), lty = c(2, 1, 2), col = "#d1495b")
  ok <- is.finite(m) & is.finite(r)
  if (sum(ok) > 10) graphics::lines(stats::lowess(m[ok], r[ok]), col = "#2e8b57", lwd = 2)
}

plot_rootogram <- function(fit, max_y = NULL, ...) {
  y <- fit$spec$y
  if (is.null(max_y)) max_y <- min(max(y, na.rm = TRUE), stats::quantile(y, 0.99, na.rm = TRUE) * 1.5 + 5)
  ks <- 0:max_y
  obs <- tabulate(factor(pmin(y, max_y + 1), levels = 0:(max_y + 1)), max_y + 2)[seq_along(ks)]
  fd <- fitted_distribution(fit)
  expd <- vapply(ks, function(k) sum(exp(eval_lpdf(fd$family, rep(k, length(y)), recycle_pars(fd$pars, length(y))))), 0)
  top <- sqrt(expd)
  bottom <- top - sqrt(obs)
  graphics::plot(range(ks) + c(-0.5, 0.5), range(c(bottom, top, 0)), type = "n", xlab = "y",
                 ylab = "sqrt(frequency)", main = "Hanging rootogram", ...)
  graphics::rect(ks - 0.4, bottom, ks + 0.4, top, col = "#c6dbef", border = "#6baed6")
  graphics::lines(ks, top, type = "b", pch = 16, col = "#d1495b", lwd = 2)
  graphics::abline(h = 0, lty = 2)
}

plot_ppc <- function(fit, nsim = 50, ...) {
  y <- fit$spec$y
  sims <- simulate(fit, nsim)
  if (fit$family$discrete) {
    mx <- max(c(y, stats::quantile(sims, 0.99)), na.rm = TRUE)
    ks <- 0:mx
    tab <- function(v) tabulate(factor(pmin(v, mx), levels = ks), length(ks)) / length(v)
    so <- apply(sims, 2, tab)
    graphics::matplot(ks, so, type = "l", lty = 1, col = grDevices::adjustcolor("#9ecae1", 0.5),
                      xlab = "y", ylab = "proportion", main = "Predictive check", ...)
    graphics::lines(ks, tab(y), type = "b", pch = 16, col = "#d1495b", lwd = 2)
  } else {
    dens <- lapply(seq_len(ncol(sims)), function(j) stats::density(sims[, j], na.rm = TRUE))
    dy <- stats::density(y, na.rm = TRUE)
    graphics::plot(dy, type = "n", main = "Predictive check", xlab = "y",
                   ylim = c(0, max(c(dy$y, unlist(lapply(dens, `[[`, "y"))))), ...)
    for (d in dens) graphics::lines(d, col = grDevices::adjustcolor("#9ecae1", 0.5))
    graphics::lines(dy, col = "#d1495b", lwd = 2)
  }
}

plot_effects <- function(fit, variable, par = "mean", level = 0.95, n = 50, ...) {
  data <- model_data(fit)
  if (is.null(variable)) {
    vars <- unique(unlist(lapply(fit$spec$preds, function(p) all.vars(p$formula))))
    vars <- intersect(vars, names(data))
    if (!length(vars)) stop_input("No covariate to plot.")
    variable <- vars[1]
  }
  nd <- reference_grid(data, variable, n)
  if (fit$family$variate_type == "categorical") {
    P <- predict(fit, nd, type = "probability")
    xv <- nd[[variable]]; xn <- if (is.numeric(xv)) xv else seq_along(xv)
    lab <- fit$spec$ylev %||% seq_len(ncol(P))
    graphics::matplot(xn, P, type = "l", lty = 1, lwd = 2, col = palette_rc(ncol(P)), ylim = c(0, 1),
                      xlab = variable, ylab = "probability", main = paste("Category probabilities by", variable), ...)
    graphics::legend("topright", legend = lab, col = palette_rc(ncol(P)), lty = 1, lwd = 2, bty = "n")
    return(invisible(NULL))
  }
  type <- if (par %in% c("mean", "median", "variance", "sd")) par else "parameter"
  pr <- tryCatch(predict(fit, nd, type = type, interval = "confidence", level = level, ndraws = 300),
                 error = function(e) NULL)
  if (is.null(pr)) { est <- predict(fit, nd, type = type); lo <- hi <- NULL } else if (is.data.frame(pr)) {
    est <- pr$estimate; lo <- pr$lower; hi <- pr$upper
  } else { est <- pr$estimate[, par]; lo <- pr$lower[, par]; hi <- pr$upper[, par] }
  if (type == "parameter" && is.matrix(est)) est <- est[, par]
  xv <- nd[[variable]]
  if (is.factor(xv) || is.character(xv)) {
    xi <- seq_along(xv)
    graphics::plot(xi, est, xaxt = "n", pch = 16, col = "#1b6ca8", xlab = variable, ylab = par,
                   ylim = range(c(est, lo, hi), na.rm = TRUE), main = paste("Effect of", variable), ...)
    graphics::axis(1, xi, as.character(xv))
    if (!is.null(lo)) graphics::arrows(xi, lo, xi, hi, angle = 90, code = 3, length = 0.05, col = "#1b6ca8")
  } else {
    graphics::plot(xv, est, type = "n", xlab = variable, ylab = par, ylim = range(c(est, lo, hi), na.rm = TRUE),
                   main = paste("Effect of", variable), ...)
    if (!is.null(lo)) graphics::polygon(c(xv, rev(xv)), c(lo, rev(hi)), col = grDevices::adjustcolor("#1b6ca8", 0.2), border = NA)
    graphics::lines(xv, est, col = "#1b6ca8", lwd = 2)
  }
}

plot_trace <- function(fit, pars = NULL, ...) {
  if (fit$inference != "bayesian") stop_unsupported("Trace plots for a frequentist fit")
  dr <- fit$draws
  if (is.null(pars)) pars <- colnames(dr)[seq_len(min(6, ncol(dr)))]
  op <- graphics::par(mfrow = c(ceiling(length(pars) / 2), min(2, length(pars))), mar = c(3, 3, 2, 1))
  on.exit(graphics::par(op))
  for (p in pars) graphics::plot(dr[, p], type = "l", col = "#1b6ca8", main = p, xlab = "", ylab = "", ...)
}

model_data <- function(fit) {
  d <- tryCatch(eval(fit$call$data, environment(fit$spec$formulas[[1]])), error = function(e) NULL)
  if (is.null(d)) stop_input("The data used for fitting are not available; give newdata explicitly.")
  as.data.frame(d)
}

reference_grid <- function(data, variable, n = 50) {
  ref <- lapply(data, function(v) {
    if (is.numeric(v)) mean(v, na.rm = TRUE)
    else if (is.factor(v)) factor(names(which.max(table(v))), levels = levels(v))
    else names(which.max(table(v)))
  })
  v <- data[[variable]]
  vals <- if (is.numeric(v)) seq(min(v, na.rm = TRUE), max(v, na.rm = TRUE), length.out = n)
          else if (is.factor(v)) factor(levels(v), levels = levels(v)) else unique(v)
  nd <- as.data.frame(lapply(ref, function(r) rep(r, length(vals))))
  for (nm in names(data)) if (is.factor(data[[nm]])) nd[[nm]] <- factor(nd[[nm]], levels = levels(data[[nm]]))
  nd[[variable]] <- vals
  nd
}

plot_ranef <- function(fit, ...) {
  spec <- fit$spec
  if (!spec$Q) stop_input("The model has no random effects.")
  if (fit$inference == "bayesian") {
    est <- colMeans(fit$u_draws); lo <- apply(fit$u_draws, 2, stats::quantile, 0.025); hi <- apply(fit$u_draws, 2, stats::quantile, 0.975)
  } else {
    inner <- laplace_inner(spec, unname(fit$coefficients), unname(fit$theta), fit$u, regcomp_control())
    sdv <- sqrt(pmax(Matrix::diag(Matrix::solve(inner$negH)), 0))
    est <- fit$u; lo <- est - 1.96 * sdv; hi <- est + 1.96 * sdv
  }
  for (j in seq_along(spec$blocks)) {
    b <- spec$blocks[[j]]; if (b$type == "smooth") next
    for (e in seq_len(b$q)) {
      idx <- spec$u_idx[[j]][(seq_len(b$L) - 1) * b$q + e]
      o <- order(est[idx]); k <- seq_along(idx)
      graphics::plot(est[idx][o], k, xlim = range(c(lo[idx], hi[idx])), pch = 16, col = "#1b6ca8", yaxt = "n",
        xlab = "random effect", ylab = b$group, main = paste(b$par, ":", b$label, b$terms[e]), ...)
      graphics::segments(lo[idx][o], k, hi[idx][o], k, col = "#1b6ca8")
      graphics::abline(v = 0, lty = 2, col = "grey50")
      if (length(idx) <= 40) graphics::axis(2, k, b$levels[o], las = 1, cex.axis = 0.6)
    }
  }
  invisible(fit)
}
