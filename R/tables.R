# Tidy summaries, regression tables and export of tables and plots.

# Split "mu:x" into the distribution parameter and the term.
split_name <- function(nm) {
  comp <- sub(":.*$", "", nm)
  list(component = comp, term = substring(nm, nchar(comp) + 2))
}

#' Tidy summaries of a fitted model
#'
#' Methods for the generics `tidy()`, `glance()` and `augment()` of the
#' generics package (also used by broom, modelsummary and gtsummary).
#' `tidy()` returns one row per coefficient with the distribution parameter
#' it belongs to (`component`), and one row per random effect variance
#' parameter (`effect = "ran_pars"`). Random effect standard deviations are
#' reported on their natural scale (`sd(...)`). For Bayesian fits the
#' estimate is the posterior mean, `std.error` the posterior standard
#' deviation, the interval a credible interval, and `rhat` and `ess` are
#' added.
#'
#' @param x A fitted model from [regcomp()].
#' @param conf.int Add confidence (credible) intervals.
#' @param conf.level Interval level.
#' @param exponentiate Exponentiate coefficients and their intervals (odds
#'   ratios, rate ratios); random effect parameters are not changed.
#' @param vcov `"model"`, `"sandwich"` or `"quasi"` standard errors (see
#'   [summary.regcomp_fit()]).
#' @param standardize Multiply each coefficient (and its interval) by the
#'   standard deviation of its design column: the change of the linear
#'   predictor per standard deviation of the covariate (partially
#'   standardized coefficients; meaningful for any family, unlike scaling
#'   the response). Applied before `exponentiate`.
#' @param data,newdata For `augment()`: data to add columns to; defaults to
#'   the analysis rows of the fit.
#' @param ... Unused.
#' @return A data frame.
#' @name tidy.regcomp_fit
#' @examples
#' d <- data.frame(x = rnorm(100)); d$y <- rpois(100, exp(0.5 + 0.3 * d$x))
#' f <- regcomp(y ~ x, family = "poisson", data = d)
#' if (requireNamespace("generics", quietly = TRUE)) {
#'   print(generics::tidy(f, exponentiate = TRUE))
#'   print(generics::glance(f))
#' }
NULL

#' @rdname tidy.regcomp_fit
#' @exportS3Method generics::tidy
tidy.regcomp_fit <- function(x, conf.int = TRUE, conf.level = 0.95, exponentiate = FALSE, vcov = "model", standardize = FALSE, ...) {
  s <- summary(x, level = conf.level, vcov = vcov, ...)$table
  bayes <- x$inference == "bayesian"
  nm <- rownames(s); sp <- split_name(nm)
  ran <- nm %in% names(x$theta)
  out <- data.frame(
    term = sp$term, component = sp$component, effect = ifelse(ran, "ran_pars", "fixed"),
    estimate = if (bayes) s$mean else s$estimate,
    std.error = if (bayes) s$sd else s$std_error,
    statistic = if (bayes) NA_real_ else s$z,
    p.value = if (bayes) NA_real_ else s$p_value,
    conf.low = s$lower, conf.high = s$upper, stringsAsFactors = FALSE)
  if (bayes) { out$rhat <- s$rhat; out$ess <- s$ess_bulk }
  # Random effect standard deviations on their natural scale.
  lsd <- ran & grepl("log_sd\\(", out$term)
  if (any(lsd)) {
    out$std.error[lsd] <- exp(out$estimate[lsd]) * out$std.error[lsd]
    out[lsd, c("estimate", "conf.low", "conf.high")] <- exp(out[lsd, c("estimate", "conf.low", "conf.high")])
    # "1 | g:log_sd((Intercept))" becomes "sd((Intercept) | g)".
    blk <- sub(":log_sd\\(.*$", "", out$term[lsd]); co <- sub("^.*:log_sd\\((.*)\\)$", "\\1", out$term[lsd])
    grp <- ifelse(grepl("\\|", blk) & !grepl("\\(", blk), trimws(sub("^.*\\|", "", blk)), blk)
    out$term[lsd] <- sprintf("sd(%s | %s)", co, grp)
    out$statistic[lsd] <- NA; out$p.value[lsd] <- NA
  }
  # Cholesky parameters of unstructured blocks: report the correlations
  # (plug in estimates, without intervals).
  ch <- which(ran & grepl(":chol\\(", out$term))
  if (length(ch)) for (j in seq_along(x$spec$blocks)) {
    b <- x$spec$blocks[[j]]; ti <- x$spec$theta_idx[[j]]
    k <- match(names(x$theta)[ti], nm); k <- k[k %in% ch]
    if (!length(k) || b$q < 2) next
    R <- stats::cov2cor(block_cov(b, x$theta[ti]))
    pr <- which(lower.tri(R), arr.ind = TRUE)
    lab <- sub("^.*:chol\\((.*)\\)$", "\\1", out$term[k])
    grp <- sub(":chol\\(.*$", "", out$term[k]); grp <- ifelse(grepl("\\|", grp) & !grepl("\\(", grp), trimws(sub("^.*\\|", "", grp)), grp)
    out$estimate[k] <- R[pr][seq_along(k)]
    out$term[k] <- sprintf("cor(%s | %s)", sub(",", ", ", lab, fixed = TRUE), grp)
    out[k, c("std.error", "statistic", "p.value", "conf.low", "conf.high")] <- NA
  }
  if (standardize) {
    # Coefficient times the standard deviation of its design column: the
    # change in the linear predictor per standard deviation of the covariate.
    for (k in x$spec$rpars) { X <- x$spec$preds[[k]]$X; j <- match(paste0(k, ":", colnames(X)), nm)
      s <- apply(X, 2, stats::sd); s[colnames(X) == "(Intercept)" | s == 0] <- 1
      out[j, c("estimate", "std.error", "conf.low", "conf.high")] <- out[j, c("estimate", "std.error", "conf.low", "conf.high")] * s }
  }
  if (exponentiate) {
    fx <- !ran
    out$std.error[fx] <- exp(out$estimate[fx]) * out$std.error[fx]
    out[fx, c("estimate", "conf.low", "conf.high")] <- exp(out[fx, c("estimate", "conf.low", "conf.high")])
  }
  if (!conf.int) out$conf.low <- out$conf.high <- NULL
  rownames(out) <- NULL
  out
}

#' @rdname tidy.regcomp_fit
#' @exportS3Method generics::glance
glance.regcomp_fit <- function(x, ...) {
  base <- data.frame(family = x$family$display_name, engine = x$engine, nobs = x$nobs, stringsAsFactors = FALSE)
  if (x$inference == "bayesian")
    return(cbind(base, draws = nrow(x$draws), max_rhat = x$diagnostics$max_rhat,
                 min_ess = x$diagnostics$min_ess_bulk, divergences = x$diagnostics$divergences %||% NA))
  cbind(base, df = x$df, logLik = x$logLik, AIC = -2 * x$logLik + 2 * x$df,
        BIC = -2 * x$logLik + log(x$nobs) * x$df, converged = isTRUE(x$diagnostics$converged))
}

#' @rdname tidy.regcomp_fit
#' @exportS3Method generics::augment
augment.regcomp_fit <- function(x, data = NULL, newdata = NULL, ...) {
  if (!is.null(newdata)) {
    out <- newdata; out$.fitted <- stats::predict(x, newdata, type = "mean"); return(out)
  }
  out <- data %||% x$model_data
  if (is.null(out)) stop_unsupported("augment() of a fit without stored analysis rows (pass data)")
  out$.fitted <- stats::fitted(x)
  out$.resid <- stats::residuals(x)
  out
}

#' Regression tables
#'
#' Builds a table of estimates for one or several fits, with rows grouped by
#' distribution parameter, random effect standard deviations and fit
#' statistics. It prints in the console and renders as a formatted table in
#' R Markdown and Quarto documents (HTML, PDF and Word output) through
#' `knitr::kable()`; [export_table()] writes it to a file.
#'
#' @param ... Fitted models, optionally named (the names become column
#'   headers).
#' @param conf.level Interval level.
#' @param exponentiate Exponentiate coefficients (all of them, or a logical
#'   vector with one entry per model).
#' @param show `"ci"` for estimates with intervals, `"se"` for estimates with
#'   standard errors.
#' @param stars Add significance stars (`* p < 0.05`, `** p < 0.01`,
#'   `*** p < 0.001`) to Wald estimates.
#' @param digits Number of decimal places.
#' @param stats Fit statistics to add: any of `"nobs"`, `"logLik"`, `"AIC"`,
#'   `"BIC"`; Bayesian fits report draws and the largest R-hat instead.
#' @param vcov Standard errors for frequentist fits: `"model"`, `"sandwich"`
#'   or `"quasi"`.
#' @param standardize Standardized coefficients (see [tidy.regcomp_fit()]).
#' @return A `regcomp_table`: a data frame with columns `parameter`, `term`
#'   and one column per model.
#' @export
#' @examples
#' d <- data.frame(x = rnorm(200)); d$y <- rnbinom(200, mu = exp(1 + 0.4 * d$x), size = 2)
#' tab <- regcomp_table(Poisson = regcomp(y ~ x, family = "poisson", data = d),
#'                      NB = regcomp(y ~ x, family = "nb2", data = d), exponentiate = TRUE)
#' tab
#' export_table(tab, tempfile(fileext = ".html"))
regcomp_table <- function(..., conf.level = 0.95, exponentiate = FALSE, show = c("ci", "se"),
                          stars = FALSE, digits = 3, stats = c("nobs", "logLik", "AIC", "BIC"), vcov = "model",
                          standardize = FALSE) {
  show <- match.arg(show)
  fits <- list(...)
  if (!length(fits) || !all(vapply(fits, inherits, TRUE, "regcomp_fit"))) stop_input("regcomp_table() needs fitted models from regcomp().")
  mn <- names(fits) %||% rep("", length(fits))
  mn[mn == ""] <- paste0("(", which(mn == ""), ")")
  ex <- rep_len(exponentiate, length(fits))
  fmt <- function(v) ifelse(is.na(v), "", formatC(v, digits = digits, format = "f"))
  cells <- lapply(seq_along(fits), function(i) {
    t <- tidy.regcomp_fit(fits[[i]], conf.level = conf.level, exponentiate = ex[i], vcov = vcov, standardize = standardize)
    st <- if (stars && !all(is.na(t$p.value))) as.character(cut(t$p.value, c(-Inf, 0.001, 0.01, 0.05, Inf), c("***", "**", "*", ""))) else ""
    st[is.na(st)] <- ""
    val <- paste0(fmt(t$estimate), st, if (show == "ci") ifelse(is.na(t$conf.low), "", paste0(" [", fmt(t$conf.low), ", ", fmt(t$conf.high), "]"))
                  else ifelse(is.na(t$std.error), "", paste0(" (", fmt(t$std.error), ")")))
    par <- ifelse(t$effect == "ran_pars", paste(t$component, "(random effects)"), t$component)
    data.frame(parameter = par, term = t$term, effect = t$effect, value = val, stringsAsFactors = FALSE)
  })
  keys <- unique(do.call(rbind, lapply(cells, `[`, c("parameter", "term", "effect"))))
  keys <- keys[order(keys$effect == "ran_pars", match(keys$parameter, unique(keys$parameter))), ]
  out <- keys[, c("parameter", "term")]
  for (i in seq_along(fits)) {
    k <- match(paste(keys$parameter, keys$term), paste(cells[[i]]$parameter, cells[[i]]$term))
    out[[mn[i]]] <- ifelse(is.na(k), "", cells[[i]]$value[k])
  }
  sr <- lapply(fits, function(f) {
    if (f$inference == "bayesian") return(c(nobs = format(f$nobs), draws = format(nrow(f$draws)),
                                            max_rhat = formatC(f$diagnostics$max_rhat, digits = 3, format = "f")))
    c(nobs = format(f$nobs), logLik = fmt(f$logLik), AIC = fmt(-2 * f$logLik + 2 * f$df),
      BIC = fmt(-2 * f$logLik + log(f$nobs) * f$df))
  })
  snames <- unique(unlist(lapply(sr, names)))
  snames <- snames[snames %in% c(stats, "draws", "max_rhat")]
  if (length(snames)) {
    srow <- data.frame(parameter = "", term = snames, stringsAsFactors = FALSE)
    for (i in seq_along(fits)) srow[[mn[i]]] <- ifelse(is.na(sr[[i]][snames]), "", sr[[i]][snames])
    out <- rbind(out, srow)
  }
  rownames(out) <- NULL
  kinds <- unique(vapply(fits, `[[`, "", "inference"))
  attr(out, "note") <- paste0(if (show == "ci") sprintf("%g%% %s intervals in brackets", 100 * conf.level,
                                                      if (identical(kinds, "bayesian")) "credible" else "confidence")
                              else "Standard errors in parentheses",
                              if (standardize) "; coefficients per standard deviation of the covariate" else "",
                              if (any(ex)) "; coefficients exponentiated" else "",
                              if (stars) "; * p < 0.05, ** p < 0.01, *** p < 0.001" else "", ".")
  class(out) <- c("regcomp_table", "data.frame")
  out
}

#' @export
print.regcomp_table <- function(x, ...) {
  d <- as.data.frame(x)
  d$parameter[c(FALSE, d$parameter[-1] == d$parameter[-nrow(d)])] <- ""
  print.data.frame(d, row.names = FALSE, right = FALSE)
  cat(attr(x, "note"), "\n")
  invisible(x)
}

#' @export
as.data.frame.regcomp_table <- function(x, ...) { attr(x, "note") <- NULL; class(x) <- "data.frame"; x }

# The kable format for the current knitr output (HTML, LaTeX or Markdown,
# which pandoc turns into Word tables).
kable_format <- function() {
  if (knitr::is_latex_output()) "latex" else if (knitr::is_html_output()) "html" else "pipe"
}

#' @exportS3Method knitr::knit_print
knit_print.regcomp_table <- function(x, ...) {
  d <- as.data.frame(x)
  d$parameter[c(FALSE, d$parameter[-1] == d$parameter[-nrow(d)])] <- ""
  fm <- kable_format()
  k <- knitr::kable(d, format = fm, col.names = c("Parameter", "Term", names(d)[-(1:2)]), align = c("l", "l", rep("r", ncol(d) - 2)))
  note <- attr(x, "note")
  knitr::asis_output(paste0(paste(k, collapse = "\n"), "\n\n", if (fm == "latex") note else paste0("*", note, "*"), "\n"))
}

#' Export tables and plots
#'
#' `export_table()` writes a [regcomp_table()] (or the table of one or
#' several fits) to a file whose format follows the extension: `.csv`,
#' `.html`, `.tex` (LaTeX), `.md` (Markdown) or `.docx` (Word, with the
#' flextable package). `export_plot()` draws a diagnostic plot of a fit
#' (see [plot.regcomp_fit()]) into a `.png`, `.jpeg`, `.tiff`, `.pdf` or
#' `.svg` file.
#'
#' @param x A `regcomp_table`, a fitted model, or a list of fitted models.
#' @param file Output file.
#' @param ... For `export_table()`, arguments to [regcomp_table()] when `x`
#'   holds fits; for `export_plot()`, arguments to [plot.regcomp_fit()].
#' @param fit A fitted model.
#' @param which The plot to draw.
#' @param width,height Size in inches.
#' @param res Resolution in pixels per inch for bitmap formats.
#' @return The file name, invisibly.
#' @export
#' @examples
#' d <- data.frame(x = rnorm(100)); d$y <- rpois(100, exp(0.5 + 0.3 * d$x))
#' f <- regcomp(y ~ x, family = "poisson", data = d)
#' export_table(f, tempfile(fileext = ".csv"))
#' export_plot(f, tempfile(fileext = ".png"), which = "qq")
export_table <- function(x, file, ...) {
  if (inherits(x, "regcomp_fit")) x <- regcomp_table(x, ...)
  else if (is.list(x) && !inherits(x, "regcomp_table") && all(vapply(x, inherits, TRUE, "regcomp_fit"))) x <- do.call(regcomp_table, c(x, list(...)))
  if (!inherits(x, "regcomp_table")) stop_input("export_table() needs a regcomp_table or fitted models.")
  ext <- tolower(tools::file_ext(file))
  d <- as.data.frame(x); note <- attr(x, "note")
  cn <- c("Parameter", "Term", names(d)[-(1:2)])
  d$parameter[c(FALSE, d$parameter[-1] == d$parameter[-nrow(d)])] <- ""
  kab <- function(fm) {
    if (!requireNamespace("knitr", quietly = TRUE)) stop_unsupported(sprintf("Writing a .%s table without the 'knitr' package; it", ext))
    paste(knitr::kable(d, format = fm, col.names = cn), collapse = "\n")
  }
  switch(ext,
    csv = utils::write.csv(stats::setNames(d, cn), file, row.names = FALSE),
    html = writeLines(c("<!DOCTYPE html>", "<html><head><meta charset=\"utf-8\"></head><body>", kab("html"),
                        paste0("<p><em>", note, "</em></p>"), "</body></html>"), file),
    tex = writeLines(c(kab("latex"), note), file),
    md = writeLines(c(kab("pipe"), "", paste0("*", note, "*")), file),
    docx = {
      if (!requireNamespace("flextable", quietly = TRUE)) stop_unsupported("Writing a .docx table without the 'flextable' package; it")
      ft <- flextable::add_footer_lines(flextable::autofit(flextable::flextable(stats::setNames(d, cn))), note)
      flextable::save_as_docx(ft, path = file)
    },
    stop_input(sprintf("Unknown table format '.%s': use .csv, .html, .tex, .md or .docx.", ext)))
  invisible(file)
}

#' @rdname export_table
#' @export
export_plot <- function(fit, file, which = "qq", width = 7, height = 5, res = 300, ...) {
  ext <- tolower(tools::file_ext(file))
  switch(ext,
    png = grDevices::png(file, width = width, height = height, units = "in", res = res),
    jpg = , jpeg = grDevices::jpeg(file, width = width, height = height, units = "in", res = res),
    tif = , tiff = grDevices::tiff(file, width = width, height = height, units = "in", res = res),
    pdf = grDevices::pdf(file, width = width, height = height),
    svg = grDevices::svg(file, width = width, height = height),
    stop_input(sprintf("Unknown plot format '.%s': use .png, .jpeg, .tiff, .pdf or .svg.", ext)))
  on.exit(grDevices::dev.off())
  graphics::plot(fit, which = which, ...)
  invisible(file)
}
