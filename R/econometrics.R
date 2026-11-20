# Linear econometric estimators: instrumental variables (two stage least
# squares), the fixed effects within estimator and the Hausman test.

new_regcomp_lm <- function(b, V, res, df, type, call, extra = list()) {
  structure(c(list(coefficients = b, vcov = V, residuals = res, df.residual = df,
                   nobs = length(res), type = type, call = call), extra), class = "regcomp_lm")
}

lm_table <- function(object) {
  b <- object$coefficients; se <- sqrt(diag(object$vcov)); t <- b / se
  data.frame(estimate = b, std_error = se, t = t, p_value = 2 * stats::pt(-abs(t), object$df.residual),
             row.names = names(b))
}

#' @export
print.regcomp_lm <- function(x, digits = 4, ...) {
  cat(sprintf("<regcomp_lm> %s, n = %d, residual df = %d\n", x$type, x$nobs, x$df.residual))
  print(round(lm_table(x), digits))
  if (!is.null(x$diagnostics)) { cat("\nDiagnostic tests:\n"); print(round(x$diagnostics, digits)) }
  invisible(x)
}
#' @export
coef.regcomp_lm <- function(object, ...) object$coefficients
#' @export
vcov.regcomp_lm <- function(object, ...) object$vcov
#' @export
residuals.regcomp_lm <- function(object, ...) object$residuals
#' @export
nobs.regcomp_lm <- function(object, ...) object$nobs
#' @exportS3Method generics::tidy
tidy.regcomp_lm <- function(x, conf.int = TRUE, conf.level = 0.95, ...) {
  tb <- lm_table(x); q <- stats::qt(1 - (1 - conf.level) / 2, x$df.residual)
  out <- data.frame(term = rownames(tb), estimate = tb$estimate, std.error = tb$std_error,
                    statistic = tb$t, p.value = tb$p_value)
  if (conf.int) { out$conf.low <- out$estimate - q * out$std.error; out$conf.high <- out$estimate + q * out$std.error }
  out
}

ols <- function(X, y) {
  q <- qr(X)
  if (q$rank < ncol(X)) stop_input(sprintf("The design is rank deficient; aliased columns: %s.",
                                           paste(colnames(X)[q$pivot[-seq_len(q$rank)]], collapse = ", ")))
  b <- qr.coef(q, y); list(b = b, res = drop(y - X %*% b), XtXi = chol2inv(qr.R(q))[order(q$pivot), order(q$pivot)])
}

f_test <- function(rss0, rss1, df1, df2) {
  f <- ((rss0 - rss1) / df1) / (rss1 / df2)
  c(df1 = df1, df2 = df2, statistic = f, p_value = stats::pf(f, df1, df2, lower.tail = FALSE))
}

#' Instrumental variable regression
#'
#' Two stage least squares for a linear model with endogenous regressors,
#' with the usual diagnostics: a weak instrument F test for each
#' endogenous regressor (the first stage F test of the excluded
#' instruments), the Wu-Hausman test of endogeneity and the Sargan test of
#' overidentifying restrictions. The results equal `ivreg::ivreg()` (and
#' `AER::ivreg()`).
#'
#' @param formula A three part formula `y ~ regressors | instruments`, as in
#'   ivreg: exogenous regressors appear on both sides, endogenous regressors
#'   only on the left of `|`.
#' @param data A data frame.
#' @return An object of class `regcomp_lm` with coefficients, covariance
#'   (`sigma^2 (X'P_Z X)^-1` with `sigma^2 = RSS/(n - p)`), residuals and a
#'   `diagnostics` table.
#' @export
#' @examples
#' set.seed(1)
#' n <- 500; z1 <- rnorm(n); z2 <- rnorm(n); u <- rnorm(n)
#' x <- z1 + 0.5 * z2 + u + rnorm(n)
#' d <- data.frame(y = 1 + 2 * x + u, x, z1, z2)
#' iv_regcomp(y ~ x | z1 + z2, data = d)
iv_regcomp <- function(formula, data) {
  cl <- match.call()
  f <- formula
  if (length(f) != 3 || !is.call(f[[3]]) || !identical(f[[3]][[1]], as.name("|")))
    stop_input("formula must look like y ~ regressors | instruments.")
  env <- environment(f)
  fx <- stats::as.formula(call("~", f[[2]], f[[3]][[2]]), env = env)
  fz <- stats::as.formula(call("~", f[[3]][[3]]), env = env)
  mf <- stats::model.frame(stats::as.formula(call("~", f[[2]], call("+", f[[3]][[2]], f[[3]][[3]])), env = env), data)
  y <- stats::model.response(stats::model.frame(fx, mf))
  X <- stats::model.matrix(fx, mf); Z <- stats::model.matrix(fz, mf)
  n <- length(y); p <- ncol(X)
  if (ncol(Z) < p) stop_input("Fewer instruments than regressors: the model is not identified.")
  Xhat <- Z %*% qr.coef(qr(Z), X); colnames(Xhat) <- colnames(X)
  st <- ols(Xhat, y)
  b <- st$b; res <- drop(y - X %*% b)
  df <- n - p; V <- sum(res^2) / df * st$XtXi
  dimnames(V) <- list(colnames(X), colnames(X)); names(b) <- colnames(X)
  endo <- setdiff(colnames(X), colnames(Z)); excl <- setdiff(colnames(Z), colnames(X))
  dg <- NULL
  if (length(endo)) {
    zin <- intersect(colnames(Z), colnames(X))
    weak <- lapply(endo, function(e) {
      r1 <- ols(Z, X[, e])$res; r0 <- if (length(zin)) ols(Z[, zin, drop = FALSE], X[, e])$res else X[, e]
      f_test(sum(r0^2), sum(r1^2), length(excl), n - ncol(Z)) })
    V1 <- X[, endo, drop = FALSE] - Xhat[, endo, drop = FALSE]
    r0 <- ols(X, y)$res; r1 <- ols(cbind(X, V1), y)$res
    wh <- f_test(sum(r0^2), sum(r1^2), length(endo), n - p - length(endo))
    dg <- rbind(do.call(rbind, weak), wh)
    rownames(dg) <- c(paste0("Weak instruments (", endo, ")"), "Wu-Hausman")
    if (ncol(Z) > p) {
      e <- ols(Z, res)$res; s <- n * (1 - sum(e^2) / sum((res - mean(res))^2))
      dg <- rbind(dg, Sargan = c(df1 = ncol(Z) - p, df2 = NA, statistic = s,
                                 p_value = stats::pchisq(s, ncol(Z) - p, lower.tail = FALSE)))
    }
    dg <- as.data.frame(dg)
  }
  new_regcomp_lm(b, V, res, df, "instrumental variables (2SLS)", cl, list(diagnostics = dg))
}

#' Fixed effects (within) estimator for panel data
#'
#' Least squares after subtracting the unit means (and, for two way
#' effects, after also projecting out time dummies), as `plm(model =
#' "within")`. Equal to least squares with unit (and time) dummies, without
#' building them.
#'
#' @param formula A formula `y ~ x1 + x2` without an intercept term
#'   (absorbed by the unit effects).
#' @param data A data frame.
#' @param id,time Names of the unit and time variables in `data`.
#' @param effect `"individual"` or `"twoways"`.
#' @return An object of class `regcomp_lm`.
#' @export
#' @examples
#' set.seed(1)
#' d <- data.frame(id = rep(1:50, each = 5), t = rep(1:5, 50), x = rnorm(250))
#' a <- rnorm(50); d$x <- d$x + a[d$id]; d$y <- 1 + 0.5 * d$x + a[d$id] + rnorm(250)
#' panel_within(y ~ x, data = d, id = "id")
panel_within <- function(formula, data, id, time = NULL, effect = c("individual", "twoways")) {
  effect <- match.arg(effect); cl <- match.call()
  if (effect == "twoways" && is.null(time)) stop_input("Two way effects need time.")
  mf <- stats::model.frame(formula, data)
  keep <- as.integer(rownames(mf)); if (anyNA(keep)) keep <- seq_len(nrow(mf))
  g <- factor(data[[id]][keep])
  y <- stats::model.response(mf)
  X <- stats::model.matrix(formula, mf); X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  if (effect == "twoways") {
    tt <- factor(data[[time]][keep]); Tm <- stats::model.matrix(~ tt)[, -1, drop = FALSE]
    X <- cbind(X, Tm)
  }
  dm <- function(v) v - stats::ave(v, g)
  Xw <- apply(X, 2, dm); if (!is.matrix(Xw)) Xw <- matrix(Xw, ncol = ncol(X), dimnames = list(NULL, colnames(X)))
  st <- ols(Xw, dm(y))
  n <- length(y); df <- n - nlevels(g) - ncol(X)
  V <- sum(st$res^2) / df * st$XtXi; dimnames(V) <- list(colnames(X), colnames(X)); b <- st$b; names(b) <- colnames(X)
  k <- if (effect == "twoways") setdiff(colnames(X), colnames(Tm)) else colnames(X)
  new_regcomp_lm(b[k], V[k, k, drop = FALSE], st$res, df, sprintf("within estimator (%s effects)", effect), cl)
}

#' Hausman test
#'
#' Compares a consistent estimator (fixed effects, instrumental variables)
#' with an estimator that is efficient under the null hypothesis (random
#' effects, least squares) over their common coefficients:
#' `(b1 - b0)' (V1 - V0)^- (b1 - b0)`, chi squared with the rank of
#' `V1 - V0` degrees of freedom. Equal to `plm::phtest()` given the same
#' fits. Works with regcompR fits (the first parameter's coefficients, by
#' term), [panel_within()] and [iv_regcomp()] results and other fits with
#' `coef()` and `vcov()` methods.
#'
#' @param consistent,efficient Fitted models.
#' @return A one row data frame with statistic, df and p value.
#' @export
#' @examples
#' set.seed(1)
#' d <- data.frame(id = factor(rep(1:50, each = 5)), x = rnorm(250))
#' d$y <- 1 + 0.5 * d$x + rnorm(50)[d$id] + rnorm(250)
#' fe <- panel_within(y ~ x, data = d, id = "id")
#' re <- regcomp(y ~ x + (1 | id), family = "normal", data = d)
#' hausman_test(fe, re)
hausman_test <- function(consistent, efficient) {
  get <- function(o) {
    if (inherits(o, "regcomp_fit")) {
      k <- primary_par(o$spec$fam, regression_pars(o$spec$fam)); i <- o$spec$beta_idx[[k]]
      nm <- colnames(o$spec$preds[[k]]$X)
      return(list(b = stats::setNames(unname(o$coefficients[i]), nm),
                  V = matrix(stats::vcov(o)[i, i], length(i), dimnames = list(nm, nm))))
    }
    list(b = stats::coef(o), V = as.matrix(stats::vcov(o)))
  }
  a <- get(consistent); e <- get(efficient)
  k <- setdiff(intersect(names(a$b), names(e$b)), "(Intercept)")
  if (!length(k)) stop_input("The two fits share no coefficients.")
  d <- a$b[k] - e$b[k]; D <- a$V[k, k, drop = FALSE] - e$V[k, k, drop = FALSE]
  ev <- eigen((D + t(D)) / 2, symmetric = TRUE); keep <- ev$values > 1e-10 * max(abs(ev$values))
  s <- sum((crossprod(ev$vectors[, keep, drop = FALSE], d))^2 / ev$values[keep])
  data.frame(statistic = s, df = sum(keep), p_value = stats::pchisq(s, sum(keep), lower.tail = FALSE))
}
