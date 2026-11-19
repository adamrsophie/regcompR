#' Residuals of a fitted model
#'
#' `type = "quantile"` gives randomized quantile residuals (Dunn and Smyth
#' 1996, <doi:10.1080/10618600.1996.10474708>): the PIT value of each
#' observation under its fitted distribution, randomized uniformly over the
#' probability mass of discrete or censored observations, mapped to the
#' standard normal scale. Under a correct model they are independent
#' standard normal. `type = "pit"` returns the PIT values themselves,
#' `"response"` the difference between the observation and the fitted mean,
#' and `"pearson"` that difference divided by the fitted standard deviation.
#'
#' @param object A fitted model.
#' @param type Residual type.
#' @param seed Optional seed for the randomization.
#' @param ... Unused.
#' @return Numeric vector.
#' @export
residuals.regcomp_fit <- function(object, type = c("quantile", "pit", "response", "pearson"),
                                  seed = NULL, ...) {
  type <- match.arg(type)
  if (!is.null(seed)) set.seed(seed)
  spec <- object$spec
  fd <- fitted_distribution(object)
  y <- spec$y
  if (type %in% c("response", "pearson")) {
    if (!is.null(spec$cmat) && any(spec$cmat[, "lower"] != spec$cmat[, "upper"]))
      stop_unsupported("Response and Pearson residuals for censored observations")
    r <- y - suppressWarnings(mean(fd))
    if (type == "pearson") r <- r / sqrt(suppressWarnings(variance(fd)))
    return(r)
  }
  u <- pit_values(spec$fam, y, fd$pars, spec$cmat)
  if (type == "pit") u else stats::qnorm(u)
}

pit_values <- function(fam, y, pars, cmat = NULL) {
  n <- length(y)
  pars <- recycle_pars(pars, n)
  if (fam$variate_type != "scalar") stop_unsupported("PIT values for a vector valued response", fam$id)
  lo <- y; up <- y
  if (!is.null(cmat)) { lo <- cmat[, "lower"]; up <- cmat[, "upper"] }
  Fup <- ifelse(is.infinite(up) & up > 0, 1, family_cdf(fam, ifelse(is.finite(up), up, 0), pars))
  lo_q <- if (fam$discrete) lo - 1 else lo
  Flo <- ifelse(is.infinite(lo) & lo < 0, 0, family_cdf(fam, ifelse(is.finite(lo_q), lo_q, 0), pars))
  exact_cont <- !fam$discrete & lo == up
  Flo[exact_cont] <- Fup[exact_cont]
  if (!is.null(cmat)) {
    tl <- cmat[, "tlower"]; tu <- cmat[, "tupper"]
    Ftl <- ifelse(is.finite(tl), family_cdf(fam, ifelse(is.finite(tl), if (fam$discrete) tl - 1 else tl, 0), pars), 0)
    Ftu <- ifelse(is.finite(tu), family_cdf(fam, ifelse(is.finite(tu), tu, 0), pars), 1)
    Fup <- (Fup - Ftl) / (Ftu - Ftl); Flo <- (Flo - Ftl) / (Ftu - Ftl)
  }
  u <- Flo + stats::runif(n) * (Fup - Flo)
  pmin(pmax(u, 1e-12), 1 - 1e-12)
}
