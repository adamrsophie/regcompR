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
  if (!is.null(seed)) local_seed(seed)
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
  if (!fam$variate_type %in% c("scalar", "categorical")) stop_unsupported("PIT values for a vector valued response", fam$id)
  if (fam$variate_type == "scalar") pars <- recycle_pars(pars, n)
  lo <- y; up <- y
  if (!is.null(cmat)) {
    lo <- pmax(cmat[, "lower"], cmat[, "tlower"]); up <- pmin(cmat[, "upper"], cmat[, "tupper"])
  }
  Fup <- ifelse(is.infinite(up) & up > 0, 1, family_cdf(fam, ifelse(is.finite(up), up, 0), pars))
  # Lower end: the left limit F(lo-) for a discrete law or an exact value
  # (so an atom, and only an atom, is randomized over); F(lo) for the lower
  # bound of a continuous censoring interval. A continuous law without
  # atoms has F(y-) = F(y) and uses no random numbers.
  left <- fam$discrete | lo == up
  lq <- ifelse(is.finite(lo), lo, 0)
  Flo <- ifelse(is.infinite(lo) & lo < 0, 0, ifelse(left, family_cdf_left(fam, lq, pars), family_cdf(fam, lq, pars)))
  if (!is.null(cmat)) {
    tl <- cmat[, "tlower"]; tu <- cmat[, "tupper"]
    tq <- ifelse(is.finite(tl), tl, 0)
    Ftl <- ifelse(is.finite(tl), if (fam$discrete) family_cdf_left(fam, tq, pars) else family_cdf(fam, tq, pars), 0)
    Ftu <- ifelse(is.finite(tu), family_cdf(fam, ifelse(is.finite(tu), tu, 0), pars), 1)
    Fup <- (Fup - Ftl) / (Ftu - Ftl); Flo <- (Flo - Ftl) / (Ftu - Ftl)
  }
  # Draw uniforms only where there is a jump to randomize over.
  jmp <- which(!is.na(Fup) & !is.na(Flo) & Fup > Flo)
  u <- Fup
  if (length(jmp)) u[jmp] <- Flo[jmp] + stats::runif(length(jmp)) * (Fup[jmp] - Flo[jmp])
  pmin(pmax(u, 1e-12), 1 - 1e-12)
}
