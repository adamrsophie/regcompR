#' Capability queries
#'
#' `capabilities()` reports, for a family or distribution, how each operation
#' is provided: `"analytic"` (closed form or exact algorithm), `"numeric"`
#' (declared numerical fallback with tolerance checks) or `"unavailable"`.
#' It also reports whether the family can be used as a regression response
#' and which fitting engines support it. `supports()` answers a single
#' question with `TRUE` or `FALSE`.
#'
#' @param x A family, distribution, identifier or fitted model.
#' @param what Operation name: `"log_prob"`, `"cdf"`, `"quantile"`, `"rng"`,
#'   `"mean"`, `"variance"`, `"skewness"`, `"kurtosis"`, `"score"`,
#'   `"hessian"` or `"regression"`.
#' @param backend Engine name: `"ml"`, `"tmb"`, `"stan"` (also `"bayes"`)
#'   or `"bugs"` (also `"jags"`, `"openbugs"`).
#' @return `capabilities()`: a named list. `supports()`: a logical scalar.
#' @export
#' @examples
#' capabilities("conway_maxwell_poisson")
#' supports("cauchy", "mean")
#' supports("poisson", backend = "stan")
capabilities <- function(x) {
  fam <- capability_family(x)
  scalar <- fam$variate_type == "scalar" || scalar_vecpar(fam)
  op <- function(slot, fallback) {
    if (!is.null(fam[[slot]])) "analytic" else if (fallback) "numeric" else "unavailable"
  }
  # A numerical CDF needs a support bounded below for discrete laws.
  lb <- tryCatch({ dp <- lapply(fam$pars, function(q) q$default %||% 1); is.finite(fam$support(dp)$lower[1]) }, error = function(e) TRUE)
  cdf_ok <- !is.null(fam$cdf) || (scalar && isTRUE(fam$numeric_cdf) && (!fam$discrete || isTRUE(lb)))
  list(
    log_prob = "analytic",
    cdf = if (!scalar) "unavailable" else op("cdf", cdf_ok),
    quantile = if (!scalar) "unavailable" else op("quantile", cdf_ok),
    rng = op("rng", scalar && cdf_ok),
    mean = op("mean", scalar && isTRUE(fam$numeric_moments)),
    variance = op("variance", scalar && isTRUE(fam$numeric_moments)),
    skewness = op("skewness", scalar && isTRUE(fam$numeric_moments)),
    kurtosis = op("kurtosis", scalar && isTRUE(fam$numeric_moments)),
    score = if (scalar) "numeric" else "unavailable",
    hessian = if (scalar) "numeric" else "unavailable",
    regression = length(regression_pars(fam)) > 0 &&
      !isTRUE(fam$extra$no_regression),
    ml = length(regression_pars(fam)) > 0 && !isTRUE(fam$extra$no_regression),
    stan = if (is.null(fam$stan) || !fam$variate_type %in% c("scalar", "categorical") && is.null(fam$stan$vector)) "unavailable"
           else if (isTRUE(fam$stan$native)) "native" else "generated",
    tmb = if (is.null(tmb_family_key(fam))) "unavailable" else "available",
    bugs = if (is.null(fam$bugs)) "unavailable"
           else if (!is.null(fam$bugs$native)) "native" else "zeros_trick",
    status = fam$status
  )
}

capability_family <- function(x) {
  if (inherits(x, "regcomp_dist")) return(x$family)
  if (inherits(x, "regcomp_fit")) return(x$family)
  get_family(x)
}

#' @rdname capabilities
#' @export
supports <- function(x, what = NULL, backend = NULL) {
  cap <- capabilities(x)
  if (!is.null(backend)) {
    backend <- match.arg(backend, c("ml", "tmb", "stan", "bayes", "bugs", "jags", "openbugs"))
    if (backend %in% c("jags", "openbugs")) backend <- "bugs"
    if (backend == "bayes") backend <- "stan"
    v <- cap[[backend]]
    return(if (is.logical(v)) v else v != "unavailable")
  }
  if (is.null(what)) stop_input("Give either 'what' or 'backend'.")
  v <- cap[[what]]
  if (is.null(v)) stop_input(sprintf("Unknown capability '%s'.", what))
  if (is.logical(v)) v else v != "unavailable"
}
