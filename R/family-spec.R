# Specification objects: ParameterSpec, SupportSpec, CapabilitySpec and the
# family (DistributionSpec) constructor used by every registered law.

domain_bounds <- list(
  real = c(-Inf, Inf), positive = c(0, Inf), nonneg = c(0, Inf),
  unit = c(0, 1), circular = c(-pi, pi), greater_than_one = c(1, Inf),
  greater_than_two = c(2, Inf), integer = c(-Inf, Inf),
  positive_integer = c(1, Inf), nonneg_integer = c(0, Inf), unit_closed = c(0, 1), symmetric_unit = c(-1, 1),
  any = c(-Inf, Inf)
)

#' Parameter specification
#'
#' Describes one parameter of a distribution: its domain, links, role and
#' how it enters regression. Used when defining distributions with
#' `custom_family()`.
#'
#' @param name Parameter name used in formulas and `distribution()`.
#' @param domain One of `"real"`, `"positive"`, `"nonneg"`, `"unit"`,
#'   `"circular"`, `"greater_than_one"`, `"greater_than_two"`,
#'   `"positive_integer"`, `"nonneg_integer"`, `"integer"`, `"unit_closed"` or
#'   `"any"`.
#' @param role Free text role such as `"location"`, `"scale"`, `"shape"`,
#'   `"dispersion"`, `"rate"`, `"skewness"`, `"tail"`, `"inflation"` or
#'   `"mixing"`.
#' @param link Default link name.
#' @param links Allowed link names.
#' @param symbol Mathematical symbol for documentation.
#' @param is_mean `TRUE` only when the parameter equals the expectation of the
#'   response for every parameter value. Never inferred from the name.
#' @param offset_allowed Whether an offset may enter this parameter's
#'   predictor.
#' @param regression Whether the parameter may have a predictor. Fixed data
#'   arguments such as binomial trials use `FALSE` together with `data = TRUE`.
#' @param data `TRUE` for data arguments (for example the number of trials).
#' @param interpretation Text shown in documentation.
#' @param lower,upper Optional numeric bounds overriding the domain.
#' @param default Optional value used by `distribution()` when omitted.
#' @param estimate For a non regression parameter, `"ordered"` or `"free"` to
#'   estimate it as a vector whose length comes from the response (ordinal
#'   thresholds); `NULL` otherwise.
#' @return An object of class `regcomp_par`.
#' @export
#' @examples
#' par_spec("mu", "positive", role = "location", is_mean = TRUE)
par_spec <- function(name, domain = "real", role = "shape", link = NULL,
                     links = NULL, symbol = name, is_mean = FALSE,
                     offset_allowed = FALSE, regression = TRUE, data = FALSE,
                     interpretation = "", lower = NULL, upper = NULL,
                     default = NULL, estimate = NULL) {
  if (!domain %in% names(domain_bounds))
    stop_input(sprintf("Unknown parameter domain '%s'.", domain))
  if (is.null(link)) link <- default_link_for(domain)
  if (is.null(links)) links <- unique(c(link, allowed_links_for(domain)))
  b <- domain_bounds[[domain]]
  structure(list(
    name = name, symbol = symbol, domain = domain, role = role, link = link,
    links = links, is_mean = is_mean, offset_allowed = offset_allowed,
    regression = regression && !data, data = data,
    interpretation = interpretation,
    lower = if (is.null(lower)) b[1] else lower,
    upper = if (is.null(upper)) b[2] else upper,
    default = default, estimate = estimate
  ), class = "regcomp_par")
}

#' @export
print.regcomp_par <- function(x, ...) {
  cat(sprintf("<parameter %s> domain %s, role %s, link %s%s\n", x$name,
              x$domain, x$role, x$link, if (x$is_mean) ", equals the mean" else ""))
  invisible(x)
}

# Check parameter values against their domains; returns a logical vector of
# valid positions (after recycling).
check_par_domain <- function(spec, value) {
  ok <- !is.na(value)
  d <- spec$domain
  # Only the "any" domain admits infinite values.
  if (d != "any") ok <- ok & is.finite(value)
  lo <- spec$lower %||% -Inf; hi <- spec$upper %||% Inf
  closed <- d %in% c("nonneg", "unit_closed", "positive_integer", "nonneg_integer", "integer")
  if (d == "circular") return(ok)
  ok <- ok & (if (closed) value >= lo & value <= hi else (value > lo | lo == -Inf) & (value < hi | hi == Inf))
  if (d %in% c("positive_integer", "nonneg_integer", "integer")) ok <- ok & abs(value - round(value)) < 1e-8
  ok & !is.na(ok)
}

valid_statuses <- c("metadata_only", "experimental", "implemented", "validated",
                    "production", "roadmap")

# Constructor for a registered family. Function slots take `pars`, a named
# list of recycled numeric vectors (matrices for vector valued families).
family_def <- function(id, display_name, category, pars,
                       lpdf, cdf = NULL, lccdf = NULL, lcdf = NULL, quantile = NULL,
                       rng = NULL, mean = NULL, variance = NULL,
                       skewness = NULL, kurtosis = NULL,
                       support = NULL, discrete = FALSE,
                       variate_type = "scalar", aliases = character(),
                       parameterization = "default",
                       formula = "", mean_formula = "", variance_formula = "",
                       moment_conditions = "", references = character(),
                       reference = NULL, stan = NULL, start = NULL,
                       status = "implemented", numerical_notes = character(),
                       identifiability_notes = character(),
                       regression_notes = character(),
                       offset_semantics = "", limits = NULL,
                       check = NULL, operators_allowed = NULL,
                       pars_fun = NULL, response_check = NULL,
                       numeric_cdf = TRUE, numeric_moments = TRUE,
                       tmb = NULL, bugs = NULL, extra = list()) {
  stopifnot(status %in% valid_statuses)
  names(pars) <- vapply(pars, `[[`, "", "name")
  if (is.null(support)) {
    support <- if (discrete) function(p) list(lower = 0, upper = Inf)
               else function(p) list(lower = -Inf, upper = Inf)
  }
  if (is.null(operators_allowed)) {
    operators_allowed <- c("truncate", "censor", "mixture", "transform")
    if (discrete) operators_allowed <- c(operators_allowed, "inflate", "hurdle",
                                         "convolution", "compound")
    else operators_allowed <- c(operators_allowed, "inflate", "hurdle",
                                "discretize", "round", "convolution", "splice")
  }
  structure(list(
    id = id, display_name = display_name, aliases = aliases,
    category = category, variate_type = variate_type,
    parameterization = parameterization, discrete = discrete,
    pars = pars, lpdf = lpdf, cdf = cdf, lccdf = lccdf, lcdf = lcdf, quantile = quantile,
    rng = rng, mean = mean, variance = variance, skewness = skewness,
    kurtosis = kurtosis, support = support, formula = formula,
    mean_formula = mean_formula, variance_formula = variance_formula,
    moment_conditions = moment_conditions, references = references,
    reference = reference, stan = stan, start = start, status = status,
    numerical_notes = numerical_notes,
    identifiability_notes = identifiability_notes,
    regression_notes = regression_notes, offset_semantics = offset_semantics,
    limits = limits, check = check, operators_allowed = operators_allowed,
    pars_fun = pars_fun, response_check = response_check,
    numeric_cdf = numeric_cdf, numeric_moments = numeric_moments,
    tmb = tmb, bugs = bugs, extra = extra, operator = NULL
  ), class = "regcomp_family")
}

#' @export
print.regcomp_family <- function(x, ...) {
  cat(sprintf("<regcomp_family %s> %s\n", x$id, x$display_name))
  cat(sprintf("  category: %s, variate: %s, %s, status: %s\n", x$category,
              x$variate_type, if (x$discrete) "discrete" else "continuous", x$status))
  if (nzchar(x$parameterization) && x$parameterization != "default")
    cat("  parameterization:", x$parameterization, "\n")
  cat("  parameters:\n")
  for (p in x$pars) cat(sprintf("    %-8s %-17s link %-9s %s\n", p$name, p$domain,
                                 if (p$data) "(data)" else p$link,
                                 if (p$is_mean) "[mean]" else ""))
  if (nzchar(x$formula)) cat("  formula:", x$formula, "\n")
  if (!is.null(x$operator)) cat("  built by operator:", x$operator$type, "\n")
  invisible(x)
}

# Names of parameters that may carry a predictor.
regression_pars <- function(fam) {
  names(Filter(function(p) isTRUE(p$regression), fam$pars))
}

data_pars <- function(fam) names(Filter(function(p) isTRUE(p$data), fam$pars))
