# The registry holds every family definition. Families are added at build
# time by the family-*.R files.

.registry <- new.env(parent = emptyenv())
.aliases <- new.env(parent = emptyenv())

register_family <- function(fam) {
  if (!is.null(fam$stan) && is.null(fam$stan$constraint) && !is.null(stan_joint_constraints[[fam$id]]))
    fam$stan$constraint <- stan_joint_constraints[[fam$id]]
  assign(fam$id, fam, envir = .registry)
  for (a in fam$aliases) assign(a, fam$id, envir = .aliases)
  invisible(fam)
}

resolve_id <- function(id) {
  if (exists(id, envir = .registry, inherits = FALSE)) return(id)
  if (exists(id, envir = .aliases, inherits = FALSE)) return(get(id, envir = .aliases))
  NULL
}

#' Look up a registered family
#'
#' `get_family()` returns the family specification for a registered
#' distribution, by identifier or alias. `families()` lists identifiers,
#' optionally filtered by category or status.
#'
#' @param id Identifier or alias, for example `"negative_binomial_nb2"` or
#'   `"nb2"`.
#' @param parameterization Optional parameterization name. When given, the
#'   identifier `paste(id, parameterization, sep = "_")` is tried as well, so
#'   `get_family("negative_binomial", "nb2")` works.
#' @param category Optional category filter.
#' @param status Optional status filter.
#' @return A `regcomp_family` object, or a character vector for `families()`.
#' @export
#' @examples
#' get_family("poisson")
#' get_family("negative_binomial", "nb2")
#' head(families(category = "count"))
get_family <- function(id, parameterization = NULL) {
  if (inherits(id, "regcomp_family")) return(id)
  if (!is.character(id) || length(id) != 1L)
    stop_input("A family identifier must be a single string.")
  cand <- c(if (!is.null(parameterization)) paste(id, parameterization, sep = "_"),
            id)
  for (cc in cand) {
    rid <- resolve_id(cc)
    if (!is.null(rid)) return(get(rid, envir = .registry))
  }
  stop_input(sprintf("Unknown distribution '%s'. See families().", id))
}

#' @rdname get_family
#' @export
families <- function(category = NULL, status = NULL) {
  ids <- sort(ls(.registry))
  keep <- vapply(ids, function(i) {
    f <- get(i, envir = .registry)
    (is.null(category) || f$category %in% category) &&
      (is.null(status) || f$status %in% status)
  }, TRUE)
  ids[keep]
}

#' Registry table
#'
#' One row per registered distribution with its category, variate type,
#' parameters, validation status and the methods available in R, the
#' maximum likelihood engines and Stan. Capabilities are computed from the
#' specification, never asserted separately.
#'
#' @param category Optional category filter.
#' @return A data frame.
#' @export
#' @examples
#' tab <- registry_table()
#' table(tab$category)
#' head(tab[, c("id", "status", "cdf", "quantile", "stan")])
registry_table <- function(category = NULL) {
  ids <- families(category)
  rows <- lapply(ids, function(i) {
    f <- get_family(i)
    cap <- capabilities(f)
    data.frame(
      id = f$id, display_name = f$display_name, category = f$category,
      variate_type = f$variate_type, discrete = f$discrete,
      parameters = paste(names(f$pars), collapse = ", "),
      status = f$status, log_prob = cap$log_prob, cdf = cap$cdf,
      quantile = cap$quantile, rng = cap$rng, mean = cap$mean,
      variance = cap$variance, regression = cap$regression,
      stan = cap$stan, tmb = cap$tmb,
      reference = if (is.null(f$reference)) NA_character_ else f$reference$package,
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
