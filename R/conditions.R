# Typed conditions. Every failure the package raises on purpose has a class
# so callers can catch it precisely with tryCatch().

regcomp_condition <- function(message, class, call = NULL, data = list()) {
  structure(
    class = c(class, "regcomp_condition", "error", "condition"),
    c(list(message = message, call = call), data)
  )
}

stop_unsupported <- function(what, family = NULL, call = sys.call(-1)) {
  msg <- if (is.null(family)) what else
    sprintf("%s is not available for the '%s' distribution.", what, family)
  stop(regcomp_condition(msg, "regcomp_unsupported", call,
                         list(operation = what, family = family)))
}

stop_domain <- function(message, call = sys.call(-1)) {
  stop(regcomp_condition(message, "regcomp_domain_error", call))
}

stop_numerical <- function(message, call = sys.call(-1), data = list()) {
  stop(regcomp_condition(message, "regcomp_numerical_failure", call, data))
}

stop_input <- function(message, call = sys.call(-1)) {
  stop(regcomp_condition(message, "regcomp_input_error", call))
}

warn_undefined_moment <- function(what, n) {
  w <- structure(
    class = c("regcomp_undefined_moment", "warning", "condition"),
    list(message = sprintf(
      "The %s does not exist for %d parameter value%s; NA returned.",
      what, n, if (n == 1) "" else "s"), call = NULL))
  warning(w)
}

warn_identifiability <- function(message) {
  w <- structure(class = c("regcomp_identifiability", "warning", "condition"),
                 list(message = message, call = NULL))
  warning(w)
}

#' Conditions raised by regcompR
#'
#' All errors raised on purpose by the package inherit from
#' `regcomp_condition` and one of the classes below, so they can be caught
#' selectively with [tryCatch()].
#'
#' * `regcomp_unsupported`: an operation (for example a quantile function or a
#'   backend) is not available for a distribution or model. The package never
#'   substitutes an unjustified approximation.
#' * `regcomp_domain_error`: a parameter or response value lies outside its
#'   domain or support.
#' * `regcomp_numerical_failure`: a series, root finder or integral did not
#'   reach its declared tolerance.
#' * `regcomp_input_error`: malformed arguments.
#'
#' Warnings use the classes `regcomp_undefined_moment` (a moment that does
#' not exist for the given parameters; `NA` is returned) and
#' `regcomp_identifiability`.
#'
#' @name regcomp-conditions
#' @examples
#' res <- tryCatch(mean(distribution("cauchy", mu = 0, sigma = 1)),
#'                 regcomp_undefined_moment = function(w) "undefined")
#' res
NULL
