# Link functions map a distribution parameter onto the real line, where its
# predictor lives. Each link carries its inverse, the derivative of the
# inverse, the parameter domain it is valid for and its Stan inverse.

new_link <- function(name, linkfun, linkinv, mu_eta, valid, stan_inv,
                     domain) {
  structure(list(name = name, linkfun = linkfun, linkinv = linkinv,
                 mu_eta = mu_eta, valid = valid, stan_inv = stan_inv,
                 domain = domain), class = "regcomp_link")
}

link_table <- list(
  identity = new_link("identity", function(m) m, function(e) e,
                      function(e) rep(1, length(e)), function(m) is.finite(m),
                      function(x) x, "real"),
  log = new_link("log", log, function(e) pmax(exp(e), .Machine$double.xmin),
                 function(e) pmax(exp(e), .Machine$double.xmin), function(m) m > 0,
                 function(x) sprintf("exp(%s)", x), "positive"),
  softplus = new_link("softplus", function(m) ifelse(m > 30, m, log(expm1(m))),
                      function(e) pmax(log1pexp(e), .Machine$double.xmin),
                      stats::plogis, function(m) m > 0,
                      function(x) sprintf("log1p_exp(%s)", x), "positive"),
  sqrt = new_link("sqrt", sqrt, function(e) e^2, function(e) 2 * e,
                  function(m) m > 0, function(x) sprintf("square(%s)", x), "positive"),
  inverse = new_link("inverse", function(m) 1 / m, function(e) 1 / e,
                     function(e) -1 / e^2, function(m) m > 0,
                     function(x) sprintf("inv(%s)", x), "positive"),
  logit = new_link("logit", stats::qlogis,
                   function(e) pmin(pmax(stats::plogis(e), 1e-300), 1 - 1e-16),
                   stats::dlogis, function(m) m > 0 & m < 1,
                   function(x) sprintf("inv_logit(%s)", x), "unit"),
  probit = new_link("probit", stats::qnorm,
                    function(e) pmin(pmax(stats::pnorm(e), 1e-300), 1 - 1e-16),
                    stats::dnorm, function(m) m > 0 & m < 1,
                    function(x) sprintf("Phi(%s)", x), "unit"),
  cloglog = new_link("cloglog", function(m) log(-log1p(-m)),
                     function(e) pmin(pmax(-expm1(-exp(e)), 1e-300), 1 - 1e-16),
                     function(e) exp(e - exp(e)), function(m) m > 0 & m < 1,
                     function(x) sprintf("inv_cloglog(%s)", x), "unit"),
  cauchit = new_link("cauchit", stats::qcauchy,
                     function(e) pmin(pmax(stats::pcauchy(e), 1e-300), 1 - 1e-16),
                     stats::dcauchy, function(m) m > 0 & m < 1,
                     function(x) sprintf("(atan(%s) / pi() + 0.5)", x), "unit"),
  loglog = new_link("loglog", function(m) -log(-log(m)),
                    function(e) pmin(pmax(exp(-exp(-e)), 1e-300), 1 - 1e-16),
                    function(e) exp(-e - exp(-e)), function(m) m > 0 & m < 1,
                    function(x) sprintf("exp(-exp(-(%s)))", x), "unit"),
  tan_half = new_link("tan_half", function(m) tan(m / 2), function(e) 2 * atan(e),
                      function(e) 2 / (1 + e^2), function(m) m > -pi & m < pi,
                      function(x) sprintf("(2 * atan(%s))", x), "circular"),
  fisherz = new_link("fisherz", atanh, function(e) pmin(pmax(tanh(e), -1 + 1e-15), 1 - 1e-15),
                     function(e) 1 - tanh(e)^2, function(m) m > -1 & m < 1,
                     function(x) sprintf("tanh(%s)", x), "symmetric_unit"),
  logshift1 = new_link("logshift1", function(m) log(m - 1), function(e) 1 + exp(e),
                       function(e) exp(e), function(m) m > 1,
                       function(x) sprintf("(1 + exp(%s))", x), "greater_than_one"),
  logshift2 = new_link("logshift2", function(m) log(m - 2), function(e) 2 + exp(e),
                       function(e) exp(e), function(m) m > 2,
                       function(x) sprintf("(2 + exp(%s))", x), "greater_than_two")
)

#' Link functions
#'
#' Returns a link object by name. Links available: `identity`, `log`,
#' `softplus`, `sqrt`, `inverse`, `logit`, `probit`, `cloglog`, `cauchit`,
#' `loglog`, `tan_half` (for circular locations), `logshift1` and
#' `logshift2` (parameters constrained above 1 or 2).
#'
#' @param name Link name.
#' @return An object of class `regcomp_link` with elements `linkfun`,
#'   `linkinv`, `mu_eta` and `valid`.
#' @export
#' @examples
#' l <- link("logit")
#' l$linkinv(0)
#' links()
link <- function(name) {
  if (inherits(name, "regcomp_link")) return(name)
  if (!is.character(name) || length(name) != 1L || !name %in% names(link_table))
    stop_input(sprintf("Unknown link '%s'. Available: %s.", format(name),
                       paste(names(link_table), collapse = ", ")))
  link_table[[name]]
}

#' @rdname link
#' @export
links <- function() names(link_table)

#' @export
print.regcomp_link <- function(x, ...) {
  cat("<regcomp_link>", x$name, "for", x$domain, "parameters\n")
  invisible(x)
}

# Default links by parameter domain.
default_link_for <- function(domain) {
  switch(domain, real = "identity", positive = "log", unit = "logit",
         circular = "tan_half", greater_than_one = "logshift1",
         greater_than_two = "logshift2", nonneg = "log", symmetric_unit = "fisherz", "identity")
}

allowed_links_for <- function(domain) {
  switch(domain,
         real = c("identity"),
         positive = c("log", "softplus", "sqrt", "inverse", "identity"),
         nonneg = c("log", "softplus", "identity"),
         unit = c("logit", "probit", "cloglog", "cauchit", "loglog", "identity"),
         circular = c("tan_half", "identity"),
         greater_than_one = c("logshift1", "identity"),
         greater_than_two = c("logshift2", "identity"),
         symmetric_unit = c("fisherz", "identity"),
         "identity")
}
