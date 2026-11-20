# Huber, stochastic frontier, independent multinomial probit and partial
# proportional odds laws.

huber_const <- function(k) sqrt(2 * pi) * (2 * stats::pnorm(k) - 1) + 2 * exp(-k^2 / 2) / k
huber_cdf_z <- function(z, k) {
  C <- huber_const(k); A <- exp(-k^2 / 2) / (k * C)
  ifelse(z < -k, exp(k * z + k^2 / 2) / (k * C),
         ifelse(z > k, 1 - exp(-k * z + k^2 / 2) / (k * C), A + sqrt(2 * pi) * (stats::pnorm(z) - stats::pnorm(-k)) / C))
}
huber_q_z <- function(pr, k) {
  C <- huber_const(k); A <- exp(-k^2 / 2) / (k * C)
  ifelse(pr < A, (log(pr * k * C) - k^2 / 2) / k,
         ifelse(pr > 1 - A, -(log((1 - pr) * k * C) - k^2 / 2) / k,
                stats::qnorm(pmin(pmax(stats::pnorm(-k) + (pr - A) * C / sqrt(2 * pi), 0), 1))))
}

register_family(family_def(
  "huber", "Huber (least favorable)", rl, aliases = c("huber_least_favorable"),
  pars = list(loc_par("mu", TRUE), scale_par("sigma"),
              P("k", "positive", role = "tuning", regression = FALSE, default = 1.345,
                interpretation = "Tuning constant: normal center within k scale units, exponential tails beyond")),
  formula = "f(y) = exp(-rho_k(z)) / (sigma C(k)), z = (y - mu)/sigma, rho_k(z) = z^2/2 for |z| <= k, k |z| - k^2/2 otherwise, C(k) = sqrt(2 pi)(2 Phi(k) - 1) + 2 exp(-k^2/2)/k",
  mean_formula = "mu", variance_formula = "sigma^2 [sqrt(2 pi)(2 Phi(k) - 1 - 2 k phi(k)) + 2 exp(-k^2/2)(k + 2/k + 2/k^3)] / C(k)",
  regression_notes = "With sigma fixed, the maximum likelihood estimate of mu is the Huber M estimate (MASS::rlm with the same scale); k is fixed (default 1.345, 95% efficiency at the normal) and can be changed with fixed = list(k = ).",
  lpdf = function(y, p) { z <- abs(y - p$mu) / p$sigma
    ifelse(z <= p$k, -z^2 / 2, p$k^2 / 2 - p$k * z) - log(p$sigma) - log(huber_const(p$k)) },
  cdf = function(q, p) huber_cdf_z((q - p$mu) / p$sigma, p$k),
  quantile = function(pr, p) p$mu + p$sigma * huber_q_z(pr, p$k),
  rng = function(n, p) p$mu + p$sigma * huber_q_z(stats::runif(n), p$k),
  mean = function(p) p$mu,
  variance = function(p) { k <- p$k
    p$sigma^2 * (sqrt(2 * pi) * (2 * stats::pnorm(k) - 1 - 2 * k * stats::dnorm(k)) + 2 * exp(-k^2 / 2) * (k + 2 / k + 2 / k^3)) / huber_const(k) },
  references = c("Huber (1964), Robust estimation of a location parameter. Annals of Mathematical Statistics, 35, 73-101.",
                 "Huber and Ronchetti (2009), Robust Statistics, 2nd ed., sec. 4.5. Wiley."),
  reference = ref_impl("definition", function(y, p) { z <- abs(y - p$mu) / p$sigma; k <- p$k
    log(ifelse(z <= k, exp(-z^2 / 2), exp(k^2 / 2 - k * z)) / (p$sigma * stats::integrate(function(u)
      ifelse(abs(u) <= k[1], exp(-u^2 / 2), exp(k[1]^2 / 2 - k[1] * abs(u))), -Inf, Inf)$value)) }),
  stan = sfun("huber", c("mu", "sigma", "k"),
    "    real z = abs(y - mu) / sigma;\n    real c = sqrt(2 * pi()) * (2 * Phi(k) - 1) + 2 * exp(-k^2 / 2) / k;\n    return (z <= k ? -0.5 * z^2 : k^2 / 2 - k * z) - log(sigma) - log(c);"),
  start = function(y) list(mu = stats::median(y), sigma = max(stats::mad(y), 1e-3)),
  status = "validated"
))

register_family(family_def(
  "stochastic_frontier", "Stochastic frontier (normal minus half normal)", rl,
  aliases = c("normal_half_normal_frontier", "sfa"),
  pars = list(loc_par("mu"), P("sigma_v", "positive", role = "scale", interpretation = "SD of the symmetric noise v"),
              P("sigma_u", "positive", role = "scale", interpretation = "Scale of the half normal inefficiency u")),
  formula = "y = mu + v - u, v ~ N(0, sigma_v^2), u ~ |N(0, sigma_u^2)|: f(y) = 2/s phi(e/s) Phi(-e lambda/s), e = y - mu, s^2 = sigma_v^2 + sigma_u^2, lambda = sigma_u/sigma_v",
  mean_formula = "mu - sigma_u sqrt(2/pi)", variance_formula = "sigma_v^2 + sigma_u^2 (1 - 2/pi)",
  regression_notes = "A production frontier (Aigner, Lovell and Schmidt 1977): mu is the frontier and u the inefficiency; for a cost frontier model -y. Equal to frontier::sfa() with the same parameters (sigmaSq = s^2, gamma = sigma_u^2 / s^2).",
  lpdf = function(y, p) { s <- sqrt(p$sigma_v^2 + p$sigma_u^2); e <- y - p$mu
    log(2) - log(s) + stats::dnorm(e / s, log = TRUE) + stats::pnorm(-e * p$sigma_u / (p$sigma_v * s), log.p = TRUE) },
  rng = function(n, p) p$mu + stats::rnorm(n, 0, p$sigma_v) - abs(stats::rnorm(n, 0, p$sigma_u)),
  mean = function(p) p$mu - p$sigma_u * sqrt(2 / pi),
  variance = function(p) p$sigma_v^2 + p$sigma_u^2 * (1 - 2 / pi),
  references = c("Aigner, Lovell and Schmidt (1977), Formulation and estimation of stochastic frontier production function models. Journal of Econometrics, 6, 21-37.",
                 "Kumbhakar and Lovell (2000), Stochastic Frontier Analysis. Cambridge University Press."),
  reference = ref_impl("sn", function(y, p) { s <- sqrt(p$sigma_v^2 + p$sigma_u^2)
    sn::dsn(y, xi = p$mu, omega = s, alpha = -p$sigma_u / p$sigma_v, log = TRUE) }),
  stan = sfun("stochastic_frontier", c("mu", "sigma_v", "sigma_u"),
    "    real s = sqrt(sigma_v^2 + sigma_u^2);\n    real e = y - mu;\n    return log2() - log(s) + std_normal_lpdf(e / s) + std_normal_lcdf(-e * sigma_u / (sigma_v * s));"),
  start = function(y) list(mu = stats::quantile(y, 0.75, names = FALSE), sigma_v = max(stats::sd(y) / 2, 1e-3), sigma_u = max(stats::sd(y) / 2, 1e-3)),
  status = "validated"
))

# Gauss-Hermite rule (Golub and Welsch) for integrals against exp(-x^2).
gauss_hermite <- function(n) {
  J <- matrix(0, n, n); off <- sqrt(seq_len(n - 1) / 2)
  J[cbind(1:(n - 1), 2:n)] <- off; J[cbind(2:n, 1:(n - 1))] <- off
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = sqrt(pi) * e$vectors[1, ]^2)
}
gh40 <- gauss_hermite(40)

# Independent multinomial probit: utilities eta_k + e_k with independent
# standard normal errors; P(Y = k) = E[prod_{j != k} Phi(U + eta_k - eta_j)].
mprobit_lpmf <- function(E, y) {
  n <- nrow(E); Ey <- E[cbind(seq_len(n), y)]
  lt <- vapply(seq_along(gh40$x), function(q) {
    L <- stats::pnorm(sqrt(2) * gh40$x[q] + Ey - E, log.p = TRUE)
    L[cbind(seq_len(n), y)] <- 0
    rowSums(L) + log(gh40$w[q] / sqrt(pi)) }, numeric(n))
  lt <- matrix(lt, n); m <- apply(lt, 1, max)
  m + log(rowSums(exp(lt - m)))
}
mprobit_family <- function(K) {
  nm <- paste0("eta", 2:K)
  family_def("multinomial_probit", sprintf("Multinomial probit, independent errors (%d categories)", K), "categorical",
    lapply(nm, function(k) P(k, "real", role = sprintf("utility of category %s versus 1", sub("eta", "", k)), link = "identity")),
    discrete = TRUE, variate_type = "categorical", numeric_moments = FALSE,
    support = function(p) list(lower = 1, upper = K),
    lpdf = function(y, p) mprobit_lpmf(cbind(0, do.call(cbind, lapply(p[nm], rep_len, length(y)))), y),
    extra = list(share_formula = TRUE), status = "validated")
}
register_family(family_def("multinomial_probit", "Multinomial probit, independent errors (nominal)", "categorical",
  list(P("eta", "real", role = "utilities versus category 1 (K - 1 values)", regression = FALSE)),
  discrete = TRUE, variate_type = "categorical", numeric_moments = FALSE,
  formula = "U_k = eta_k + e_k, e_k independent N(0, 1), eta_1 = 0; Y = argmax U_k",
  lpdf = function(y, p) mprobit_lpmf(matrix(c(0, p$eta), length(y), length(p$eta) + 1, byrow = TRUE), y),
  regression_notes = "Independent (identity covariance) errors: the probabilities are one dimensional integrals, computed with 40 point Gauss-Hermite quadrature. Correlated error covariances (MNP) are not available. In regcomp() the family expands to parameters eta2, ..., etaK that share the main formula unless given their own.",
  extra = list(expand = function(K) mprobit_family(K)),
  references = "Hausman and Wise (1978), A conditional probit model for qualitative choice. Econometrica, 46, 403-426.",
  reference = ref_impl("mvtnorm", function(y, p) { e <- c(0, p$eta); K <- length(e)
    vapply(y, function(k) log(mvtnorm::pmvnorm(upper = e[k] - e[-k], sigma = diag(K - 1) + 1, algorithm = mvtnorm::Miwa())[1]), 0) }),
  status = "validated"))

#' Partial proportional odds family
#'
#' A cumulative link model whose thresholds depend on covariates:
#' `P(Y <= k) = F(theta_k + z' gamma_k - mu)`, where `mu` is the
#' proportional (parallel) part of the main formula and `nonparallel` gives
#' the covariates `z` with category specific effects. This is the nominal
#' effects model of `ordinal::clm(nominal = )` and the partial proportional
#' odds model of Peterson and Harrell (1990). In `regcomp()` the family
#' adds parameters `d1`, ..., `d(K-1)` (the threshold shifts `z' gamma_k`),
#' each with the `nonparallel` formula unless given its own.
#'
#' @param nonparallel One sided formula of the covariates with category
#'   specific effects. A covariate should not also appear in the main
#'   formula.
#' @param link `"logit"`, `"probit"`, `"cloglog"`, `"loglog"` or `"cauchit"`.
#' @return A family object for `regcomp()`.
#' @export
#' @examples
#' set.seed(1)
#' d <- data.frame(x = rnorm(400), z = rnorm(400))
#' eta <- 0.8 * d$x
#' u <- runif(400)
#' d$y <- factor(1 + (u > plogis(-1 - eta)) + (u > plogis(1 + 0.5 * d$z - eta)))
#' f <- regcomp(y ~ x, family = partial_odds(~ z), data = d)
#' coef(f)
partial_odds <- function(nonparallel, link = c("logit", "probit", "cloglog", "loglog", "cauchit")) {
  link <- match.arg(link)
  if (!inherits(nonparallel, "formula") || length(nonparallel) != 2) stop_input("nonparallel must be a one sided formula, as in ~ z.")
  base <- get_family(paste0("cumulative_", link))
  base$id <- paste0("partial_odds_", link)
  base$display_name <- sprintf("Partial proportional odds (cumulative %s)", link)
  base$stan <- NULL; base$bugs <- NULL; base$tmb <- NULL
  base$extra$expand <- function(K) {
    dn <- paste0("d", seq_len(K - 1))
    f <- base
    f$pars <- c(f$pars, stats::setNames(lapply(dn, function(k) P(k, "real", role = sprintf("shift of threshold %s", sub("d", "", k)), link = "identity")), dn))
    f$extra <- list(no_intercept = c("mu", "scale", dn), fixed_unless_formula = "scale",
                    default_formula = stats::setNames(rep(list(nonparallel), length(dn)), dn))
    lF <- ord_logF[[link]]
    f$check <- NULL
    f$lpdf <- function(y, p) { n <- length(y); s <- rep_len(p$scale %||% 1, n)
      C <- matrix(p$cuts, n, K - 1, byrow = TRUE) + do.call(cbind, lapply(p[dn], rep_len, n))
      C <- cbind(-Inf, C, Inf)
      b <- (C[cbind(seq_len(n), y + 1)] - p$mu) / s; a <- (C[cbind(seq_len(n), y)] - p$mu) / s
      out <- suppressWarnings(ifelse(a + b > 0, log_diff_exp(lF(a, FALSE), lF(b, FALSE)), log_diff_exp(lF(b), lF(a))))
      out[!(b > a) | !is.finite(out)] <- -Inf
      out }
    f
  }
  base
}
