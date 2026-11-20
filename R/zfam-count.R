# Count distributions on the non-negative (or positive) integers.

cnt <- "count"
log_mean_par <- function(name = "mu", ...) P(name, "positive", role = "location",
  is_mean = TRUE, offset_allowed = TRUE, ...)

# Poisson ---------------------------------------------------------------
register_family(family_def(
  "poisson", "Poisson", cnt, discrete = TRUE, aliases = c("PO", "pois"),
  pars = list(P("lambda", "positive", role = "rate", is_mean = TRUE,
                offset_allowed = TRUE, symbol = "\\lambda",
                interpretation = "Mean and variance")),
  formula = "p(y) = exp(-lambda) lambda^y / y!",
  mean_formula = "lambda", variance_formula = "lambda",
  lpdf = function(y, p) stats::dpois(y, p$lambda, log = TRUE),
  cdf = function(q, p) stats::ppois(q, p$lambda),
  lccdf = function(q, p) stats::ppois(q, p$lambda, lower.tail = FALSE, log.p = TRUE),
  lcdf = function(q, p) stats::ppois(q, p$lambda, log.p = TRUE),
  quantile = function(pr, p) stats::qpois(pr, p$lambda),
  rng = function(n, p) stats::rpois(n, p$lambda),
  mean = function(p) p$lambda, variance = function(p) p$lambda,
  skewness = function(p) p$lambda^-0.5, kurtosis = function(p) 1 / p$lambda,
  offset_semantics = "log(lambda) = eta + log(exposure); exposure must be positive.",
  references = "Johnson, Kemp and Kotz (2005), Univariate Discrete Distributions, 3rd ed., ch. 4.",
  reference = ref_impl("gamlss.dist", function(y, p) gamlss.dist::dPO(y, mu = p$lambda, log = TRUE),
                       function(q, p) gamlss.dist::pPO(q, mu = p$lambda)),
  stan = stan_native("poisson_lpmf({y} | {lambda})", "poisson_lcdf({y} | {lambda})",
                     "poisson_lccdf({y} | {lambda})", "poisson_rng({lambda})"),
  tmb = "poisson",
  bugs = list(native = "dpois({lambda})",
              loglik = "{y} * log({lambda}) - {lambda} - logfact({y})"),
  start = function(y) list(lambda = max(mean(y), 0.1)),
  status = "validated"
))

# Geometric variants ---------------------------------------------------------
register_family(family_def(
  "geometric_failures", "Geometric (failures before first success)", cnt,
  discrete = TRUE, aliases = c("GEOMo", "geometric"),
  pars = list(P("prob", "unit", role = "probability", symbol = "p",
                interpretation = "Success probability")),
  formula = "p(y) = p (1 - p)^y, y = 0, 1, ...",
  mean_formula = "(1 - p) / p", variance_formula = "(1 - p) / p^2",
  lpdf = function(y, p) stats::dgeom(y, p$prob, log = TRUE),
  cdf = function(q, p) stats::pgeom(q, p$prob),
  quantile = function(pr, p) stats::qgeom(pr, p$prob),
  rng = function(n, p) stats::rgeom(n, p$prob),
  mean = function(p) (1 - p$prob) / p$prob,
  variance = function(p) (1 - p$prob) / p$prob^2,
  references = "Johnson, Kemp and Kotz (2005), ch. 5.2.",
  reference = ref_impl("gamlss.dist", function(y, p) gamlss.dist::dGEOMo(y, mu = p$prob, log = TRUE)),
  stan = stan_generated("
  real rc_geom_lpmf(int y, real p) { return log(p) + y * log1m(p); }
  real rc_geom_lcdf(int y, real p) { return log1m_exp((y + 1) * log1m(p)); }
  real rc_geom_lccdf(int y, real p) { return (y + 1) * log1m(p); }",
                        "rc_geom_lpmf({y} | {prob})", "rc_geom_lcdf({y} | {prob})",
                        "rc_geom_lccdf({y} | {prob})",
                        "neg_binomial_rng(1, {prob} / (1 - {prob}))"),
  start = function(y) list(prob = 1 / (1 + mean(y))),
  status = "validated"
))

register_family(family_def(
  "geometric_trials", "Geometric (trials until first success)", cnt,
  discrete = TRUE,
  pars = list(P("prob", "unit", role = "probability", symbol = "p")),
  support = function(p) list(lower = 1, upper = Inf),
  formula = "p(y) = p (1 - p)^(y - 1), y = 1, 2, ...",
  mean_formula = "1 / p", variance_formula = "(1 - p) / p^2",
  lpdf = function(y, p) stats::dgeom(y - 1, p$prob, log = TRUE),
  cdf = function(q, p) stats::pgeom(q - 1, p$prob),
  quantile = function(pr, p) stats::qgeom(pr, p$prob) + 1,
  rng = function(n, p) stats::rgeom(n, p$prob) + 1,
  mean = function(p) 1 / p$prob, variance = function(p) (1 - p$prob) / p$prob^2,
  references = "Johnson, Kemp and Kotz (2005), ch. 5.2.",
  reference = ref_impl("stats", function(y, p) stats::dgeom(y - 1, p$prob, log = TRUE)),
  stan = stan_generated("
  real rc_geomt_lpmf(int y, real p) { return log(p) + (y - 1) * log1m(p); }
  real rc_geomt_lcdf(int y, real p) { return log1m_exp(y * log1m(p)); }
  real rc_geomt_lccdf(int y, real p) { return y * log1m(p); }",
                        "rc_geomt_lpmf({y} | {prob})", "rc_geomt_lcdf({y} | {prob})",
                        "rc_geomt_lccdf({y} | {prob})",
                        "(1 + neg_binomial_rng(1, {prob} / (1 - {prob})))"),
  start = function(y) list(prob = 1 / max(mean(y), 1.01)),
  status = "validated"
))

register_family(family_def(
  "geometric_mean", "Geometric (mean parameterization, GAMLSS GEOM)", cnt,
  discrete = TRUE, aliases = "GEOM",
  pars = list(log_mean_par()),
  formula = "p(y) = mu^y / (1 + mu)^(y + 1)",
  mean_formula = "mu", variance_formula = "mu + mu^2",
  lpdf = function(y, p) stats::dgeom(y, 1 / (1 + p$mu), log = TRUE),
  cdf = function(q, p) stats::pgeom(q, 1 / (1 + p$mu)),
  quantile = function(pr, p) stats::qgeom(pr, 1 / (1 + p$mu)),
  rng = function(n, p) stats::rgeom(n, 1 / (1 + p$mu)),
  mean = function(p) p$mu, variance = function(p) p$mu + p$mu^2,
  references = "Rigby, Stasinopoulos, Heller and De Bastiani (2019), Distributions for Modeling Location, Scale, and Shape, sec. 22.",
  reference = ref_impl("gamlss.dist", function(y, p) gamlss.dist::dGEOM(y, mu = p$mu, log = TRUE)),
  stan = stan_native("neg_binomial_2_lpmf({y} | {mu}, 1)", "neg_binomial_2_lcdf({y} | {mu}, 1)",
                     "neg_binomial_2_lccdf({y} | {mu}, 1)", "neg_binomial_2_rng({mu}, 1)"),
  start = function(y) list(mu = max(mean(y), 0.1)),
  status = "validated"
))

# Negative binomial ------------------------------------------------------------
nb_size_family <- function(id, name, aliases, sizefun, dispname, formula,
                           varformula, varfun, refun, stan_size, interp,
                           tmb = NULL) {
  register_family(family_def(
    id, name, cnt, discrete = TRUE, aliases = aliases,
    pars = list(log_mean_par(), P(dispname, "positive", role = "dispersion",
                                  interpretation = interp)),
    formula = formula, mean_formula = "mu", variance_formula = varformula,
    lpdf = function(y, p) stats::dnbinom(y, size = sizefun(p), mu = p$mu, log = TRUE),
    cdf = function(q, p) stats::pnbinom(q, size = sizefun(p), mu = p$mu),
    lccdf = function(q, p) stats::pnbinom(q, size = sizefun(p), mu = p$mu,
                                         lower.tail = FALSE, log.p = TRUE),
    lcdf = function(q, p) stats::pnbinom(q, size = sizefun(p), mu = p$mu, log.p = TRUE),
    quantile = function(pr, p) stats::qnbinom(pr, size = sizefun(p), mu = p$mu),
    rng = function(n, p) stats::rnbinom(n, size = sizefun(p), mu = p$mu),
    mean = function(p) p$mu, variance = varfun,
    references = c("Hilbe (2011), Negative Binomial Regression, 2nd ed.",
                   "Brooks et al. (2017), glmmTMB, The R Journal 9(2)."),
    reference = refun,
    offset_semantics = "log(mu) = eta + log(exposure). The dispersion predictor takes no offset.",
    stan = stan_native(sprintf("neg_binomial_2_lpmf({y} | {mu}, %s)", stan_size),
                       sprintf("neg_binomial_2_lcdf({y} | {mu}, %s)", stan_size),
                       sprintf("neg_binomial_2_lccdf({y} | {mu}, %s)", stan_size),
                       sprintf("neg_binomial_2_rng({mu}, %s)", stan_size)),
    limits = "Poisson as the size tends to infinity.",
    bugs = list(loglik = sprintf("loggam({y} + %1$s) - loggam(%1$s) - logfact({y}) + %1$s * log(%1$s / (%1$s + {mu})) + {y} * log({mu} / (%1$s + {mu}))",
                                 paste0("(", gsub("inv\\(\\{sigma\\}\\)", "(1 / {sigma})", stan_size), ")"))),
    tmb = tmb,
    start = function(y) {
      m <- max(mean(y), 0.1); v <- max(stats::var(y), m * 1.1)
      s <- list(mu = m)
      s[[dispname]] <- switch(id,
        negative_binomial_nb2 = m^2 / (v - m),
        negative_binomial_nb1 = (v - m) / m,
        negative_binomial_gamlss_nbi = (v - m) / m^2,
        negative_binomial_gamlss_nbii = (v - m) / m, 1)
      s
    },
    status = "validated"
  ))
}

nb_size_family("negative_binomial_nb2", "Negative binomial (NB2)",
  c("nb2", "negative_binomial", "negbin"), function(p) p$phi, "phi",
  "p(y) = C(y + phi - 1, y) (mu / (mu + phi))^y (phi / (mu + phi))^phi",
  "mu + mu^2 / phi", function(p) p$mu + p$mu^2 / p$phi,
  ref_impl("stats", function(y, p) stats::dnbinom(y, size = p$phi, mu = p$mu, log = TRUE)),
  "{phi}", "Size (inverse dispersion); larger means closer to Poisson", tmb = "nb2")

nb_size_family("negative_binomial_nb1", "Negative binomial (NB1)",
  c("nb1"), function(p) p$mu / p$phi, "phi",
  "NB with size mu / phi (glmmTMB nbinom1)",
  "mu (1 + phi)", function(p) p$mu * (1 + p$phi),
  ref_impl("glmmTMB", function(y, p) stats::dnbinom(y, size = p$mu / p$phi, mu = p$mu, log = TRUE),
           note = "glmmTMB nbinom1 uses Var = mu (1 + phi), i.e. size = mu / phi."),
  "{mu} / {phi}", "Variance to mean ratio minus one", tmb = "nb1")

nb_size_family("negative_binomial_gamlss_nbi", "Negative binomial type I (GAMLSS NBI)",
  c("NBI"), function(p) 1 / p$sigma, "sigma",
  "NB with size 1 / sigma (GAMLSS NBI)",
  "mu + sigma mu^2", function(p) p$mu + p$sigma * p$mu^2,
  ref_impl("gamlss.dist", function(y, p) gamlss.dist::dNBI(y, mu = p$mu, sigma = p$sigma, log = TRUE)),
  "inv({sigma})", "Dispersion; Var = mu + sigma mu^2")

nb_size_family("negative_binomial_gamlss_nbii", "Negative binomial type II (GAMLSS NBII)",
  c("NBII"), function(p) p$mu / p$sigma, "sigma",
  "NB with size mu / sigma (GAMLSS NBII)",
  "mu (1 + sigma)", function(p) p$mu * (1 + p$sigma),
  ref_impl("gamlss.dist", function(y, p) gamlss.dist::dNBII(y, mu = p$mu, sigma = p$sigma, log = TRUE)),
  "{mu} / {sigma}", "Dispersion; Var = mu (1 + sigma)")

register_family(family_def(
  "negative_binomial_family", "Negative binomial family (GAMLSS NBF)", cnt,
  discrete = TRUE, aliases = "NBF",
  pars = list(log_mean_par(), P("sigma", "positive", role = "dispersion"),
              P("nu", "real", role = "power", interpretation = "Variance power")),
  formula = "NB with size mu^(2 - nu) / sigma",
  mean_formula = "mu", variance_formula = "mu + sigma mu^nu",
  lpdf = function(y, p) stats::dnbinom(y, size = p$mu^(2 - p$nu) / p$sigma, mu = p$mu, log = TRUE),
  cdf = function(q, p) stats::pnbinom(q, size = p$mu^(2 - p$nu) / p$sigma, mu = p$mu),
  quantile = function(pr, p) stats::qnbinom(pr, size = p$mu^(2 - p$nu) / p$sigma, mu = p$mu),
  rng = function(n, p) stats::rnbinom(n, size = p$mu^(2 - p$nu) / p$sigma, mu = p$mu),
  mean = function(p) p$mu, variance = function(p) p$mu + p$sigma * p$mu^p$nu,
  references = "Rigby et al. (2019), sec. 22.4.",
  reference = ref_impl("gamlss.dist", function(y, p) gamlss.dist::dNBF(y, mu = p$mu, sigma = p$sigma, nu = p$nu, log = TRUE)),
  stan = stan_native("neg_binomial_2_lpmf({y} | {mu}, pow({mu}, 2 - {nu}) / {sigma})",
                     "neg_binomial_2_lcdf({y} | {mu}, pow({mu}, 2 - {nu}) / {sigma})",
                     "neg_binomial_2_lccdf({y} | {mu}, pow({mu}, 2 - {nu}) / {sigma})",
                     "neg_binomial_2_rng({mu}, pow({mu}, 2 - {nu}) / {sigma})"),
  identifiability_notes = "nu is weakly identified unless mu varies widely across observations.",
  start = function(y) list(mu = max(mean(y), 0.1), sigma = 0.5, nu = 2),
  status = "validated"
))

register_family(family_def(
  "beta_negative_binomial", "Beta negative binomial (GAMLSS BNB)", cnt,
  discrete = TRUE, aliases = "BNB",
  pars = list(log_mean_par(), P("sigma", "positive", role = "dispersion"),
              P("nu", "positive", role = "dispersion")),
  formula = "BNB with n = mu nu / sigma, m = 1 / sigma + 1, k = 1 / nu: log p = lbeta(y + n, m + k) - lbeta(n, m) - lgamma(y + 1) - lgamma(k) + lgamma(y + k)",
  mean_formula = "mu", variance_formula = "k n (m + k - 1) (m + n - 1) / ((m - 2) (m - 1)^2)",
  moment_conditions = "The mean is always finite (m > 1); the variance is finite iff sigma < 1 (m > 2).",
  lpdf = function(y, p) {
    m <- 1 / p$sigma + 1; n <- p$mu * p$nu / p$sigma; k <- 1 / p$nu
    lbeta(y + n, m + k) - lbeta(n, m) - lgamma(y + 1) - lgamma(k) + lgamma(y + k)
  },
  mean = function(p) p$mu + 0 * p$sigma,
  variance = function(p) { m <- 1 / p$sigma + 1; n <- p$mu * p$nu / p$sigma; k <- 1 / p$nu
    ifelse(m > 2, k * n * (m + k - 1) * (m + n - 1) / ((m - 2) * (m - 1)^2), Inf) },
  rng = function(n, p) {
    m <- 1 / p$sigma + 1; nn <- p$mu * p$nu / p$sigma; k <- 1 / p$nu
    pr <- stats::rbeta(n, m, nn)
    stats::rnbinom(n, size = k, prob = pr)
  },
  references = c("Rigby et al. (2019), sec. 22.5.", "Stan Functions Reference, beta_neg_binomial."),
  reference = ref_impl("gamlss.dist", function(y, p) gamlss.dist::dBNB(y, mu = p$mu, sigma = p$sigma, nu = p$nu, log = TRUE)),
  stan = stan_generated("
  real rc_bnb_lpmf(int y, real mu, real sigma, real nu) {
    real m = inv(sigma) + 1;
    real n = mu * nu / sigma;
    real k = inv(nu);
    return lbeta(y + n, m + k) - lbeta(n, m) - lgamma(y + 1) - lgamma(k) + lgamma(y + k);
  }
  int rc_bnb_rng(real mu, real sigma, real nu) {
    real pr = beta_rng(inv(sigma) + 1, mu * nu / sigma);
    return neg_binomial_rng(inv(nu), pr / (1 - pr));
  }", "rc_bnb_lpmf({y} | {mu}, {sigma}, {nu})",
                        rng = "rc_bnb_rng({mu}, {sigma}, {nu})"),
  start = function(y) list(mu = max(mean(y), 0.1), sigma = 0.2, nu = 0.5),
  status = "validated"
))
