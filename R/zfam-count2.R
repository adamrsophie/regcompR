# Further count distributions with infinite support.

GSD <- "Rigby, Stasinopoulos, Heller and De Bastiani (2019), Distributions for Modeling Location, Scale, and Shape. Chapman and Hall/CRC."
gdd <- gd
gammasgn <- function(x) ifelse(x > 0, 1, (-1)^(floor(-x) + 1))
cnt_par <- function(name, domain = "positive", role = "shape", ...) P(name, domain, role = role, ...)

# Conway Maxwell Poisson ---------------------------------------------------
cmp_lt <- function(k, p) k * log(p$lambda) - p$nu * lgamma(k + 1)
cmp_nc <- normalized_count(cmp_lt)
register_family(family_def(
  "conway_maxwell_poisson", "Conway-Maxwell-Poisson", cnt, discrete = TRUE, aliases = c("cmp", "com_poisson"),
  pars = list(cnt_par("lambda", role = "rate", offset_allowed = TRUE), cnt_par("nu", role = "dispersion",
              interpretation = "nu < 1 overdispersion, nu > 1 underdispersion, nu = 1 Poisson")),
  formula = "p(y) = lambda^y / (y!)^nu / Z(lambda, nu)",
  lpdf = cmp_nc$lpmf, cdf = cmp_nc$cdf, mean = cmp_nc$mean, variance = cmp_nc$variance,
  regression_notes = "lambda is not the mean; use conway_maxwell_poisson_mean for a mean parameterization.",
  references = c("Conway and Maxwell (1962), J. Industrial Engineering, 12, 132-136.",
                 "Shmueli, Minka, Kadane, Borle and Boatwright (2005), JRSS C, 54, 127-142."),
  reference = ref_impl("COMPoissonReg", function(y, p) COMPoissonReg::dcmp(y, p$lambda, p$nu, log = TRUE)),
  stan = stan_generated(stan_series_count("rc_cmp", c("lambda", "nu"), "    return k * log(lambda) - nu * lgamma(k + 1);"),
                        "rc_cmp_lpmf({y} | {lambda}, {nu})"),
  start = function(y) list(lambda = max(mean(y), 0.1), nu = 1), status = "validated"))

cmp_lambda_for_mean <- function(mu, nu) vapply(seq_along(mu), function(i) {
  f <- function(ll) cmp_nc$mean(list(lambda = exp(ll), nu = nu[i])) - mu[i]
  a <- nu[i] * log(max(mu[i], 1e-8)); lo <- a - 5; hi <- a + 5
  while (f(lo) > 0) lo <- lo - 5
  while (f(hi) < 0) hi <- hi + 5
  exp(stats::uniroot(f, c(lo, hi), tol = 1e-12)$root) }, 0)
register_family(family_def(
  "conway_maxwell_poisson_mean", "Conway-Maxwell-Poisson (mean parameterization)", cnt, discrete = TRUE, aliases = "cmp_mu",
  pars = list(log_mean_par("mu"), cnt_par("nu", role = "dispersion")),
  formula = "CMP(lambda, nu) with lambda solving E(Y) = mu",
  lpdf = function(y, p) cmp_nc$lpmf(y, list(lambda = cmp_lambda_for_mean(p$mu, p$nu), nu = p$nu)),
  cdf = function(q, p) cmp_nc$cdf(q, list(lambda = cmp_lambda_for_mean(p$mu, p$nu), nu = p$nu)),
  mean = function(p) p$mu,
  variance = function(p) cmp_nc$variance(list(lambda = cmp_lambda_for_mean(p$mu, p$nu), nu = p$nu)),
  numerical_notes = "lambda is found by root finding on the series mean; slower than the rate form.",
  references = "Huang (2017), Mean-parametrized Conway-Maxwell-Poisson regression models for dispersed counts. Statistical Modelling, 17, 359-380.",
  reference = ref_impl("COMPoissonReg", function(y, p) COMPoissonReg::dcmp(y, cmp_lambda_for_mean(p$mu, p$nu), p$nu, log = TRUE)),
  start = function(y) list(mu = max(mean(y), 0.1), nu = 1), status = "validated"))

# Generalized Poisson (Consul) --------------------------------------------
gp_lt <- function(k, p) {
  a <- p$theta + p$lambda * k
  ifelse(a > 0, log(p$theta) + (k - 1) * log(a) - a - lgamma(k + 1), -Inf)
}
gp_nc <- normalized_count(gp_lt)
gp_check <- function(p) if (any(p$lambda <= -1 | p$lambda >= 1 | p$lambda < -p$theta / 4, na.rm = TRUE))
  "generalized_poisson needs max(-1, -theta/4) < lambda < 1."
register_family(family_def(
  "generalized_poisson", "Generalized Poisson (Consul)", cnt, discrete = TRUE, aliases = "GP",
  pars = list(cnt_par("theta", role = "rate"), P("lambda", "real", role = "dispersion", link = "identity",
              interpretation = "lambda > 0 overdispersion, < 0 underdispersion")),
  check = gp_check,
  formula = "p(y) = theta (theta + lambda y)^(y - 1) exp(-theta - lambda y) / y!; for lambda < 0 truncated where theta + lambda y > 0 and renormalized",
  mean_formula = "theta / (1 - lambda) for lambda >= 0", variance_formula = "theta / (1 - lambda)^3 for lambda >= 0",
  lpdf = function(y, p) { out <- gp_lt(y, p); neg <- p$lambda < 0
    if (any(neg)) out[neg] <- out[neg] - gp_nc$logZ(lapply(p, function(v) rep_len(v, length(y))[neg]))
    out },
  mean = function(p) ifelse(p$lambda >= 0, p$theta / (1 - p$lambda), gp_nc$mean(p)),
  variance = function(p) ifelse(p$lambda >= 0, p$theta / (1 - p$lambda)^3, gp_nc$variance(p)),
  references = "Consul and Jain (1973), A generalization of the Poisson distribution. Technometrics, 15, 791-799.",
  reference = ref_impl("definition", function(y, p) gp_nc$lpmf(y, p), note = "Exact series normalization."),
  stan = stan_generated(stan_series_count("rc_gpois", c("theta", "lambda"),
    "    real a = theta + lambda * k;\n    if (a <= 0) return negative_infinity();\n    return log(theta) + (k - 1) * log(a) - a - lgamma(k + 1);"),
    "rc_gpois_lpmf({y} | {theta}, {lambda})"),
  start = function(y) list(theta = max(mean(y), 0.1), lambda = 0), status = "validated"))

register_family(family_def(
  "generalized_poisson_mean", "Generalized Poisson (mean and dispersion, GP-1)", cnt, discrete = TRUE, aliases = "GP1",
  pars = list(log_mean_par("mu"), P("lambda", "real", role = "dispersion", link = "identity")),
  check = function(p) if (any(p$lambda <= -1 | p$lambda >= 1, na.rm = TRUE)) "lambda must lie in (-1, 1).",
  formula = "Consul generalized Poisson with theta = mu (1 - lambda)",
  mean_formula = "mu (for lambda >= 0)", variance_formula = "mu / (1 - lambda)^2",
  lpdf = function(y, p) { q <- list(theta = p$mu * (1 - p$lambda), lambda = p$lambda); out <- gp_lt(y, q); neg <- p$lambda < 0
    if (any(neg)) out[neg] <- out[neg] - gp_nc$logZ(lapply(q, function(v) rep_len(v, length(y))[neg])); out },
  mean = function(p) ifelse(p$lambda >= 0, p$mu, gp_nc$mean(list(theta = p$mu * (1 - p$lambda), lambda = p$lambda))),
  variance = function(p) ifelse(p$lambda >= 0, p$mu / (1 - p$lambda)^2, gp_nc$variance(list(theta = p$mu * (1 - p$lambda), lambda = p$lambda))),
  references = "Famoye (1993), Restricted generalized Poisson regression model. Comm. Stat. Theory Methods, 22, 1335-1354.",
  reference = ref_impl("definition", function(y, p) gp_nc$lpmf(y, list(theta = p$mu * (1 - p$lambda), lambda = p$lambda))),
  stan = stan_generated(stan_series_count("rc_gpois", c("theta", "lambda"),
    "    real a = theta + lambda * k;\n    if (a <= 0) return negative_infinity();\n    return log(theta) + (k - 1) * log(a) - a - lgamma(k + 1);"),
    "rc_gpois_lpmf({y} | {mu} * (1 - {lambda}), {lambda})"),
  start = function(y) list(mu = max(mean(y), 0.1), lambda = 0), status = "validated"))

# Double Poisson (Efron) ---------------------------------------------------
dpo_lt <- function(k, p) {
  xl <- ifelse(k > 0, k * log(pmax(k, 1)), 0)
  -0.5 * log(p$sigma) - p$mu / p$sigma - lgamma(k + 1) + xl - k + k * log(p$mu) / p$sigma + k / p$sigma - xl / p$sigma
}
dpo_nc <- normalized_count(dpo_lt)
register_family(family_def(
  "double_poisson", "Double Poisson (Efron)", cnt, discrete = TRUE, aliases = "DPO",
  pars = list(P("mu", "positive", role = "location", offset_allowed = TRUE), cnt_par("sigma", role = "dispersion")),
  formula = "p(y) proportional to sigma^(-1/2) exp(-mu/sigma) (e^-y y^y / y!) (e mu / y)^(y / sigma)",
  regression_notes = "mu is approximately the mean and sigma approximately the variance to mean ratio.",
  lpdf = dpo_nc$lpmf, cdf = dpo_nc$cdf, mean = dpo_nc$mean, variance = dpo_nc$variance,
  references = c("Efron (1986), Double exponential families and their use in generalized linear regression. JASA, 81, 709-721.", GSD),
  reference = ref_impl("gamlss.dist", gdd("dDPO")),
  stan = stan_generated(stan_series_count("rc_dpo", c("mu", "sigma"),
    "    real xl = k > 0 ? k * log(k) : 0;\n    return -0.5 * log(sigma) - mu / sigma - lgamma(k + 1) + xl - k + k * log(mu) / sigma + k / sigma - xl / sigma;"),
    "rc_dpo_lpmf({y} | {mu}, {sigma})"),
  start = function(y) list(mu = max(mean(y), 0.1), sigma = 1), status = "validated"))

# Hyper-Poisson -------------------------------------------------------------
hp_lt <- function(k, p) k * log(p$lambda) - (lgamma(p$gamma + k) - lgamma(p$gamma))
hp_nc <- normalized_count(hp_lt)
register_family(family_def(
  "hyper_poisson", "Hyper-Poisson (Bardwell and Crow)", cnt, discrete = TRUE,
  pars = list(cnt_par("lambda", role = "rate"), cnt_par("gamma", role = "dispersion", interpretation = "gamma > 1 overdispersion, < 1 underdispersion")),
  formula = "p(y) = lambda^y / ((gamma)_y 1F1(1; gamma; lambda))",
  lpdf = hp_nc$lpmf, cdf = hp_nc$cdf, mean = hp_nc$mean, variance = hp_nc$variance,
  references = "Bardwell and Crow (1964), A two-parameter family of hyper-Poisson distributions. JASA, 59, 133-141.",
  reference = ref_impl("definition", function(y, p) hp_lt(y, p) - log(vapply(seq_along(y), function(i)
    sum(exp(hp_lt(0:3000, list(lambda = p$lambda[i], gamma = p$gamma[i])))), 0))),
  stan = stan_generated(stan_series_count("rc_hyperpois", c("lambda", "gamma"),
    "    return k * log(lambda) - (lgamma(gamma + k) - lgamma(gamma));"), "rc_hyperpois_lpmf({y} | {lambda}, {gamma})"),
  start = function(y) list(lambda = max(mean(y), 0.1), gamma = 1), status = "validated"))

# Gamma count (Winkelmann) -----------------------------------------------
gc_cdfg <- function(k, p) ifelse(k <= 0, 1, stats::pgamma(p$lambda, shape = k * p$alpha))
register_family(family_def(
  "gamma_count", "Gamma count (Winkelmann)", cnt, discrete = TRUE,
  pars = list(cnt_par("lambda", role = "rate", interpretation = "Event rate times exposure"),
              cnt_par("alpha", role = "dispersion", interpretation = "alpha > 1 underdispersion, < 1 overdispersion")),
  formula = "P(Y = y) = G(alpha y, lambda) - G(alpha (y + 1), lambda), G the regularized lower incomplete gamma (G(0, .) = 1)",
  lpdf = function(y, p) log(pmax(gc_cdfg(y, p) - gc_cdfg(y + 1, p), 0)),
  cdf = function(q, p) 1 - gc_cdfg(floor(q) + 1, p),
  references = "Winkelmann (1995), Duration dependence and dispersion in count-data models. JBES, 13, 467-474.",
  reference = ref_impl("definition", function(y, p) log(stats::pgamma(p$lambda, y * p$alpha) - stats::pgamma(p$lambda, (y + 1) * p$alpha)) ),
  stan = sfun("gamma_count", c("lambda", "alpha"),
    "    real a = y == 0 ? 1 : gamma_p(y * alpha, lambda);\n    return log(a - gamma_p((y + 1) * alpha, lambda));", discrete = TRUE),
  start = function(y) list(lambda = max(mean(y), 0.1), alpha = 1), status = "validated"))

# Poisson inverse Gaussian and Sichel --------------------------------------
lbk <- function(nu, x) log_besselK(x, nu)
register_family(family_def(
  "poisson_inverse_gaussian", "Poisson inverse Gaussian", cnt, discrete = TRUE, aliases = c("PIG", "pig"),
  pars = list(log_mean_par("mu"), cnt_par("sigma", role = "dispersion")),
  formula = "p(y) = sqrt(2 a / pi) exp(1/sigma) mu^y K_{y - 1/2}(a) / ((a sigma)^y y!), a^2 = 1/sigma^2 + 2 mu / sigma",
  mean_formula = "mu", variance_formula = "mu + sigma mu^2",
  lpdf = function(y, p) { a <- sqrt(1 / p$sigma^2 + 2 * p$mu / p$sigma)
    0.5 * log(2 * a / pi) + 1 / p$sigma + y * log(p$mu) + lbk(y - 0.5, a) - y * log(a * p$sigma) - lgamma(y + 1) },
  mean = function(p) p$mu, variance = function(p) p$mu + p$sigma * p$mu^2,
  references = c("Holla (1966), On a Poisson-inverse Gaussian distribution. Metrika, 11, 115-121.", GSD),
  reference = ref_impl("gamlss.dist", gdd("dPIG")),
  stan = sfun("pig", c("mu", "sigma"),
    "    real a = sqrt(1 / square(sigma) + 2 * mu / sigma);\n    return 0.5 * log(2 * a / pi()) + 1 / sigma + y * log(mu) + rc_log_besselK(y - 0.5, a) - y * log(a * sigma) - lgamma(y + 1);",
    discrete = TRUE, needs = "log_besselK"),
  start = function(y) list(mu = max(mean(y), 0.1), sigma = 0.5), status = "validated"))

register_family(family_def(
  "sichel", "Sichel (Poisson generalized inverse Gaussian, mean parameterization)", cnt, discrete = TRUE, aliases = c("SICHEL", "poisson_gig"),
  pars = list(log_mean_par("mu"), cnt_par("sigma", role = "dispersion"), P("nu", "real", role = "shape")),
  formula = "p(y) = (mu/c)^y K_{y + nu}(a) / (y! (a sigma)^(y + nu) K_nu(1/sigma)), c = K_{nu+1}(1/sigma) / K_nu(1/sigma), a^2 = 1/sigma^2 + 2 mu / (c sigma)",
  mean_formula = "mu",
  lpdf = function(y, p) { cc <- exp(lbk(p$nu + 1, 1 / p$sigma) - lbk(p$nu, 1 / p$sigma))
    a <- sqrt(1 / p$sigma^2 + 2 * p$mu / (cc * p$sigma))
    y * log(p$mu / cc) + lbk(y + p$nu, a) - lgamma(y + 1) - (y + p$nu) * log(a * p$sigma) - lbk(p$nu, 1 / p$sigma) },
  mean = function(p) p$mu,
  references = c("Sichel (1971), On a family of discrete distributions particularly suited to represent long-tailed frequency data.", GSD),
  reference = ref_impl("gamlss.dist", gdd("dSICHEL")),
  stan = sfun("sichel", c("mu", "sigma", "nu"),
    "    real c = exp(rc_log_besselK(nu + 1, 1 / sigma) - rc_log_besselK(nu, 1 / sigma));\n    real a = sqrt(1 / square(sigma) + 2 * mu / (c * sigma));\n    return y * log(mu / c) + rc_log_besselK(y + nu, a) - lgamma(y + 1) - (y + nu) * log(a * sigma) - rc_log_besselK(nu, 1 / sigma);",
    discrete = TRUE, needs = "log_besselK"),
  start = function(y) list(mu = max(mean(y), 0.1), sigma = 0.5, nu = -0.5), status = "validated"))

# Delaporte -------------------------------------------------------------
del_lpmf <- function(y, p) vapply(seq_along(y), function(i) {
  j <- 0:y[i]; lam <- p$mu[i] * p$nu[i]; m2 <- p$mu[i] * (1 - p$nu[i])
  log_sum_exp(stats::dpois(y[i] - j, lam, log = TRUE) + stats::dnbinom(j, size = 1 / p$sigma[i], mu = m2, log = TRUE)) }, 0)
register_family(family_def(
  "delaporte", "Delaporte", cnt, discrete = TRUE, aliases = "DEL",
  pars = list(log_mean_par("mu"), cnt_par("sigma", role = "dispersion"), P("nu", "unit", role = "shape", interpretation = "Poisson share of the mean")),
  formula = "Y = Poisson(mu nu) + negative binomial(mean mu (1 - nu), size 1/sigma)",
  mean_formula = "mu", variance_formula = "mu + sigma mu^2 (1 - nu)^2",
  lpdf = function(y, p) del_lpmf(y, p),
  rng = function(n, p) stats::rpois(n, p$mu * p$nu) + stats::rnbinom(n, size = 1 / p$sigma, mu = p$mu * (1 - p$nu)),
  mean = function(p) p$mu, variance = function(p) p$mu + p$sigma * p$mu^2 * (1 - p$nu)^2,
  references = c("Delaporte (1959); Ruohonen (1988), Scandinavian Actuarial Journal.", GSD),
  reference = ref_impl("gamlss.dist", gdd("dDEL")),
  stan = sfun("delaporte", c("mu", "sigma", "nu"),
    "    vector[y + 1] l;\n    for (j in 0:y) l[j + 1] = poisson_lpmf(y - j | mu * nu) + neg_binomial_2_lpmf(j | mu * (1 - nu), 1 / sigma);\n    return log_sum_exp(l);",
    discrete = TRUE),
  start = function(y) list(mu = max(mean(y), 0.1), sigma = 0.5, nu = 0.5), status = "validated"))

# Logarithmic families -----------------------------------------------------
register_family(family_def(
  "logarithmic", "Logarithmic (log-series)", cnt, discrete = TRUE, aliases = c("log_series", "LG"),
  pars = list(P("theta", "unit", role = "shape")),
  formula = "p(y) = -theta^y / (y log(1 - theta)), y = 1, 2, ...",
  mean_formula = "-theta / ((1 - theta) log(1 - theta))",
  lpdf = function(y, p) y * log(p$theta) - log(y) - log(-log1p(-p$theta)),
  mean = function(p) -p$theta / ((1 - p$theta) * log1p(-p$theta)),
  variance = function(p) { L <- log1p(-p$theta); -p$theta * (p$theta + L) / ((1 - p$theta)^2 * L^2) },
  support = function(p) list(lower = 1, upper = Inf),
  references = c("Fisher, Corbet and Williams (1943), J. Animal Ecology, 12, 42-58.", JKK),
  reference = ref_impl("extraDistr", function(y, p) extraDistr::dlgser(y, p$theta, log = TRUE)),
  stan = sfun("logarithmic", "theta", "    return y * log(theta) - log(y) - log(-log1m(theta));", discrete = TRUE),
  bugs = list(loglik = "{y} * log({theta}) - log({y}) - log(-log(1 - {theta}))"),
  start = function(y) list(theta = 0.5), status = "validated"))

register_family(family_def(
  "generalized_log_series", "Generalized logarithmic series (Jain and Gupta)", cnt, discrete = TRUE,
  pars = list(P("theta", "unit", role = "shape"), P("beta", "greater_than_one", role = "shape")),
  check = function(p) if (any(p$theta * p$beta >= 1, na.rm = TRUE)) "Need theta < 1 / beta.",
  formula = "p(y) = -1/log(1 - theta) Gamma(beta y) / (Gamma(y + 1) Gamma(beta y - y + 1)) theta^y (1 - theta)^(beta y - y)",
  lpdf = function(y, p) -log(-log1p(-p$theta)) + lgamma(p$beta * y) - lgamma(y + 1) - lgamma(p$beta * y - y + 1) +
    y * log(p$theta) + (p$beta * y - y) * log1p(-p$theta),
  support = function(p) list(lower = 1, upper = Inf),
  references = "Jain and Gupta (1973), A logarithmic series type distribution. Trabajos de Estadistica, 24, 99-105.",
  reference = ref_impl("logarithmic at beta = 1", function(y, p) -log(-log1p(-p$theta)) + lgamma(p$beta * y) - lgamma(y + 1) -
    lgamma(p$beta * y - y + 1) + y * log(p$theta) + (p$beta * y - y) * log(1 - p$theta), note = "Definition; reduces to the logarithmic law at beta = 1."),
  stan = sfun("glsd", c("theta", "beta"),
    "    return -log(-log1m(theta)) + lgamma(beta * y) - lgamma(y + 1) - lgamma(beta * y - y + 1) + y * log(theta) + (beta * y - y) * log1m(theta);", discrete = TRUE),
  start = function(y) list(theta = 0.3, beta = 1.2), status = "implemented"))

# Power law counts -------------------------------------------------------
register_family(family_def(
  "zeta", "Zeta (discrete Pareto, Zipf on all positive integers)", cnt, discrete = TRUE, aliases = c("ZIPF", "zipf"),
  pars = list(P("s", "greater_than_one", role = "exponent")),
  formula = "p(y) = y^(-s) / zeta(s), y = 1, 2, ...",
  moment_conditions = "Mean exists iff s > 2, variance iff s > 3.",
  lpdf = function(y, p) -p$s * log(y) - log(zeta_fn(p$s)),
  mean = function(p) ifelse(p$s > 2, zeta_fn(p$s - 1) / zeta_fn(p$s), NA_real_),
  variance = function(p) ifelse(p$s > 3, zeta_fn(p$s - 2) / zeta_fn(p$s) - (zeta_fn(p$s - 1) / zeta_fn(p$s))^2, NA_real_),
  support = function(p) list(lower = 1, upper = Inf),
  references = JKK, reference = ref_impl("VGAM", function(y, p) VGAM::dzeta(y, p$s - 1, log = TRUE)),
  stan = sfun("zeta", "s", "    return -s * log(y) - log(rc_zeta(s));", discrete = TRUE, needs = "zeta"),
  start = function(y) list(s = 2.5), status = "validated"))

register_family(family_def(
  "yule_simon", "Yule-Simon", cnt, discrete = TRUE, aliases = "YULE",
  pars = list(cnt_par("rho", role = "shape")),
  formula = "p(y) = rho B(y, rho + 1), y = 1, 2, ...",
  mean_formula = "rho / (rho - 1) for rho > 1", variance_formula = "rho^2 / ((rho - 1)^2 (rho - 2)) for rho > 2",
  lpdf = function(y, p) log(p$rho) + lbeta(y, p$rho + 1),
  cdf = function(q, p) 1 - floor(q) * beta(floor(q), p$rho + 1),
  mean = function(p) ifelse(p$rho > 1, p$rho / (p$rho - 1), NA_real_),
  variance = function(p) ifelse(p$rho > 2, p$rho^2 / ((p$rho - 1)^2 * (p$rho - 2)), NA_real_),
  support = function(p) list(lower = 1, upper = Inf),
  references = "Simon (1955), On a class of skew distribution functions. Biometrika, 42, 425-440.",
  reference = ref_impl("VGAM", function(y, p) VGAM::dyules(y, p$rho, log = TRUE)),
  stan = sfun("yule_simon", "rho", "    return log(rho) + lbeta(y, rho + 1);", discrete = TRUE),
  bugs = list(loglik = "log({rho}) + loggam({y}) + loggam({rho} + 1) - loggam({y} + {rho} + 1)"),
  start = function(y) list(rho = 2), status = "validated"))

register_family(family_def(
  "waring", "Waring (GAMLSS mean parameterization)", cnt, discrete = TRUE, aliases = "WARING",
  pars = list(log_mean_par("mu"), cnt_par("sigma", role = "dispersion")),
  formula = "p(y) = B(y + mu/sigma, 1/sigma + 2) / B(mu/sigma, 1/sigma + 1)",
  mean_formula = "mu", variance_formula = "mu (mu + 1) (1 + sigma) / (1 - sigma) for sigma < 1",
  lpdf = function(y, p) lbeta(y + p$mu / p$sigma, 1 / p$sigma + 2) - lbeta(p$mu / p$sigma, 1 / p$sigma + 1),
  mean = function(p) p$mu,
  references = c("Irwin (1975), JRSS A, 138, 18-31.", GSD), reference = ref_impl("gamlss.dist", gdd("dWARING")),
  stan = sfun("waring", c("mu", "sigma"), "    return lbeta(y + mu / sigma, 1 / sigma + 2) - lbeta(mu / sigma, 1 / sigma + 1);", discrete = TRUE),
  start = function(y) list(mu = max(mean(y), 0.1), sigma = 0.5), status = "validated"))

register_family(family_def(
  "borel", "Borel", cnt, discrete = TRUE,
  pars = list(P("mu", "unit", role = "shape", interpretation = "Mean offspring number of the branching process")),
  formula = "p(y) = exp(-mu y) (mu y)^(y - 1) / y!, y = 1, 2, ...",
  mean_formula = "1 / (1 - mu)", variance_formula = "mu / (1 - mu)^3",
  lpdf = function(y, p) -p$mu * y + (y - 1) * log(p$mu * y) - lgamma(y + 1),
  mean = function(p) 1 / (1 - p$mu), variance = function(p) p$mu / (1 - p$mu)^3,
  support = function(p) list(lower = 1, upper = Inf),
  references = "Borel (1942), Sur l'emploi du theoreme de Bernoulli. Comptes Rendus, 214, 452-456.",
  reference = ref_impl("VGAM", function(y, p) VGAM::dbort(y, Qsize = 1, a = p$mu, log = TRUE), note = "Borel-Tanner with Q = 1."),
  stan = sfun("borel", "mu", "    return -mu * y + (y - 1) * log(mu * y) - lgamma(y + 1);", discrete = TRUE),
  start = function(y) list(mu = 0.5), status = "validated"))

register_family(family_def(
  "borel_tanner", "Borel-Tanner", cnt, discrete = TRUE,
  pars = list(P("mu", "unit", role = "shape"), P("n", "positive_integer", role = "initial", data = TRUE)),
  formula = "p(y) = n / y exp(-mu y) (mu y)^(y - n) / (y - n)!, y = n, n + 1, ...",
  mean_formula = "n / (1 - mu)", variance_formula = "n mu / (1 - mu)^3",
  lpdf = function(y, p) log(p$n) - log(y) - p$mu * y + (y - p$n) * log(p$mu * y) - lgamma(y - p$n + 1),
  mean = function(p) p$n / (1 - p$mu), variance = function(p) p$n * p$mu / (1 - p$mu)^3,
  support = function(p) list(lower = p$n, upper = Inf),
  references = "Haight and Breuer (1960), The Borel-Tanner distribution. Biometrika, 47, 143-150.",
  reference = ref_impl("VGAM", function(y, p) mapply(function(x, Q, a) VGAM::dbort(x, Qsize = Q, a = a, log = TRUE), y, p$n, p$mu)),
  stan = sfun("borel_tanner", c("mu", "n"), "    return log(n) - log(y) - mu * y + (y - n) * log(mu * y) - lgamma(y - n + 1);", discrete = TRUE),
  status = "validated"))

register_family(family_def(
  "gauss_kuzmin", "Gauss-Kuzmin", cnt, discrete = TRUE,
  pars = list(P("scale_unused", "positive", role = "none", regression = FALSE, default = 1,
                interpretation = "Placeholder; the law has no parameters")),
  formula = "p(k) = -log2(1 - 1/(k + 1)^2), k = 1, 2, ...",
  moment_conditions = "The mean is infinite.",
  lpdf = function(y, p) log(-log2(1 - 1 / (y + 1)^2)),
  cdf = function(q, p) ifelse(q < 1, 0, 1 - log2((floor(q) + 2) / (floor(q) + 1))),
  quantile = function(pr, p) pmax(ceiling(1 / (2^(1 - pr) - 1) - 2 - 1e-9), 1),
  mean = function(p) rep(Inf, length(p[[1]])),
  support = function(p) list(lower = 1, upper = Inf), numeric_moments = FALSE,
  references = "Kuzmin (1928); Khinchin (1964), Continued Fractions.",
  reference = ref_impl("closed form", function(y, p) log(log2(1 + 1 / (y * (y + 2))))),
  status = "validated"))

# Hermite and Polya-Aeppli -------------------------------------------------
register_family(family_def(
  "hermite", "Hermite (Poisson plus twice Poisson)", cnt, discrete = TRUE,
  pars = list(cnt_par("a1", role = "rate"), cnt_par("a2", role = "rate")),
  formula = "Y = X1 + 2 X2, X1 ~ Poisson(a1), X2 ~ Poisson(a2); p(y) = exp(-a1 - a2) sum_j a1^(y - 2j) a2^j / ((y - 2j)! j!)",
  mean_formula = "a1 + 2 a2", variance_formula = "a1 + 4 a2",
  lpdf = function(y, p) vapply(seq_along(y), function(i) { j <- 0:floor(y[i] / 2)
    -p$a1[i] - p$a2[i] + log_sum_exp((y[i] - 2 * j) * log(p$a1[i]) + j * log(p$a2[i]) - lgamma(y[i] - 2 * j + 1) - lgamma(j + 1)) }, 0),
  rng = function(n, p) stats::rpois(n, p$a1) + 2 * stats::rpois(n, p$a2),
  mean = function(p) p$a1 + 2 * p$a2, variance = function(p) p$a1 + 4 * p$a2,
  references = "Kemp and Kemp (1965), Some properties of the Hermite distribution. Biometrika, 52, 381-394.",
  reference = ref_impl("hermite", function(y, p) log(mapply(function(x, a, b) hermite::dhermite(x, a, b, m = 2), y, p$a1, p$a2))),
  stan = sfun("hermite", c("a1", "a2"),
    "    int J = y / 2;\n    vector[J + 1] l;\n    for (j in 0:J) l[j + 1] = (y - 2 * j) * log(a1) + j * log(a2) - lgamma(y - 2 * j + 1) - lgamma(j + 1);\n    return -a1 - a2 + log_sum_exp(l);",
    discrete = TRUE),
  start = function(y) list(a1 = max(mean(y) / 2, 0.1), a2 = max(mean(y) / 4, 0.1)), status = "validated"))

register_family(family_def(
  "polya_aeppli", "Polya-Aeppli (geometric Poisson)", cnt, discrete = TRUE,
  pars = list(cnt_par("lambda", role = "rate"), P("prob", "unit", role = "shape")),
  formula = "Poisson(lambda) number of clusters, cluster sizes geometric on 1, 2, ... with P(size = k) = (1 - prob) prob^(k - 1)",
  mean_formula = "lambda / (1 - prob)", variance_formula = "lambda (1 + prob) / (1 - prob)^2",
  lpdf = function(y, p) vapply(seq_along(y), function(i) {
    if (y[i] == 0) return(-p$lambda[i]); k <- 1:y[i]
    -p$lambda[i] + log_sum_exp(lchoose(y[i] - 1, k - 1) + k * log(p$lambda[i] * (1 - p$prob[i]) / p$prob[i]) - lgamma(k + 1)) + y[i] * log(p$prob[i]) }, 0),
  mean = function(p) p$lambda / (1 - p$prob), variance = function(p) p$lambda * (1 + p$prob) / (1 - p$prob)^2,
  references = "Johnson, Kemp and Kotz (2005), sec. 9.7.",
  reference = ref_impl("polyaAeppli", function(y, p) log(mapply(polyaAeppli::dPolyaAeppli, y, p$lambda, p$prob))),
  stan = sfun("polya_aeppli", c("lambda", "prob"),
    "    if (y == 0) return -lambda;\n    vector[y] l;\n    for (k in 1:y) l[k] = lchoose(y - 1, k - 1) + k * log(lambda * (1 - prob) / prob) - lgamma(k + 1);\n    return -lambda + log_sum_exp(l) + y * log(prob);",
    discrete = TRUE),
  start = function(y) list(lambda = max(mean(y) / 2, 0.1), prob = 0.3), status = "validated"))

# Poisson lognormal -------------------------------------------------------
# Integrate over the latent normal z in log space, centered at the mode of
# g(z) = y (mu + sigma z) - exp(mu + sigma z) - z^2/2 with a window of 40
# curvature standard deviations, so narrow peaks at large y are not missed.
pln_lpmf <- function(y, p) vapply(seq_along(y), function(i) {
  yy <- y[i]; mu <- p$mu[i]; s <- p$sigma[i]
  g <- function(z) yy * (mu + s * z) - exp(mu + s * z) - z^2 / 2
  zs <- stats::uniroot(function(z) yy * s - s * exp(mu + s * z) - z, c(-40, 40), extendInt = "yes", tol = 1e-12)$root
  sd <- 1 / sqrt(s^2 * exp(mu + s * zs) + 1); g0 <- g(zs)
  v <- stats::integrate(function(z) exp(g(z) - g0), zs - 40 * sd, zs + 40 * sd, rel.tol = 1e-12, subdivisions = 1000L)$value
  log(v) + g0 - lgamma(yy + 1) - 0.5 * log(2 * pi) }, 0)
register_family(family_def(
  "poisson_lognormal", "Poisson lognormal", cnt, discrete = TRUE, aliases = "poilog",
  pars = list(P("mu", "real", role = "location", offset_allowed = TRUE, interpretation = "Mean of log lambda"),
              cnt_par("sigma", role = "dispersion", interpretation = "Standard deviation of log lambda")),
  formula = "p(y) = int Poisson(y | exp(mu + sigma z)) phi(z) dz",
  mean_formula = "exp(mu + sigma^2 / 2)", variance_formula = "m + m^2 (exp(sigma^2) - 1), m = exp(mu + sigma^2/2)",
  lpdf = pln_lpmf,
  rng = function(n, p) stats::rpois(n, exp(p$mu + p$sigma * stats::rnorm(n))),
  mean = function(p) exp(p$mu + p$sigma^2 / 2),
  variance = function(p) { m <- exp(p$mu + p$sigma^2 / 2); m + m^2 * (exp(p$sigma^2) - 1) },
  numerical_notes = "Adaptive quadrature over the latent normal; Stan code uses an observation level random effect instead (see regression_notes).",
  regression_notes = "Equivalent to a Poisson model with an observation level normal random effect (1 | obs).",
  references = "Bulmer (1974), On fitting the Poisson lognormal distribution to species-abundance data. Biometrics, 30, 101-110.",
  reference = ref_impl("poilog", function(y, p) log(mapply(poilog::dpoilog, y, p$mu, p$sigma))),
  start = function(y) list(mu = log(max(mean(y), 0.1)), sigma = 0.5), status = "validated"))

# Skellam -------------------------------------------------------------------
register_family(family_def(
  "skellam", "Skellam (difference of two Poisson variables)", cnt, discrete = TRUE,
  pars = list(cnt_par("mu1", role = "rate"), cnt_par("mu2", role = "rate")),
  formula = "p(k) = exp(-mu1 - mu2) (mu1 / mu2)^(k/2) I_|k|(2 sqrt(mu1 mu2))",
  mean_formula = "mu1 - mu2", variance_formula = "mu1 + mu2",
  lpdf = function(y, p) { z <- 2 * sqrt(p$mu1 * p$mu2)
    -p$mu1 - p$mu2 + y / 2 * log(p$mu1 / p$mu2) + log(besselI(z, abs(y), expon.scaled = TRUE)) + z },
  rng = function(n, p) stats::rpois(n, p$mu1) - stats::rpois(n, p$mu2),
  mean = function(p) p$mu1 - p$mu2, variance = function(p) p$mu1 + p$mu2,
  support = function(p) list(lower = -Inf, upper = Inf),
  references = "Skellam (1946), The frequency distribution of the difference between two Poisson variates. JRSS, 109, 296.",
  reference = ref_impl("skellam", function(y, p) skellam::dskellam(y, p$mu1, p$mu2, log = TRUE)),
  stan = sfun("skellam", c("mu1", "mu2"),
    "    real z = 2 * sqrt(mu1 * mu2);\n    return -mu1 - mu2 + y * 0.5 * log(mu1 / mu2) + log_modified_bessel_first_kind(abs(y), z);",
    discrete = TRUE),
  start = function(y) list(mu1 = max(mean(y[y > 0]), 1, na.rm = TRUE), mu2 = 1), status = "validated"))

# Discrete Weibull (Nakagawa and Osaki) ---------------------------------------
register_family(family_def(
  "discrete_weibull", "Discrete Weibull (type I, Nakagawa and Osaki)", cnt, discrete = TRUE, aliases = "DW",
  pars = list(P("q", "unit", role = "shape"), cnt_par("beta", role = "shape")),
  formula = "P(Y >= y) = q^(y^beta); p(y) = q^(y^beta) - q^((y + 1)^beta), y = 0, 1, ...",
  lpdf = function(y, p) log(p$q^(y^p$beta) - p$q^((y + 1)^p$beta)),
  cdf = function(q, p) 1 - p$q^((floor(q) + 1)^p$beta),
  quantile = function(pr, p) pmax(ceiling((log(1 - pr) / log(p$q))^(1 / p$beta) - 1 - 1e-10), 0),
  references = "Nakagawa and Osaki (1975), The discrete Weibull distribution. IEEE Trans. Reliability, 24, 300-301.",
  reference = ref_impl("DiscreteWeibull", function(y, p) log(mapply(DiscreteWeibull::ddweibull, y, p$q, p$beta, MoreArgs = list(zero = TRUE)))),
  stan = sfun("discrete_weibull", c("q", "beta"),
    "    return log_diff_exp(pow(y, beta) * log(q), pow(y + 1, beta) * log(q));", discrete = TRUE,
    lcdf = NULL),
  start = function(y) list(q = 0.7, beta = 1), status = "validated"))

# Bell ------------------------------------------------------------------------
bell_numbers <- function(nmax) {
  B <- numeric(nmax + 1); B[1] <- 0
  row <- 0
  for (n in seq_len(nmax)) {
    nr <- numeric(n + 1); nr[1] <- row[length(row)]
    for (j in seq_len(n)) nr[j + 1] <- log_add_exp(nr[j], row[j])
    B[n + 1] <- nr[1]; row <- nr
  }
  B[2] <- 0
  B
}
lbell <- function(y) { b <- bell_numbers(max(c(y, 1))); b[y + 1] }
register_family(family_def(
  "bell", "Bell", cnt, discrete = TRUE,
  pars = list(cnt_par("theta", role = "shape")),
  formula = "p(y) = theta^y exp(1 - e^theta) B_y / y!, B_y the Bell numbers",
  mean_formula = "theta e^theta", variance_formula = "theta (1 + theta) e^theta",
  lpdf = function(y, p) y * log(p$theta) + 1 - exp(p$theta) + lbell(y) - lgamma(y + 1),
  mean = function(p) p$theta * exp(p$theta), variance = function(p) p$theta * (1 + p$theta) * exp(p$theta),
  regression_notes = "Bell regression usually models the mean mu = theta exp(theta); use links on theta or transform with W0(mu).",
  references = "Castellares, Ferrari and Lemonte (2018), On the Bell distribution and its associated regression model for count data. Applied Mathematical Modelling, 56, 172-185.",
  reference = ref_impl("definition", function(y, p) y * log(p$theta) + 1 - exp(p$theta) +
    log(vapply(y, function(k) { if (k == 0) return(1); s <- 0; for (j in 0:k) s <- s + sum((-1)^(0:j) * choose(j, 0:j) * (j - 0:j)^k) / factorial(j); s }, 0)) - lgamma(y + 1)),
  start = function(y) list(theta = 1), status = "validated"))

register_family(family_def(
  "poisson_lindley", "Poisson-Lindley (Sankaran)", cnt, discrete = TRUE,
  pars = list(cnt_par("theta", role = "shape")),
  formula = "p(y) = theta^2 (y + theta + 2) / (theta + 1)^(y + 3)",
  mean_formula = "(theta + 2) / (theta (theta + 1))",
  lpdf = function(y, p) 2 * log(p$theta) + log(y + p$theta + 2) - (y + 3) * log1p(p$theta),
  mean = function(p) (p$theta + 2) / (p$theta * (p$theta + 1)),
  variance = function(p) (p$theta^3 + 4 * p$theta^2 + 6 * p$theta + 2) / (p$theta^2 * (p$theta + 1)^2),
  references = "Sankaran (1970), The discrete Poisson-Lindley distribution. Biometrics, 26, 145-149.",
  reference = ref_impl("definition", function(y, p) log(vapply(seq_along(y), function(i) stats::integrate(function(l)
    stats::dpois(y[i], l) * p$theta[i]^2 / (p$theta[i] + 1) * (1 + l) * exp(-p$theta[i] * l), 0, Inf, rel.tol = 1e-12)$value, 0)),
    note = "Poisson mixture over the Lindley density."),
  stan = sfun("poisson_lindley", "theta", "    return 2 * log(theta) + log(y + theta + 2) - (y + 3) * log1p(theta);", discrete = TRUE),
  bugs = list(loglik = "2 * log({theta}) + log({y} + {theta} + 2) - ({y} + 3) * log(1 + {theta})"),
  start = function(y) list(theta = 1), status = "validated"))

register_family(family_def(
  "extended_negative_binomial", "Extended negative binomial (Engen; Hess, Lewis and Liebscher)", cnt, discrete = TRUE,
  pars = list(P("m", "positive_integer", role = "truncation", data = TRUE),
              P("r", "real", role = "shape", link = "identity", interpretation = "r > -m; r in (-1, 0) for m = 1 gives Engen's law"),
              P("prob", "unit", role = "shape")),
  formula = "f(k) = choose(k + r - 1, k) p^k / ((1 - p)^(-r) - sum_{j < m} choose(j + r - 1, j) p^j), k >= m",
  lpdf = function(y, p) vapply(seq_along(y), function(i) {
    r <- p$r[i]; pr <- p$prob[i]; m <- p$m[i]
    # signed log of choose(k + r - 1, k) = Gamma(k + r) / (k! Gamma(r))
    lsc <- function(k) list(l = lgamma(k + r) - lgamma(k + 1) - lgamma(r),
                            s = gammasgn(k + r) * gammasgn(r))
    j <- seq_len(m) - 1
    tj <- vapply(j, function(k) { z <- lsc(k); z$s * exp(z$l + k * log(pr)) }, 0)
    den <- (1 - pr)^(-r) - sum(tj)
    z <- lsc(y[i])
    z$l + y[i] * log(pr) - log(abs(den)) + log(z$s * sign(den)) }, 0),
  support = function(p) list(lower = p$m, upper = Inf),
  numerical_notes = "Uses gamma ratios with sign for negative r; moderate y only.",
  references = "Hess, Lewis and Liebscher (2002), Extended negative binomial distribution. ASTIN Bulletin, 32, 283-303.",
  reference = ref_impl("definition (r > 0 reduces to a zero truncated negative binomial for m = 1)",
    function(y, p) stats::dnbinom(y, p$r, 1 - p$prob, log = TRUE) - log1p(-stats::dnbinom(0, p$r, 1 - p$prob))),
  status = "implemented"))

pf_lt <- function(k, p) -p$b * log(k) - p$c * log(k)^2
pf_nc <- normalized_count(pf_lt, k0 = 1)
register_family(family_def(
  "parabolic_fractal", "Parabolic fractal", cnt, discrete = TRUE,
  pars = list(P("b", "real", role = "shape", link = "identity"), cnt_par("c", role = "curvature")),
  formula = "f(n) proportional to n^(-b) exp(-c (log n)^2), n = 1, 2, ...",
  lpdf = pf_nc$lpmf, cdf = pf_nc$cdf, mean = pf_nc$mean, variance = pf_nc$variance,
  support = function(p) list(lower = 1, upper = Inf),
  references = "Laherrere and Deheuvels (1996), Comptes Rendus de l'Academie des Sciences, 322, 535-541.",
  reference = ref_impl("definition", function(y, p) pf_lt(y, p) - log(vapply(seq_along(y), function(i)
    sum(exp(pf_lt(1:200000, list(b = p$b[i], c = p$c[i])))), 0))),
  stan = stan_generated(stan_series_count("rc_parafrac", c("b", "c"), "    return -b * log(k) - c * square(log(k));", k0 = 1),
                        "rc_parafrac_lpmf({y} | {b}, {c})"),
  start = function(y) list(b = 1, c = 0.5), status = "validated"))

# Discrete compound Poisson by Panjer recursion -----------------------------
dcp_pmf <- function(N, lambda, alpha) {
  K <- length(alpha); f <- numeric(N + 1); f[1] <- exp(-lambda)
  if (N > 0) for (n in 1:N) { k <- 1:min(n, K); f[n + 1] <- lambda / n * sum(k * alpha[k] * f[n - k + 1]) }
  f
}
register_family(family_def(
  "discrete_compound_poisson", "Discrete compound Poisson (stuttering Poisson)", cnt, discrete = TRUE,
  variate_type = "scalar_vector_parameter",
  pars = list(cnt_par("lambda", role = "rate"), P("alpha", "unit", role = "jump probabilities", regression = FALSE,
              interpretation = "alpha_k = P(jump size k), summing to one")),
  formula = "pgf exp(lambda sum_k alpha_k (z^k - 1)); pmf by Panjer recursion f(n) = lambda/n sum_k k alpha_k f(n - k)",
  lpdf = function(y, p) { f <- dcp_pmf(max(c(y, 0)), p$lambda, p$alpha); out <- rep(-Inf, length(y))
    ok <- y >= 0 & y == round(y); out[ok] <- log(f[y[ok] + 1]); out },
  rng = function(n, p) vapply(stats::rpois(n, p$lambda), function(m) sum(sample.int(length(p$alpha), m, TRUE, p$alpha)), 0),
  references = "Panjer (1981), Recursive evaluation of a family of compound distributions. ASTIN Bulletin, 12, 22-26.",
  reference = ref_impl("convolution", function(y, p) log(vapply(y, function(k) {
    tot <- 0; conv <- 1
    for (n in 0:(k + 60)) { val <- if (k + 1 <= length(conv)) conv[k + 1] else 0
      tot <- tot + stats::dpois(n, p$lambda) * val; conv <- stats::convolve(conv, rev(c(0, p$alpha)), type = "open") }
    tot }, 0)), note = "Direct sum over the Poisson count of convolution powers."),
  status = "validated"))

# Hardy (golf hole scores) ------------------------------------------------------
hardy_pmf1 <- function(n, p, q, m) {
  js <- if (m %% 2 == 1) ((m + 1) / 2):m else (m / 2):m
  r <- 1 - p - q
  tot <- 0
  for (j in js) {
    if (n - j < 0) next
    A <- if (2 * j - m - 1 >= 0) choose(j - 1, 2 * j - m - 1) * p^(m - j + 1) * r^(2 * j - m - 1) else 0
    B <- choose(j, 2 * j - m) * p^(m - j) * r^(2 * j - m)
    tot <- tot + choose(n - 1, n - j) * q^(n - j) * (A + B)
  }
  tot
}
register_family(family_def(
  "hardy", "Hardy (golf hole score)", cnt, discrete = TRUE,
  pars = list(P("p", "unit", role = "good shot probability"), P("q", "unit", role = "bad shot probability"),
              P("m", "positive_integer", role = "par", data = TRUE)),
  check = function(pp) if (any(pp$p + pp$q >= 1, na.rm = TRUE)) "Need p + q < 1.",
  formula = "Hardy (1945) model of hole scores: shots are good (p), bad (q) or ordinary (1 - p - q); see Wikipedia 'Hardy distribution'",
  lpdf = function(y, pp) log(mapply(hardy_pmf1, y, pp$p, pp$q, pp$m)),
  mean = function(pp) mapply(function(p, q, m) -sum((m + 1 - 1:m) * p^(0:(m - 1)) / (q - 1)^(1:m)), pp$p, pp$q, pp$m),
  support = function(pp) list(lower = 1, upper = Inf),
  references = "Hardy (1945), A mathematical study of golf hole scores. Mathematical Gazette.",
  reference = ref_impl("definition", function(y, pp) log(mapply(hardy_pmf1, y, pp$p, pp$q, pp$m))),
  status = "implemented"))

register_family(family_def(
  "boltzmann", "Boltzmann (Gibbs) distribution over discrete states", "categorical", discrete = TRUE,
  variate_type = "scalar_vector_parameter",
  pars = list(P("energies", "real", role = "state energies", regression = FALSE),
              cnt_par("kT", role = "temperature", interpretation = "Boltzmann constant times temperature")),
  formula = "P(state i) = exp(-e_i / kT) / sum_j exp(-e_j / kT)",
  lpdf = function(y, p) { l <- -p$energies / p$kT; l <- l - log_sum_exp(l); out <- rep(-Inf, length(y))
    ok <- y >= 1 & y <= length(l) & y == round(y); out[ok] <- l[y[ok]]; out },
  rng = function(n, p) sample.int(length(p$energies), n, TRUE, exp(-p$energies / p$kT)),
  references = "Landau and Lifshitz (1980), Statistical Physics, 3rd ed.",
  reference = ref_impl("categorical", function(y, p) log(exp(-p$energies / p$kT) / sum(exp(-p$energies / p$kT)))[y]),
  status = "validated"))

# Discrete phase type -----------------------------------------------------------
register_family(family_def(
  "discrete_phase_type", "Discrete phase-type", cnt, discrete = TRUE, variate_type = "scalar_matrix_parameter",
  pars = list(P("alpha", "unit_closed", role = "initial distribution", regression = FALSE),
              P("Tmat", "unit_closed", role = "substochastic transition matrix", regression = FALSE)),
  formula = "P(Y = k) = alpha T^(k - 1) t, t = (I - T) 1, k = 1, 2, ...",
  lpdf = function(y, p) { Tm <- as.matrix(p$Tmat); t0 <- 1 - rowSums(Tm); a <- p$alpha
    vapply(y, function(k) { if (k < 1 || k != round(k)) return(-Inf); v <- a
      if (k > 1) for (i in seq_len(k - 1)) v <- v %*% Tm
      log(sum(v * t0)) }, 0) },
  mean = function(p) { Tm <- as.matrix(p$Tmat); sum(p$alpha %*% solve(diag(nrow(Tm)) - Tm)) },
  rng = function(n, p) { Tm <- as.matrix(p$Tmat); K <- nrow(Tm); P0 <- cbind(Tm, 1 - rowSums(Tm))
    vapply(seq_len(n), function(i) { s <- sample.int(K, 1, prob = p$alpha); k <- 1
      repeat { s <- sample.int(K + 1, 1, prob = P0[s, ]); if (s == K + 1) return(k); k <- k + 1 } }, 0) },
  references = "Neuts (1981), Matrix-Geometric Solutions in Stochastic Models. Johns Hopkins.",
  reference = ref_impl("geometric when 1 x 1", function(y, p) log(p$alpha * p$Tmat^(y - 1) * (1 - p$Tmat))),
  status = "validated"))

# Zero truncated laws registered by name ------------------------------------------
register_family(modify_family(truncate_dist("poisson", lower = 1),
  id = "zero_truncated_poisson", display_name = "Zero-truncated Poisson", aliases = c("ZTP", "positive_poisson"),
  reference = ref_impl("extraDistr", function(y, p) extraDistr::dtpois(y, p$lambda, a = 0, log = TRUE)),
  status = "validated"))
register_family(modify_family(truncate_dist("negative_binomial_nb2", lower = 1),
  id = "zero_truncated_negative_binomial", display_name = "Zero-truncated negative binomial (NB2)", aliases = "ZTNB",
  reference = ref_impl("countreg form", function(y, p) stats::dnbinom(y, size = p$phi, mu = p$mu, log = TRUE) -
    log1p(-stats::dnbinom(0, size = p$phi, mu = p$mu))),
  status = "validated"))
