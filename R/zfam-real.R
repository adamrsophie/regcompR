# Continuous distributions on the real line.

rl <- "real"
std_pcdf <- function(z) stats::pnorm(z)

register_family(family_def(
  "normal", "Normal", rl, aliases = c("gaussian", "NO"),
  pars = list(loc_par("mu", TRUE), scale_par("sigma")),
  formula = "f(y) = exp(-(y - mu)^2 / (2 sigma^2)) / (sigma sqrt(2 pi))",
  mean_formula = "mu", variance_formula = "sigma^2",
  lpdf = function(y, p) stats::dnorm(y, p$mu, p$sigma, log = TRUE),
  cdf = function(q, p) stats::pnorm(q, p$mu, p$sigma),
  lccdf = function(q, p) stats::pnorm(q, p$mu, p$sigma, lower.tail = FALSE, log.p = TRUE),
  quantile = function(pr, p) stats::qnorm(pr, p$mu, p$sigma),
  rng = function(n, p) stats::rnorm(n, p$mu, p$sigma),
  mean = function(p) p$mu, variance = function(p) p$sigma^2,
  skewness = function(p) 0 * p$mu, kurtosis = function(p) 0 * p$mu,
  references = "Johnson, Kotz and Balakrishnan (1994), Continuous Univariate Distributions, vol. 1, ch. 13.",
  reference = ref_impl("gamlss.dist", gd("dNO"), gp("pNO")),
  stan = stan_native("normal_lpdf({y} | {mu}, {sigma})", "normal_lcdf({y} | {mu}, {sigma})",
                     "normal_lccdf({y} | {mu}, {sigma})", "normal_rng({mu}, {sigma})"),
  bugs = list(native = "dnorm({mu}, 1 / pow({sigma}, 2))",
              loglik = "-0.5 * log(2 * 3.141592653589793) - log({sigma}) - 0.5 * pow(({y} - {mu}) / {sigma}, 2)"),
  tmb = "normal", status = "validated"
))

register_family(family_def(
  "student_t", "Student t", rl, aliases = c("t", "TF"),
  pars = list(loc_par("mu"), scale_par("sigma"), P("nu", "positive", role = "tail", symbol = "nu",
                                                    interpretation = "Degrees of freedom")),
  formula = "f(y) = Gamma((nu + 1)/2) / (Gamma(nu/2) sqrt(nu pi) sigma) (1 + z^2/nu)^(-(nu + 1)/2), z = (y - mu)/sigma",
  mean_formula = "mu for nu > 1", variance_formula = "sigma^2 nu / (nu - 2) for nu > 2",
  moment_conditions = "Mean exists iff nu > 1; variance iff nu > 2; skewness iff nu > 3; kurtosis iff nu > 4.",
  lpdf = function(y, p) stats::dt((y - p$mu) / p$sigma, p$nu, log = TRUE) - log(p$sigma),
  cdf = function(q, p) stats::pt((q - p$mu) / p$sigma, p$nu),
  quantile = function(pr, p) p$mu + p$sigma * stats::qt(pr, p$nu),
  rng = function(n, p) p$mu + p$sigma * stats::rt(n, p$nu),
  mean = function(p) ifelse(p$nu > 1, p$mu, NA_real_),
  variance = function(p) ifelse(p$nu > 2, p$sigma^2 * p$nu / (p$nu - 2), NA_real_),
  skewness = function(p) ifelse(p$nu > 3, 0, NA_real_),
  kurtosis = function(p) ifelse(p$nu > 4, 6 / (p$nu - 4), NA_real_),
  references = "Johnson, Kotz and Balakrishnan (1995), Continuous Univariate Distributions, vol. 2, ch. 28.",
  reference = ref_impl("gamlss.dist", gd("dTF"), gp("pTF")),
  stan = stan_native("student_t_lpdf({y} | {nu}, {mu}, {sigma})", "student_t_lcdf({y} | {nu}, {mu}, {sigma})",
                     "student_t_lccdf({y} | {nu}, {mu}, {sigma})", "student_t_rng({nu}, {mu}, {sigma})"),
  bugs = list(native = "dt({mu}, 1 / pow({sigma}, 2), {nu})",
              loglik = "loggam(({nu} + 1) / 2) - loggam({nu} / 2) - 0.5 * log({nu} * 3.141592653589793) - log({sigma}) - ({nu} + 1) / 2 * log(1 + pow(({y} - {mu}) / {sigma}, 2) / {nu})"),
  identifiability_notes = "nu is weakly identified when it is large; consider a prior or a fixed value.",
  tmb = "student_t",
  start = function(y) list(mu = stats::median(y), sigma = max(stats::mad(y), 1e-3), nu = 5),
  status = "validated"
))

register_family(family_def(
  "cauchy", "Cauchy", rl,
  pars = list(loc_par("mu"), scale_par("sigma")),
  formula = "f(y) = 1 / (pi sigma (1 + z^2)), z = (y - mu)/sigma",
  mean_formula = "does not exist", variance_formula = "does not exist",
  moment_conditions = "No moments of order 1 or higher exist.",
  lpdf = function(y, p) stats::dcauchy(y, p$mu, p$sigma, log = TRUE),
  cdf = function(q, p) stats::pcauchy(q, p$mu, p$sigma),
  quantile = function(pr, p) stats::qcauchy(pr, p$mu, p$sigma),
  rng = function(n, p) stats::rcauchy(n, p$mu, p$sigma),
  mean = function(p) rep(NA_real_, length(p$mu)), variance = function(p) rep(NA_real_, length(p$mu)),
  skewness = function(p) rep(NA_real_, length(p$mu)), kurtosis = function(p) rep(NA_real_, length(p$mu)),
  references = "Johnson, Kotz and Balakrishnan (1994), vol. 1, ch. 16.",
  reference = ref_impl("stats", function(y, p) stats::dcauchy(y, p$mu, p$sigma, log = TRUE)),
  stan = stan_native("cauchy_lpdf({y} | {mu}, {sigma})", "cauchy_lcdf({y} | {mu}, {sigma})",
                     "cauchy_lccdf({y} | {mu}, {sigma})", "cauchy_rng({mu}, {sigma})"),
  bugs = list(native = "dt({mu}, 1 / pow({sigma}, 2), 1)",
              loglik = "-log(3.141592653589793 * {sigma}) - log(1 + pow(({y} - {mu}) / {sigma}, 2))"),
  start = function(y) list(mu = stats::median(y), sigma = max(stats::IQR(y) / 2, 1e-3)),
  status = "validated"
))

register_family(family_def(
  "logistic", "Logistic", rl, aliases = "LO",
  pars = list(loc_par("mu", TRUE), scale_par("sigma")),
  formula = "f(y) = exp(-z) / (sigma (1 + exp(-z))^2), z = (y - mu)/sigma",
  mean_formula = "mu", variance_formula = "sigma^2 pi^2 / 3",
  lpdf = function(y, p) stats::dlogis(y, p$mu, p$sigma, log = TRUE),
  cdf = function(q, p) stats::plogis(q, p$mu, p$sigma),
  quantile = function(pr, p) stats::qlogis(pr, p$mu, p$sigma),
  rng = function(n, p) stats::rlogis(n, p$mu, p$sigma),
  mean = function(p) p$mu, variance = function(p) p$sigma^2 * pi^2 / 3,
  skewness = function(p) 0 * p$mu, kurtosis = function(p) 1.2 + 0 * p$mu,
  references = "Johnson, Kotz and Balakrishnan (1995), vol. 2, ch. 23.",
  reference = ref_impl("gamlss.dist", gd("dLO"), gp("pLO")),
  stan = stan_native("logistic_lpdf({y} | {mu}, {sigma})", "logistic_lcdf({y} | {mu}, {sigma})",
                     "logistic_lccdf({y} | {mu}, {sigma})", "logistic_rng({mu}, {sigma})"),
  bugs = list(native = "dlogis({mu}, 1 / {sigma})",
              loglik = "-({y} - {mu}) / {sigma} - log({sigma}) - 2 * log(1 + exp(-({y} - {mu}) / {sigma}))"),
  status = "validated"
))

register_family(family_def(
  "laplace", "Laplace (double exponential)", rl, aliases = c("double_exponential"),
  pars = list(loc_par("mu", TRUE), scale_par("sigma")),
  formula = "f(y) = exp(-|y - mu| / sigma) / (2 sigma)",
  mean_formula = "mu", variance_formula = "2 sigma^2",
  lpdf = function(y, p) -abs(y - p$mu) / p$sigma - log(2 * p$sigma),
  cdf = function(q, p) ifelse(q < p$mu, 0.5 * exp((q - p$mu) / p$sigma), 1 - 0.5 * exp(-(q - p$mu) / p$sigma)),
  quantile = function(pr, p) p$mu - p$sigma * sign(pr - 0.5) * log(1 - 2 * abs(pr - 0.5)),
  rng = function(n, p) p$mu + p$sigma * (stats::rexp(n) - stats::rexp(n)),
  mean = function(p) p$mu, variance = function(p) 2 * p$sigma^2,
  skewness = function(p) 0 * p$mu, kurtosis = function(p) 3 + 0 * p$mu,
  numerical_notes = "The density is not differentiable at mu; finite difference derivatives use one sided steps there.",
  references = "Kotz, Kozubowski and Podgorski (2001), The Laplace Distribution and Generalizations. Birkhauser.",
  reference = ref_impl("extraDistr", function(y, p) extraDistr::dlaplace(y, p$mu, p$sigma, log = TRUE)),
  stan = stan_native("double_exponential_lpdf({y} | {mu}, {sigma})", "double_exponential_lcdf({y} | {mu}, {sigma})",
                     "double_exponential_lccdf({y} | {mu}, {sigma})", "double_exponential_rng({mu}, {sigma})"),
  bugs = list(native = "ddexp({mu}, 1 / {sigma})", loglik = "-abs({y} - {mu}) / {sigma} - log(2 * {sigma})"),
  status = "validated"
))

register_family(family_def(
  "asymmetric_laplace", "Asymmetric Laplace (quantile parameterization)", rl, aliases = "ald",
  pars = list(loc_par("mu"), scale_par("sigma"),
              P("tau", "unit", role = "quantile", interpretation = "Quantile level; mu is the tau quantile")),
  formula = "f(y) = tau (1 - tau) / sigma exp(-rho_tau((y - mu)/sigma)), rho_tau(u) = u (tau - I(u < 0))",
  mean_formula = "mu + sigma (1 - 2 tau) / (tau (1 - tau))",
  variance_formula = "sigma^2 (1 - 2 tau + 2 tau^2) / (tau^2 (1 - tau)^2)",
  regression_notes = "Commonly used as a working likelihood for quantile regression with tau fixed; likelihood based intervals then depend on the asymmetric Laplace assumption.",
  lpdf = function(y, p) {
    u <- (y - p$mu) / p$sigma
    log(p$tau) + log1p(-p$tau) - log(p$sigma) - u * (p$tau - (u < 0))
  },
  cdf = function(q, p) ifelse(q < p$mu, p$tau * exp((1 - p$tau) * (q - p$mu) / p$sigma),
                              1 - (1 - p$tau) * exp(-p$tau * (q - p$mu) / p$sigma)),
  quantile = function(pr, p) ifelse(pr < p$tau, p$mu + p$sigma / (1 - p$tau) * log(pr / p$tau),
                                    p$mu - p$sigma / p$tau * log((1 - pr) / (1 - p$tau))),
  mean = function(p) p$mu + p$sigma * (1 - 2 * p$tau) / (p$tau * (1 - p$tau)),
  variance = function(p) p$sigma^2 * (1 - 2 * p$tau + 2 * p$tau^2) / (p$tau^2 * (1 - p$tau)^2),
  references = "Yu and Zhang (2005), A three-parameter asymmetric Laplace distribution and its extension. Communications in Statistics, 34, 1867-1879.",
  reference = ref_impl("ald", function(y, p) log(mapply(ald::dALD, y, p$mu, p$sigma, p$tau))),
  stan = sfun("asymmetric_laplace", c("mu", "sigma", "tau"),
    "    real u = (y - mu) / sigma;\n    return log(tau) + log1m(tau) - log(sigma) - u * (tau - (u < 0));",
    lcdf = "    if (y < mu) return log(tau) + (1 - tau) * (y - mu) / sigma;\n    return log1m((1 - tau) * exp(-tau * (y - mu) / sigma));"),
  bugs = list(loglik = "log({tau}) + log(1 - {tau}) - log({sigma}) - ({y} - {mu}) / {sigma} * ({tau} - step({mu} - {y}))"),
  start = function(y) list(mu = stats::median(y), sigma = max(stats::sd(y) / 4, 1e-3), tau = 0.5),
  status = "validated"
))

gn_cdf <- function(z, b) 0.5 + sign(z) * 0.5 * stats::pgamma(abs(z)^b, 1 / b)
gn_q <- function(pr, b) sign(pr - 0.5) * stats::qgamma(2 * abs(pr - 0.5), 1 / b)^(1 / b)

register_family(family_def(
  "generalized_normal", "Generalized normal (exponential power, Nadarajah form)", rl,
  pars = list(loc_par("mu", TRUE), P("alpha", "positive", role = "scale"),
              P("beta", "positive", role = "tail", interpretation = "Shape; 2 gives the normal, 1 the Laplace")),
  formula = "f(y) = beta / (2 alpha Gamma(1/beta)) exp(-(|y - mu| / alpha)^beta)",
  mean_formula = "mu", variance_formula = "alpha^2 Gamma(3/beta) / Gamma(1/beta)",
  lpdf = function(y, p) log(p$beta) - log(2 * p$alpha) - lgamma(1 / p$beta) - (abs(y - p$mu) / p$alpha)^p$beta,
  cdf = function(q, p) gn_cdf((q - p$mu) / p$alpha, p$beta),
  quantile = function(pr, p) p$mu + p$alpha * gn_q(pr, p$beta),
  mean = function(p) p$mu, variance = function(p) p$alpha^2 * exp(lgamma(3 / p$beta) - lgamma(1 / p$beta)),
  kurtosis = function(p) exp(lgamma(5 / p$beta) + lgamma(1 / p$beta) - 2 * lgamma(3 / p$beta)) - 3,
  skewness = function(p) 0 * p$mu,
  references = "Nadarajah (2005), A generalized normal distribution. Journal of Applied Statistics, 32, 685-694.",
  reference = ref_impl("gamlss.dist", function(y, p) gamlss.dist::dPE2(y, p$mu, p$alpha, p$beta, log = TRUE)),
  stan = sfun("generalized_normal", c("mu", "alpha", "beta"),
    "    return log(beta) - log(2 * alpha) - lgamma(1 / beta) - pow(abs(y - mu) / alpha, beta);",
    lcdf = "    real z = (y - mu) / alpha;\n    return log(0.5 + (z > 0 ? 1 : -1) * 0.5 * gamma_p(1 / beta, pow(abs(z), beta)));"),
  bugs = list(loglik = "log({beta}) - log(2 * {alpha}) - loggam(1 / {beta}) - pow(abs({y} - {mu}) / {alpha}, {beta})"),
  start = function(y) list(mu = mean(y), alpha = stats::sd(y), beta = 2),
  status = "validated"
))

register_family(family_def(
  "power_exponential_type2", "Power exponential type 2 (GAMLSS PE2)", rl, aliases = "PE2",
  pars = list(loc_par("mu", TRUE), scale_par("sigma"), P("nu", "positive", role = "tail")),
  formula = "f(y) = nu / (2 sigma Gamma(1/nu)) exp(-|z|^nu), z = (y - mu)/sigma (same law as generalized_normal)",
  mean_formula = "mu", variance_formula = "sigma^2 Gamma(3/nu) / Gamma(1/nu)",
  lpdf = function(y, p) log(p$nu) - log(2 * p$sigma) - lgamma(1 / p$nu) - abs((y - p$mu) / p$sigma)^p$nu,
  cdf = function(q, p) gn_cdf((q - p$mu) / p$sigma, p$nu),
  quantile = function(pr, p) p$mu + p$sigma * gn_q(pr, p$nu),
  mean = function(p) p$mu, variance = function(p) p$sigma^2 * exp(lgamma(3 / p$nu) - lgamma(1 / p$nu)),
  references = GS, reference = ref_impl("gamlss.dist", gd("dPE2"), gp("pPE2")),
  stan = sfun("pe2", c("mu", "sigma", "nu"),
    "    return log(nu) - log(2 * sigma) - lgamma(1 / nu) - pow(abs(y - mu) / sigma, nu);",
    lcdf = "    real z = (y - mu) / sigma;\n    return log(0.5 + (z > 0 ? 1 : -1) * 0.5 * gamma_p(1 / nu, pow(abs(z), nu)));"),
  bugs = list(loglik = "log({nu}) - log(2 * {sigma}) - loggam(1 / {nu}) - pow(abs({y} - {mu}) / {sigma}, {nu})"),
  start = function(y) list(mu = mean(y), sigma = stats::sd(y), nu = 2),
  status = "validated"
))

pe_logc <- function(nu) 0.5 * (-(2 / nu) * log(2) + lgamma(1 / nu) - lgamma(3 / nu))
register_family(family_def(
  "power_exponential", "Power exponential (GAMLSS PE, sigma is the standard deviation)", rl, aliases = "PE",
  pars = list(loc_par("mu", TRUE), P("sigma", "positive", role = "scale", interpretation = "Standard deviation"),
              P("nu", "positive", role = "tail")),
  formula = "f(y) = nu exp(-0.5 |z / c|^nu) / (sigma c 2^(1 + 1/nu) Gamma(1/nu)), c^2 = 2^(-2/nu) Gamma(1/nu) / Gamma(3/nu)",
  mean_formula = "mu", variance_formula = "sigma^2",
  lpdf = function(y, p) {
    lc <- pe_logc(p$nu); z <- (y - p$mu) / p$sigma
    -log(p$sigma) + log(p$nu) - lc - 0.5 * abs(z / exp(lc))^p$nu - (1 + 1 / p$nu) * log(2) - lgamma(1 / p$nu)
  },
  cdf = function(q, p) {
    z <- (q - p$mu) / p$sigma; cc <- exp(pe_logc(p$nu))
    0.5 * (1 + stats::pgamma(0.5 * abs(z / cc)^p$nu, 1 / p$nu) * sign(z))
  },
  mean = function(p) p$mu, variance = function(p) p$sigma^2,
  references = GS, reference = ref_impl("gamlss.dist", gd("dPE"), gp("pPE")),
  stan = sfun("pe", c("mu", "sigma", "nu"),
    "    real lc = 0.5 * (-(2 / nu) * log(2) + lgamma(1 / nu) - lgamma(3 / nu));\n    real z = (y - mu) / sigma;\n    return -log(sigma) + log(nu) - lc - 0.5 * pow(abs(z / exp(lc)), nu) - (1 + 1 / nu) * log(2) - lgamma(1 / nu);",
    lcdf = "    real lc = 0.5 * (-(2 / nu) * log(2) + lgamma(1 / nu) - lgamma(3 / nu));\n    real z = (y - mu) / sigma;\n    return log(0.5 * (1 + (z > 0 ? 1 : -1) * gamma_p(1 / nu, 0.5 * pow(abs(z / exp(lc)), nu))));"),
  start = function(y) list(mu = mean(y), sigma = stats::sd(y), nu = 2),
  status = "validated"
))

register_family(family_def(
  "generalized_t", "Generalized t (GAMLSS GT)", rl, aliases = "GT",
  pars = list(loc_par("mu"), scale_par("sigma"), P("nu", "positive", role = "tail"), P("tau", "positive", role = "tail")),
  formula = "f(y) = tau / (2 sigma nu^(1/tau) B(1/tau, nu)) (1 + |z|^tau / nu)^(-(nu + 1/tau))",
  mean_formula = "mu for nu tau > 1", variance_formula = "sigma^2 nu^(2/tau) B(3/tau, nu - 2/tau) / B(1/tau, nu) for nu tau > 2",
  moment_conditions = "Mean exists iff nu tau > 1; variance iff nu tau > 2.",
  lpdf = function(y, p) {
    zt <- abs((y - p$mu) / p$sigma)^p$tau
    log(p$tau) - log(2 * p$sigma) - (1 / p$tau) * log(p$nu) - lgamma(1 / p$tau) - lgamma(p$nu) +
      lgamma(p$nu + 1 / p$tau) - (p$nu + 1 / p$tau) * log1p(zt / p$nu)
  },
  mean = function(p) ifelse(p$nu * p$tau > 1, p$mu, NA_real_),
  variance = function(p) ifelse(p$nu * p$tau > 2, p$sigma^2 * p$nu^(2 / p$tau) *
    exp(lbeta(3 / p$tau, p$nu - 2 / p$tau) - lbeta(1 / p$tau, p$nu)), NA_real_),
  references = "McDonald and Newey (1988), Partially adaptive estimation of regression models via the generalized t distribution. Econometric Theory, 4, 428-457.",
  reference = ref_impl("gamlss.dist", gd("dGT"), gp("pGT")),
  stan = sfun("gt", c("mu", "sigma", "nu", "tau"),
    "    real zt = pow(abs((y - mu) / sigma), tau);\n    return log(tau) - log(2 * sigma) - (1 / tau) * log(nu) - lgamma(1 / tau) - lgamma(nu) + lgamma(nu + 1 / tau) - (nu + 1 / tau) * log1p(zt / nu);"),
  bugs = list(loglik = "log({tau}) - log(2 * {sigma}) - (1 / {tau}) * log({nu}) - loggam(1 / {tau}) - loggam({nu}) + loggam({nu} + 1 / {tau}) - ({nu} + 1 / {tau}) * log(1 + pow(abs(({y} - {mu}) / {sigma}), {tau}) / {nu})"),
  start = function(y) list(mu = stats::median(y), sigma = stats::sd(y), nu = 5, tau = 2),
  status = "validated"
))

owens_t_r <- function(h, a) {
  vapply(seq_along(h), function(i) {
    if (!is.finite(h[i])) return(0)
    stats::integrate(function(x) exp(-0.5 * h[i]^2 * (1 + x^2)) / (1 + x^2), 0, a[i],
                     rel.tol = 1e-12)$value / (2 * pi)
  }, 0)
}

register_family(family_def(
  "skew_normal", "Skew normal (Azzalini)", rl, aliases = "SN",
  pars = list(loc_par("mu"), scale_par("sigma"), P("alpha", "real", role = "skewness")),
  formula = "f(y) = 2 / sigma phi(z) Phi(alpha z), z = (y - mu)/sigma",
  mean_formula = "mu + sigma delta sqrt(2/pi), delta = alpha / sqrt(1 + alpha^2)",
  variance_formula = "sigma^2 (1 - 2 delta^2 / pi)",
  regression_notes = "mu is a location parameter, not the mean, unless alpha = 0.",
  lpdf = function(y, p) {
    z <- (y - p$mu) / p$sigma
    log(2) - log(p$sigma) + stats::dnorm(z, log = TRUE) + stats::pnorm(p$alpha * z, log.p = TRUE)
  },
  cdf = function(q, p) {
    z <- (q - p$mu) / p$sigma
    pmin(pmax(stats::pnorm(z) - 2 * owens_t_r(z, rep_len(p$alpha, length(z))), 0), 1)
  },
  rng = function(n, p) {
    d <- p$alpha / sqrt(1 + p$alpha^2)
    u0 <- stats::rnorm(n); v <- stats::rnorm(n)
    p$mu + p$sigma * (d * abs(u0) + sqrt(1 - d^2) * v)
  },
  mean = function(p) p$mu + p$sigma * p$alpha / sqrt(1 + p$alpha^2) * sqrt(2 / pi),
  variance = function(p) p$sigma^2 * (1 - 2 * p$alpha^2 / (1 + p$alpha^2) / pi),
  references = "Azzalini (1985), A class of distributions which includes the normal ones. Scandinavian Journal of Statistics, 12, 171-178.",
  reference = ref_impl("sn", function(y, p) sn::dsn(y, p$mu, p$sigma, p$alpha, log = TRUE),
                       function(q, p) sn::psn(q, p$mu, p$sigma, p$alpha)),
  stan = stan_native("skew_normal_lpdf({y} | {mu}, {sigma}, {alpha})", "skew_normal_lcdf({y} | {mu}, {sigma}, {alpha})",
                     "skew_normal_lccdf({y} | {mu}, {sigma}, {alpha})", "skew_normal_rng({mu}, {sigma}, {alpha})"),
  bugs = list(loglik = "log(2) - log({sigma}) - 0.5 * log(2 * 3.141592653589793) - 0.5 * pow(({y} - {mu}) / {sigma}, 2) + log(phi({alpha} * ({y} - {mu}) / {sigma}))"),
  start = function(y) list(mu = mean(y), sigma = stats::sd(y), alpha = 0),
  status = "validated"
))

st_par <- function() list(loc_par("mu"), scale_par("sigma"), P("nu", "real", role = "skewness"),
                          P("tau", "positive", role = "tail"))
register_family(family_def(
  "skew_t_type1", "Skew t type 1 (GAMLSS ST1, Azzalini type)", rl, aliases = "ST1",
  pars = st_par(), formula = "f(y) = 2/sigma f_t(z; tau) F_t(nu z; tau)",
  lpdf = function(y, p) { z <- (y - p$mu) / p$sigma
    stats::pt(p$nu * z, p$tau, log.p = TRUE) + stats::dt(z, p$tau, log = TRUE) + log(2) - log(p$sigma) },
  references = GS, reference = ref_impl("gamlss.dist", gd("dST1")),
  stan = sfun("st1", c("mu", "sigma", "nu", "tau"),
    "    real z = (y - mu) / sigma;\n    return student_t_lcdf(nu * z | tau, 0, 1) + student_t_lpdf(z | tau, 0, 1) + log(2) - log(sigma);"),
  start = function(y) list(mu = mean(y), sigma = stats::sd(y), nu = 0, tau = 10),
  status = "validated"
))

register_family(family_def(
  "skew_t_type2", "Skew t type 2 (GAMLSS ST2, Azzalini and Capitanio)", rl, aliases = "ST2",
  pars = st_par(), formula = "f(y) = 2/sigma f_t(z; tau) F_t(nu sqrt((tau + 1)/(tau + z^2)) z; tau + 1)",
  lpdf = function(y, p) { z <- (y - p$mu) / p$sigma; lam <- (p$tau + 1) / (p$tau + z^2)
    stats::pt(p$nu * sqrt(lam) * z, p$tau + 1, log.p = TRUE) + stats::dt(z, p$tau, log = TRUE) + log(2) - log(p$sigma) },
  references = c(GS, "Azzalini and Capitanio (2003), JRSS B, 65, 367-389."),
  reference = ref_impl("gamlss.dist", gd("dST2")),
  stan = sfun("st2", c("mu", "sigma", "nu", "tau"),
    "    real z = (y - mu) / sigma;\n    real lam = (tau + 1) / (tau + z * z);\n    return student_t_lcdf(nu * sqrt(lam) * z | tau + 1, 0, 1) + student_t_lpdf(z | tau, 0, 1) + log(2) - log(sigma);"),
  start = function(y) list(mu = mean(y), sigma = stats::sd(y), nu = 0, tau = 10),
  status = "validated"
))

register_family(family_def(
  "skew_t_type3", "Skew t type 3 (GAMLSS ST3, Fernandez and Steel)", rl, aliases = "ST3",
  pars = list(loc_par("mu"), scale_par("sigma"), P("nu", "positive", role = "skewness"), P("tau", "positive", role = "tail")),
  formula = "f(y) = 2 nu / (sigma (1 + nu^2)) f_t(nu z) for y < mu, f_t(z / nu) otherwise",
  lpdf = function(y, p) { z <- (y - p$mu) / p$sigma
    ifelse(y < p$mu, stats::dt(p$nu * z, p$tau, log = TRUE), stats::dt(z / p$nu, p$tau, log = TRUE)) +
      log(2 * p$nu / (1 + p$nu^2)) - log(p$sigma) },
  cdf = function(q, p) { z <- (q - p$mu) / p$sigma
    ifelse(q < p$mu, 2 / (1 + p$nu^2) * stats::pt(p$nu * z, p$tau),
           1 / (1 + p$nu^2) + 2 * p$nu^2 / (1 + p$nu^2) * (stats::pt(z / p$nu, p$tau) - 0.5)) },
  references = c(GS, "Fernandez and Steel (1998), JASA, 93, 359-371."),
  reference = ref_impl("gamlss.dist", gd("dST3"), gp("pST3")),
  stan = sfun("st3", c("mu", "sigma", "nu", "tau"),
    "    real z = (y - mu) / sigma;\n    real lp = y < mu ? student_t_lpdf(nu * z | tau, 0, 1) : student_t_lpdf(z / nu | tau, 0, 1);\n    return lp + log(2 * nu / (1 + square(nu))) - log(sigma);",
    lcdf = "    real z = (y - mu) / sigma;\n    if (y < mu) return log(2 / (1 + square(nu))) + student_t_lcdf(nu * z | tau, 0, 1);\n    return log(1 / (1 + square(nu)) + 2 * square(nu) / (1 + square(nu)) * (student_t_cdf(z / nu | tau, 0, 1) - 0.5));"),
  start = function(y) list(mu = mean(y), sigma = stats::sd(y), nu = 1, tau = 10),
  status = "validated"
))

register_family(family_def(
  "skew_t_type4", "Skew t type 4 (GAMLSS ST4, spliced t)", rl, aliases = "ST4",
  pars = list(loc_par("mu"), scale_par("sigma"), P("nu", "positive", role = "tail"), P("tau", "positive", role = "tail")),
  formula = "Spliced Student t: df nu below mu, df tau above mu, with continuity at mu",
  lpdf = function(y, p) {
    lk1 <- lgamma(0.5) + lgamma(p$nu / 2) - lgamma((p$nu + 1) / 2) + 0.5 * log(p$nu)
    lk2 <- lgamma(0.5) + lgamma(p$tau / 2) - lgamma((p$tau + 1) / 2) + 0.5 * log(p$tau)
    lk <- lk2 - lk1; k <- exp(lk); z <- (y - p$mu) / p$sigma
    ifelse(y < p$mu, stats::dt(z, p$nu, log = TRUE), stats::dt(z, p$tau, log = TRUE) + lk) + log(2 / (1 + k)) - log(p$sigma)
  },
  references = GS, reference = ref_impl("gamlss.dist", gd("dST4")),
  stan = sfun("st4", c("mu", "sigma", "nu", "tau"),
    "    real lk1 = lgamma(0.5) + lgamma(nu / 2) - lgamma((nu + 1) / 2) + 0.5 * log(nu);\n    real lk2 = lgamma(0.5) + lgamma(tau / 2) - lgamma((tau + 1) / 2) + 0.5 * log(tau);\n    real lk = lk2 - lk1;\n    real z = (y - mu) / sigma;\n    real lp = y < mu ? student_t_lpdf(z | nu, 0, 1) : student_t_lpdf(z | tau, 0, 1) + lk;\n    return lp + log(2 / (1 + exp(lk))) - log(sigma);"),
  start = function(y) list(mu = mean(y), sigma = stats::sd(y), nu = 10, tau = 10),
  status = "validated"
))

register_family(family_def(
  "skew_t_type5", "Skew t type 5 (GAMLSS ST5, Jones and Faddy)", rl, aliases = "ST5",
  pars = st_par(), formula = "Jones and Faddy (2003) skew t with a = (2/tau + lambda)/2, b = (2/tau - lambda)/2",
  lpdf = function(y, p) {
    z <- (y - p$mu) / p$sigma; ve <- 2 / p$tau; lam <- 2 * p$nu / (p$tau * sqrt(2 * p$tau + p$nu^2))
    a <- (ve + lam) / 2; b <- (ve - lam) / 2; s <- sqrt(a + b + z^2)
    (a + 0.5) * log1p(z / s) + (b + 0.5) * log1p(-z / s) - (a + b - 1) * log(2) - 0.5 * log(a + b) - lbeta(a, b) - log(p$sigma)
  },
  references = c(GS, "Jones and Faddy (2003), JRSS B, 65, 159-174."),
  reference = ref_impl("gamlss.dist", gd("dST5")),
  stan = sfun("st5", c("mu", "sigma", "nu", "tau"),
    "    real z = (y - mu) / sigma;\n    real ve = 2 / tau;\n    real lam = 2 * nu / (tau * sqrt(2 * tau + square(nu)));\n    real a = (ve + lam) / 2;\n    real b = (ve - lam) / 2;\n    real s = sqrt(a + b + square(z));\n    return (a + 0.5) * log1p(z / s) + (b + 0.5) * log1p(-z / s) - (a + b - 1) * log(2) - 0.5 * log(a + b) - lbeta(a, b) - log(sigma);"),
  start = function(y) list(mu = mean(y), sigma = stats::sd(y), nu = 0, tau = 1),
  status = "validated"
))

register_family(family_def(
  "sinh_arcsinh", "Sinh arcsinh (Jones and Pewsey, GAMLSS SHASHo)", rl, aliases = c("SHASHo", "SHASH"),
  pars = list(loc_par("mu"), scale_par("sigma"), P("nu", "real", role = "skewness"), P("tau", "positive", role = "tail")),
  formula = "f(y) = tau c / (sigma sqrt(2 pi (1 + z^2))) exp(-r^2/2), r = sinh(tau asinh(z) - nu), c = cosh(tau asinh(z) - nu)",
  lpdf = function(y, p) {
    z <- (y - p$mu) / p$sigma; a <- p$tau * asinh(z) - p$nu
    -log(p$sigma) + log(p$tau) - 0.5 * log(2 * pi) - 0.5 * log1p(z^2) + log(cosh(a)) - 0.5 * sinh(a)^2
  },
  cdf = function(q, p) stats::pnorm(sinh(p$tau * asinh((q - p$mu) / p$sigma) - p$nu)),
  quantile = function(pr, p) p$mu + p$sigma * sinh((asinh(stats::qnorm(pr)) + p$nu) / p$tau),
  references = "Jones and Pewsey (2009), Sinh-arcsinh distributions. Biometrika, 96, 761-780.",
  reference = ref_impl("gamlss.dist", gd("dSHASHo"), gp("pSHASHo")),
  stan = sfun("shasho", c("mu", "sigma", "nu", "tau"),
    "    real z = (y - mu) / sigma;\n    real a = tau * asinh(z) - nu;\n    return -log(sigma) + log(tau) - 0.5 * log(2 * pi()) - 0.5 * log1p(square(z)) + log(cosh(a)) - 0.5 * square(sinh(a));",
    lcdf = "    return std_normal_lcdf(sinh(tau * asinh((y - mu) / sigma) - nu));",
    lccdf = "    return std_normal_lccdf(sinh(tau * asinh((y - mu) / sigma) - nu));"),
  bugs = list(loglik = "-log({sigma}) + log({tau}) - 0.5 * log(2 * 3.141592653589793) - 0.5 * log(1 + pow(({y} - {mu}) / {sigma}, 2)) + log(cosh({tau} * log(({y} - {mu}) / {sigma} + sqrt(1 + pow(({y} - {mu}) / {sigma}, 2))) - {nu})) - 0.5 * pow(sinh({tau} * log(({y} - {mu}) / {sigma} + sqrt(1 + pow(({y} - {mu}) / {sigma}, 2))) - {nu}), 2)"),
  start = function(y) list(mu = mean(y), sigma = stats::sd(y), nu = 0, tau = 1),
  status = "validated"
))

jsu_aux <- function(p) {
  rtau <- 1 / p$tau; w <- ifelse(rtau < 1e-7, 1, exp(rtau^2)); omega <- -p$nu * rtau
  cc <- (0.5 * (w - 1) * (w * cosh(2 * omega) + 1))^(-0.5)
  list(rtau = rtau, w = w, omega = omega, c = cc)
}
register_family(family_def(
  "johnson_su", "Johnson SU (GAMLSS JSU, mean and standard deviation parameterization)", rl, aliases = "JSU",
  pars = list(loc_par("mu", TRUE), P("sigma", "positive", role = "scale", interpretation = "Standard deviation"),
              P("nu", "real", role = "skewness"), P("tau", "positive", role = "tail")),
  formula = "Johnson SU reparameterized so that mu is the mean and sigma the standard deviation",
  mean_formula = "mu", variance_formula = "sigma^2",
  lpdf = function(y, p) {
    a <- jsu_aux(p)
    z <- (y - (p$mu + a$c * p$sigma * sqrt(a$w) * sinh(a$omega))) / (a$c * p$sigma)
    r <- -p$nu + asinh(z) / a$rtau
    -log(p$sigma) - log(a$c) - log(a$rtau) - 0.5 * log1p(z^2) - 0.5 * log(2 * pi) - 0.5 * r^2
  },
  cdf = function(q, p) { a <- jsu_aux(p)
    z <- (q - (p$mu + a$c * p$sigma * sqrt(a$w) * sinh(a$omega))) / (a$c * p$sigma)
    stats::pnorm(-p$nu + asinh(z) / a$rtau) },
  quantile = function(pr, p) { a <- jsu_aux(p)
    z <- sinh(a$rtau * (stats::qnorm(pr) + p$nu))
    p$mu + a$c * p$sigma * sqrt(a$w) * sinh(a$omega) + a$c * p$sigma * z },
  mean = function(p) p$mu, variance = function(p) p$sigma^2,
  references = c(GS, "Johnson (1949), Biometrika, 36, 149-176."),
  reference = ref_impl("gamlss.dist", gd("dJSU"), gp("pJSU")),
  stan = sfun("jsu", c("mu", "sigma", "nu", "tau"),
    "    real rtau = 1 / tau;\n    real w = exp(square(rtau));\n    real omega = -nu * rtau;\n    real c = pow(0.5 * (w - 1) * (w * cosh(2 * omega) + 1), -0.5);\n    real z = (y - (mu + c * sigma * sqrt(w) * sinh(omega))) / (c * sigma);\n    real r = -nu + asinh(z) / rtau;\n    return -log(sigma) - log(c) - log(rtau) - 0.5 * log1p(square(z)) - 0.5 * log(2 * pi()) - 0.5 * square(r);",
    lcdf = "    real rtau = 1 / tau;\n    real w = exp(square(rtau));\n    real omega = -nu * rtau;\n    real c = pow(0.5 * (w - 1) * (w * cosh(2 * omega) + 1), -0.5);\n    real z = (y - (mu + c * sigma * sqrt(w) * sinh(omega))) / (c * sigma);\n    return std_normal_lcdf(-nu + asinh(z) / rtau);"),
  start = function(y) list(mu = mean(y), sigma = stats::sd(y), nu = 0, tau = 1),
  status = "validated"
))

register_family(family_def(
  "ex_gaussian", "Exponentially modified Gaussian (ex-Gaussian, GAMLSS exGAUS)", rl, aliases = c("exGAUS", "exgaussian"),
  pars = list(loc_par("mu"), scale_par("sigma"), P("nu", "positive", role = "scale", interpretation = "Mean of the exponential component")),
  formula = "Y = N(mu, sigma^2) + Exp(mean nu): f(y) = exp(-(y - mu)/nu + sigma^2/(2 nu^2)) Phi((y - mu)/sigma - sigma/nu) / nu",
  mean_formula = "mu + nu", variance_formula = "sigma^2 + nu^2",
  lpdf = function(y, p) {
    z <- y - p$mu - p$sigma^2 / p$nu
    -log(p$nu) - (z + p$sigma^2 / (2 * p$nu)) / p$nu + stats::pnorm(z / p$sigma, log.p = TRUE)
  },
  cdf = function(q, p) stats::pnorm((q - p$mu) / p$sigma) -
    exp(-(q - p$mu) / p$nu + p$sigma^2 / (2 * p$nu^2) + stats::pnorm((q - p$mu) / p$sigma - p$sigma / p$nu, log.p = TRUE)),
  rng = function(n, p) stats::rnorm(n, p$mu, p$sigma) + stats::rexp(n, 1 / p$nu),
  mean = function(p) p$mu + p$nu, variance = function(p) p$sigma^2 + p$nu^2,
  references = c(GS, "Stan Functions Reference, exp_mod_normal."),
  reference = ref_impl("gamlss.dist", gd("dexGAUS"), gp("pexGAUS")),
  stan = stan_native("exp_mod_normal_lpdf({y} | {mu}, {sigma}, inv({nu}))", "exp_mod_normal_lcdf({y} | {mu}, {sigma}, inv({nu}))",
                     "exp_mod_normal_lccdf({y} | {mu}, {sigma}, inv({nu}))", "exp_mod_normal_rng({mu}, {sigma}, inv({nu}))"),
  bugs = list(loglik = "-log({nu}) - ({y} - {mu} - pow({sigma}, 2) / {nu} + pow({sigma}, 2) / (2 * {nu})) / {nu} + log(phi(({y} - {mu} - pow({sigma}, 2) / {nu}) / {sigma}))"),
  start = function(y) list(mu = mean(y) - stats::sd(y) / 2, sigma = stats::sd(y) / 2, nu = stats::sd(y) / 2),
  status = "validated"
))

mills <- function(z) exp(stats::pnorm(z, lower.tail = FALSE, log.p = TRUE) - stats::dnorm(z, log = TRUE))
register_family(family_def(
  "normal_laplace", "Normal Laplace (Reed and Jorgensen)", rl,
  pars = list(loc_par("mu"), scale_par("sigma"), P("alpha", "positive", role = "tail"), P("beta", "positive", role = "tail")),
  formula = "Y = N(mu, sigma^2) + asymmetric Laplace(alpha, beta): f(y) = alpha beta/(alpha + beta) phi(z) [R(alpha sigma - z) + R(beta sigma + z)], R the Mills ratio",
  mean_formula = "mu + 1/alpha - 1/beta", variance_formula = "sigma^2 + 1/alpha^2 + 1/beta^2",
  lpdf = function(y, p) {
    z <- (y - p$mu) / p$sigma
    l1 <- stats::pnorm(p$alpha * p$sigma - z, lower.tail = FALSE, log.p = TRUE)
    l2 <- stats::pnorm(p$beta * p$sigma + z, lower.tail = FALSE, log.p = TRUE)
    log(p$alpha) + log(p$beta) - log(p$alpha + p$beta) +
      log_add_exp(l1 + (p$alpha * p$sigma - z)^2 / 2 - z^2 / 2, l2 + (p$beta * p$sigma + z)^2 / 2 - z^2 / 2)
  },
  cdf = function(q, p) {
    z <- (q - p$mu) / p$sigma
    a <- p$alpha; b <- p$beta; s <- p$sigma
    stats::pnorm(z) - stats::dnorm(z) * (b * mills(a * s - z) - a * mills(b * s + z)) / (a + b)
  },
  rng = function(n, p) stats::rnorm(n, p$mu, p$sigma) + stats::rexp(n, p$alpha) - stats::rexp(n, p$beta),
  mean = function(p) p$mu + 1 / p$alpha - 1 / p$beta,
  variance = function(p) p$sigma^2 + 1 / p$alpha^2 + 1 / p$beta^2,
  references = "Reed and Jorgensen (2004), The double Pareto-lognormal distribution. Communications in Statistics, 33, 1733-1753.",
  reference = ref_impl("convolution", function(y, p) vapply(seq_along(y), function(i) log(stats::integrate(
    function(e) stats::dnorm(y[i] - e, p$mu[i], p$sigma[i]) * ifelse(e >= 0, p$alpha[i] * p$beta[i] / (p$alpha[i] + p$beta[i]) * exp(-p$alpha[i] * e),
                                                                   p$alpha[i] * p$beta[i] / (p$alpha[i] + p$beta[i]) * exp(p$beta[i] * e)), -Inf, Inf, rel.tol = 1e-12)$value), 0),
    note = "Numerical convolution of the normal and asymmetric Laplace densities."),
  stan = sfun("normal_laplace", c("mu", "sigma", "alpha", "beta"),
    "    real z = (y - mu) / sigma;\n    real l1 = std_normal_lccdf(alpha * sigma - z) + square(alpha * sigma - z) / 2;\n    real l2 = std_normal_lccdf(beta * sigma + z) + square(beta * sigma + z) / 2;\n    return log(alpha) + log(beta) - log(alpha + beta) + log_sum_exp(l1, l2) - square(z) / 2;"),
  start = function(y) list(mu = mean(y), sigma = stats::sd(y) / 2, alpha = 2 / stats::sd(y), beta = 2 / stats::sd(y)),
  status = "validated"
))
