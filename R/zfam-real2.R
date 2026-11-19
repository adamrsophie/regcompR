# Further continuous distributions on the real line (and variable support).

rfam <- function(id, display, pars, ...) register_family(family_def(id, display, "real", pars, ...))

# Quantile defined laws: density -1/Q'(u) at u = F(y) found by root finding.
qdef_cdf <- function(Q, lo = 1e-12) function(q, p) pervec(q, p, function(y, pp) {
  a <- Q(lo, pp); b <- Q(1 - lo, pp)
  if (y <= a) return(0); if (y >= b) return(1)
  stats::uniroot(function(u) Q(u, pp) - y, c(lo, 1 - lo), tol = 1e-13)$root })
qdef_lpdf <- function(Q, dQ) { cdfq <- qdef_cdf(Q)
  function(y, p) { u <- cdfq(y, p); out <- -log(dQ(pmin(pmax(u, 1e-12), 1 - 1e-12), p)); out[u <= 0 | u >= 1] <- -Inf; out } }

# Extreme value laws -------------------------------------------------------------
rfam("gumbel", "Gumbel (maximum)", list(loc_par("mu"), scale_par("beta")), aliases = c("extreme_value_type1", "GU_max"),
  formula = "F(y) = exp(-exp(-(y - mu)/beta))", mean_formula = "mu + 0.5772 beta",
  lpdf = function(y, p) { z <- (y - p$mu) / p$beta; -log(p$beta) - z - exp(-z) },
  cdf = function(q, p) exp(-exp(-(q - p$mu) / p$beta)), quantile = function(pr, p) p$mu - p$beta * log(-log(pr)),
  mean = function(p) p$mu + 0.5772156649015329 * p$beta, variance = function(p) pi^2 / 6 * p$beta^2,
  references = JKB, reference = ref_impl("extraDistr", function(y, p) extraDistr::dgumbel(y, p$mu, p$beta, log = TRUE)),
  stan = stan_native("gumbel_lpdf({y} | {mu}, {beta})", "gumbel_lcdf({y} | {mu}, {beta})", "gumbel_lccdf({y} | {mu}, {beta})", "gumbel_rng({mu}, {beta})"),
  status = "validated")

rfam("gumbel_min", "Gumbel (minimum, GAMLSS GU)", list(loc_par("mu"), scale_par("sigma")), aliases = "GU",
  lpdf = function(y, p) { z <- (y - p$mu) / p$sigma; -log(p$sigma) + z - exp(z) },
  cdf = function(q, p) 1 - exp(-exp((q - p$mu) / p$sigma)), quantile = function(pr, p) p$mu + p$sigma * log(-log1p(-pr)),
  mean = function(p) p$mu - 0.5772156649015329 * p$sigma,
  references = GS, reference = ref_impl("evd (reflected maximum Gumbel)", function(y, p) evd::dgumbel(-y, -p$mu, p$sigma, log = TRUE)), status = "validated")

gev_t <- function(y, p) { z <- (y - p$mu) / p$sigma
  ifelse(abs(p$xi) < 1e-9, exp(-z), pmax(1 + p$xi * z, 0)^(-1 / p$xi)) }
gev_support <- function(p) { b <- p$mu - p$sigma / p$xi
  list(lower = ifelse(p$xi > 0, b, -Inf), upper = ifelse(p$xi < 0, b, Inf)) }
rfam("gev", "Generalized extreme value", list(loc_par("mu"), scale_par("sigma"), P("xi", "real", role = "tail", link = "identity")),
  aliases = c("generalized_extreme_value", "fisher_tippett"), support = gev_support,
  formula = "F(y) = exp(-t(y)), t(y) = (1 + xi (y - mu)/sigma)^(-1/xi)",
  lpdf = function(y, p) { t <- gev_t(y, p); ok <- t > 0 & is.finite(t)
    ifelse(ok, -log(p$sigma) + (p$xi + 1) * log(t) - t, -Inf) },
  cdf = function(q, p) exp(-gev_t(q, p)),
  quantile = function(pr, p) ifelse(abs(p$xi) < 1e-9, p$mu - p$sigma * log(-log(pr)), p$mu + p$sigma * ((-log(pr))^(-p$xi) - 1) / p$xi),
  mean = function(p) ifelse(p$xi < 1, p$mu + p$sigma * (gamma(1 - p$xi) - 1) / p$xi, Inf),
  references = "Coles (2001), An Introduction to Statistical Modeling of Extreme Values. Springer.",
  reference = ref_impl("evd", function(y, p) mapply(function(x, a, b, c) evd::dgev(x, a, b, c, log = TRUE), y, p$mu, p$sigma, p$xi)),
  stan = sfun("gev", c("mu", "sigma", "xi"),
    "    real z = (y - mu) / sigma;\n    real t;\n    if (abs(xi) < 1e-9) t = exp(-z); else { if (1 + xi * z <= 0) return negative_infinity(); t = pow(1 + xi * z, -1 / xi); }\n    return -log(sigma) + (xi + 1) * log(t) - t;"),
  start = function(y) list(mu = mean(y), sigma = stats::sd(y), xi = 0.1), status = "validated")

rfam("generalized_pareto", "Generalized Pareto", list(P("mu", "real", role = "threshold", link = "identity", data = TRUE),
  scale_par("sigma"), P("xi", "real", role = "tail", link = "identity")), aliases = c("gpd", "GPD"),
  support = function(p) list(lower = p$mu, upper = ifelse(p$xi < 0, p$mu - p$sigma / p$xi, Inf)),
  formula = "F(y) = 1 - (1 + xi (y - mu)/sigma)^(-1/xi)",
  lpdf = function(y, p) { z <- (y - p$mu) / p$sigma
    ifelse(z < 0 | 1 + p$xi * z <= 0, -Inf, ifelse(abs(p$xi) < 1e-9, -log(p$sigma) - z, -log(p$sigma) - (1 / p$xi + 1) * log1p(p$xi * z))) },
  cdf = function(q, p) { z <- pmax((q - p$mu) / p$sigma, 0); ifelse(abs(p$xi) < 1e-9, 1 - exp(-z), 1 - pmax(1 + p$xi * z, 0)^(-1 / p$xi)) },
  quantile = function(pr, p) ifelse(abs(p$xi) < 1e-9, p$mu - p$sigma * log1p(-pr), p$mu + p$sigma * ((1 - pr)^(-p$xi) - 1) / p$xi),
  mean = function(p) ifelse(p$xi < 1, p$mu + p$sigma / (1 - p$xi), Inf),
  references = "Pickands (1975), Statistical inference using extreme order statistics. Ann. Stat., 3, 119-131.",
  reference = ref_impl("evd", function(y, p) mapply(function(x, a, b, c) evd::dgpd(x, a, b, c, log = TRUE), y, p$mu, p$sigma, p$xi)),
  stan = sfun("gpd", c("mu", "sigma", "xi"),
    "    real z = (y - mu) / sigma;\n    if (z < 0 || 1 + xi * z <= 0) return negative_infinity();\n    if (abs(xi) < 1e-9) return -log(sigma) - z;\n    return -log(sigma) - (1 / xi + 1) * log1p(xi * z);"),
  status = "validated")

rfam("extended_gpd", "Extended generalized Pareto (Naveau et al., model 1)", list(scale_par("sigma"), P("xi", "real", role = "tail", link = "identity"), shp("kappa")),
  aliases = "egpd", support = function(p) list(lower = 0, upper = Inf),
  formula = "F(y) = H(y)^kappa with H the GPD(0, sigma, xi) CDF",
  lpdf = function(y, p) { H <- 1 - pmax(1 + p$xi * y / p$sigma, 0)^(-1 / p$xi)
    log(p$kappa) + (p$kappa - 1) * log(H) - log(p$sigma) - (1 / p$xi + 1) * log1p(p$xi * y / p$sigma) },
  cdf = function(q, p) (1 - pmax(1 + p$xi * pmax(q, 0) / p$sigma, 0)^(-1 / p$xi))^p$kappa,
  quantile = function(pr, p) p$sigma * ((1 - pr^(1 / p$kappa))^(-p$xi) - 1) / p$xi,
  references = "Naveau, Huser, Ribereau and Hannart (2016), Modeling jointly low, moderate, and heavy rainfall intensities without a threshold selection. Water Resources Research, 52, 2753-2769.",
  reference = ref_impl("mev", function(y, p) log(p$kappa) + (p$kappa - 1) * log(evd::pgpd(y, 0, p$sigma, p$xi)) + mapply(function(x, b, c) evd::dgpd(x, 0, b, c, log = TRUE), y, p$sigma, p$xi)),
  status = "validated")

# Logistic type laws ----------------------------------------------------------------
rfam("generalized_logistic", "Generalized logistic (type IV, log-F)", list(loc_par("mu"), scale_par("sigma"), shp("p"), shp("q")),
  aliases = c("log_F", "type_iv_logistic"),
  formula = "f(y) = exp(-q z) / (sigma B(p, q) (1 + exp(-z))^(p + q)), z = (y - mu)/sigma",
  lpdf = function(y, r) { z <- (y - r$mu) / r$sigma; -r$q * z - (r$p + r$q) * log1pexp(-z) - log(r$sigma) - lbeta(r$p, r$q) },
  cdf = function(x, r) stats::pbeta(stats::plogis((x - r$mu) / r$sigma), r$p, r$q),
  quantile = function(pr, r) r$mu + r$sigma * stats::qlogis(stats::qbeta(pr, r$p, r$q)),
  references = "Prentice (1976), A generalization of the probit and logit methods for dose response curves. Biometrics, 32, 761-768.",
  reference = ref_impl("glogis", function(y, r) stats::dbeta(stats::plogis((y - r$mu) / r$sigma), r$p, r$q, log = TRUE) + stats::dlogis(y, r$mu, r$sigma, log = TRUE)),
  stan = sfun("glogis4", c("mu", "sigma", "p", "q"),
    "    real z = (y - mu) / sigma;\n    return -q * z - (p + q) * log1p_exp(-z) - log(sigma) - lbeta(p, q);"),
  status = "validated")

rfam("hyperbolic_secant", "Hyperbolic secant", list(loc_par("mu"), scale_par("sigma")),
  lpdf = function(y, p) { z <- pi * (y - p$mu) / (2 * p$sigma); -log(2 * p$sigma) - log(cosh(z)) },
  cdf = function(q, p) 2 / pi * atan(exp(pi * (q - p$mu) / (2 * p$sigma))),
  quantile = function(pr, p) p$mu + 2 * p$sigma / pi * log(tan(pi * pr / 2)),
  mean = function(p) p$mu, variance = function(p) p$sigma^2,
  references = "Baten (1934); Talacko (1956), Trabajos de Estadistica, 7, 159-174.",
  reference = ref_impl("definition", function(y, p) log(1 / (2 * p$sigma) / cosh(pi / 2 * (y - p$mu) / p$sigma))), status = "validated")

rfam("moyal", "Moyal", list(loc_par("mu"), scale_par("sigma")),
  lpdf = function(y, p) { z <- (y - p$mu) / p$sigma; -0.5 * (z + exp(-z)) - 0.5 * log(2 * pi) - log(p$sigma) },
  cdf = function(q, p) 2 * stats::pnorm(exp(-(q - p$mu) / (2 * p$sigma)), lower.tail = FALSE),
  quantile = function(pr, p) p$mu - 2 * p$sigma * log(stats::qnorm(1 - pr / 2)),
  mean = function(p) p$mu + p$sigma * (0.5772156649015329 + log(2)), variance = function(p) pi^2 / 2 * p$sigma^2,
  references = "Moyal (1955), Theory of ionization fluctuations. Phil. Mag., 46, 263-280.",
  reference = ref_impl("chi-squared transform", function(y, p) { z <- (y - p$mu) / p$sigma; stats::dchisq(exp(-z), 1, log = TRUE) - z - log(p$sigma) }),
  status = "validated")

rfam("log_gamma", "Generalized log-gamma (Prentice / Lawless)", list(loc_par("mu"), scale_par("sigma"), P("lambda", "real", role = "shape", link = "identity")),
  aliases = "generalized_log_gamma",
  formula = "W = (Y - mu)/sigma; lambda W + log(1/lambda^2) is log Gamma(1/lambda^2)",
  lpdf = function(y, p) { w <- (y - p$mu) / p$sigma; l <- p$lambda
    ifelse(abs(l) < 1e-6, stats::dnorm(w, log = TRUE),
      log(abs(l)) + (1 / l^2) * log(1 / l^2) + (1 / l^2) * (l * w - exp(l * w)) - lgamma(1 / l^2)) - log(p$sigma) },
  cdf = function(q, p) { w <- (q - p$mu) / p$sigma; l <- p$lambda; k <- 1 / l^2
    ifelse(abs(l) < 1e-6, stats::pnorm(w), ifelse(l > 0, stats::pgamma(k * exp(l * w), k), stats::pgamma(k * exp(l * w), k, lower.tail = FALSE))) },
  references = "Lawless (2003), Statistical Models and Methods for Lifetime Data, 2nd ed., sec. 5.3.",
  reference = ref_impl("flexsurv gengamma (log scale)", function(y, p) flexsurv::dgengamma(exp(y), p$mu, p$sigma, p$lambda, log = TRUE) + y),
  status = "validated")

# Noncentral and generalized t ----------------------------------------------------
rfam("noncentral_t", "Noncentral t", list(P("df", "positive"), P("ncp", "real", role = "noncentrality", link = "identity")),
  lpdf = function(y, p) stats::dt(y, p$df, p$ncp, log = TRUE), cdf = function(q, p) stats::pt(q, p$df, p$ncp),
  quantile = function(pr, p) stats::qt(pr, p$df, p$ncp),
  references = JKB, reference = ref_impl("integral", function(y, p) log(vapply(seq_along(y), function(i) stats::integrate(function(v)
    stats::dnorm(y[i] * sqrt(v / p$df[i]), p$ncp[i]) * sqrt(v / p$df[i]) * stats::dchisq(v, p$df[i]), 0, Inf, rel.tol = 1e-12)$value, 0))),
  status = "validated")

rfam("skewed_generalized_t", "Skewed generalized t (Theodossiou)", list(loc_par("mu"), scale_par("sigma"), P("lambda", "symmetric_unit", role = "skewness"),
  shp("p"), shp("q")), aliases = "SGT",
  formula = "f(y) = p / (2 sigma q^(1/p) B(1/p, q) (|y - mu|^p / (q sigma^p (1 + lambda sign(y - mu))^p) + 1)^(1/p + q))",
  lpdf = function(y, r) { d <- y - r$mu
    log(r$p) - log(2 * r$sigma) - log(r$q) / r$p - lbeta(1 / r$p, r$q) -
      (1 / r$p + r$q) * log1p(abs(d)^r$p / (r$q * r$sigma^r$p * (1 + r$lambda * sign(d))^r$p)) },
  references = "Theodossiou (1998), Financial data and the skewed generalized t distribution. Management Science, 44, 1650-1661.",
  reference = ref_impl("sgt", function(y, r) sgt::dsgt(y, r$mu, r$sigma, r$lambda, r$p, r$q, mean.cent = FALSE, var.adj = FALSE, log = TRUE)),
  status = "validated")

# Normal variance mean mixtures ----------------------------------------------------
rfam("nig", "Normal inverse Gaussian", list(loc_par("mu"), shp("alpha"), P("beta", "real", role = "skewness", link = "identity"), scale_par("delta")),
  aliases = "normal_inverse_gaussian", check = function(p) if (any(abs(p$beta) >= p$alpha, na.rm = TRUE)) "Need |beta| < alpha.",
  formula = "f(y) = alpha delta K1(alpha s) / (pi s) exp(delta gamma + beta (y - mu)), s = sqrt(delta^2 + (y - mu)^2), gamma = sqrt(alpha^2 - beta^2)",
  lpdf = function(y, p) { s <- sqrt(p$delta^2 + (y - p$mu)^2); g <- sqrt(p$alpha^2 - p$beta^2)
    log(p$alpha * p$delta) + log_besselK(p$alpha * s, 1) - log(pi * s) + p$delta * g + p$beta * (y - p$mu) },
  mean = function(p) p$mu + p$delta * p$beta / sqrt(p$alpha^2 - p$beta^2),
  variance = function(p) p$delta * p$alpha^2 / (p$alpha^2 - p$beta^2)^1.5,
  references = "Barndorff-Nielsen (1997), Normal inverse Gaussian distributions and stochastic volatility modelling. Scand. J. Stat., 24, 1-13.",
  reference = ref_impl("GeneralizedHyperbolic", function(y, p) log(mapply(function(x, m, d, a, b) GeneralizedHyperbolic::dnig(x, m, d, a, b), y, p$mu, p$delta, p$alpha, p$beta))),
  stan = sfun("nig", c("mu", "alpha", "beta", "delta"),
    "    real s = sqrt(square(delta) + square(y - mu));\n    return log(alpha * delta) + log_modified_bessel_second_kind_frac(1, alpha * s) - log(pi() * s) + delta * sqrt(square(alpha) - square(beta)) + beta * (y - mu);"),
  status = "validated")

rfam("generalized_hyperbolic", "Generalized hyperbolic", list(P("lambda", "real", link = "identity"), shp("alpha"),
  P("beta", "real", link = "identity"), scale_par("delta"), loc_par("mu")), aliases = "GH",
  check = function(p) if (any(abs(p$beta) >= p$alpha, na.rm = TRUE)) "Need |beta| < alpha.",
  lpdf = function(y, p) { g <- sqrt(p$alpha^2 - p$beta^2); s <- sqrt(p$delta^2 + (y - p$mu)^2)
    p$lambda * log(g / p$delta) - 0.5 * log(2 * pi) - log_besselK(p$delta * g, p$lambda) +
      (p$lambda - 0.5) * (log(s) - log(p$alpha)) + log_besselK(p$alpha * s, p$lambda - 0.5) + p$beta * (y - p$mu) },
  references = "Barndorff-Nielsen (1977), Exponentially decreasing distributions for the logarithm of particle size. Proc. Roy. Soc. A, 353, 401-419.",
  reference = ref_impl("GeneralizedHyperbolic", function(y, p) log(GeneralizedHyperbolic::dghyp(y, p$mu[1], p$delta[1], p$alpha[1], p$beta[1], p$lambda[1]))),
  status = "validated")

rfam("hyperbolic", "Hyperbolic", list(shp("alpha"), P("beta", "real", link = "identity"), scale_par("delta"), loc_par("mu")),
  check = function(p) if (any(abs(p$beta) >= p$alpha, na.rm = TRUE)) "Need |beta| < alpha.",
  lpdf = function(y, p) { g <- sqrt(p$alpha^2 - p$beta^2)
    log(g) - log(2 * p$alpha * p$delta) - log_besselK(p$delta * g, 1) - p$alpha * sqrt(p$delta^2 + (y - p$mu)^2) + p$beta * (y - p$mu) },
  references = "Barndorff-Nielsen (1977).",
  reference = ref_impl("GeneralizedHyperbolic", function(y, p) log(GeneralizedHyperbolic::dhyperb(y, p$mu[1], p$delta[1], p$alpha[1], p$beta[1]))),
  status = "validated")

pIV_logk <- function(p) vapply(seq_along(p$m), function(i) { pp <- lapply(p, `[`, i)
  -log(stats::integrate(function(z) exp(-pp$m * log1p(z^2) - pp$nu * atan(z)), -Inf, Inf, rel.tol = 1e-12)$value) }, 0)
rfam("pearson_iv", "Pearson type IV", list(P("m", "positive", role = "shape", lower = 0.5), P("nu", "real", link = "identity"), scale_par("a"), loc_par("lambda")),
  formula = "f(y) = k (1 + z^2)^(-m) exp(-nu atan(z)) / a, z = (y - lambda)/a",
  lpdf = function(y, p) { z <- (y - p$lambda) / p$a; pIV_logk(p) - p$m * log1p(z^2) - p$nu * atan(z) - log(p$a) },
  numerical_notes = "Normalizing constant by quadrature.",
  references = "Pearson (1895); Heinrich (2004), A guide to the Pearson type IV distribution. CDF/MEMO/STATISTICS/PUBLIC/6820.",
  reference = ref_impl("PearsonDS", function(y, p) PearsonDS::dpearsonIV(y, p$m[1], p$nu[1], p$lambda[1], p$a[1], log = TRUE)),
  status = "validated")

# Stable and related laws ------------------------------------------------------------
# Stable density in Nolan's S1 parameterization via the Zolotarev integral
# (Nolan 1997, Theorem 1): standardized Z with f(z) = f(-z; alpha, -beta).
dstable1 <- function(z, a, b) {
  if (abs(a - 1) < 1e-8) {
    if (abs(b) < 1e-10) return(stats::dcauchy(z))
    if (b < 0) return(dstable1(-z, a, -b))
    x <- z
    V <- function(t) 2 / pi * ((pi / 2 + b * t) / cos(t)) * exp((pi / 2 + b * t) * tan(t) / b)
    e <- exp(-pi * x / (2 * b))
    return(e / (2 * b) * stats::integrate(function(t) { v <- V(t); ifelse(is.finite(v), v * exp(-e * v), 0) },
                                          -pi / 2, pi / 2, subdivisions = 2000L, rel.tol = 1e-10, stop.on.error = FALSE)$value)
  }
  if (z < 0) return(dstable1(-z, a, -b))
  zeta <- -b * tan(pi * a / 2); th0 <- atan(b * tan(pi * a / 2)) / a
  if (z < 1e-10) return(gamma(1 + 1 / a) * cos(th0) / (pi * (1 + zeta^2)^(1 / (2 * a))))
  V <- function(t) (cos(a * th0))^(1 / (a - 1)) * (cos(t) / sin(a * (th0 + t)))^(a / (a - 1)) * cos(a * th0 + (a - 1) * t) / cos(t)
  k <- z^(a / (a - 1))
  a * z^(1 / (a - 1)) / (pi * abs(a - 1)) * stats::integrate(function(t) { v <- V(t); ifelse(is.finite(v) & v > 0, v * exp(-k * v), 0) },
    -th0, pi / 2, subdivisions = 2000L, rel.tol = 1e-10, stop.on.error = FALSE)$value
}
# S1 with alpha = 1 adds beta (2/pi) gamma log(gamma) to the location (Nolan 2020, def. 1.7).
s1_shift <- function(a, b, g) if (abs(a - 1) < 1e-8) b * 2 / pi * g * log(g) else 0
dstable_s1 <- function(x, a, b, g, d) dstable1((x - d - s1_shift(a, b, g)) / g, a, b) / g
# Distribution function by the same integral representation (Nolan 1997, Theorem 1).
pstable1 <- function(z, a, b) {
  if (abs(a - 1) < 1e-8) {
    if (abs(b) < 1e-10) return(stats::pcauchy(z))
    if (b < 0) return(1 - pstable1(-z, a, -b))
    V <- function(t) 2 / pi * ((pi / 2 + b * t) / cos(t)) * exp((pi / 2 + b * t) * tan(t) / b)
    e <- exp(-pi * z / (2 * b))
    return(stats::integrate(function(t) { v <- V(t); ifelse(is.finite(v), exp(-e * v), 0) }, -pi / 2, pi / 2,
                            subdivisions = 2000L, rel.tol = 1e-10, stop.on.error = FALSE)$value / pi)
  }
  if (z < 0) return(1 - pstable1(-z, a, -b))
  th0 <- atan(b * tan(pi * a / 2)) / a
  if (z < 1e-12) return((pi / 2 - th0) / pi)
  V <- function(t) (cos(a * th0))^(1 / (a - 1)) * (cos(t) / sin(a * (th0 + t)))^(a / (a - 1)) * cos(a * th0 + (a - 1) * t) / cos(t)
  k <- z^(a / (a - 1))
  I <- stats::integrate(function(t) { v <- V(t); ifelse(is.finite(v) & v > 0, exp(-k * v), 0) }, -th0, pi / 2,
                        subdivisions = 2000L, rel.tol = 1e-10, stop.on.error = FALSE)$value / pi
  if (a < 1) (pi / 2 - th0) / pi + I else 1 - I
}
pstable_s1 <- function(x, a, b, g, d) pstable1((x - d - s1_shift(a, b, g)) / g, a, b)
rfam("stable", "Stable (alpha stable, S1 parameterization)", list(P("alpha", "unit", role = "index", interpretation = "Use (0, 2]; scaled to (0, 2) via 2 * alpha"),
  P("beta", "symmetric_unit", role = "skewness"), scale_par("gamma"), loc_par("delta")), aliases = c("levy_stable", "alpha_stable"),
  formula = "log phi(t) = -(gamma |t|)^A (1 - i beta sign(t) tan(pi A / 2)) + i delta t, A = 2 alpha",
  lpdf = function(y, p) log(pervec(y, p, function(x, pp) dstable_s1(x, 2 * pp$alpha, pp$beta, pp$gamma, pp$delta))),
  cdf = function(q, p) pervec(q, p, function(x, pp) pstable_s1(x, 2 * pp$alpha, pp$beta, pp$gamma, pp$delta)),
  numeric_moments = FALSE,
  numerical_notes = "Index A = 2 alpha so that the unit link covers (0, 2). Density by the Zolotarev integral (Nolan 1997).",
  references = c("Nolan (2020), Univariate Stable Distributions. Springer.", "Zolotarev (1986), One-dimensional Stable Distributions. AMS."),
  reference = ref_impl("stabledist", function(y, p) log(mapply(function(x, a, b, g, d) stabledist::dstable(x, 2 * a, b, g, d, pm = 1), y, p$alpha, p$beta, p$gamma, p$delta))),
  status = "validated")

rfam("landau", "Landau", list(loc_par("mu"), scale_par("c")),
  formula = "f(x) = (1/pi) int_0^inf exp(-t log t - z t) sin(pi t) dt / c, z = (x - mu)/c",
  lpdf = function(y, p) log(pervec(y, p, function(x, pp) dstable_s1(x, 1, 1, pi / 2 * pp$c, pp$mu))),
  cdf = function(q, p) pervec(q, p, function(x, pp) pstable_s1(x, 1, 1, pi / 2 * pp$c, pp$mu)),
  numeric_moments = FALSE,
  references = "Landau (1944), On the energy loss of fast particles by ionization. J. Phys. USSR, 8, 201-205.",
  reference = ref_impl("Landau integral", function(y, p) log(vapply((y - p$mu) / p$c, function(z) stats::integrate(function(t)
    exp(-t * log(t) - z * t) * sin(pi * t) / pi, 0, Inf, subdivisions = 2000L, rel.tol = 1e-10, stop.on.error = FALSE)$value, 0)) - log(p$c),
    note = "Standard Landau equals the S1 stable law with alpha = 1, beta = 1, gamma = pi/2 (Nolan 2020)."),
  status = "validated")

rfam("holtsmark", "Holtsmark (symmetric stable with index 3/2)", list(loc_par("mu"), scale_par("c")),
  lpdf = function(y, p) log(pervec(y, p, function(x, pp) dstable_s1(x, 1.5, 0, pp$c, pp$mu))),
  cdf = function(q, p) pervec(q, p, function(x, pp) pstable_s1(x, 1.5, 0, pp$c, pp$mu)),
  numeric_moments = FALSE, references = "Holtsmark (1919), Ann. Phys., 363, 577-630.",
  reference = ref_impl("stabledist", function(y, p) stabledist::dstable(y, 1.5, 0, p$c, p$mu, pm = 1, log = TRUE)), status = "validated")

# Geometric stable law as Y = mu E + E^(1/A) S, E ~ Exp(1), S ~ S1(A, beta, sigma, 0)
# (Kozubowski and Rachev 1999, A != 1).
gstable_mix <- function(x, A, b, s, mu) stats::integrate(function(e) vapply(e, function(ee)
  dstable_s1(x - mu * ee, A, b, s * ee^(1 / A), 0), 0) * exp(-e), 0, Inf, rel.tol = 1e-8, stop.on.error = FALSE)$value

rfam("linnik", "Linnik (symmetric geometric stable)", list(P("alpha", "unit", role = "index", interpretation = "Index A = 2 alpha in (0, 2)"), scale_par("sigma")),
  formula = "phi(t) = 1 / (1 + |sigma t|^A)",
  lpdf = function(y, p) log(pervec(y, p, function(x, pp) gstable_mix(x, 2 * pp$alpha, 0, pp$sigma, 0))),
  numeric_moments = FALSE, references = "Linnik (1963); Kotz, Kozubowski and Podgorski (2001), The Laplace Distribution and Generalizations.",
  reference = ref_impl("characteristic function inversion", function(y, p) log(pervec(y, p, function(x, pp)
    cf_density(x, function(t) 1 / (1 + (pp$sigma * t)^(2 * pp$alpha)))))),
  status = "validated")

rfam("geometric_stable", "Geometric stable", list(P("alpha", "unit", role = "index"), P("beta", "symmetric_unit", role = "skewness"), scale_par("sigma"), loc_par("mu")),
  formula = "phi(t) = 1 / (1 + (sigma |t|)^A (1 - i beta sign(t) tan(pi A/2)) - i mu t), A = 2 alpha",
  lpdf = function(y, p) log(pervec(y, p, function(x, pp) gstable_mix(x, 2 * pp$alpha, pp$beta, pp$sigma, pp$mu))),
  numeric_moments = FALSE, references = "Kozubowski and Rachev (1999), Univariate geometric stable laws. J. Comput. Anal. Appl., 1, 177-217.",
  reference = ref_impl("characteristic function inversion", function(y, p) log(pervec(y, p, function(x, pp) { A <- 2 * pp$alpha
    cf_density(x, function(t) 1 / (1 + (pp$sigma * t)^A * (1 - 1i * pp$beta * tan(pi * A / 2)) - 1i * pp$mu * t)) }))),
  status = "validated")

rfam("voigt", "Voigt profile", list(loc_par("mu"), scale_par("sigma"), scale_par("gamma")),
  formula = "Convolution of Normal(mu, sigma) and Cauchy(0, gamma)",
  lpdf = function(y, p) log(pervec(y, p, function(x, pp) cf_density(x - pp$mu, function(t) exp(-pp$sigma^2 * t^2 / 2 - pp$gamma * t)))),
  rng = function(n, p) stats::rnorm(n, p$mu, p$sigma) + stats::rcauchy(n, 0, p$gamma), numeric_moments = FALSE,
  references = "Olivero and Longbothum (1977), J. Quant. Spectrosc. Radiat. Transfer, 17, 233-236.",
  reference = ref_impl("convolution", function(y, p) log(pervec(y, p, function(x, pp) stats::integrate(function(t)
    stats::dnorm(t, pp$mu, pp$sigma) * stats::dcauchy(x - t, 0, pp$gamma), -Inf, Inf, rel.tol = 1e-11)$value))), status = "validated")

rfam("gamma_difference", "Gamma difference", list(shp("a1"), P("b1", "positive", role = "rate"), shp("a2"), P("b2", "positive", role = "rate")),
  formula = "Y = X1 - X2, Xj ~ Gamma(aj, rate bj)",
  lpdf = function(y, p) log(pervec(y, p, function(x, pp) cf_density(x, function(t) (1 - 1i * t / pp$b1)^(-pp$a1) * (1 + 1i * t / pp$b2)^(-pp$a2)))),
  rng = function(n, p) stats::rgamma(n, p$a1, p$b1) - stats::rgamma(n, p$a2, p$b2),
  mean = function(p) p$a1 / p$b1 - p$a2 / p$b2, variance = function(p) p$a1 / p$b1^2 + p$a2 / p$b2^2,
  references = "Klar (2015), A note on gamma difference distributions. J. Stat. Comput. Simul., 85, 3708-3715.",
  reference = ref_impl("convolution", function(y, p) log(pervec(y, p, function(x, pp) stats::integrate(function(t)
    stats::dgamma(t, pp$a1, pp$b1) * stats::dgamma(t - x, pp$a2, pp$b2), max(x, 0), Inf, rel.tol = 1e-11)$value))), status = "validated")

rfam("fisher_z", "Fisher z", list(P("d1", "positive"), P("d2", "positive")),
  formula = "z = log(F)/2, F ~ F(d1, d2)",
  lpdf = function(y, p) log(2) + p$d1 / 2 * log(p$d1) + p$d2 / 2 * log(p$d2) + p$d1 * y - lbeta(p$d1 / 2, p$d2 / 2) - (p$d1 + p$d2) / 2 * log(p$d1 * exp(2 * y) + p$d2),
  cdf = function(q, p) stats::pf(exp(2 * q), p$d1, p$d2), quantile = function(pr, p) log(stats::qf(pr, p$d1, p$d2)) / 2,
  references = "Fisher (1924), Proc. Int. Cong. Math., Toronto, 2, 805-813.",
  reference = ref_impl("F transform", function(y, p) stats::df(exp(2 * y), p$d1, p$d2, log = TRUE) + log(2) + 2 * y), status = "validated")

rfam("behrens_fisher", "Behrens-Fisher", list(P("nu1", "positive"), P("nu2", "positive"), P("theta", "unit", role = "angle", interpretation = "Angle in (0, pi/2) scaled to (0, 1)")),
  formula = "T = T1 cos(phi) - T2 sin(phi), Tj ~ t(nu_j), phi = pi theta / 2",
  lpdf = function(y, p) log(pervec(y, p, function(x, pp) { ph <- pi * pp$theta / 2; c1 <- cos(ph); s1 <- sin(ph)
    stats::integrate(function(t) stats::dt(t, pp$nu1) * stats::dt((t * c1 - x) / s1, pp$nu2) / s1, -Inf, Inf, rel.tol = 1e-10)$value })),
  rng = function(n, p) stats::rt(n, p$nu1) * cos(pi * p$theta / 2) - stats::rt(n, p$nu2) * sin(pi * p$theta / 2),
  references = "Fisher (1935); Sukhatme (1938), Sankhya, 4, 39-48.",
  reference = ref_impl("simulation free check: theta -> 0 gives t(nu1)", function(y, p) log(pervec(y, p, function(x, pp) { ph <- pi * pp$theta / 2
    stats::integrate(function(t) stats::dt(t, pp$nu2) * stats::dt((x + t * sin(ph)) / cos(ph), pp$nu1) / cos(ph), -Inf, Inf, rel.tol = 1e-10)$value }))),
  status = "validated")

champ_lognorm <- function(p) vapply(seq_along(p$alpha), function(i) { pp <- lapply(p, `[`, i)
  log(stats::integrate(function(z) 1 / (cosh(pp$alpha * z) + pp$lambda), -Inf, Inf, rel.tol = 1e-12)$value) }, 0)
rfam("champernowne", "Champernowne", list(shp("alpha"), P("lambda", "real", role = "shape", link = "identity", lower = -1), loc_par("y0")),
  formula = "f(y) = n / (cosh(alpha (y - y0)) + lambda)",
  lpdf = function(y, p) -log(cosh(p$alpha * (y - p$y0)) + p$lambda) - champ_lognorm(p),
  references = "Champernowne (1953), A model of income distribution. Economic J., 63, 318-351.",
  reference = ref_impl("closed form for lambda = 0", function(y, p) log(p$alpha / pi / cosh(p$alpha * (y - p$y0)))), status = "validated")

rfam("tracy_widom", "Tracy-Widom", list(P("beta", "positive_integer", role = "symmetry class (1, 2 or 4)", data = TRUE)),
  lpdf = function(y, p) { if (!requireNamespace("RMTstat", quietly = TRUE)) stop_unsupported("Tracy-Widom without RMTstat", "tracy_widom")
    RMTstat::dtw(y, p$beta[1], log = TRUE) },
  cdf = function(q, p) RMTstat::ptw(q, p$beta[1]), numeric_moments = FALSE,
  numerical_notes = "Uses the tabulated Painleve II solution in RMTstat (Suggests).",
  references = "Tracy and Widom (1994), Level-spacing distributions and the Airy kernel. Comm. Math. Phys., 159, 151-174.",
  reference = ref_impl("RMTstat", function(y, p) RMTstat::dtw(y, p$beta[1], log = TRUE)), status = "implemented")

# Quantile defined laws --------------------------------------------------------------
tl_Q <- function(u, p) ifelse(abs(p$lambda) < 1e-9, log(u / (1 - u)), (u^p$lambda - (1 - u)^p$lambda) / p$lambda)
tl_dQ <- function(u, p) u^(p$lambda - 1) + (1 - u)^(p$lambda - 1)
rfam("tukey_lambda", "Tukey lambda", list(P("lambda", "real", role = "shape", link = "identity")),
  support = function(p) list(lower = ifelse(p$lambda > 0, -1 / p$lambda, -Inf), upper = ifelse(p$lambda > 0, 1 / p$lambda, Inf)),
  formula = "Q(u) = (u^lambda - (1 - u)^lambda) / lambda",
  lpdf = qdef_lpdf(tl_Q, tl_dQ), cdf = qdef_cdf(tl_Q), quantile = tl_Q,
  references = "Tukey (1962); Joiner and Rosenblatt (1971), JASA, 66, 394-399.",
  reference = ref_impl("gld (RS with lambda3 = lambda4)", function(y, p) log(gld::dgl(y, c(0, p$lambda[1], p$lambda[1], p$lambda[1]), param = "rs"))),
  status = "validated")

gld_Q <- function(u, p) p$l1 + (ifelse(abs(p$l3) < 1e-9, log(u), (u^p$l3 - 1) / p$l3) - ifelse(abs(p$l4) < 1e-9, log1p(-u), ((1 - u)^p$l4 - 1) / p$l4)) / p$l2
gld_dQ <- function(u, p) (u^(p$l3 - 1) + (1 - u)^(p$l4 - 1)) / p$l2
rfam("generalized_lambda", "Generalized lambda (FKML)", list(loc_par("l1"), scale_par("l2"), P("l3", "real", link = "identity"), P("l4", "real", link = "identity")),
  aliases = "gld", formula = "Q(u) = l1 + ((u^l3 - 1)/l3 - ((1 - u)^l4 - 1)/l4) / l2",
  support = function(p) list(lower = gld_Q(0, p), upper = gld_Q(1, p)),
  lpdf = qdef_lpdf(gld_Q, gld_dQ), cdf = qdef_cdf(gld_Q), quantile = gld_Q,
  references = "Freimer, Kollia, Mudholkar and Lin (1988), A study of the generalized Tukey lambda family. Comm. Stat. Theory Methods, 17, 3547-3567.",
  reference = ref_impl("gld", function(y, p) log(gld::dgl(y, c(p$l1[1], p$l2[1], p$l3[1], p$l4[1]), param = "fkml"))), status = "validated")

wk_Q <- function(u, p) p$xi + p$alpha / p$beta * (1 - (1 - u)^p$beta) - p$gamma / p$delta * (1 - (1 - u)^(-p$delta))
wk_dQ <- function(u, p) p$alpha * (1 - u)^(p$beta - 1) + p$gamma * (1 - u)^(-p$delta - 1)
rfam("wakeby", "Wakeby", list(loc_par("xi"), P("alpha", "nonneg"), P("beta", "real", link = "identity"), P("gamma", "nonneg"), P("delta", "real", link = "identity")),
  support = function(p) list(lower = p$xi, upper = ifelse(p$delta < 0 | p$gamma == 0, p$xi + p$alpha / p$beta - p$gamma / p$delta, Inf)),
  formula = "Q(u) = xi + alpha/beta (1 - (1 - u)^beta) - gamma/delta (1 - (1 - u)^(-delta))",
  lpdf = qdef_lpdf(wk_Q, wk_dQ), cdf = qdef_cdf(wk_Q), quantile = wk_Q,
  references = "Houghton (1978), Birth of a parent: The Wakeby distribution for modeling flood flows. Water Resources Research, 14, 1105-1109.",
  reference = ref_impl("lmomco", function(y, p) log(lmomco::pdfwak(y, lmomco::vec2par(c(p$xi[1], p$alpha[1], p$beta[1], p$gamma[1], p$delta[1]), type = "wak")))),
  status = "validated")
