# Continuous distributions on the positive half line.

pos <- "positive"
pos_support <- function(p) list(lower = 0, upper = Inf)
posfam <- function(id, display, pars, ..., support = pos_support)
  register_family(family_def(id, display, pos, pars, support = support, ...))
# Exponential family core ---------------------------------------------------
posfam("exponential", "Exponential", list(P("lambda", "positive", role = "rate", offset_allowed = TRUE)),
  aliases = c("exp"), formula = "f(y) = lambda exp(-lambda y)",
  mean_formula = "1 / lambda", variance_formula = "1 / lambda^2",
  lpdf = function(y, p) stats::dexp(y, p$lambda, log = TRUE), cdf = function(q, p) stats::pexp(q, p$lambda),
  quantile = function(pr, p) stats::qexp(pr, p$lambda), rng = function(n, p) stats::rexp(n, p$lambda),
  mean = function(p) 1 / p$lambda, variance = function(p) 1 / p$lambda^2,
  references = JKB, reference = ref_impl("gamlss.dist", function(y, p) gamlss.dist::dEXP(y, 1 / p$lambda, log = TRUE)),
  stan = stan_native("exponential_lpdf({y} | {lambda})", "exponential_lcdf({y} | {lambda})",
                     "exponential_lccdf({y} | {lambda})", "exponential_rng({lambda})"),
  bugs = list(native = "dexp({lambda})", loglik = "log({lambda}) - {lambda} * {y}"),
  start = function(y) list(lambda = 1 / mean(y)), status = "validated")

posfam("exponential_mean", "Exponential (mean parameterization)", list(mean_par()), aliases = "EXP",
  formula = "f(y) = exp(-y / mu) / mu", mean_formula = "mu", variance_formula = "mu^2",
  lpdf = function(y, p) stats::dexp(y, 1 / p$mu, log = TRUE), cdf = function(q, p) stats::pexp(q, 1 / p$mu),
  quantile = function(pr, p) stats::qexp(pr, 1 / p$mu), rng = function(n, p) stats::rexp(n, 1 / p$mu),
  mean = function(p) p$mu, variance = function(p) p$mu^2,
  references = JKB, reference = ref_impl("gamlss.dist", function(y, p) gamlss.dist::dEXP(y, p$mu, log = TRUE)),
  stan = stan_native("exponential_lpdf({y} | inv({mu}))", "exponential_lcdf({y} | inv({mu}))",
                     "exponential_lccdf({y} | inv({mu}))", "exponential_rng(inv({mu}))"),
  bugs = list(native = "dexp(1 / {mu})", loglik = "-log({mu}) - {y} / {mu}"),
  start = function(y) list(mu = mean(y)), status = "validated")

posfam("gamma", "Gamma (mean and coefficient of variation, GAMLSS GA)", list(mean_par(), P("sigma", "positive", role = "dispersion",
  interpretation = "Coefficient of variation")), aliases = c("GA"),
  formula = "shape = 1 / sigma^2, scale = mu sigma^2", mean_formula = "mu", variance_formula = "sigma^2 mu^2",
  lpdf = function(y, p) stats::dgamma(y, 1 / p$sigma^2, scale = p$mu * p$sigma^2, log = TRUE),
  cdf = function(q, p) stats::pgamma(q, 1 / p$sigma^2, scale = p$mu * p$sigma^2),
  quantile = function(pr, p) stats::qgamma(pr, 1 / p$sigma^2, scale = p$mu * p$sigma^2),
  rng = function(n, p) stats::rgamma(n, 1 / p$sigma^2, scale = p$mu * p$sigma^2),
  mean = function(p) p$mu, variance = function(p) p$sigma^2 * p$mu^2,
  references = c(JKB, GS), reference = ref_impl("gamlss.dist", gd("dGA"), gp("pGA")),
  regression_notes = "The GLM gamma model with log link; the GLM dispersion equals sigma^2.",
  stan = stan_native("gamma_lpdf({y} | inv_square({sigma}), inv_square({sigma}) / {mu})",
                     "gamma_lcdf({y} | inv_square({sigma}), inv_square({sigma}) / {mu})",
                     "gamma_lccdf({y} | inv_square({sigma}), inv_square({sigma}) / {mu})",
                     "gamma_rng(inv_square({sigma}), inv_square({sigma}) / {mu})"),
  bugs = list(native = "dgamma(1 / pow({sigma}, 2), 1 / (pow({sigma}, 2) * {mu}))"),
  start = function(y) list(mu = mean(y), sigma = stats::sd(y) / mean(y)), status = "validated")

posfam("gamma_rate", "Gamma (shape and rate)", list(shp("shape"), P("rate", "positive", role = "rate", offset_allowed = TRUE)),
  aliases = c("erlang_general"), formula = "f(y) = rate^shape y^(shape - 1) exp(-rate y) / Gamma(shape)",
  mean_formula = "shape / rate", variance_formula = "shape / rate^2",
  lpdf = function(y, p) stats::dgamma(y, p$shape, p$rate, log = TRUE), cdf = function(q, p) stats::pgamma(q, p$shape, p$rate),
  quantile = function(pr, p) stats::qgamma(pr, p$shape, p$rate), rng = function(n, p) stats::rgamma(n, p$shape, p$rate),
  mean = function(p) p$shape / p$rate, variance = function(p) p$shape / p$rate^2,
  references = JKB, reference = ref_impl("gamlss.dist", function(y, p) gamlss.dist::dGA(y, p$shape / p$rate, 1 / sqrt(p$shape), log = TRUE)),
  stan = stan_native("gamma_lpdf({y} | {shape}, {rate})", "gamma_lcdf({y} | {shape}, {rate})",
                     "gamma_lccdf({y} | {shape}, {rate})", "gamma_rng({shape}, {rate})"),
  bugs = list(native = "dgamma({shape}, {rate})"), status = "validated")

posfam("erlang", "Erlang", list(P("k", "positive_integer", role = "shape", data = TRUE), P("rate", "positive", role = "rate")),
  formula = "Gamma with integer shape k", mean_formula = "k / rate",
  lpdf = function(y, p) stats::dgamma(y, p$k, p$rate, log = TRUE), cdf = function(q, p) stats::pgamma(q, p$k, p$rate),
  quantile = function(pr, p) stats::qgamma(pr, p$k, p$rate), rng = function(n, p) stats::rgamma(n, p$k, p$rate),
  mean = function(p) p$k / p$rate, variance = function(p) p$k / p$rate^2,
  references = JKB, reference = ref_impl("definition", function(y, p) p$k * log(p$rate) + (p$k - 1) * log(y) - p$rate * y - lfactorial(p$k - 1)),
  stan = stan_native("gamma_lpdf({y} | {k}, {rate})", "gamma_lcdf({y} | {k}, {rate})", "gamma_lccdf({y} | {k}, {rate})"),
  status = "validated")

posfam("inverse_gamma", "Inverse gamma", list(shp("shape"), P("scale", "positive", role = "scale")),
  formula = "f(y) = scale^shape y^(-shape - 1) exp(-scale / y) / Gamma(shape)",
  mean_formula = "scale / (shape - 1) for shape > 1",
  lpdf = function(y, p) p$shape * log(p$scale) - lgamma(p$shape) - (p$shape + 1) * log(y) - p$scale / y,
  cdf = function(q, p) stats::pgamma(p$scale / q, p$shape, lower.tail = FALSE),
  quantile = function(pr, p) p$scale / stats::qgamma(pr, p$shape, lower.tail = FALSE),
  mean = function(p) ifelse(p$shape > 1, p$scale / (p$shape - 1), NA),
  variance = function(p) ifelse(p$shape > 2, p$scale^2 / ((p$shape - 1)^2 * (p$shape - 2)), NA),
  references = JKB, reference = ref_impl("actuar", function(y, p) actuar::dinvgamma(y, p$shape, scale = p$scale, log = TRUE)),
  stan = stan_native("inv_gamma_lpdf({y} | {shape}, {scale})", "inv_gamma_lcdf({y} | {shape}, {scale})",
                     "inv_gamma_lccdf({y} | {shape}, {scale})", "inv_gamma_rng({shape}, {scale})"),
  status = "validated")

posfam("inverse_gaussian", "Inverse Gaussian (Wald)", list(mean_par(), P("sigma", "positive", role = "dispersion")),
  aliases = c("IG", "wald"), formula = "f(y) = exp(-(y - mu)^2 / (2 mu^2 sigma^2 y)) / sqrt(2 pi sigma^2 y^3)",
  mean_formula = "mu", variance_formula = "sigma^2 mu^3",
  lpdf = function(y, p) -0.5 * log(2 * pi * p$sigma^2 * y^3) - (y - p$mu)^2 / (2 * p$mu^2 * p$sigma^2 * y),
  cdf = function(q, p) { l <- 1 / p$sigma^2; a <- sqrt(l / q)
    stats::pnorm(a * (q / p$mu - 1)) + exp(2 * l / p$mu + stats::pnorm(-a * (q / p$mu + 1), log.p = TRUE)) },
  rng = function(n, p) { nu <- stats::rnorm(n)^2; l <- 1 / p$sigma^2
    x <- p$mu + p$mu^2 * nu / (2 * l) - p$mu / (2 * l) * sqrt(4 * p$mu * l * nu + p$mu^2 * nu^2)
    ifelse(stats::runif(n) <= p$mu / (p$mu + x), x, p$mu^2 / x) },
  mean = function(p) p$mu, variance = function(p) p$sigma^2 * p$mu^3,
  references = c(JKB, GS), reference = ref_impl("gamlss.dist", gd("dIG"), gp("pIG")),
  stan = sfun("inverse_gaussian", c("mu", "sigma"),
    "    return -0.5 * log(2 * pi() * square(sigma) * pow(y, 3)) - square(y - mu) / (2 * square(mu) * square(sigma) * y);"),
  bugs = list(loglik = "-0.5 * log(2 * 3.141592653589793 * pow({sigma}, 2) * pow({y}, 3)) - pow({y} - {mu}, 2) / (2 * pow({mu}, 2) * pow({sigma}, 2) * {y})"),
  start = function(y) list(mu = mean(y), sigma = sqrt(stats::var(y) / mean(y)^3)), status = "validated")

posfam("lognormal", "Lognormal", list(loc_par("mu"), scale_par("sigma")), aliases = c("LOGNO", "lnorm"),
  formula = "log Y ~ Normal(mu, sigma)", mean_formula = "exp(mu + sigma^2 / 2)",
  variance_formula = "(exp(sigma^2) - 1) exp(2 mu + sigma^2)",
  lpdf = function(y, p) stats::dlnorm(y, p$mu, p$sigma, log = TRUE), cdf = function(q, p) stats::plnorm(q, p$mu, p$sigma),
  quantile = function(pr, p) stats::qlnorm(pr, p$mu, p$sigma), rng = function(n, p) stats::rlnorm(n, p$mu, p$sigma),
  mean = function(p) exp(p$mu + p$sigma^2 / 2), variance = function(p) (exp(p$sigma^2) - 1) * exp(2 * p$mu + p$sigma^2),
  references = c(JKB, GS), reference = ref_impl("gamlss.dist", gd("dLOGNO"), gp("pLOGNO")),
  stan = stan_native("lognormal_lpdf({y} | {mu}, {sigma})", "lognormal_lcdf({y} | {mu}, {sigma})",
                     "lognormal_lccdf({y} | {mu}, {sigma})", "lognormal_rng({mu}, {sigma})"),
  bugs = list(native = "dlnorm({mu}, 1 / pow({sigma}, 2))"),
  start = function(y) list(mu = mean(log(y)), sigma = stats::sd(log(y))), status = "validated")

posfam("weibull", "Weibull", list(shp("shape"), P("scale", "positive", role = "scale", offset_allowed = TRUE)),
  aliases = c("WEI"), formula = "F(y) = 1 - exp(-(y / scale)^shape)",
  mean_formula = "scale Gamma(1 + 1/shape)", variance_formula = "scale^2 (Gamma(1 + 2/shape) - Gamma(1 + 1/shape)^2)",
  lpdf = function(y, p) stats::dweibull(y, p$shape, p$scale, log = TRUE), cdf = function(q, p) stats::pweibull(q, p$shape, p$scale),
  lccdf = function(q, p) -(q / p$scale)^p$shape,
  lcdf = function(q, p) stats::pweibull(q, p$shape, p$scale, log.p = TRUE),
  quantile = function(pr, p) stats::qweibull(pr, p$shape, p$scale), rng = function(n, p) stats::rweibull(n, p$shape, p$scale),
  mean = function(p) p$scale * gamma(1 + 1 / p$shape),
  variance = function(p) p$scale^2 * (gamma(1 + 2 / p$shape) - gamma(1 + 1 / p$shape)^2),
  references = c(JKB, GS), reference = ref_impl("gamlss.dist", function(y, p) gamlss.dist::dWEI(y, p$scale, p$shape, log = TRUE)),
  stan = stan_native("weibull_lpdf({y} | {shape}, {scale})", "weibull_lcdf({y} | {shape}, {scale})",
                     "weibull_lccdf({y} | {shape}, {scale})", "weibull_rng({shape}, {scale})"),
  bugs = list(native = "dweib({shape}, pow({scale}, -{shape}))"),
  start = function(y) list(shape = 1.2, scale = mean(y)), status = "validated")

posfam("weibull_mean", "Weibull (mean parameterization, GAMLSS WEI3)", list(mean_par(), shp("sigma")), aliases = "WEI3",
  formula = "Weibull with shape sigma and scale mu / Gamma(1 + 1/sigma)", mean_formula = "mu",
  lpdf = function(y, p) stats::dweibull(y, p$sigma, p$mu / gamma(1 + 1 / p$sigma), log = TRUE),
  cdf = function(q, p) stats::pweibull(q, p$sigma, p$mu / gamma(1 + 1 / p$sigma)),
  quantile = function(pr, p) stats::qweibull(pr, p$sigma, p$mu / gamma(1 + 1 / p$sigma)),
  rng = function(n, p) stats::rweibull(n, p$sigma, p$mu / gamma(1 + 1 / p$sigma)), mean = function(p) p$mu,
  references = GS, reference = ref_impl("gamlss.dist", gd("dWEI3"), gp("pWEI3")),
  stan = stan_native("weibull_lpdf({y} | {sigma}, {mu} / tgamma(1 + inv({sigma})))",
                     "weibull_lcdf({y} | {sigma}, {mu} / tgamma(1 + inv({sigma})))",
                     "weibull_lccdf({y} | {sigma}, {mu} / tgamma(1 + inv({sigma})))"),
  start = function(y) list(mu = mean(y), sigma = 1.2), status = "validated")

gg_lpdf <- function(y, p) { th <- 1 / (p$sigma^2 * p$nu^2); z <- (y / p$mu)^p$nu
  log(abs(p$nu)) + th * log(th) + th * log(z) - th * z - lgamma(th) - log(y) }
posfam("generalized_gamma", "Generalized gamma (GAMLSS GG, Lopatatzidis and Green)",
  list(P("mu", "positive", role = "location", offset_allowed = TRUE), P("sigma", "positive", role = "scale"),
       P("nu", "real", role = "shape", link = "identity")), aliases = c("GG", "stacy"),
  formula = "theta = 1/(sigma^2 nu^2), z = (y/mu)^nu, Z ~ Gamma(theta, rate theta)",
  lpdf = function(y, p) ifelse(abs(p$nu) < 1e-6, stats::dlnorm(y, log(p$mu), p$sigma, log = TRUE), gg_lpdf(y, p)),
  cdf = function(q, p) { th <- 1 / (p$sigma^2 * p$nu^2); z <- (q / p$mu)^p$nu
    ifelse(p$nu > 0, stats::pgamma(z * th, th), stats::pgamma(z * th, th, lower.tail = FALSE)) },
  quantile = function(pr, p) { th <- 1 / (p$sigma^2 * p$nu^2)
    p$mu * (stats::qgamma(ifelse(p$nu > 0, pr, 1 - pr), th) / th)^(1 / p$nu) },
  references = c("Stacy (1962), A generalization of the gamma distribution. Ann. Math. Stat., 33, 1187-1192.", GS),
  reference = ref_impl("gamlss.dist", gd("dGG"), gp("pGG")),
  stan = sfun("gengamma", c("mu", "sigma", "nu"),
    "    if (abs(nu) < 1e-6) return lognormal_lpdf(y | log(mu), sigma);\n    real th = 1 / square(sigma * nu);\n    real z = pow(y / mu, nu);\n    return log(abs(nu)) + th * log(th) + th * log(z) - th * z - lgamma(th) - log(y);"),
  start = function(y) list(mu = stats::median(y), sigma = 0.5, nu = 1), status = "validated")

posfam("generalized_gamma_stacy", "Generalized gamma (Stacy: scale, d, p)", list(P("a", "positive", role = "scale"), shp("d"), shp("p")),
  formula = "f(y) = (p / a^d) y^(d - 1) exp(-(y/a)^p) / Gamma(d / p)",
  lpdf = function(y, q) log(q$p) - q$d * log(q$a) + (q$d - 1) * log(y) - (y / q$a)^q$p - lgamma(q$d / q$p),
  cdf = function(x, q) stats::pgamma((x / q$a)^q$p, q$d / q$p),
  quantile = function(pr, q) q$a * stats::qgamma(pr, q$d / q$p)^(1 / q$p),
  mean = function(q) q$a * gamma((q$d + 1) / q$p) / gamma(q$d / q$p),
  references = "Stacy (1962).", reference = ref_impl("flexsurv", function(y, q) flexsurv::dgengamma.orig(y, shape = q$p, scale = q$a, k = q$d / q$p, log = TRUE)),
  stan = sfun("gengamma_stacy", c("a", "d", "p"),
    "    return log(p) - d * log(a) + (d - 1) * log(y) - pow(y / a, p) - lgamma(d / p);"),
  status = "validated")

# Log-location-scale and heavy tailed laws -----------------------------------
posfam("log_logistic", "Log-logistic (Fisk)", list(P("scale", "positive", role = "scale", offset_allowed = TRUE), shp("shape")),
  aliases = c("fisk", "llogis"), formula = "F(y) = 1 / (1 + (y / scale)^(-shape))",
  mean_formula = "scale (pi/shape) / sin(pi/shape) for shape > 1",
  lpdf = function(y, p) log(p$shape) - log(p$scale) + (p$shape - 1) * log(y / p$scale) - 2 * log1p((y / p$scale)^p$shape),
  cdf = function(q, p) 1 / (1 + (q / p$scale)^(-p$shape)),
  quantile = function(pr, p) p$scale * (pr / (1 - pr))^(1 / p$shape),
  mean = function(p) ifelse(p$shape > 1, p$scale * (pi / p$shape) / sin(pi / p$shape), NA),
  references = KLUG, reference = ref_impl("actuar", function(y, p) actuar::dllogis(y, p$shape, scale = p$scale, log = TRUE)),
  stan = stan_native("loglogistic_lpdf({y} | {scale}, {shape})", "loglogistic_cdf({y} | {scale}, {shape})"),
  status = "validated")

posfam("pareto", "Pareto (type I)", list(P("xm", "positive", role = "minimum", data = TRUE), shp("alpha")),
  formula = "f(y) = alpha xm^alpha / y^(alpha + 1), y >= xm",
  mean_formula = "alpha xm / (alpha - 1) for alpha > 1",
  lpdf = function(y, p) ifelse(y >= p$xm, log(p$alpha) + p$alpha * log(p$xm) - (p$alpha + 1) * log(y), -Inf),
  cdf = function(q, p) ifelse(q < p$xm, 0, 1 - (p$xm / q)^p$alpha),
  quantile = function(pr, p) p$xm * (1 - pr)^(-1 / p$alpha),
  mean = function(p) ifelse(p$alpha > 1, p$alpha * p$xm / (p$alpha - 1), Inf),
  support = function(p) list(lower = p$xm, upper = Inf),
  references = JKB, reference = ref_impl("extraDistr", function(y, p) extraDistr::dpareto(y, p$alpha, p$xm, log = TRUE)),
  stan = stan_native("pareto_lpdf({y} | {xm}, {alpha})", "pareto_lcdf({y} | {xm}, {alpha})", "pareto_lccdf({y} | {xm}, {alpha})"),
  status = "validated")

posfam("lomax", "Lomax (Pareto type II)", list(P("scale", "positive", role = "scale"), shp("alpha")),
  aliases = c("pareto_type2", "PARETO2"), formula = "f(y) = alpha / scale (1 + y / scale)^(-(alpha + 1))",
  lpdf = function(y, p) log(p$alpha) - log(p$scale) - (p$alpha + 1) * log1p(y / p$scale),
  cdf = function(q, p) 1 - (1 + q / p$scale)^(-p$alpha),
  quantile = function(pr, p) p$scale * ((1 - pr)^(-1 / p$alpha) - 1),
  mean = function(p) ifelse(p$alpha > 1, p$scale / (p$alpha - 1), Inf),
  references = KLUG, reference = ref_impl("extraDistr", function(y, p) extraDistr::dlomax(y, 1 / p$scale, p$alpha, log = TRUE)),
  stan = stan_native("pareto_type_2_lpdf({y} | 0, {scale}, {alpha})", "pareto_type_2_lcdf({y} | 0, {scale}, {alpha})",
                     "pareto_type_2_lccdf({y} | 0, {scale}, {alpha})"),
  status = "validated")

posfam("burr", "Burr type XII (Singh-Maddala)", list(P("scale", "positive", role = "scale"), shp("c"), shp("k")),
  aliases = c("burr12", "singh_maddala"), formula = "F(y) = 1 - (1 + (y/scale)^c)^(-k)",
  lpdf = function(y, p) log(p$c) + log(p$k) - log(p$scale) + (p$c - 1) * log(y / p$scale) - (p$k + 1) * log1p((y / p$scale)^p$c),
  cdf = function(q, p) 1 - (1 + (q / p$scale)^p$c)^(-p$k),
  quantile = function(pr, p) p$scale * ((1 - pr)^(-1 / p$k) - 1)^(1 / p$c),
  references = c("Burr (1942), Cumulative frequency functions. Ann. Math. Stat., 13, 215-232.", KLUG),
  reference = ref_impl("actuar", function(y, p) actuar::dburr(y, p$k, p$c, scale = p$scale, log = TRUE)),
  stan = sfun("burr", c("scale", "c", "k"),
    "    return log(c) + log(k) - log(scale) + (c - 1) * log(y / scale) - (k + 1) * log1p(pow(y / scale, c));",
    lccdf = "    return -k * log1p(pow(y / scale, c));"),
  status = "validated")

posfam("dagum", "Dagum (Burr type III, inverse Burr)", list(P("scale", "positive", role = "scale"), shp("a"), shp("p")),
  aliases = c("burr3", "inverse_burr"), formula = "F(y) = (1 + (y/scale)^(-a))^(-p)",
  lpdf = function(y, q) log(q$a) + log(q$p) - log(y) + q$a * q$p * log(y / q$scale) - (q$p + 1) * log1p((y / q$scale)^q$a),
  cdf = function(x, q) (1 + (x / q$scale)^(-q$a))^(-q$p),
  quantile = function(pr, q) q$scale * (pr^(-1 / q$p) - 1)^(-1 / q$a),
  references = c("Dagum (1977), A new model of personal income distribution. Economie Appliquee, 30, 413-437.", KLUG),
  reference = ref_impl("actuar", function(y, q) actuar::dinvburr(y, q$p, q$a, scale = q$scale, log = TRUE)),
  stan = sfun("dagum", c("scale", "a", "p"),
    "    return log(a) + log(p) - log(y) + a * p * log(y / scale) - (p + 1) * log1p(pow(y / scale, a));"),
  status = "validated")

posfam("generalized_beta2", "Generalized beta of the second kind (GB2)",
  list(P("scale", "positive", role = "scale"), P("a", "real", role = "shape", link = "identity"), shp("p"), shp("q")),
  aliases = c("GB2", "transformed_beta"),
  formula = "f(y) = |a| y^(a p - 1) / (scale^(a p) B(p, q) (1 + (y/scale)^a)^(p + q))",
  lpdf = function(y, r) log(abs(r$a)) + (r$a * r$p - 1) * log(y) - r$a * r$p * log(r$scale) - lbeta(r$p, r$q) - (r$p + r$q) * log1p((y / r$scale)^r$a),
  cdf = function(x, r) { u <- (x / r$scale)^r$a; v <- stats::pbeta(u / (1 + u), r$p, r$q); ifelse(r$a > 0, v, 1 - v) },
  quantile = function(pr, r) { b <- stats::qbeta(ifelse(r$a > 0, pr, 1 - pr), r$p, r$q); r$scale * (b / (1 - b))^(1 / r$a) },
  references = c("McDonald (1984), Some generalized functions for the size distribution of income. Econometrica, 52, 647-663.", KLUG),
  reference = ref_impl("gamlss.dist", function(y, r) gamlss.dist::dGB2(y, r$scale, r$a, r$p, r$q, log = TRUE)),
  stan = sfun("gb2", c("scale", "a", "p", "q"),
    "    return log(abs(a)) + (a * p - 1) * log(y) - a * p * log(scale) - lbeta(p, q) - (p + q) * log1p(pow(y / scale, a));"),
  status = "validated")

posfam("beta_prime", "Beta prime (inverted beta)", list(shp("alpha"), shp("beta")), aliases = c("inverted_beta", "beta2"),
  formula = "f(y) = y^(alpha - 1) (1 + y)^(-alpha - beta) / B(alpha, beta)",
  mean_formula = "alpha / (beta - 1) for beta > 1",
  lpdf = function(y, p) (p$alpha - 1) * log(y) - (p$alpha + p$beta) * log1p(y) - lbeta(p$alpha, p$beta),
  cdf = function(q, p) stats::pbeta(q / (1 + q), p$alpha, p$beta),
  quantile = function(pr, p) { b <- stats::qbeta(pr, p$alpha, p$beta); b / (1 - b) },
  mean = function(p) ifelse(p$beta > 1, p$alpha / (p$beta - 1), Inf),
  references = JKB, reference = ref_impl("extraDistr", function(y, p) extraDistr::dbetapr(y, p$alpha, p$beta, log = TRUE)),
  stan = sfun("beta_prime", c("alpha", "beta"), "    return (alpha - 1) * log(y) - (alpha + beta) * log1p(y) - lbeta(alpha, beta);"),
  status = "validated")

posfam("beta_prime_mean", "Beta prime (mean and precision, Bourguignon et al.)", list(mean_par(), P("phi", "positive", role = "precision")),
  formula = "alpha = mu (1 + phi), beta = phi + 2", mean_formula = "mu", variance_formula = "mu (1 + mu) / phi",
  lpdf = function(y, p) { a <- p$mu * (1 + p$phi); b <- p$phi + 2; (a - 1) * log(y) - (a + b) * log1p(y) - lbeta(a, b) },
  cdf = function(q, p) stats::pbeta(q / (1 + q), p$mu * (1 + p$phi), p$phi + 2),
  mean = function(p) p$mu, variance = function(p) p$mu * (1 + p$mu) / p$phi,
  references = "Bourguignon, Santos-Neto and de Castro (2021), A new regression model for positive random variables with skewed and long tail. Metron, 79, 33-55.",
  reference = ref_impl("extraDistr", function(y, p) extraDistr::dbetapr(y, p$mu * (1 + p$phi), p$phi + 2, log = TRUE)),
  stan = sfun("beta_prime_mean", c("mu", "phi"),
    "    real a = mu * (1 + phi);\n    real b = phi + 2;\n    return (a - 1) * log(y) - (a + b) * log1p(y) - lbeta(a, b);"),
  start = function(y) list(mu = mean(y), phi = 2), status = "validated")

posfam("frechet", "Frechet (inverse Weibull)", list(shp("alpha"), P("scale", "positive", role = "scale")),
  aliases = c("inverse_weibull", "gumbel_type2"), formula = "F(y) = exp(-(y / scale)^(-alpha))",
  lpdf = function(y, p) log(p$alpha) - log(p$scale) - (1 + p$alpha) * log(y / p$scale) - (y / p$scale)^(-p$alpha),
  cdf = function(q, p) exp(-(q / p$scale)^(-p$alpha)), quantile = function(pr, p) p$scale * (-log(pr))^(-1 / p$alpha),
  mean = function(p) ifelse(p$alpha > 1, p$scale * gamma(1 - 1 / p$alpha), Inf),
  references = JKB, reference = ref_impl("extraDistr", function(y, p) extraDistr::dfrechet(y, p$alpha, 0, p$scale, log = TRUE)),
  stan = stan_native("frechet_lpdf({y} | {alpha}, {scale})", "frechet_lcdf({y} | {alpha}, {scale})", "frechet_lccdf({y} | {alpha}, {scale})"),
  status = "validated")

posfam("gompertz", "Gompertz", list(P("eta", "positive", role = "shape"), P("b", "positive", role = "rate")),
  formula = "F(y) = 1 - exp(-eta (exp(b y) - 1))",
  lpdf = function(y, p) log(p$eta) + log(p$b) + p$b * y + p$eta - p$eta * exp(p$b * y),
  cdf = function(q, p) 1 - exp(-p$eta * expm1(p$b * q)), quantile = function(pr, p) log1p(-log1p(-pr) / p$eta) / p$b,
  references = "Gompertz (1825), Phil. Trans. Roy. Soc., 115, 513-583.",
  reference = ref_impl("extraDistr", function(y, p) extraDistr::dgompertz(y, p$eta * p$b, p$b, log = TRUE)),
  stan = sfun("gompertz", c("eta", "b"), "    return log(eta) + log(b) + b * y + eta - eta * exp(b * y);",
              lccdf = "    return -eta * expm1(b * y);"),
  status = "validated")

posfam("defective_gompertz", "Defective Gompertz (cure fraction)", list(P("eta", "positive", role = "shape"), P("b", "real", role = "rate", link = "identity",
  interpretation = "b < 0 gives a cure fraction exp(eta / b)")),
  formula = "S(y) = exp(-eta (exp(b y) - 1) / b); improper when b < 0",
  lpdf = function(y, p) log(p$eta) + p$b * y - p$eta * ifelse(p$b == 0, y, expm1(p$b * y) / p$b),
  cdf = function(q, p) -expm1(-p$eta * ifelse(p$b == 0, q, expm1(p$b * q) / p$b)),
  quantile = function(pr, p) { a <- -p$b * log1p(-pr) / p$eta
    ifelse(p$b == 0, -log1p(-pr) / p$eta, ifelse(a <= -1, Inf, log1p(a) / p$b)) },
  numeric_moments = FALSE, numeric_cdf = FALSE,
  regression_notes = "Used for cure rate regression; the cure fraction is exp(eta / b) for b < 0.",
  references = "Gieser, Chang, Rao, Shuster and Pullen (1998), Modelling cure rates using the Gompertz model with covariate information. Statistics in Medicine, 17, 831-839.",
  reference = ref_impl("definition", function(y, p) log(p$eta) + p$b * y - p$eta * (exp(p$b * y) - 1) / p$b),
  status = "implemented")

posfam("shifted_gompertz", "Shifted Gompertz", list(P("b", "positive", role = "scale"), shp("eta")),
  formula = "F(y) = (1 - exp(-b y)) exp(-eta exp(-b y))",
  lpdf = function(y, p) { e <- exp(-p$b * y); log(p$b) + log(e) - p$eta * e + log1p(p$eta * (1 - e)) },
  cdf = function(q, p) (1 - exp(-p$b * q)) * exp(-p$eta * exp(-p$b * q)),
  references = "Bemmaor (1994), Modeling the diffusion of new durable goods.",
  reference = ref_impl("definition", function(y, p) log(p$b * exp(-p$b * y) * exp(-p$eta * exp(-p$b * y)) * (1 + p$eta * (1 - exp(-p$b * y))))),
  status = "validated")

posfam("rayleigh", "Rayleigh", list(P("sigma", "positive", role = "scale")),
  formula = "f(y) = y / sigma^2 exp(-y^2 / (2 sigma^2))", mean_formula = "sigma sqrt(pi/2)",
  lpdf = function(y, p) log(y) - 2 * log(p$sigma) - y^2 / (2 * p$sigma^2),
  cdf = function(q, p) 1 - exp(-q^2 / (2 * p$sigma^2)), quantile = function(pr, p) p$sigma * sqrt(-2 * log1p(-pr)),
  mean = function(p) p$sigma * sqrt(pi / 2), variance = function(p) (4 - pi) / 2 * p$sigma^2,
  references = JKB, reference = ref_impl("extraDistr", function(y, p) extraDistr::drayleigh(y, p$sigma, log = TRUE)),
  stan = stan_native("rayleigh_lpdf({y} | {sigma})", "rayleigh_lcdf({y} | {sigma})", "rayleigh_lccdf({y} | {sigma})", "rayleigh_rng({sigma})"),
  status = "validated")

posfam("rice", "Rice", list(P("nu", "nonneg", role = "location"), P("sigma", "positive", role = "scale")),
  formula = "f(y) = y / sigma^2 exp(-(y^2 + nu^2) / (2 sigma^2)) I0(y nu / sigma^2)",
  lpdf = function(y, p) { z <- y * p$nu / p$sigma^2
    log(y) - 2 * log(p$sigma) - (y - p$nu)^2 / (2 * p$sigma^2) + log(besselI(z, 0, expon.scaled = TRUE)) },
  references = "Rice (1945), Mathematical analysis of random noise. Bell System Technical Journal, 24, 46-156.",
  reference = ref_impl("VGAM", function(y, p) VGAM::drice(y, p$sigma, p$nu, log = TRUE)),
  stan = sfun("rice", c("nu", "sigma"),
    "    return log(y) - 2 * log(sigma) - (square(y) + square(nu)) / (2 * square(sigma)) + log_modified_bessel_first_kind(0, y * nu / square(sigma));"),
  status = "validated")

posfam("nakagami", "Nakagami", list(P("m", "positive", role = "shape", lower = 0.5), P("omega", "positive", role = "spread")),
  formula = "f(y) = 2 m^m / (Gamma(m) omega^m) y^(2m - 1) exp(-m y^2 / omega)",
  lpdf = function(y, p) log(2) + p$m * log(p$m) - lgamma(p$m) - p$m * log(p$omega) + (2 * p$m - 1) * log(y) - p$m * y^2 / p$omega,
  cdf = function(q, p) stats::pgamma(p$m * q^2 / p$omega, p$m),
  quantile = function(pr, p) sqrt(stats::qgamma(pr, p$m) * p$omega / p$m),
  mean = function(p) exp(lgamma(p$m + 0.5) - lgamma(p$m)) * sqrt(p$omega / p$m),
  references = "Nakagami (1960), The m-distribution.",
  reference = ref_impl("VGAM", function(y, p) VGAM::dnaka(y, scale = p$omega, shape = p$m, log = TRUE)),
  stan = sfun("nakagami", c("m", "omega"),
    "    return log2() + m * log(m) - lgamma(m) - m * log(omega) + (2 * m - 1) * log(y) - m * square(y) / omega;"),
  status = "validated")

posfam("maxwell_boltzmann", "Maxwell-Boltzmann", list(P("a", "positive", role = "scale")), aliases = "maxwell",
  formula = "f(y) = sqrt(2/pi) y^2 exp(-y^2 / (2 a^2)) / a^3",
  lpdf = function(y, p) 0.5 * log(2 / pi) + 2 * log(y) - y^2 / (2 * p$a^2) - 3 * log(p$a),
  cdf = function(q, p) stats::pchisq(q^2 / p$a^2, 3), quantile = function(pr, p) p$a * sqrt(stats::qchisq(pr, 3)),
  mean = function(p) 2 * p$a * sqrt(2 / pi),
  references = JKB, reference = ref_impl("VGAM", function(y, p) VGAM::dmaxwell(y, 1 / p$a^2, log = TRUE)),
  status = "validated")

posfam("chi_squared", "Chi-squared", list(P("df", "positive", role = "degrees of freedom")), aliases = "chisq",
  formula = "Gamma(df/2, rate 1/2)", lpdf = function(y, p) stats::dchisq(y, p$df, log = TRUE),
  cdf = function(q, p) stats::pchisq(q, p$df), quantile = function(pr, p) stats::qchisq(pr, p$df),
  mean = function(p) p$df, variance = function(p) 2 * p$df,
  references = JKB, reference = ref_impl("gamma", function(y, p) stats::dgamma(y, p$df / 2, 0.5, log = TRUE)),
  stan = stan_native("chi_square_lpdf({y} | {df})", "chi_square_lcdf({y} | {df})", "chi_square_lccdf({y} | {df})"),
  status = "validated")

posfam("noncentral_chi_squared", "Noncentral chi-squared", list(P("df", "positive"), P("ncp", "nonneg", role = "noncentrality")),
  lpdf = function(y, p) stats::dchisq(y, p$df, p$ncp, log = TRUE), cdf = function(q, p) stats::pchisq(q, p$df, p$ncp),
  quantile = function(pr, p) stats::qchisq(pr, p$df, p$ncp), mean = function(p) p$df + p$ncp,
  variance = function(p) 2 * (p$df + 2 * p$ncp), formula = "Poisson(ncp/2) mixture of central chi-squared laws",
  references = JKB, reference = ref_impl("series", function(y, p) log(vapply(seq_along(y), function(i)
    sum(stats::dpois(0:400, p$ncp[i] / 2) * stats::dchisq(y[i], p$df[i] + 2 * (0:400))), 0))),
  status = "validated")

posfam("chi", "Chi", list(P("k", "positive", role = "degrees of freedom")),
  formula = "Y = sqrt(chi-squared(k))",
  lpdf = function(y, p) (1 - p$k / 2) * log(2) + (p$k - 1) * log(y) - y^2 / 2 - lgamma(p$k / 2),
  cdf = function(q, p) stats::pchisq(q^2, p$k), quantile = function(pr, p) sqrt(stats::qchisq(pr, p$k)),
  mean = function(p) sqrt(2) * exp(lgamma((p$k + 1) / 2) - lgamma(p$k / 2)),
  references = JKB, reference = ref_impl("chisq", function(y, p) stats::dchisq(y^2, p$k, log = TRUE) + log(2 * y)),
  status = "validated")

posfam("inverse_chi_squared", "Inverse chi-squared (scaled)", list(P("df", "positive"), P("tau2", "positive", role = "scale")),
  aliases = "scaled_inverse_chi_squared", formula = "Inverse gamma(df/2, df tau2/2)",
  lpdf = function(y, p) { a <- p$df / 2; b <- p$df * p$tau2 / 2; a * log(b) - lgamma(a) - (a + 1) * log(y) - b / y },
  cdf = function(q, p) stats::pgamma(p$df * p$tau2 / (2 * q), p$df / 2, lower.tail = FALSE),
  references = "Gelman et al. (2013), Bayesian Data Analysis, 3rd ed.",
  reference = ref_impl("extraDistr", function(y, p) extraDistr::dinvchisq(y, p$df, p$tau2, log = TRUE)),
  stan = stan_native("scaled_inv_chi_square_lpdf({y} | {df}, sqrt({tau2}))"), status = "validated")

posfam("f", "F (Fisher-Snedecor)", list(P("df1", "positive"), P("df2", "positive")), aliases = c("fisher_snedecor", "F"),
  lpdf = function(y, p) stats::df(y, p$df1, p$df2, log = TRUE), cdf = function(q, p) stats::pf(q, p$df1, p$df2),
  quantile = function(pr, p) stats::qf(pr, p$df1, p$df2),
  mean = function(p) ifelse(p$df2 > 2, p$df2 / (p$df2 - 2), NA),
  references = JKB, reference = ref_impl("beta prime", function(y, p) { a <- p$df1 / 2; b <- p$df2 / 2; s <- p$df1 / p$df2
    log(s) + (a - 1) * log(s * y) - (a + b) * log1p(s * y) - lbeta(a, b) }),
  stan = sfun("fdist", c("df1", "df2"),
    "    real a = df1 / 2;\n    real b = df2 / 2;\n    real s = df1 / df2;\n    return log(s) + (a - 1) * log(s * y) - (a + b) * log1p(s * y) - lbeta(a, b);"),
  status = "validated")

posfam("noncentral_f", "Noncentral F", list(P("df1", "positive"), P("df2", "positive"), P("ncp", "nonneg")),
  lpdf = function(y, p) stats::df(y, p$df1, p$df2, p$ncp, log = TRUE), cdf = function(q, p) stats::pf(q, p$df1, p$df2, p$ncp),
  references = JKB, reference = ref_impl("series", function(y, p) log(vapply(seq_along(y), function(i) { j <- 0:400
    s <- p$df1[i] / p$df2[i]; sum(stats::dpois(j, p$ncp[i] / 2) * stats::dbeta(s * y[i] / (1 + s * y[i]), p$df1[i] / 2 + j, p$df2[i] / 2) * s / (1 + s * y[i])^2) }, 0))),
  status = "validated")

posfam("hotelling_t2", "Hotelling T-squared", list(P("p", "positive_integer", role = "dimension", data = TRUE), P("m", "positive", role = "degrees of freedom")),
  formula = "(m - p + 1) / (p m) T^2 ~ F(p, m - p + 1)",
  lpdf = function(y, q) { s <- (q$m - q$p + 1) / (q$p * q$m); stats::df(s * y, q$p, q$m - q$p + 1, log = TRUE) + log(s) },
  cdf = function(x, q) stats::pf((q$m - q$p + 1) / (q$p * q$m) * x, q$p, q$m - q$p + 1),
  references = "Hotelling (1931), Ann. Math. Stat., 2, 360-378.",
  reference = ref_impl("F", function(y, q) { s <- (q$m - q$p + 1) / (q$p * q$m); stats::df(s * y, q$p, q$m - q$p + 1, log = TRUE) + log(s) }),
  status = "validated")

posfam("levy", "Levy", list(P("mu", "real", role = "location", link = "identity", data = TRUE), P("c", "positive", role = "scale")),
  formula = "f(y) = sqrt(c / (2 pi)) exp(-c / (2 (y - mu))) / (y - mu)^(3/2), y > mu",
  lpdf = function(y, p) ifelse(y > p$mu, 0.5 * log(p$c / (2 * pi)) - p$c / (2 * (y - p$mu)) - 1.5 * log(y - p$mu), -Inf),
  cdf = function(q, p) ifelse(q > p$mu, 2 * stats::pnorm(-sqrt(p$c / (q - p$mu))), 0),
  quantile = function(pr, p) p$mu + p$c / stats::qnorm(1 - pr / 2)^2,
  mean = function(p) rep(Inf, length(p$c)), numeric_moments = FALSE,
  support = function(p) list(lower = p$mu, upper = Inf),
  references = JKB, reference = ref_impl("rmutil", function(y, p) log(rmutil::dlevy(y, p$mu, p$c))),
  stan = sfun("levy", c("mu", "c"), "    return 0.5 * log(c / (2 * pi())) - c / (2 * (y - mu)) - 1.5 * log(y - mu);"),
  status = "validated")

posfam("half_normal", "Half-normal", list(P("sigma", "positive", role = "scale")),
  lpdf = function(y, p) log(2) + stats::dnorm(y, 0, p$sigma, log = TRUE), cdf = function(q, p) 2 * stats::pnorm(q, 0, p$sigma) - 1,
  quantile = function(pr, p) stats::qnorm((1 + pr) / 2, 0, p$sigma), mean = function(p) p$sigma * sqrt(2 / pi),
  variance = function(p) p$sigma^2 * (1 - 2 / pi),
  references = JKB, reference = ref_impl("extraDistr", function(y, p) extraDistr::dhnorm(y, p$sigma, log = TRUE)),
  stan = sfun("half_normal", "sigma", "    return log2() + normal_lpdf(y | 0, sigma);"), status = "validated")

posfam("half_cauchy", "Half-Cauchy", list(P("sigma", "positive", role = "scale")),
  lpdf = function(y, p) log(2) + stats::dcauchy(y, 0, p$sigma, log = TRUE), cdf = function(q, p) 2 * stats::pcauchy(q, 0, p$sigma) - 1,
  quantile = function(pr, p) stats::qcauchy((1 + pr) / 2, 0, p$sigma), numeric_moments = FALSE,
  references = JKB, reference = ref_impl("extraDistr", function(y, p) extraDistr::dhcauchy(y, p$sigma, log = TRUE)),
  stan = sfun("half_cauchy", "sigma", "    return log2() + cauchy_lpdf(y | 0, sigma);"), status = "validated")

posfam("half_t", "Half-t", list(P("nu", "positive"), P("sigma", "positive", role = "scale")),
  lpdf = function(y, p) log(2) + stats::dt(y / p$sigma, p$nu, log = TRUE) - log(p$sigma),
  cdf = function(q, p) 2 * stats::pt(q / p$sigma, p$nu) - 1, quantile = function(pr, p) p$sigma * stats::qt((1 + pr) / 2, p$nu),
  references = JKB, reference = ref_impl("extraDistr", function(y, p) extraDistr::dht(y, p$nu, p$sigma, log = TRUE)),
  stan = sfun("half_t", c("nu", "sigma"), "    return log2() + student_t_lpdf(y | nu, 0, sigma);"), status = "validated")

posfam("folded_normal", "Folded normal", list(P("mu", "real", role = "location", link = "identity"), P("sigma", "positive", role = "scale")),
  lpdf = function(y, p) log(stats::dnorm(y, p$mu, p$sigma) + stats::dnorm(-y, p$mu, p$sigma)),
  cdf = function(q, p) stats::pnorm(q, p$mu, p$sigma) - stats::pnorm(-q, p$mu, p$sigma),
  rng = function(n, p) abs(stats::rnorm(n, p$mu, p$sigma)),
  references = "Leone, Nelson and Nottingham (1961), Technometrics, 3, 543-550.",
  reference = ref_impl("VGAM", function(y, p) VGAM::dfoldnorm(y, p$mu, p$sigma, log = TRUE)),
  stan = sfun("folded_normal", c("mu", "sigma"), "    return log_sum_exp(normal_lpdf(y | mu, sigma), normal_lpdf(-y | mu, sigma));"),
  status = "validated")

# Lindley type ---------------------------------------------------------------------
posfam("lindley", "Lindley", list(P("theta", "positive", role = "rate")),
  formula = "f(y) = theta^2 / (1 + theta) (1 + y) exp(-theta y)", mean_formula = "(theta + 2) / (theta (theta + 1))",
  lpdf = function(y, p) 2 * log(p$theta) - log1p(p$theta) + log1p(y) - p$theta * y,
  cdf = function(q, p) 1 - (1 + p$theta + p$theta * q) / (1 + p$theta) * exp(-p$theta * q),
  mean = function(p) (p$theta + 2) / (p$theta * (p$theta + 1)),
  references = "Lindley (1958), Fiducial distributions and Bayes' theorem. JRSS B, 20, 102-107.",
  reference = ref_impl("VGAM", function(y, p) VGAM::dlind(y, p$theta, log = TRUE)),
  stan = sfun("lindley", "theta", "    return 2 * log(theta) - log1p(theta) + log1p(y) - theta * y;"), status = "validated")

posfam("xgamma", "Xgamma (Sen, Maiti and Chandra)", list(P("theta", "positive", role = "rate")),
  formula = "f(y) = theta^2 / (1 + theta) (1 + theta y^2 / 2) exp(-theta y)",
  lpdf = function(y, p) 2 * log(p$theta) - log1p(p$theta) + log1p(p$theta * y^2 / 2) - p$theta * y,
  cdf = function(q, p) 1 - (1 + p$theta + p$theta * q + p$theta^2 * q^2 / 2) / (1 + p$theta) * exp(-p$theta * q),
  references = "Sen, Maiti and Chandra (2016), The xgamma distribution. J. Stat. Theory Appl.",
  reference = ref_impl("mixture", function(y, p) log(p$theta / (1 + p$theta) * stats::dexp(y, p$theta) + 1 / (1 + p$theta) * stats::dgamma(y, 3, p$theta))),
  status = "validated")

posfam("generalized_lindley", "Generalized (exponentiated) Lindley", list(P("theta", "positive", role = "rate"), shp("alpha")),
  formula = "F(y) = (1 - (1 + theta + theta y)/(1 + theta) exp(-theta y))^alpha",
  lpdf = function(y, p) { G <- 1 - (1 + p$theta + p$theta * y) / (1 + p$theta) * exp(-p$theta * y)
    log(p$alpha) + 2 * log(p$theta) - log1p(p$theta) + log1p(y) - p$theta * y + (p$alpha - 1) * log(G) },
  cdf = function(q, p) (1 - (1 + p$theta + p$theta * q) / (1 + p$theta) * exp(-p$theta * q))^p$alpha,
  references = "Nadarajah, Bakouch and Tahmasbi (2011), A generalized Lindley distribution. Sankhya B, 73, 331-359.",
  reference = ref_impl("exponentiated Lindley", function(y, p) { G <- 1 - (1 + p$theta + p$theta * y) / (1 + p$theta) * exp(-p$theta * y)
    log(p$alpha * G^(p$alpha - 1) * p$theta^2 / (1 + p$theta) * (1 + y) * exp(-p$theta * y)) }),
  status = "validated")

# Birnbaum-Saunders family ----------------------------------------------------------
bs_lpdf <- function(y, a, b) { z <- (sqrt(y / b) - sqrt(b / y)) / a
  stats::dnorm(z, log = TRUE) + log(sqrt(y / b) + sqrt(b / y)) - log(2 * a * y) }
posfam("birnbaum_saunders", "Birnbaum-Saunders (fatigue life)", list(P("alpha", "positive", role = "shape"), P("beta", "positive", role = "scale", offset_allowed = TRUE)),
  aliases = c("fatigue_life", "BS"), formula = "Phi((sqrt(y/beta) - sqrt(beta/y)) / alpha)",
  mean_formula = "beta (1 + alpha^2 / 2)", regression_notes = "beta is the median.",
  lpdf = function(y, p) bs_lpdf(y, p$alpha, p$beta),
  cdf = function(q, p) stats::pnorm((sqrt(q / p$beta) - sqrt(p$beta / q)) / p$alpha),
  quantile = function(pr, p) { z <- p$alpha * stats::qnorm(pr) / 2; p$beta * (z + sqrt(z^2 + 1))^2 },
  mean = function(p) p$beta * (1 + p$alpha^2 / 2), variance = function(p) (p$alpha * p$beta)^2 * (1 + 5 * p$alpha^2 / 4),
  references = "Birnbaum and Saunders (1969), A new family of life distributions. J. Appl. Prob., 6, 319-327.",
  reference = ref_impl("VGAM", function(y, p) VGAM::dbisa(y, scale = p$beta, shape = p$alpha, log = TRUE)),
  stan = sfun("bisa", c("alpha", "beta"),
    "    real z = (sqrt(y / beta) - sqrt(beta / y)) / alpha;\n    return std_normal_lpdf(z) + log(sqrt(y / beta) + sqrt(beta / y)) - log(2 * alpha * y);"),
  status = "validated")

posfam("birnbaum_saunders_mean", "Birnbaum-Saunders (mean and precision, Santos-Neto et al.)", list(mean_par(), P("delta", "positive", role = "precision")),
  formula = "alpha = sqrt(2 / delta), beta = delta mu / (delta + 1)", mean_formula = "mu",
  variance_formula = "mu^2 (2 delta + 5) / (delta + 1)^2",
  lpdf = function(y, p) bs_lpdf(y, sqrt(2 / p$delta), p$delta * p$mu / (p$delta + 1)),
  cdf = function(q, p) { a <- sqrt(2 / p$delta); b <- p$delta * p$mu / (p$delta + 1); stats::pnorm((sqrt(q / b) - sqrt(b / q)) / a) },
  mean = function(p) p$mu, variance = function(p) p$mu^2 * (2 * p$delta + 5) / (p$delta + 1)^2,
  references = "Santos-Neto, Cysneiros, Leiva and Barros (2016), Reparameterized Birnbaum-Saunders regression models with varying precision. Electronic J. Stat., 10, 2825-2855.",
  reference = ref_impl("VGAM", function(y, p) VGAM::dbisa(y, scale = p$delta * p$mu / (p$delta + 1), shape = sqrt(2 / p$delta), log = TRUE)),
  start = function(y) list(mu = mean(y), delta = 5), status = "validated")

posfam("generalized_birnbaum_saunders", "Generalized Birnbaum-Saunders (power, GBS)", list(mean_par(), scale_par("sigma"), P("nu", "unit", role = "shape")),
  formula = "F(y) = Phi((y^(1 - nu)/sqrt(mu) - sqrt(mu)/y^nu) / sigma)",
  lpdf = function(y, p) -0.5 * log(2 * pi) - log(p$sigma) - 0.5 * log(p$mu) - p$nu * log(y) + log(1 - p$nu + p$nu * p$mu / y) -
    (y - p$mu)^2 / (2 * p$sigma^2 * p$mu * y^(2 * p$nu)),
  cdf = function(q, p) stats::pnorm((q^(1 - p$nu) / sqrt(p$mu) - sqrt(p$mu) / q^p$nu) / p$sigma),
  references = "arXiv 2609.21063 (2026), A new generalized Birnbaum-Saunders regression model.",
  reference = ref_impl("BS at nu = 1/2", function(y, p) log(stats::dnorm((y^(1 - p$nu) / sqrt(p$mu) - sqrt(p$mu) / y^p$nu) / p$sigma) *
    ((1 - p$nu) * y^(-p$nu) / sqrt(p$mu) + p$nu * sqrt(p$mu) * y^(-p$nu - 1)) / p$sigma)),
  status = "validated")

# Others on the Wikipedia list --------------------------------------------------------
posfam("log_cauchy", "Log-Cauchy", list(loc_par("mu"), scale_par("sigma")),
  lpdf = function(y, p) stats::dcauchy(log(y), p$mu, p$sigma, log = TRUE) - log(y),
  cdf = function(q, p) stats::pcauchy(log(q), p$mu, p$sigma), quantile = function(pr, p) exp(stats::qcauchy(pr, p$mu, p$sigma)),
  numeric_moments = FALSE, references = JKB,
  reference = ref_impl("transform", function(y, p) log(stats::dcauchy(log(y), p$mu, p$sigma) / y)), status = "validated")

posfam("log_laplace", "Log-Laplace", list(loc_par("mu"), scale_par("b")),
  lpdf = function(y, p) -log(2 * p$b) - abs(log(y) - p$mu) / p$b - log(y),
  cdf = function(q, p) { z <- (log(q) - p$mu) / p$b; ifelse(z < 0, exp(z) / 2, 1 - exp(-z) / 2) },
  references = "Kozubowski and Podgorski (2003).",
  reference = ref_impl("extraDistr", function(y, p) extraDistr::dlaplace(log(y), p$mu, p$b, log = TRUE) - log(y)), status = "validated")

posfam("exponential_logarithmic", "Exponential-logarithmic", list(P("p", "unit", role = "shape"), P("beta", "positive", role = "rate")),
  formula = "f(y) = -beta (1 - p) exp(-beta y) / (log(p) (1 - (1 - p) exp(-beta y)))",
  lpdf = function(y, q) log(q$beta) + log1p(-q$p) - q$beta * y - log(-log(q$p)) - log1p(-(1 - q$p) * exp(-q$beta * y)),
  cdf = function(x, q) 1 - log1p(-(1 - q$p) * exp(-q$beta * x)) / log(q$p),
  references = "Tahmasbi and Rezaei (2008), A two-parameter lifetime distribution with decreasing failure rate. CSDA, 52, 3889-3901.",
  reference = ref_impl("definition", function(y, q) log(-q$beta * (1 - q$p) * exp(-q$beta * y) / (log(q$p) * (1 - (1 - q$p) * exp(-q$beta * y))))),
  status = "validated")

posfam("generalized_exponential", "Generalized (exponentiated) exponential", list(shp("alpha"), P("lambda", "positive", role = "rate")),
  formula = "F(y) = (1 - exp(-lambda y))^alpha",
  lpdf = function(y, p) log(p$alpha) + log(p$lambda) - p$lambda * y + (p$alpha - 1) * log1p(-exp(-p$lambda * y)),
  cdf = function(q, p) (-expm1(-p$lambda * q))^p$alpha, quantile = function(pr, p) -log1p(-pr^(1 / p$alpha)) / p$lambda,
  references = "Gupta and Kundu (1999), Generalized exponential distributions. Aust. N. Z. J. Stat., 41, 173-188.",
  reference = ref_impl("definition", function(y, p) log(p$alpha * p$lambda * exp(-p$lambda * y) * (1 - exp(-p$lambda * y))^(p$alpha - 1))),
  stan = sfun("gen_exp", c("alpha", "lambda"), "    return log(alpha) + log(lambda) - lambda * y + (alpha - 1) * log1m_exp(-lambda * y);"),
  status = "validated")

posfam("exponentiated_weibull", "Exponentiated Weibull", list(shp("shape"), P("scale", "positive", role = "scale"), shp("alpha")),
  formula = "F(y) = (1 - exp(-(y/scale)^shape))^alpha",
  lpdf = function(y, p) log(p$alpha) + stats::dweibull(y, p$shape, p$scale, log = TRUE) + (p$alpha - 1) * stats::pweibull(y, p$shape, p$scale, log.p = TRUE),
  cdf = function(q, p) stats::pweibull(q, p$shape, p$scale)^p$alpha,
  quantile = function(pr, p) stats::qweibull(pr^(1 / p$alpha), p$shape, p$scale),
  references = "Mudholkar and Srivastava (1993), Exponentiated Weibull family for analyzing bathtub failure-rate data. IEEE Trans. Reliability, 42, 299-302.",
  reference = ref_impl("reliaR", function(y, p) log(p$alpha * stats::dweibull(y, p$shape, p$scale) * stats::pweibull(y, p$shape, p$scale)^(p$alpha - 1))),
  stan = sfun("exp_weibull", c("shape", "scale", "alpha"),
    "    return log(alpha) + weibull_lpdf(y | shape, scale) + (alpha - 1) * weibull_lcdf(y | shape, scale);"),
  status = "validated")

posfam("chen", "Chen (bathtub)", list(shp("lambda"), shp("beta")),
  formula = "F(y) = 1 - exp(lambda (1 - exp(y^beta)))",
  lpdf = function(y, p) log(p$lambda) + log(p$beta) + (p$beta - 1) * log(y) + y^p$beta + p$lambda * (1 - exp(y^p$beta)),
  cdf = function(q, p) 1 - exp(p$lambda * (1 - exp(q^p$beta))), quantile = function(pr, p) log1p(-log1p(-pr) / p$lambda)^(1 / p$beta),
  references = "Chen (2000), A new two-parameter lifetime distribution with bathtub shape or increasing failure rate function. Stat. Prob. Letters, 49, 155-161.",
  reference = ref_impl("definition", function(y, p) log(p$lambda * p$beta * y^(p$beta - 1) * exp(y^p$beta) * exp(p$lambda * (1 - exp(y^p$beta))))),
  status = "validated")

posfam("power_generalized_weibull", "Power generalized Weibull (Nikulin and Haghighi)", list(P("scale", "positive", role = "scale"), shp("nu"), shp("gamma")),
  formula = "S(y) = exp(1 - (1 + (y/scale)^nu)^(1/gamma))",
  lpdf = function(y, p) { z <- (y / p$scale)^p$nu; w <- (1 + z)^(1 / p$gamma)
    log(p$nu / p$gamma) - log(y) + log(z) + (1 / p$gamma - 1) * log1p(z) + 1 - w },
  cdf = function(q, p) 1 - exp(1 - (1 + (q / p$scale)^p$nu)^(1 / p$gamma)),
  quantile = function(pr, p) p$scale * ((1 - log1p(-pr))^p$gamma - 1)^(1 / p$nu),
  references = "Nikulin and Haghighi (2006), A chi-squared test for the generalized power Weibull family. J. Math. Sci., 133, 1333-1341.",
  reference = ref_impl("definition", function(y, p) log(p$nu / (p$gamma * p$scale^p$nu) * y^(p$nu - 1) * (1 + (y / p$scale)^p$nu)^(1 / p$gamma - 1) * exp(1 - (1 + (y / p$scale)^p$nu)^(1 / p$gamma)))),
  status = "validated")

posfam("generalized_inverse_gaussian", "Generalized inverse Gaussian", list(P("lambda", "real", role = "shape", link = "identity"), P("chi", "positive"), P("psi", "positive")),
  aliases = "GIG", formula = "f(y) = (psi/chi)^(lambda/2) / (2 K_lambda(sqrt(chi psi))) y^(lambda - 1) exp(-(chi/y + psi y)/2)",
  lpdf = function(y, p) p$lambda / 2 * log(p$psi / p$chi) - log(2) - log_besselK(sqrt(p$chi * p$psi), p$lambda) +
    (p$lambda - 1) * log(y) - (p$chi / y + p$psi * y) / 2,
  mean = function(p) sqrt(p$chi / p$psi) * exp(log_besselK(sqrt(p$chi * p$psi), p$lambda + 1) - log_besselK(sqrt(p$chi * p$psi), p$lambda)),
  references = "Jorgensen (1982), Statistical Properties of the Generalized Inverse Gaussian Distribution. Springer.",
  reference = ref_impl("GIGrvg", function(y, p) log(mapply(function(x, l, c, s) GIGrvg::dgig(x, l, c, s), y, p$lambda, p$chi, p$psi))),
  stan = sfun("gig", c("lambda", "chi", "psi"),
    "    return lambda / 2 * log(psi / chi) - log2() - rc_log_besselK(lambda, sqrt(chi * psi)) + (lambda - 1) * log(y) - (chi / y + psi * y) / 2;",
    needs = "log_besselK"),
  status = "validated")

posfam("hypoexponential", "Hypoexponential (two distinct rates)", list(P("lambda1", "positive", role = "rate"), P("lambda2", "positive", role = "rate")),
  formula = "f(y) = l1 l2 / (l2 - l1) (exp(-l1 y) - exp(-l2 y))",
  lpdf = function(y, p) log(p$lambda1 * p$lambda2 / abs(p$lambda2 - p$lambda1) * abs(exp(-p$lambda1 * y) - exp(-p$lambda2 * y))),
  cdf = function(q, p) 1 - (p$lambda2 * exp(-p$lambda1 * q) - p$lambda1 * exp(-p$lambda2 * q)) / (p$lambda2 - p$lambda1),
  rng = function(n, p) stats::rexp(n, p$lambda1) + stats::rexp(n, p$lambda2),
  mean = function(p) 1 / p$lambda1 + 1 / p$lambda2,
  references = "Ross (2014), Introduction to Probability Models.",
  reference = ref_impl("convolution", function(y, p) log(vapply(seq_along(y), function(i) stats::integrate(function(t)
    stats::dexp(t, p$lambda1[i]) * stats::dexp(y[i] - t, p$lambda2[i]), 0, y[i], rel.tol = 1e-12)$value, 0))),
  status = "validated")

posfam("hyperexponential", "Hyperexponential (two phase)", list(P("prob", "unit", role = "weight"), P("lambda1", "positive"), P("lambda2", "positive")),
  formula = "f(y) = p l1 exp(-l1 y) + (1 - p) l2 exp(-l2 y)",
  lpdf = function(y, p) log(p$prob * stats::dexp(y, p$lambda1) + (1 - p$prob) * stats::dexp(y, p$lambda2)),
  cdf = function(q, p) p$prob * stats::pexp(q, p$lambda1) + (1 - p$prob) * stats::pexp(q, p$lambda2),
  mean = function(p) p$prob / p$lambda1 + (1 - p$prob) / p$lambda2,
  references = "Ross (2014).", reference = ref_impl("mixture", function(y, p) log(p$prob * p$lambda1 * exp(-p$lambda1 * y) + (1 - p$prob) * p$lambda2 * exp(-p$lambda2 * y))),
  status = "validated")

posfam("bimodal_gamma", "Bimodal gamma (Vila et al.)", list(P("alpha", "positive"), P("beta", "positive"), P("delta", "real", role = "shape", link = "identity")),
  formula = "f(y) = (1 + (1 - delta y)^2) gamma(y | alpha, beta) / Z, Z = 2 - 2 delta alpha/beta + delta^2 alpha (alpha + 1)/beta^2",
  lpdf = function(y, p) { Z <- 2 - 2 * p$delta * p$alpha / p$beta + p$delta^2 * p$alpha * (p$alpha + 1) / p$beta^2
    log1p((1 - p$delta * y)^2) - log(Z) + stats::dgamma(y, p$alpha, p$beta, log = TRUE) },
  references = "Vila, Ferreira, Saulo, Prataviera and Ortega (2020), A bimodal gamma distribution: Properties, regression model and applications. Statistics, 54, 469-493.",
  reference = ref_impl("definition", function(y, p) { Z <- vapply(seq_along(y), function(i) stats::integrate(function(t)
    (1 + (1 - p$delta[i] * t)^2) * stats::dgamma(t, p$alpha[i], p$beta[i]), 0, Inf)$value, 0)
    log((1 + (1 - p$delta * y)^2) * stats::dgamma(y, p$alpha, p$beta) / Z) }),
  status = "validated")

posfam("mittag_leffler", "Mittag-Leffler (Pillai, heavy tailed)", list(P("alpha", "unit", role = "index"), P("scale", "positive", role = "scale")),
  formula = "F(y) = 1 - E_alpha(-(y/scale)^alpha), E_alpha the Mittag-Leffler function",
  cdf = function(q, p) 1 - ml_fun(-(q / p$scale)^p$alpha, p$alpha),
  lpdf = function(y, p) log(ml_density(y / p$scale, p$alpha)) - log(p$scale),
  rng = function(n, p) { u <- stats::runif(n); v <- stats::runif(n); a <- p$alpha
    -p$scale * log(u) * (sin(a * pi) / tan(a * pi * v) - cos(a * pi))^(1 / a) },
  numeric_moments = FALSE,
  references = "Pillai (1990), On Mittag-Leffler functions and related distributions. Ann. Inst. Stat. Math., 42, 157-161.",
  reference = ref_impl("power series", function(y, p) log(vapply(seq_along(y), function(i) { k <- 1:150; z <- y[i] / p$scale[i]; a <- p$alpha[i]
    -sum((-1)^k * k * a * z^(a * k - 1) / gamma(a * k + 1)) / p$scale[i] }, 0)), note = "Series for the density, valid for moderate y."),
  status = "implemented")

# Mittag-Leffler helpers: E_a(z) for z <= 0 and the density, via the
# integral representation for 0 < a < 1 (Gorenflo, Loutchko and Luchko 2002):
# P(Y > y) = (sin(pi a)/pi) int_0^inf r^(a-1) exp(-r y) / (r^(2a) + 2 r^a cos(pi a) + 1) dr.
ml_fun <- function(z, a) {
  r <- recycle_pars(list(z = z, a = a))
  vapply(seq_along(r$z), function(i) { x <- (-r$z[i])^(1 / r$a[i]); al <- r$a[i]
    if (x == 0) return(1)
    stats::integrate(function(s) sin(pi * al) / pi * s^(al - 1) * exp(-s * x) / (s^(2 * al) + 2 * s^al * cos(pi * al) + 1),
                     0, Inf, rel.tol = 1e-10)$value }, 0)
}
ml_density <- function(y, a) {
  r <- recycle_pars(list(y = y, a = a))
  vapply(seq_along(r$y), function(i) { al <- r$a[i]; x <- r$y[i]
    stats::integrate(function(s) sin(pi * al) / pi * s^al * exp(-s * x) / (s^(2 * al) + 2 * s^al * cos(pi * al) + 1),
                     0, Inf, rel.tol = 1e-10)$value }, 0)
}
