# Continuous distributions on bounded intervals.

bnd <- "bounded"
unit_support <- function(p) list(lower = 0, upper = 1)
bfam <- function(id, display, pars, ..., support = unit_support)
  register_family(family_def(id, display, bnd, pars, support = support, ...))
unit_mean <- function(name = "mu") P(name, "unit", role = "location", is_mean = TRUE)
ab_support <- function(p) list(lower = p$a, upper = p$b)
ab_pars <- function() list(P("a", "real", role = "lower", link = "identity", data = TRUE),
                           P("b", "real", role = "upper", link = "identity", data = TRUE))

# Beta family -----------------------------------------------------------------------
bfam("beta", "Beta (mean and precision, Ferrari and Cribari-Neto)", list(unit_mean(), P("phi", "positive", role = "precision")),
  aliases = c("beta_proportion", "beta_regression"),
  formula = "Beta(mu phi, (1 - mu) phi)", mean_formula = "mu", variance_formula = "mu (1 - mu) / (1 + phi)",
  lpdf = function(y, p) stats::dbeta(y, p$mu * p$phi, (1 - p$mu) * p$phi, log = TRUE),
  cdf = function(q, p) stats::pbeta(q, p$mu * p$phi, (1 - p$mu) * p$phi),
  quantile = function(pr, p) stats::qbeta(pr, p$mu * p$phi, (1 - p$mu) * p$phi),
  rng = function(n, p) stats::rbeta(n, p$mu * p$phi, (1 - p$mu) * p$phi),
  mean = function(p) p$mu, variance = function(p) p$mu * (1 - p$mu) / (1 + p$phi),
  references = "Ferrari and Cribari-Neto (2004), Beta regression for modelling rates and proportions. J. Appl. Stat., 31, 799-815.",
  reference = ref_impl("gamlss.dist", function(y, p) gamlss.dist::dBE(y, p$mu, 1 / sqrt(p$phi + 1), log = TRUE)),
  stan = stan_native("beta_proportion_lpdf({y} | {mu}, {phi})", "beta_proportion_lcdf({y} | {mu}, {phi})",
                     "beta_proportion_lccdf({y} | {mu}, {phi})", "beta_proportion_rng({mu}, {phi})"),
  bugs = list(native = "dbeta({mu} * {phi}, (1 - {mu}) * {phi})"),
  start = function(y) list(mu = mean(y), phi = max(mean(y) * (1 - mean(y)) / stats::var(y) - 1, 1)), status = "validated")

bfam("beta_sigma", "Beta (GAMLSS BE: mean and sigma)", list(unit_mean(), P("sigma", "unit", role = "dispersion")), aliases = "BE",
  formula = "phi = (1 - sigma^2) / sigma^2", variance_formula = "sigma^2 mu (1 - mu)",
  lpdf = function(y, p) { f <- (1 - p$sigma^2) / p$sigma^2; stats::dbeta(y, p$mu * f, (1 - p$mu) * f, log = TRUE) },
  cdf = function(q, p) { f <- (1 - p$sigma^2) / p$sigma^2; stats::pbeta(q, p$mu * f, (1 - p$mu) * f) },
  quantile = function(pr, p) { f <- (1 - p$sigma^2) / p$sigma^2; stats::qbeta(pr, p$mu * f, (1 - p$mu) * f) },
  mean = function(p) p$mu, variance = function(p) p$sigma^2 * p$mu * (1 - p$mu),
  references = GS, reference = ref_impl("gamlss.dist", gd("dBE"), gp("pBE")),
  stan = stan_native("beta_proportion_lpdf({y} | {mu}, (1 - square({sigma})) / square({sigma}))"), status = "validated")

bfam("beta_shape", "Beta (two shapes)", list(shp("alpha"), shp("beta")),
  lpdf = function(y, p) stats::dbeta(y, p$alpha, p$beta, log = TRUE), cdf = function(q, p) stats::pbeta(q, p$alpha, p$beta),
  quantile = function(pr, p) stats::qbeta(pr, p$alpha, p$beta), rng = function(n, p) stats::rbeta(n, p$alpha, p$beta),
  mean = function(p) p$alpha / (p$alpha + p$beta),
  references = JKB, reference = ref_impl("definition", function(y, p) (p$alpha - 1) * log(y) + (p$beta - 1) * log1p(-y) - lbeta(p$alpha, p$beta)),
  stan = stan_native("beta_lpdf({y} | {alpha}, {beta})", "beta_lcdf({y} | {alpha}, {beta})", "beta_lccdf({y} | {alpha}, {beta})", "beta_rng({alpha}, {beta})"),
  bugs = list(native = "dbeta({alpha}, {beta})"), status = "validated")

bfam("beta4", "Four-parameter beta", c(list(shp("alpha"), shp("beta")), ab_pars()), support = ab_support,
  lpdf = function(y, p) stats::dbeta((y - p$a) / (p$b - p$a), p$alpha, p$beta, log = TRUE) - log(p$b - p$a),
  cdf = function(q, p) stats::pbeta((q - p$a) / (p$b - p$a), p$alpha, p$beta),
  quantile = function(pr, p) p$a + (p$b - p$a) * stats::qbeta(pr, p$alpha, p$beta),
  mean = function(p) p$a + (p$b - p$a) * p$alpha / (p$alpha + p$beta),
  references = JKB, reference = ref_impl("ExtDist", function(y, p) stats::dbeta((y - p$a) / (p$b - p$a), p$alpha, p$beta, log = TRUE) - log(p$b - p$a), note = "Location scale of the beta."),
  status = "validated")

bfam("pert", "PERT (beta PERT)", list(P("a", "real", role = "minimum", link = "identity", data = TRUE),
  P("m", "real", role = "mode", link = "identity"), P("b", "real", role = "maximum", link = "identity", data = TRUE)),
  support = ab_support, formula = "Beta4 with alpha = 1 + 4 (m - a)/(b - a), beta = 1 + 4 (b - m)/(b - a)",
  lpdf = function(y, p) { al <- 1 + 4 * (p$m - p$a) / (p$b - p$a); be <- 1 + 4 * (p$b - p$m) / (p$b - p$a)
    stats::dbeta((y - p$a) / (p$b - p$a), al, be, log = TRUE) - log(p$b - p$a) },
  cdf = function(q, p) stats::pbeta((q - p$a) / (p$b - p$a), 1 + 4 * (p$m - p$a) / (p$b - p$a), 1 + 4 * (p$b - p$m) / (p$b - p$a)),
  mean = function(p) (p$a + 4 * p$m + p$b) / 6,
  check = function(p) if (any(p$m < p$a | p$m > p$b, na.rm = TRUE)) "The mode m must lie in [a, b].",
  references = "Malcolm, Roseboom, Clark and Fazar (1959), Operations Research, 7, 646-669.",
  reference = ref_impl("beta4", function(y, p) { al <- 1 + 4 * (p$m - p$a) / (p$b - p$a); be <- 1 + 4 * (p$b - p$m) / (p$b - p$a)
    (al - 1) * log(y - p$a) + (be - 1) * log(p$b - y) - lbeta(al, be) - (al + be - 1) * log(p$b - p$a) }), status = "validated")

bfam("kumaraswamy", "Kumaraswamy", list(shp("a"), shp("b")), aliases = "kw",
  formula = "F(y) = 1 - (1 - y^a)^b",
  lpdf = function(y, p) log(p$a) + log(p$b) + (p$a - 1) * log(y) + (p$b - 1) * log1p(-y^p$a),
  cdf = function(q, p) 1 - (1 - q^p$a)^p$b, quantile = function(pr, p) (1 - (1 - pr)^(1 / p$b))^(1 / p$a),
  mean = function(p) p$b * beta(1 + 1 / p$a, p$b),
  references = "Kumaraswamy (1980), J. Hydrology, 46, 79-88.",
  reference = ref_impl("extraDistr", function(y, p) extraDistr::dkumar(y, p$a, p$b, log = TRUE)),
  stan = sfun("kumaraswamy", c("a", "b"), "    return log(a) + log(b) + (a - 1) * log(y) + (b - 1) * log1m(pow(y, a));",
              lccdf = "    return b * log1m(pow(y, a));"),
  status = "validated")

bfam("kumaraswamy_median", "Kumaraswamy (median parameterization, Mitnik and Baek)", list(P("mu", "unit", role = "median"), shp("phi")),
  formula = "a = phi, b = log(1/2) / log(1 - mu^phi)", regression_notes = "mu is the median.",
  lpdf = function(y, p) { a <- p$phi; b <- log(0.5) / log1p(-p$mu^a); log(a) + log(b) + (a - 1) * log(y) + (b - 1) * log1p(-y^a) },
  cdf = function(q, p) { a <- p$phi; b <- log(0.5) / log1p(-p$mu^a); 1 - (1 - q^a)^b },
  quantile = function(pr, p) { a <- p$phi; b <- log(0.5) / log1p(-p$mu^a); (1 - (1 - pr)^(1 / b))^(1 / a) },
  references = "Mitnik and Baek (2013), The Kumaraswamy distribution: median-dispersion re-parameterizations for regression modeling. Statistical Papers, 54, 177-192.",
  reference = ref_impl("extraDistr", function(y, p) extraDistr::dkumar(y, p$phi, log(0.5) / log1p(-p$mu^p$phi), log = TRUE)),
  start = function(y) list(mu = stats::median(y), phi = 2), status = "validated")

bfam("uniform", "Continuous uniform", ab_pars(), aliases = c("unif", "rectangular"), support = ab_support,
  lpdf = function(y, p) ifelse(y >= p$a & y <= p$b, -log(p$b - p$a), -Inf), cdf = function(q, p) stats::punif(q, p$a, p$b),
  quantile = function(pr, p) stats::qunif(pr, p$a, p$b), rng = function(n, p) stats::runif(n, p$a, p$b),
  mean = function(p) (p$a + p$b) / 2, variance = function(p) (p$b - p$a)^2 / 12,
  references = JKB, reference = ref_impl("stats", function(y, p) stats::dunif(y, p$a, p$b, log = TRUE)),
  stan = stan_native("uniform_lpdf({y} | {a}, {b})", "uniform_lcdf({y} | {a}, {b})", "uniform_lccdf({y} | {a}, {b})"), status = "validated")

bfam("arcsine", "Arcsine", ab_pars(), support = ab_support,
  lpdf = function(y, p) -log(pi) - 0.5 * log((y - p$a) * (p$b - y)),
  cdf = function(q, p) 2 / pi * asin(sqrt((q - p$a) / (p$b - p$a))),
  quantile = function(pr, p) p$a + (p$b - p$a) * sin(pi * pr / 2)^2,
  mean = function(p) (p$a + p$b) / 2, variance = function(p) (p$b - p$a)^2 / 8,
  references = JKB, reference = ref_impl("beta(1/2, 1/2)", function(y, p) stats::dbeta((y - p$a) / (p$b - p$a), 0.5, 0.5, log = TRUE) - log(p$b - p$a)),
  status = "validated")

bfam("triangular", "Triangular", list(P("a", "real", role = "minimum", link = "identity", data = TRUE),
  P("c", "real", role = "mode", link = "identity"), P("b", "real", role = "maximum", link = "identity", data = TRUE)),
  support = ab_support,
  lpdf = function(y, p) log(ifelse(y < p$a | y > p$b, 0, ifelse(y < p$c, 2 * (y - p$a) / ((p$b - p$a) * (p$c - p$a)), 2 * (p$b - y) / ((p$b - p$a) * (p$b - p$c))))),
  cdf = function(q, p) ifelse(q <= p$a, 0, ifelse(q >= p$b, 1, ifelse(q < p$c, (q - p$a)^2 / ((p$b - p$a) * (p$c - p$a)), 1 - (p$b - q)^2 / ((p$b - p$a) * (p$b - p$c))))),
  mean = function(p) (p$a + p$b + p$c) / 3,
  references = JKB, reference = ref_impl("extraDistr", function(y, p) extraDistr::dtriang(y, p$a, p$b, p$c, log = TRUE)), status = "validated")

bfam("trapezoidal", "Trapezoidal", list(P("a", "real", link = "identity", data = TRUE), P("b", "real", link = "identity"),
  P("c", "real", link = "identity"), P("d", "real", link = "identity", data = TRUE)),
  support = function(p) list(lower = p$a, upper = p$d),
  lpdf = function(y, p) { h <- 2 / (p$d + p$c - p$a - p$b)
    log(ifelse(y < p$a | y > p$d, 0, ifelse(y < p$b, h * (y - p$a) / (p$b - p$a), ifelse(y <= p$c, h, h * (p$d - y) / (p$d - p$c))))) },
  references = "Kacker and Lawrence (2007), Metrologia, 44, 117-127.",
  reference = ref_impl("trapezoid", function(y, p) log(trapezoid::dtrapezoid(y, p$a, p$b, p$c, p$d))), status = "validated")

# Irwin-Hall density and CDF by the alternating sum, evaluated in the lower
# half of the support (symmetry f(x) = f(n - x), F(x) = 1 - F(n - x)) where
# the cancellation is mild.
# ponytail: still loses digits for n above about 40; switch to log scale
# compensated sums if larger n is needed.
ih_cdf <- function(x, n) vapply(seq_along(x), function(i) { if (x[i] <= 0) return(0); if (x[i] >= n[i]) return(1)
  up <- x[i] > n[i] / 2; z <- if (up) n[i] - x[i] else x[i]; k <- 0:floor(z)
  v <- sum((-1)^k * choose(n[i], k) * (z - k)^n[i]) / factorial(n[i]); if (up) 1 - v else v }, 0)
ih_dens <- function(x, n) vapply(seq_along(x), function(i) { if (x[i] <= 0 || x[i] >= n[i]) return(0)
  z <- min(x[i], n[i] - x[i]); k <- 0:floor(z)
  max(0, sum((-1)^k * choose(n[i], k) * (z - k)^(n[i] - 1)) / factorial(n[i] - 1)) }, 0)
bfam("irwin_hall", "Irwin-Hall (sum of uniforms)", list(P("n", "positive_integer", role = "count", data = TRUE)),
  support = function(p) list(lower = 0, upper = p$n),
  lpdf = function(y, p) log(ih_dens(y, rep_len(p$n, length(y)))),
  cdf = function(q, p) ih_cdf(q, rep_len(p$n, length(q))),
  rng = function(n, p) vapply(seq_len(n), function(i) sum(stats::runif(rep_len(p$n, n)[i])), 0),
  mean = function(p) p$n / 2, variance = function(p) p$n / 12,
  references = "Irwin (1927); Hall (1927), Biometrika, 19.",
  reference = ref_impl("convolution", function(y, p) log(vapply(seq_along(y), function(i) { n <- p$n[i]
    if (n == 1) return(1); stats::integrate(function(t) vapply(t, function(tt) { k <- 0:floor(tt); if (tt <= 0 || tt >= n - 1) return(0)
      sum((-1)^k * choose(n - 1, k) * (tt - k)^(n - 2)) / factorial(n - 2) }, 0), max(0, y[i] - 1), min(n - 1, y[i]))$value }, 0))),
  status = "validated")

bfam("bates", "Bates (mean of uniforms)", list(P("n", "positive_integer", role = "count", data = TRUE)),
  lpdf = function(y, p) { n <- rep_len(p$n, length(y)); log(n * ih_dens(n * y, n)) },
  cdf = function(q, p) ih_cdf(q * p$n, rep_len(p$n, length(q))),
  mean = function(p) 0.5 + 0 * p$n, variance = function(p) 1 / (12 * p$n),
  references = "Bates (1955), Ann. Math. Stat., 26, 705-720.",
  reference = ref_impl("Irwin-Hall", function(y, p) log(p$n) + log(vapply(seq_along(y), function(i) { n <- p$n[i]; x <- n * y[i]; k <- 0:floor(x)
    sum((-1)^k * choose(n, k) * (x - k)^(n - 1)) / factorial(n - 1) }, 0))), status = "validated")

bfam("logit_normal", "Logit-normal", list(loc_par("mu"), scale_par("sigma")), aliases = "LOGITNO",
  lpdf = function(y, p) stats::dnorm(stats::qlogis(y), p$mu, p$sigma, log = TRUE) - log(y) - log1p(-y),
  cdf = function(q, p) stats::pnorm(stats::qlogis(q), p$mu, p$sigma), quantile = function(pr, p) stats::plogis(stats::qnorm(pr, p$mu, p$sigma)),
  references = "Atchison and Shen (1980), Logistic-normal distributions. Biometrika, 67, 261-272.",
  reference = ref_impl("logitnorm", function(y, p) logitnorm::dlogitnorm(y, p$mu, p$sigma, log = TRUE)),
  stan = sfun("logit_normal", c("mu", "sigma"), "    return normal_lpdf(logit(y) | mu, sigma) - log(y) - log1m(y);"), status = "validated")

bfam("u_quadratic", "U-quadratic", ab_pars(), support = ab_support,
  lpdf = function(y, p) { al <- 12 / (p$b - p$a)^3; be <- (p$a + p$b) / 2; log(al) + 2 * log(abs(y - be)) },
  cdf = function(q, p) { al <- 12 / (p$b - p$a)^3; be <- (p$a + p$b) / 2; al / 3 * ((q - be)^3 + (be - p$a)^3) },
  references = "Wikipedia, U-quadratic distribution.",
  reference = ref_impl("definition", function(y, p) log(12 / (p$b - p$a)^3 * (y - (p$a + p$b) / 2)^2)), status = "validated")

bfam("wigner_semicircle", "Wigner semicircle", list(P("R", "positive", role = "radius")), support = function(p) list(lower = -p$R, upper = p$R),
  lpdf = function(y, p) log(2 / (pi * p$R^2)) + 0.5 * log(pmax(p$R^2 - y^2, 0)),
  cdf = function(q, p) { x <- pmin(pmax(q / p$R, -1), 1); 0.5 + (x * sqrt(1 - x^2) + asin(x)) / pi },
  mean = function(p) 0 * p$R, variance = function(p) p$R^2 / 4,
  references = "Wigner (1955), Ann. Math., 62, 548-564.",
  reference = ref_impl("scaled beta(3/2, 3/2)", function(y, p) stats::dbeta((y + p$R) / (2 * p$R), 1.5, 1.5, log = TRUE) - log(2 * p$R)), status = "validated")

bfam("raised_cosine", "Raised cosine", list(loc_par("mu"), P("s", "positive", role = "half width")),
  support = function(p) list(lower = p$mu - p$s, upper = p$mu + p$s),
  lpdf = function(y, p) log((1 + cos(pi * (y - p$mu) / p$s)) / (2 * p$s)),
  cdf = function(q, p) { z <- pmin(pmax((q - p$mu) / p$s, -1), 1); 0.5 * (1 + z + sin(pi * z) / pi) },
  references = "Wikipedia, Raised cosine distribution.",
  reference = ref_impl("definition", function(y, p) log(1 / (2 * p$s) * (1 + cos((y - p$mu) / p$s * pi)))), status = "validated")

bfam("reciprocal", "Reciprocal (log-uniform)", list(P("a", "positive", role = "lower", data = TRUE), P("b", "positive", role = "upper", data = TRUE)),
  aliases = "log_uniform", support = ab_support,
  lpdf = function(y, p) -log(y) - log(log(p$b / p$a)), cdf = function(q, p) log(q / p$a) / log(p$b / p$a),
  quantile = function(pr, p) p$a * (p$b / p$a)^pr,
  references = "Hamming (1970), On the distribution of numbers. Bell System Technical Journal, 49, 1609-1625.",
  reference = ref_impl("transform", function(y, p) stats::dunif(log(y), log(p$a), log(p$b), log = TRUE) - log(y)), status = "validated")

cb_logC <- function(lam) ifelse(abs(lam - 0.5) < 1e-6, log(2), log(2 * atanh(1 - 2 * lam) / (1 - 2 * lam)))
bfam("continuous_bernoulli", "Continuous Bernoulli", list(P("lambda", "unit", role = "shape")),
  formula = "f(y) = C(lambda) lambda^y (1 - lambda)^(1 - y)",
  lpdf = function(y, p) cb_logC(p$lambda) + y * log(p$lambda) + (1 - y) * log1p(-p$lambda),
  cdf = function(q, p) ifelse(abs(p$lambda - 0.5) < 1e-6, q, (p$lambda^q * (1 - p$lambda)^(1 - q) + p$lambda - 1) / (2 * p$lambda - 1)),
  references = "Loaiza-Ganem and Cunningham (2019), The continuous Bernoulli. NeurIPS.",
  reference = ref_impl("definition", function(y, p) log(2 * atanh(1 - 2 * p$lambda) / (1 - 2 * p$lambda) * p$lambda^y * (1 - p$lambda)^(1 - y))),
  status = "validated")

# Cobin: exponential dispersion family on [0, 1] with B(theta) = log((e^theta - 1)/theta)
cobin_B <- function(th) ifelse(abs(th) < 1e-8, th / 2, log(abs(expm1(th) / th)))
bfam("cobin", "Continuous binomial (cobin, Lee, Kim and Ji)", list(P("theta", "real", role = "natural", link = "identity"),
  P("lambda", "positive_integer", role = "dispersion", data = TRUE, interpretation = "Integer dispersion; mean of lambda continuous Bernoullis")),
  formula = "f(y) = h(y, lambda) exp(lambda (theta y - B(theta))), B(theta) = log((e^theta - 1)/theta), h from the Irwin-Hall density",
  mean_formula = "B'(theta) = 1/(1 - e^-theta) - 1/theta",
  lpdf = function(y, p) log(p$lambda) + log(ih_dens(p$lambda * y, rep_len(p$lambda, length(y)))) + p$lambda * (p$theta * y - cobin_B(p$theta)),
  mean = function(p) ifelse(abs(p$theta) < 1e-8, 0.5, 1 / (1 - exp(-p$theta)) - 1 / p$theta),
  regression_notes = "Canonical link: theta = x'beta, i.e. mu = B'(x'beta) (the cobit link).",
  references = "Lee, Kim and Ji (2025), Scalable and robust regression models for continuous proportional data. arXiv 2504.15269.",
  reference = ref_impl("definition (lambda = 1 continuous Bernoulli)", function(y, p) log(p$lambda) + log(ih_dens(p$lambda * y, rep_len(p$lambda, length(y)))) + p$lambda * (p$theta * y - log((exp(p$theta) - 1) / p$theta))),
  status = "validated")

# Unit distributions for proportions -----------------------------------------------
bfam("simplex", "Simplex (Barndorff-Nielsen and Jorgensen)", list(unit_mean(), scale_par("sigma")), aliases = "SIMPLEX",
  formula = "f(y) = [2 pi sigma^2 (y(1 - y))^3]^(-1/2) exp(-d(y; mu) / (2 sigma^2)), d = (y - mu)^2 / (y (1 - y) mu^2 (1 - mu)^2)",
  lpdf = function(y, p) -0.5 * log(2 * pi * p$sigma^2 * (y * (1 - y))^3) - (y - p$mu)^2 / (y * (1 - y) * p$mu^2 * (1 - p$mu)^2) / (2 * p$sigma^2),
  mean = function(p) p$mu,
  references = c("Barndorff-Nielsen and Jorgensen (1991), Some parametric models on the simplex. J. Multivariate Anal., 39, 106-116.", GS),
  reference = ref_impl("gamlss.dist", gd("dSIMPLEX")),
  stan = sfun("simplex", c("mu", "sigma"),
    "    return -0.5 * log(2 * pi() * square(sigma) * pow(y * (1 - y), 3)) - square(y - mu) / (y * (1 - y) * square(mu) * square(1 - mu)) / (2 * square(sigma));"),
  status = "validated")

bfam("unit_lindley", "Unit-Lindley (Mazucheli et al.)", list(unit_mean()),
  formula = "f(y) = (1 - mu)^2 / (mu (1 - y)^3) exp(-y (1 - mu) / (mu (1 - y))), mean mu",
  lpdf = function(y, p) 2 * log1p(-p$mu) - log(p$mu) - 3 * log1p(-y) - y * (1 - p$mu) / (p$mu * (1 - y)),
  cdf = function(q, p) { a <- (1 - p$mu) / p$mu; x <- pmin(pmax(q, 0), 1 - 1e-12)
    ifelse(q >= 1, 1, 1 - (1 - x * a / ((1 + a) * (x - 1))) * exp(-a * x / (1 - x))) },
  mean = function(p) p$mu,
  references = "Mazucheli, Menezes and Chakraborty (2019), On the one parameter unit-Lindley distribution and its associated regression model for proportion data. J. Appl. Stat., 46, 700-714.",
  reference = ref_impl("transform of Lindley", function(y, p) { th <- (1 - p$mu) / p$mu; x <- y / (1 - y)
    2 * log(th) - log1p(th) + log1p(x) - th * x - 2 * log1p(-y) }),
  stan = sfun("unit_lindley", "mu", "    return 2 * log1m(mu) - log(mu) - 3 * log1m(y) - y * (1 - mu) / (mu * (1 - y));"),
  start = function(y) list(mu = mean(y)), status = "validated")

bfam("unit_weibull", "Unit-Weibull (quantile parameterization)", list(P("mu", "unit", role = "quantile"), shp("beta"),
  P("tau", "unit", role = "quantile level", regression = FALSE, default = 0.5)),
  formula = "F(y) = exp(log(tau) (log y / log mu)^beta)",
  lpdf = function(y, p) { r <- log(y) / log(p$mu); log(p$beta) - log(y) - log(-log(p$mu)) + log(-log(p$tau)) + (p$beta - 1) * log(r) + log(p$tau) * r^p$beta },
  cdf = function(q, p) exp(log(p$tau) * (log(q) / log(p$mu))^p$beta),
  quantile = function(pr, p) exp(log(p$mu) * (log(pr) / log(p$tau))^(1 / p$beta)),
  references = "Mazucheli, Menezes, Fernandes, de Oliveira and Ghitany (2020), The unit-Weibull distribution as an alternative to the Kumaraswamy distribution for the modeling of quantiles conditional on covariates. J. Appl. Stat., 47, 954-974.",
  reference = ref_impl("transform of Weibull", function(y, p) { a <- -log(p$tau) / (-log(p$mu))^p$beta; x <- -log(y)
    log(a * p$beta) + (p$beta - 1) * log(x) - a * x^p$beta - log(y) }), status = "validated")

bfam("unit_gamma", "Unit-gamma (mean parameterization)", list(unit_mean(), shp("phi")),
  formula = "Y = exp(-X), X ~ Gamma(phi, rate theta), theta = mu^(1/phi) / (1 - mu^(1/phi))", mean_formula = "mu",
  lpdf = function(y, p) { th <- p$mu^(1 / p$phi) / (1 - p$mu^(1 / p$phi)); stats::dgamma(-log(y), p$phi, th, log = TRUE) - log(y) },
  cdf = function(q, p) { th <- p$mu^(1 / p$phi) / (1 - p$mu^(1 / p$phi)); stats::pgamma(-log(q), p$phi, th, lower.tail = FALSE) },
  mean = function(p) p$mu,
  references = "Mousa, El-Sheikh and Abdel-Fattah (2016), A gamma regression for bounded continuous variables. Advances and Applications in Statistics, 49, 305-326.",
  reference = ref_impl("definition", function(y, p) { th <- p$mu^(1 / p$phi) / (1 - p$mu^(1 / p$phi))
    p$phi * log(th) - lgamma(p$phi) + (th - 1) * log(y) + (p$phi - 1) * log(-log(y)) }), status = "validated")

bfam("topp_leone", "Topp-Leone", list(shp("nu")), formula = "F(y) = (y (2 - y))^nu",
  lpdf = function(y, p) log(2 * p$nu) + (p$nu - 1) * log(y * (2 - y)) + log1p(-y),
  cdf = function(q, p) (q * (2 - q))^p$nu, quantile = function(pr, p) 1 - sqrt(1 - pr^(1 / p$nu)),
  references = "Topp and Leone (1955), JASA, 50, 209-219.",
  reference = ref_impl("definition", function(y, p) log(2 * p$nu * (1 - y) * (y * (2 - y))^(p$nu - 1))), status = "validated")

bfam("johnson_sb", "Johnson SB", list(P("gamma", "real", link = "identity"), shp("delta"), loc_par("xi"), scale_par("lambda")),
  support = function(p) list(lower = p$xi, upper = p$xi + p$lambda),
  lpdf = function(y, p) { z <- (y - p$xi) / p$lambda
    log(p$delta) - log(p$lambda) - log(z * (1 - z)) + stats::dnorm(p$gamma + p$delta * log(z / (1 - z)), log = TRUE) },
  cdf = function(q, p) { z <- (q - p$xi) / p$lambda; stats::pnorm(p$gamma + p$delta * stats::qlogis(z)) },
  quantile = function(pr, p) p$xi + p$lambda * stats::plogis((stats::qnorm(pr) - p$gamma) / p$delta),
  references = "Johnson (1949), Systems of frequency curves generated by methods of translation. Biometrika, 36, 149-176.",
  reference = ref_impl("SuppDists", function(y, p) log(SuppDists::dJohnson(y, list(gamma = p$gamma[1], delta = p$delta[1], xi = p$xi[1], lambda = p$lambda[1], type = "SB")))),
  status = "validated")

bfam("marchenko_pastur", "Marchenko-Pastur (ratio lambda <= 1)", list(P("lambda", "unit", role = "ratio"), P("s2", "positive", role = "scale")),
  support = function(p) list(lower = p$s2 * (1 - sqrt(p$lambda))^2, upper = p$s2 * (1 + sqrt(p$lambda))^2),
  lpdf = function(y, p) { lo <- p$s2 * (1 - sqrt(p$lambda))^2; hi <- p$s2 * (1 + sqrt(p$lambda))^2
    log(sqrt(pmax((hi - y) * (y - lo), 0)) / (2 * pi * p$s2 * p$lambda * y)) },
  mean = function(p) p$s2,
  references = "Marchenko and Pastur (1967), Mat. Sbornik, 72, 507-536.",
  reference = ref_impl("RMTstat", function(y, p) RMTstat::dmp(y, svr = 1 / p$lambda, var = p$s2, log = TRUE)), status = "validated")

bfam("ordered_beta", "Ordered beta (Kubinec)", list(unit_mean(), P("phi", "positive", role = "precision"),
  P("k1", "real", role = "lower cutpoint", link = "identity"), P("k2", "real", role = "upper cutpoint", link = "identity")),
  support = function(p) list(lower = 0, upper = 1),
  formula = "P(0) = 1 - plogis(eta - k1), P(1) = plogis(eta - k2), (0,1): (plogis(eta - k1) - plogis(eta - k2)) Beta(mu phi, (1 - mu) phi), eta = logit(mu)",
  lpdf = function(y, p) { e <- stats::qlogis(p$mu)
    ifelse(y == 0, stats::plogis(e - p$k1, lower.tail = FALSE, log.p = TRUE),
    ifelse(y == 1, stats::plogis(e - p$k2, log.p = TRUE),
      log(stats::plogis(e - p$k1) - stats::plogis(e - p$k2)) + stats::dbeta(y, p$mu * p$phi, (1 - p$mu) * p$phi, log = TRUE))) },
  extra = list(atom_points = c(0, 1)),
  atom = function(y, p) { e <- stats::qlogis(p$mu)
    ifelse(y == 0, stats::plogis(e - p$k1, lower.tail = FALSE), ifelse(y == 1, stats::plogis(e - p$k2), 0)) },
  numeric_cdf = FALSE, numeric_moments = FALSE,
  references = "Kubinec (2023), Ordered beta regression: A parsimonious, well-fitting model for continuous data with lower and upper bounds. Political Analysis, 31, 519-536.",
  reference = ref_impl("definition", function(y, p) { e <- stats::qlogis(p$mu)
    log(stats::plogis(e - p$k1) - stats::plogis(e - p$k2)) + stats::dbeta(y, p$mu * p$phi, (1 - p$mu) * p$phi, log = TRUE) }),
  stan = sfun("ordbeta", c("mu", "phi", "k1", "k2"),
    "    real e = logit(mu);\n    if (y == 0) return log1m_inv_logit(e - k1);\n    if (y == 1) return log_inv_logit(e - k2);\n    return log_diff_exp(log_inv_logit(e - k1), log_inv_logit(e - k2)) + beta_proportion_lpdf(y | mu, phi);"),
  start = function(y) list(mu = mean(y), phi = 2, k1 = -3, k2 = 3), status = "implemented")

bfam("libby_novick", "Libby-Novick generalized beta", list(shp("alpha"), shp("beta"), shp("lambda")),
  formula = "f(y) = lambda^alpha y^(alpha - 1) (1 - y)^(beta - 1) / (B(alpha, beta) (1 - (1 - lambda) y)^(alpha + beta))",
  lpdf = function(y, p) p$alpha * log(p$lambda) + (p$alpha - 1) * log(y) + (p$beta - 1) * log1p(-y) - lbeta(p$alpha, p$beta) -
    (p$alpha + p$beta) * log1p(-(1 - p$lambda) * y),
  cdf = function(q, p) stats::pbeta(p$lambda * q / (1 - (1 - p$lambda) * q), p$alpha, p$beta),
  references = "Libby and Novick (1982), Multivariate generalized beta distributions with applications to utility assessment. J. Educational Statistics, 7, 271-294.",
  reference = ref_impl("transform of beta", function(y, p) { u <- p$lambda * y / (1 - (1 - p$lambda) * y)
    stats::dbeta(u, p$alpha, p$beta, log = TRUE) + log(p$lambda) - 2 * log(1 - (1 - p$lambda) * y) }), status = "validated")

bfam("generalized_beta1", "Generalized beta of the first kind (GB1)", list(P("scale", "positive", role = "scale", data = TRUE),
  P("a", "real", role = "shape", link = "identity"), shp("p"), shp("q")),
  support = function(r) list(lower = ifelse(r$a > 0, 0, r$scale), upper = ifelse(r$a > 0, r$scale, Inf)),
  check = function(r) if (any(r$a == 0)) "GB1 needs a != 0.",
  lpdf = function(y, r) log(abs(r$a)) + (r$a * r$p - 1) * log(y) + (r$q - 1) * log1p(-(y / r$scale)^r$a) - r$a * r$p * log(r$scale) - lbeta(r$p, r$q),
  cdf = function(x, r) { v <- stats::pbeta((x / r$scale)^r$a, r$p, r$q); ifelse(r$a > 0, v, 1 - v) },
  references = "McDonald (1984), Econometrica, 52, 647-663.",
  reference = ref_impl("transform of beta", function(y, r) stats::dbeta((y / r$scale)^r$a, r$p, r$q, log = TRUE) + log(abs(r$a)) + (r$a - 1) * log(y) - r$a * log(r$scale)),
  status = "validated")

bfam("log_lindley", "Log-Lindley (Gomez-Deniz et al.)", list(shp("sigma"), P("lambda", "nonneg", role = "shape")),
  formula = "f(y) = sigma^2 / (1 + lambda sigma) (lambda - log y) y^(sigma - 1)",
  lpdf = function(y, p) 2 * log(p$sigma) - log1p(p$lambda * p$sigma) + log(p$lambda - log(y)) + (p$sigma - 1) * log(y),
  cdf = function(q, p) { x <- pmax(q, 1e-300); ifelse(q <= 0, 0, pmin(x^p$sigma * (1 + p$lambda * p$sigma - p$sigma * log(x)) / (1 + p$lambda * p$sigma), 1)) },
  references = "Gomez-Deniz, Sordo and Calderin-Ojeda (2014), The log-Lindley distribution as an alternative to the beta regression model with applications in insurance. Insurance: Mathematics and Economics, 54, 49-57.",
  reference = ref_impl("definition", function(y, p) log(p$sigma^2 / (1 + p$lambda * p$sigma) * (p$lambda - log(y)) * y^(p$sigma - 1))), status = "validated")

bfam("unit_modified_weibull", "Unit modified Weibull", list(shp("alpha"), P("gamma", "positive"), P("lambda", "nonneg")),
  formula = "F(y) = exp(-alpha (-log y)^gamma y^(-lambda))",
  lpdf = function(y, p) { L <- -log(y); log(p$alpha) + p$gamma * log(L) - (p$lambda + 1) * log(y) + log(p$lambda + p$gamma / L) - p$alpha * L^p$gamma * y^(-p$lambda) },
  cdf = function(q, p) exp(-p$alpha * (-log(q))^p$gamma * q^(-p$lambda)),
  references = "arXiv 2508.17359 (2025), Unit-modified Weibull distribution and quantile regression model.",
  reference = ref_impl("definition", function(y, p) log(p$alpha * (-log(y))^p$gamma / (log(y) * y^(p$lambda + 1)) * (p$lambda * log(y) - p$gamma) *
    exp(-p$alpha * (-log(y))^p$gamma * y^(-p$lambda)))), status = "validated")

bfam("mbur", "Median based unit Rayleigh (MBUR)", list(P("alpha", "positive", role = "scale")),
  formula = "Y = exp(-X) with X the median of three Rayleigh(alpha) draws: F(y) = 3 y^(2/alpha^2) - 2 y^(3/alpha^2)",
  lpdf = function(y, p) { k <- 1 / p$alpha^2; log(6 * k) + (2 * k - 1) * log(y) + log1p(-y^k) },
  cdf = function(q, p) { k <- 1 / p$alpha^2; 3 * q^(2 * k) - 2 * q^(3 * k) },
  references = "Attia (2024), A new one parameter unit distribution: median based unit Rayleigh (MBUR). arXiv 2410.14857.",
  reference = ref_impl("order statistic", function(y, p) { k <- 1 / p$alpha^2; u <- y^k; log(6 * u * (1 - u) * k * y^(k - 1)) }),
  status = "validated")
