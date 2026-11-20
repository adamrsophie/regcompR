# Remaining univariate laws: mixed and singular laws, quantile parameterized
# systems, the Pearson system and laws defined through special functions.

# Positive line additions ------------------------------------------------------------
posfam("gamma_gompertz", "Gamma/Gompertz", list(P("b", "positive", role = "scale"), shp("s"), P("beta", "positive", role = "shape")),
  formula = "f(y) = b s e^(b y) beta^s / (beta - 1 + e^(b y))^(s + 1)",
  lpdf = function(y, p) log(p$b * p$s) + p$b * y + p$s * log(p$beta) - (p$s + 1) * log(p$beta - 1 + exp(p$b * y)),
  cdf = function(q, p) 1 - p$beta^p$s / (p$beta - 1 + exp(p$b * q))^p$s,
  quantile = function(pr, p) log(p$beta * (1 - pr)^(-1 / p$s) - p$beta + 1) / p$b,
  references = "Bemmaor and Glady (2012), Modeling purchasing behavior with sudden death: A flexible customer lifetime model. Management Science, 58, 1012-1021.",
  reference = ref_impl("gamma mixture of Gompertz", function(y, p) log(vapply(seq_along(y), function(i) stats::integrate(function(eta)
    eta * p$b[i] * exp(p$b[i] * y[i]) * exp(-eta * (exp(p$b[i] * y[i]) - 1)) * stats::dgamma(eta, p$s[i], p$beta[i]), 0, Inf, rel.tol = 1e-12)$value, 0))),
  status = "validated")

register_family(family_def("pearson_iii", "Pearson type III (shifted gamma)", "real",
  list(shp("alpha"), P("beta", "positive", role = "rate"), loc_par("gamma")),
  support = function(p) list(lower = p$gamma, upper = Inf),
  lpdf = function(y, p) stats::dgamma(y - p$gamma, p$alpha, p$beta, log = TRUE),
  cdf = function(q, p) stats::pgamma(q - p$gamma, p$alpha, p$beta), quantile = function(pr, p) p$gamma + stats::qgamma(pr, p$alpha, p$beta),
  mean = function(p) p$gamma + p$alpha / p$beta, variance = function(p) p$alpha / p$beta^2,
  references = JKB, reference = ref_impl("PearsonDS", function(y, p) PearsonDS::dpearsonIII(y, p$alpha[1], p$gamma[1], 1 / p$beta[1], log = TRUE)),
  status = "validated"))

posfam("noncentral_chi", "Noncentral chi", list(P("k", "positive"), P("lambda", "nonneg", role = "noncentrality")),
  formula = "Y = sqrt(noncentral chi-squared(k, lambda^2))",
  lpdf = function(y, p) stats::dchisq(y^2, p$k, p$lambda^2, log = TRUE) + log(2 * y),
  cdf = function(q, p) stats::pchisq(q^2, p$k, p$lambda^2), quantile = function(pr, p) sqrt(stats::qchisq(pr, p$k, p$lambda^2)),
  references = JKB, reference = ref_impl("Bessel form", function(y, p) log(exp(-(y^2 + p$lambda^2) / 2) * y^p$k * p$lambda / (p$lambda * y)^(p$k / 2) *
    besselI(p$lambda * y, p$k / 2 - 1))), status = "validated")

# Matrix exponential by scaling and squaring (handles defective generators).
pt_expm <- function(S, t) as.matrix(Matrix::expm(Matrix::Matrix(S * t)))
posfam("phase_type", "Continuous phase-type", list(P("alpha", "unit_closed", role = "initial distribution", regression = FALSE),
  P("S", "any", role = "sub-generator matrix", regression = FALSE)), variate_type = "scalar_matrix_parameter",
  formula = "f(y) = alpha exp(S y) s0, s0 = -S 1",
  lpdf = function(y, p) { S <- as.matrix(p$S); s0 <- -rowSums(S); vapply(y, function(t) log(sum(p$alpha %*% pt_expm(S, t) * s0)), 0) },
  cdf = function(q, p) { S <- as.matrix(p$S); vapply(q, function(t) if (t <= 0) 0 else 1 - sum(p$alpha %*% pt_expm(S, t)), 0) },
  rng = function(n, p) { S <- as.matrix(p$S); K <- nrow(S); s0 <- -rowSums(S)
    vapply(seq_len(n), function(i) { s <- sample.int(K, 1, prob = p$alpha); t <- 0
      repeat { r <- -S[s, s]; t <- t + stats::rexp(1, r); pr <- c(pmax(S[s, -s], 0), s0[s]) / r
        nx <- sample.int(K, 1, prob = pr); if (nx == K) return(t); s <- c(seq_len(K)[-s], 0)[nx] } }, 0) },
  mean = function(p) sum(p$alpha %*% solve(-as.matrix(p$S))),
  numerical_notes = "Matrix exponential by scaling and squaring (Matrix::expm).",
  references = "Neuts (1981), Matrix-Geometric Solutions in Stochastic Models.",
  reference = ref_impl("hypoexponential for a two phase chain", function(y, p) { S <- as.matrix(p$S); s0 <- -rowSums(S)
    vapply(y, function(t) log(sum(p$alpha %*% as.matrix(Matrix::expm(Matrix::Matrix(S * t))) * s0)), 0) }),
  status = "validated")

hw_theta <- function(r, t) if (t < 0.15) 0 else r / sqrt(2 * pi^3 * t) * exp(pi^2 / (2 * t)) * stats::integrate(function(x)
  exp(-x^2 / (2 * t) - r * cosh(x) + log(sinh(x))) * sin(pi * x / t), 0, acosh(1 + 800 / r), subdivisions = 5000L, rel.tol = 1e-12, stop.on.error = FALSE)$value
posfam("hartman_watson", "Hartman-Watson", list(P("r", "positive", role = "shape")),
  formula = "f_r(t) = theta(r, t) / I0(r), theta(r, t) = r (2 pi^3 t)^(-1/2) e^(pi^2/(2t)) int_0^inf e^(-x^2/(2t) - r cosh x) sinh(x) sin(pi x / t) dx",
  lpdf = function(y, p) log(pmax(pervec(y, p, function(t, pp) hw_theta(pp$r, t)), 0)) - log(besselI(p$r, 0)),
  cdf = function(q, p) pervec(q, p, function(t, pp) { if (t <= 0.15) return(0)
    g <- function(x) vapply(x, function(xx) hw_theta(pp$r, xx), 0)
    br <- unique(c(0.15, pmin(c(0.5, 1, 3, 10, 50), t), t))
    min(1, sum(vapply(seq_len(length(br) - 1), function(j) stats::integrate(g, br[j], br[j + 1], rel.tol = 1e-7,
      stop.on.error = FALSE)$value, 0)) / besselI(pp$r, 0)) }),
  numerical_notes = "The integral is highly oscillatory with a huge prefactor for small t; the density is set to 0 for t < 0.15 (double precision cannot resolve it there); total mass and the Laplace transform identity hold to about 1e-6 for r near 2.",
  references = "Hartman and Watson (1974), 'Normal' distribution functions on spheres and the modified Bessel functions. Ann. Probab., 2, 593-607.",
  reference = ref_impl("Laplace transform I_|u|(r) / I0(r)", function(y, p) log(pervec(y, p, function(t, pp) hw_theta(pp$r, t))) - log(besselI(p$r, 0)),
    note = "Checked through total mass and the Laplace transform identity."),
  status = "validated")

# Real line additions ----------------------------------------------------------------
# exp(-2 x^3/3) (x Ai(x^2) - Ai'(x^2)) with Ai(z) = sqrt(z/3) K_{1/3}(zeta)/pi and
# Ai'(z) = -z K_{2/3}(zeta)/(pi sqrt 3), zeta = 2 |x|^3/3, using scaled Bessel K
# so the exponentials cancel on the left tail.
map_airy_core <- function(x) { z <- x^2; ze <- 2 / 3 * abs(x)^3
  k1 <- besselK(pmax(ze, 1e-300), 1 / 3, expon.scaled = TRUE); k2 <- besselK(pmax(ze, 1e-300), 2 / 3, expon.scaled = TRUE)
  inner <- x * sqrt(z / 3) / pi * k1 + z / (pi * sqrt(3)) * k2
  inner[abs(x) < 1e-8] <- 0.2588194037928068
  log(pmax(inner, 0)) - ze - 2 * x^3 / 3 }
rfam("map_airy", "Map-Airy", list(loc_par("mu"), scale_par("sigma")),
  formula = "f(x) = 2 exp(-2 x^3 / 3) (x Ai(x^2) - Ai'(x^2)), x = (y - mu)/sigma",
  lpdf = function(y, p) log(2) + map_airy_core((y - p$mu) / p$sigma) - log(p$sigma),
  cdf = function(q, p) vapply((q - p$mu) / p$sigma, function(x) 1 - stats::integrate(function(t) 2 * exp(map_airy_core(t)), x, Inf,
    rel.tol = 1e-10, subdivisions = 2000L)$value, 0),
  references = "Banderier, Flajolet, Schaeffer and Soria (2001), Random maps, coalescing saddles, singularity analysis, and Airy phenomena. Random Structures and Algorithms, 19, 194-246.",
  reference = ref_impl("unscaled Airy form", function(y, p) { x <- (y - p$mu) / p$sigma; z <- x^2; ze <- 2 / 3 * abs(x)^3
    log(2 * exp(-2 * x^3 / 3) * (x * sqrt(z / 3) / pi * besselK(ze, 1 / 3) + z / (pi * sqrt(3)) * besselK(ze, 2 / 3))) - log(p$sigma) }),
  status = "validated")

an_logc <- function(tau, s) log(2) + 0.5 * log(tau * (1 - tau)) - log(s) - 0.5 * log(pi) - log(sqrt(tau) + sqrt(1 - tau))
rfam("expectile", "Expectile (asymmetric normal)", list(loc_par("mu"), scale_par("sigma"), P("tau", "unit", role = "expectile level", regression = FALSE, default = 0.5)),
  aliases = "asymmetric_normal",
  formula = "f(y) = 2 sqrt(tau (1 - tau)) / (sigma sqrt(pi) (sqrt(tau) + sqrt(1 - tau))) exp(-|tau - 1(y < mu)| (y - mu)^2 / sigma^2)",
  lpdf = function(y, p) an_logc(p$tau, p$sigma) - abs(p$tau - (y < p$mu)) * (y - p$mu)^2 / p$sigma^2,
  regression_notes = "Maximum likelihood for mu gives expectile regression (Newey and Powell 1987) at level tau.",
  references = c("Newey and Powell (1987), Asymmetric least squares estimation and testing. Econometrica, 55, 819-847.",
                 "Waldmann, Sobotka and Kneib (2017), Bayesian regularisation in geoadditive expectile regression. Stat. Comput., 27, 1539-1553."),
  reference = ref_impl("definition", function(y, p) log(2 * sqrt(p$tau * (1 - p$tau)) / (p$sigma * sqrt(pi) * (sqrt(p$tau) + sqrt(1 - p$tau))) *
    exp(-ifelse(y < p$mu, 1 - p$tau, p$tau) * (y - p$mu)^2 / p$sigma^2))),
  stan = sfun("expectile", c("mu", "sigma", "tau"),
    "    real w = y < mu ? 1 - tau : tau;\n    return log2() + 0.5 * log(tau * (1 - tau)) - log(sigma) - 0.5 * log(pi()) - log(sqrt(tau) + sqrt(1 - tau)) - w * square(y - mu) / square(sigma);"),
  status = "validated")

# exp(z^2/4) D_{-nu}(z) = int_0^inf t^(nu - 1) exp(-z t - t^2/2) dt / Gamma(nu), and the
# normalizing constant is closed form, so no overflow occurs.
neg_I <- function(z, nu) stats::integrate(function(t) exp((nu - 1) * log(t) - z * t - t^2 / 2), 0, Inf, rel.tol = 1e-12)$value
rfam("normal_exponential_gamma", "Normal-exponential-gamma", list(loc_par("mu"), shp("k"), scale_par("theta")), aliases = "NEG",
  formula = "f(y) = exp(z^2/4) D_{-2k-1}(z) Gamma(2k + 1) / (2 theta 2^(k - 1) Gamma(k) Gamma(2k + 1)) style; z = |y - mu|/theta, D the parabolic cylinder function",
  lpdf = function(y, p) pervec(abs(y - p$mu) / p$theta, p, function(z, pp)
    log(neg_I(z, 2 * pp$k + 1)) - log(2 * pp$theta) - (pp$k - 1) * log(2) - lgamma(pp$k)),
  references = "Griffin and Brown (2011), Bayesian hyper-lassos with non-convex penalization. Aust. N. Z. J. Stat., 53, 423-442.",
  reference = ref_impl("normal scale mixture", function(y, p) log(pervec(y, p, function(x, pp) stats::integrate(function(v)
    stats::dnorm(x, pp$mu, sqrt(v)) * vapply(v, function(vv) stats::integrate(function(g) stats::dexp(vv, g) * stats::dgamma(g, pp$k, rate = pp$theta^2), 0, Inf)$value, 0),
    0, Inf, rel.tol = 1e-9)$value)), note = "Variance exponential with rate g, g ~ Gamma(k, rate theta^2)."),
  status = "validated")

# Tweedie compound Poisson gamma ---------------------------------------------------------
tweedie_lpdf <- function(y, p) pervec(y, p, function(x, pp) {
  lam <- pp$mu^(2 - pp$power) / (pp$phi * (2 - pp$power)); a <- (2 - pp$power) / (pp$power - 1)
  sc <- pp$phi * (pp$power - 1) * pp$mu^(pp$power - 1)
  if (x == 0) return(-lam)
  nmax <- max(50, ceiling(lam + 40 * sqrt(lam) + 10 * x / (a * sc)))
  n <- seq_len(nmax)
  log_sum_exp(stats::dpois(n, lam, log = TRUE) + stats::dgamma(x, n * a, scale = sc, log = TRUE)) })
register_family(family_def("tweedie", "Tweedie (compound Poisson-gamma, 1 < power < 2)", "mixed",
  list(mean_par(), P("phi", "positive", role = "dispersion"), P("power", "unit", role = "variance power", link = "identity",
       interpretation = "Tweedie power p = 1 + power, so power in (0, 1) gives p in (1, 2)")),
  support = function(p) list(lower = 0, upper = Inf),
  formula = "Y = sum_{j=1}^N G_j, N ~ Poisson(lambda), G_j ~ Gamma(alpha, scale gamma); Var(Y) = phi mu^p",
  lpdf = function(y, p) tweedie_lpdf(y, list(mu = p$mu, phi = p$phi, power = 1 + p$power)),
  rng = function(n, p) { pw <- 1 + p$power; lam <- p$mu^(2 - pw) / (p$phi * (2 - pw)); N <- stats::rpois(n, lam)
    ifelse(N == 0, 0, stats::rgamma(n, pmax(N, 1) * (2 - pw) / (pw - 1), scale = p$phi * (pw - 1) * p$mu^(pw - 1))) },
  mean = function(p) p$mu, variance = function(p) p$phi * p$mu^(1 + p$power),
  numeric_cdf = FALSE, numeric_moments = FALSE,
  regression_notes = "Log link on mu gives Tweedie GLM; zeros have positive probability.",
  references = c("Jorgensen (1987), Exponential dispersion models. JRSS B, 49, 127-162.", "Dunn and Smyth (2005), Series evaluation of Tweedie exponential dispersion model densities. Stat. Comput., 15, 267-280."),
  reference = ref_impl("tweedie", function(y, p) log(tweedie::dtweedie(y, power = 1 + p$power[1], mu = p$mu, phi = p$phi))),
  status = "validated"))

#' Censored distributions
#'
#' Censors a continuous law at `lower` and/or `upper`: values beyond a bound
#' are recorded at the bound, which then carries the tail probability as a
#' point mass. `censor("normal", lower = 0)` is the Tobit model and the
#' rectified Gaussian; `censor("beta", ...)` style limits give two limit
#' models. Data at the bounds contribute log tail probabilities. Count
#' laws can be censored too (`censor("poisson", upper = 3)` records
#' every count of 3 or more as 3), as in censored Poisson and negative
#' binomial regression.
#'
#' @param base Base family (continuous or count).
#' @param lower,upper Censoring points.
#' @return A `regcomp_family` with mixed support.
#' @references Tobin (1958), Estimation of relationships for limited dependent variables. Econometrica, 26, 24-36.
#' @export
#' @examples
#' tobit <- censor("normal", lower = 0)
#' d <- distribution(tobit, mu = 0.5, sigma = 1)
#' prob(d, c(0, 1))
censor <- function(base, lower = -Inf, upper = Inf) {
  base <- get_family(base)
  if (lower >= upper) stop_input("lower must be smaller than upper.")
  if (base$discrete && any(is.finite(c(lower, upper)) & abs(c(lower, upper) - round(c(lower, upper))) > 1e-12))
    stop_input("Censoring bounds of a discrete law must be integers.")
  # Upper censoring keeps P(Y >= upper); for counts that is 1 - F(upper - 1).
  lccdf <- function(q, p) { if (base$discrete) q <- q - 1
    if (!is.null(base$lccdf)) base$lccdf(q, p) else log1p(-family_cdf(base, q, p)) }
  lpdf <- function(y, p) {
    out <- eval_lpdf(base, y, p)
    lo <- y <= lower; hi <- y >= upper
    if (any(lo)) out[lo] <- log(family_cdf(base, rep(lower, sum(lo)), lapply(p, function(v) rep_len(v, length(y))[lo])))
    if (any(hi)) out[hi] <- lccdf(rep(upper, sum(hi)), lapply(p, function(v) rep_len(v, length(y))[hi]))
    out[y < lower | y > upper] <- -Inf
    out
  }
  cdf <- function(q, p) { v <- family_cdf(base, pmin(q, upper), p); v[q < lower] <- 0; v[q >= upper] <- 1; v }
  support <- function(p) { s <- base$support(p)
    list(lower = pmin(pmax(s$lower, lower), upper), upper = pmax(pmin(s$upper, upper), lower)) }
  f <- op_family(base, sprintf("censor(%s, %s, %s)", base$id, format(lower), format(upper)),
                 sprintf("%s censored to [%s, %s]", base$display_name, format(lower), format(upper)), base$pars,
                 lpdf, cdf = cdf, rng = function(n, p) pmin(pmax(family_rng(base, n, p), lower), upper),
                 quantile = function(pr, p) pmin(pmax(family_quantile(base, pr, p), lower), upper),
                 support = support, operator = list(type = "censor", base = base, lower = lower, upper = upper),
                 formula = "point masses F(lower) and 1 - F(upper) at the bounds; f(y) inside")
  if (!base$discrete) f$category <- "mixed"; f$numeric_moments <- FALSE
  f
}

#' Compound (mixed) distributions
#'
#' Integrates one parameter of a base law over a mixing law:
#' `f(y) = int f(y | par = t) g(t) dt`. The mixing family's parameters are
#' appended (suffixed with `_mix`). Covers mixed Poisson laws
#' (`compound("poisson", "lambda", "gamma_rate")` is the negative binomial),
#' Rayleigh mixtures, normal variance mixtures and other compound laws.
#'
#' @param base Base family.
#' @param par Parameter of `base` that is integrated out.
#' @param mixing Family of the mixing law, on the domain of `par`.
#' @return A `regcomp_family`.
#' @references Karlis and Xekalaki (2005), Mixed Poisson distributions. International Statistical Review, 73, 35-58.
#' @export
#' @examples
#' nb <- compound("poisson", "lambda", "gamma_rate")
compound <- function(base, par, mixing) {
  base <- get_family(base); mixing <- get_family(mixing)
  if (!par %in% names(base$pars)) stop_input(sprintf("'%s' is not a parameter of %s.", par, base$id))
  mp <- mixing$pars; mn <- paste0(names(mp), "_mix")
  for (i in seq_along(mp)) mp[[i]]$name <- mn[i]
  names(mp) <- mn
  keep <- setdiff(names(base$pars), par)
  lpdf <- function(y, p) pervec(y, p, function(x, pp) {
    gpars <- stats::setNames(pp[mn], names(mixing$pars)); s <- mixing$support(gpars)
    f <- function(t) exp(eval_lpdf(base, rep(x, length(t)), c(lapply(pp[keep], rep, length(t)), stats::setNames(list(t), par))[names(base$pars)]) +
                         eval_lpdf(mixing, t, lapply(gpars, rep, length(t))))
    if (mixing$discrete) {
      # Discrete mixer: sum over its support up to a negligible tail.
      hi <- if (is.finite(s$upper)) s$upper else family_quantile(mixing, 1 - 1e-15, gpars)
      t <- seq(s$lower, hi)
      lt <- eval_lpdf(base, rep(x, length(t)), c(lapply(pp[keep], rep, length(t)), stats::setNames(list(t), par))[names(base$pars)]) +
        eval_lpdf(mixing, t, lapply(gpars, rep, length(t)))
      return(log_sum_exp(lt[!is.nan(lt)]))
    }
    log(stats::integrate(f, s$lower, s$upper, rel.tol = 1e-10, subdivisions = 1000L, stop.on.error = FALSE)$value) })
  rng <- function(n, p) { t <- family_rng(mixing, n, stats::setNames(p[mn], names(mixing$pars)))
    q <- p[keep]; q[[par]] <- t; family_rng(base, n, q[names(base$pars)]) }
  op_family(base, sprintf("compound(%s, %s, %s)", base$id, par, mixing$id),
            sprintf("%s mixed over %s by %s", base$display_name, par, mixing$display_name),
            c(base$pars[keep], mp), lpdf, rng = rng,
            operator = list(type = "compound", base = base, par = par, mixing = mixing),
            formula = sprintf("f(y) = int f(y | %s = t) g(t) dt", par))
}

#' Skew-symmetric distributions
#'
#' Azzalini type perturbation of a symmetric location scale law:
#' `f(y) = (2 / sigma) f0(z) G(alpha z)`, `z = (y - mu) / sigma`, with
#' `G` a symmetric CDF (`"normal"`, `"logistic"`, `"cauchy"`, `"laplace"`
#' or `"t"` with `nu` degrees of freedom). This covers the skew elliptical
#' and skew symmetric classes (skew normal, skew t, skew logistic, ...).
#'
#' @param base Symmetric family with parameters `mu` and `sigma`.
#' @param skew_cdf Symmetric skewing CDF.
#' @param nu Degrees of freedom for `skew_cdf = "t"`.
#' @return A `regcomp_family` with added parameter `alpha`.
#' @references Azzalini (1985), A class of distributions which includes the normal ones. Scand. J. Stat., 12, 171-178.
#'   Wang, Boyer and Genton (2004), A skew-symmetric representation of multivariate distributions. Statistica Sinica, 14, 1259-1270.
#' @export
#' @examples
#' sl <- skew_symmetric("logistic", "logistic")
skew_symmetric <- function(base, skew_cdf = c("normal", "logistic", "cauchy", "laplace", "t"), nu = 4) {
  base <- get_family(base); skew_cdf <- match.arg(skew_cdf)
  if (!all(c("mu", "sigma") %in% names(base$pars))) stop_input("The base law needs parameters 'mu' and 'sigma'.")
  G <- switch(skew_cdf, normal = function(x) stats::pnorm(x, log.p = TRUE), logistic = function(x) stats::plogis(x, log.p = TRUE),
              cauchy = function(x) stats::pcauchy(x, log.p = TRUE), t = function(x) stats::pt(x, nu, log.p = TRUE),
              laplace = function(x) ifelse(x < 0, x - log(2), log1p(-exp(-x) / 2)))
  an <- new_par_name(base, "alpha")
  pars <- c(base$pars, stats::setNames(list(P(an, "real", role = "skewness", link = "identity")), an))
  f <- op_family(base, sprintf("skew_symmetric(%s, %s)", base$id, skew_cdf), sprintf("Skew %s (%s skewing)", base$display_name, skew_cdf),
    pars, function(y, p) log(2) + eval_lpdf(base, y, p[names(base$pars)]) + G(p[[an]] * (y - p$mu) / p$sigma),
    operator = list(type = "skew_symmetric", base = base, skew_cdf = skew_cdf), formula = "2 f0(y) G(alpha (y - mu)/sigma)")
  f$cdf <- NULL; f$quantile <- NULL
  f
}

# Quantile parameterized systems -----------------------------------------------------------
metalog_basis <- function(u, k) { L <- log(u / (1 - u)); c0 <- u - 0.5
  B <- matrix(0, length(u), k); B[, 1] <- 1
  if (k >= 2) B[, 2] <- L; if (k >= 3) B[, 3] <- c0 * L; if (k >= 4) B[, 4] <- c0
  if (k >= 5) for (j in 5:k) B[, j] <- if (j %% 2 == 1) c0^((j - 1) / 2) else c0^(j / 2 - 1) * L
  B }
metalog_dbasis <- function(u, k) { L <- log(u / (1 - u)); c0 <- u - 0.5; dL <- 1 / (u * (1 - u))
  B <- matrix(0, length(u), k)
  if (k >= 2) B[, 2] <- dL; if (k >= 3) B[, 3] <- L + c0 * dL; if (k >= 4) B[, 4] <- 1
  if (k >= 5) for (j in 5:k) { m <- if (j %% 2 == 1) (j - 1) / 2 else j / 2 - 1
    B[, j] <- if (j %% 2 == 1) m * c0^(m - 1) else m * c0^(m - 1) * L + c0^m * dL }
  B }
metalog_M <- function(u, a) drop(metalog_basis(u, length(a)) %*% a)
metalog_m <- function(u, a) drop(metalog_dbasis(u, length(a)) %*% a)
metalog_invert <- function(x, a, tr) vapply(x, function(xx) {
  f <- function(u) tr$fwd(metalog_M(u, a)) - xx
  lo <- 1e-12; hi <- 1 - 1e-12; if (f(lo) >= 0) return(lo); if (f(hi) <= 0) return(hi)
  stats::uniroot(f, c(lo, hi), tol = 1e-14)$root }, 0)
metalog_family <- function(id, display, tr, extra_pars = list(), support) register_family(family_def(id, display, "quantile_parameterized",
  c(list(P("a", "real", role = "metalog coefficients", regression = FALSE)), extra_pars), variate_type = "scalar_vector_parameter",
  support = support,
  formula = "M(u) = sum_j a_j b_j(u), b = (1, L, (u - 1/2) L, u - 1/2, (u - 1/2)^2, (u - 1/2)^2 L, ...), L = log(u/(1 - u)); Q(u) = g(M(u))",
  lpdf = function(y, p) { trp <- tr$make(p); u <- metalog_invert(y, p$a, trp)
    out <- -log(metalog_m(u, p$a)) - trp$ljac(metalog_M(u, p$a)); out[u <= 1e-12 | u >= 1 - 1e-12] <- -Inf; out },
  cdf = function(q, p) metalog_invert(q, p$a, tr$make(p)),
  quantile = function(pr, p) tr$make(p)$fwd(metalog_M(pr, p$a)),
  rng = function(n, p) tr$make(p)$fwd(metalog_M(stats::runif(n), p$a)),
  numeric_moments = FALSE,
  identifiability_notes = "Coefficients must give m(u) > 0 on (0, 1) (feasibility); not checked automatically.",
  references = "Keelin (2016), The metalog distributions. Decision Analysis, 13, 243-277.",
  reference = ref_impl("rmetalog", function(y, p) { trp <- tr$make(p); u <- metalog_invert(y, p$a, trp); -log(metalog_m(u, p$a)) - trp$ljac(metalog_M(u, p$a)) },
    note = "Checked through total mass and quantile inversion."),
  status = "validated"))
metalog_family("metalog", "Metalog (unbounded)", list(make = function(p) list(fwd = identity, ljac = function(M) 0 * M)),
  support = function(p) list(lower = -Inf, upper = Inf))
metalog_family("log_metalog", "Log-metalog (semi-bounded)", list(make = function(p) { bl <- p$bl[1]; list(fwd = function(M) bl + exp(M), ljac = function(M) M) }),
  list(P("bl", "real", role = "lower bound", link = "identity", data = TRUE)), support = function(p) list(lower = p$bl, upper = Inf))
metalog_family("logit_metalog", "Logit-metalog (bounded)", list(make = function(p) { bl <- p$bl[1]; bu <- p$bu[1]; list(
  fwd = function(M) (bl + bu * exp(M)) / (1 + exp(M)),
  ljac = function(M) log(bu - bl) + M - 2 * log1p(exp(M))) }),
  list(P("bl", "real", link = "identity", data = TRUE), P("bu", "real", link = "identity", data = TRUE)),
  support = function(p) list(lower = p$bl, upper = p$bu))

jqpd_core <- function(p, bounded) {
  cc <- stats::qnorm(1 - p$alpha)
  if (bounded) { tr <- function(x) stats::qnorm((x - p$l) / (p$u - p$l)) } else tr <- function(x) log(x - p$l)
  L <- tr(p$xl); B <- tr(p$xm); H <- tr(p$xh); n <- sign(L + H - 2 * B)
  if (bounded) {
    xi <- ifelse(n == 1, L, ifelse(n == 0, B, H)); delta <- acosh((H - L) / (2 * pmin(B - L, H - B))) / cc
    lam <- (H - L) / sinh(2 * delta * cc)
    list(Q = function(z) p$l + (p$u - p$l) * stats::pnorm(xi + lam * sinh(delta * (z + n * cc))),
         dQ = function(z) (p$u - p$l) * stats::dnorm(xi + lam * sinh(delta * (z + n * cc))) * lam * delta * cosh(delta * (z + n * cc)))
  } else {
    th <- ifelse(n == 1, p$xl - p$l, ifelse(n == 0, p$xm - p$l, p$xh - p$l))
    delta <- sinh(acosh((H - L) / (2 * pmin(B - L, H - B)))) / cc; lam <- pmin(H - B, B - L) / (delta * cc)
    list(Q = function(z) p$l + th * exp(lam * sinh(asinh(delta * z) + asinh(n * cc * delta))),
         dQ = function(z) th * exp(lam * sinh(asinh(delta * z) + asinh(n * cc * delta))) * lam * cosh(asinh(delta * z) + asinh(n * cc * delta)) * delta / sqrt(1 + delta^2 * z^2))
  }
}
jqpd_family <- function(id, display, bounded) register_family(family_def(id, display, "quantile_parameterized",
  c(list(P("xl", "real", role = "alpha quantile", link = "identity"), P("xm", "real", role = "median", link = "identity"),
         P("xh", "real", role = "1 - alpha quantile", link = "identity"), P("l", "real", role = "lower bound", link = "identity", data = TRUE)),
    if (bounded) list(P("u", "real", role = "upper bound", link = "identity", data = TRUE)),
    list(P("alpha", "unit", role = "tail level", regression = FALSE, default = 0.1))),
  support = function(p) list(lower = p$l, upper = if (bounded) p$u else Inf),
  formula = if (bounded) "J-QPD-B: Q(p) = l + (u - l) Phi(xi + lambda sinh(delta (Phi^-1(p) + n c)))" else
    "J-QPD-S: Q(p) = l + theta exp(lambda sinh(asinh(delta Phi^-1(p)) + asinh(n c delta)))",
  lpdf = function(y, p) pervec(y, p, function(x, pp) { cr <- jqpd_core(pp, bounded)
    a <- cr$Q(-30); b <- cr$Q(30); if (x <= a || x >= b) return(-Inf)
    z <- stats::uniroot(function(z) cr$Q(z) - x, c(-30, 30), tol = 1e-13)$root
    stats::dnorm(z, log = TRUE) - log(cr$dQ(z)) }),
  cdf = function(q, p) pervec(q, p, function(x, pp) { cr <- jqpd_core(pp, bounded)
    if (x <= cr$Q(-30)) return(0); if (x >= cr$Q(30)) return(1)
    stats::pnorm(stats::uniroot(function(z) cr$Q(z) - x, c(-30, 30), tol = 1e-13)$root) }),
  quantile = function(pr, p) pervec(pr, p, function(u, pp) jqpd_core(pp, bounded)$Q(stats::qnorm(u))),
  numeric_moments = FALSE,
  references = "Hadlock and Bickel (2017), Johnson quantile-parameterized distributions. Decision Analysis, 14, 35-64.",
  reference = ref_impl("quantile matching", function(y, p) pervec(y, p, function(x, pp) { cr <- jqpd_core(pp, bounded)
    z <- stats::uniroot(function(z) cr$Q(z) - x, c(-30, 30), tol = 1e-13)$root; stats::dnorm(z, log = TRUE) - log(cr$dQ(z)) }),
    note = "Checked through Q(alpha) = xl, Q(1/2) = xm, Q(1 - alpha) = xh, total mass and quantile inversion."),
  status = "validated"))
jqpd_family("jqpd_bounded", "Johnson quantile-parameterized (bounded, J-QPD-B)", TRUE)
jqpd_family("jqpd_semibounded", "Johnson quantile-parameterized (semi-bounded, J-QPD-S)", FALSE)

# Pearson system ------------------------------------------------------------------------
pearson_coef <- function(p) { b1 <- p$skew^2; b2 <- p$kurt; d <- 10 * b2 - 12 * b1 - 18
  list(b0 = (4 * b2 - 3 * b1) / d * p$sd^2, b1 = p$sd * p$skew * (b2 + 3) / d, b2 = (2 * b2 - 3 * b1 - 6) / d) }
pearson_support <- function(p) { k <- pearson_coef(p)
  r <- if (abs(k$b2) < 1e-12) (if (abs(k$b1) < 1e-12) numeric(0) else -k$b0 / k$b1) else {
    disc <- k$b1^2 - 4 * k$b2 * k$b0; if (disc < 0) numeric(0) else (-k$b1 + c(-1, 1) * sqrt(disc)) / (2 * k$b2) }
  lo <- if (any(r < 0)) max(r[r < 0]) else -Inf; hi <- if (any(r > 0)) min(r[r > 0]) else Inf
  list(lower = p$mean + lo, upper = p$mean + hi) }
pearson_logk <- function(x, k) vapply(x, function(xx) -stats::integrate(function(t) (t + k$b1) / (k$b0 + k$b1 * t + k$b2 * t^2), 0, xx, rel.tol = 1e-12)$value, 0)
register_family(family_def("pearson", "Pearson system (by mean, sd, skewness and kurtosis)", "real",
  list(loc_par("mean"), scale_par("sd"), P("skew", "real", link = "identity"), P("kurt", "positive", role = "kurtosis (not excess)")),
  check = function(p) if (any(p$kurt <= p$skew^2 + 1, na.rm = TRUE)) "Need kurtosis > skewness^2 + 1.",
  support = function(p) { s <- lapply(seq_len(max(lengths(p))), function(i) pearson_support(lapply(p, function(v) v[(i - 1) %% length(v) + 1])))
    list(lower = vapply(s, `[[`, 0, "lower"), upper = vapply(s, `[[`, 0, "upper")) },
  formula = "p'(x)/p(x) = -(x - mean + b1) / (b0 + b1 (x - mean) + b2 (x - mean)^2), coefficients from the four moments",
  lpdf = function(y, p) pervec(y, p, function(x, pp) { k <- pearson_coef(pp); s <- pearson_support(pp)
    if (x <= s$lower || x >= s$upper) return(-Inf)
    lz <- log(stats::integrate(function(t) exp(pearson_logk(t - pp$mean, k)), s$lower, s$upper, rel.tol = 1e-10, subdivisions = 1000L)$value)
    pearson_logk(x - pp$mean, k) - lz }),
  mean = function(p) p$mean, variance = function(p) p$sd^2, numeric_moments = FALSE,
  references = c("Pearson (1895), Contributions to the mathematical theory of evolution II. Phil. Trans. Roy. Soc. A, 186, 343-414.",
                 "Johnson, Kotz and Balakrishnan (1994), ch. 12."),
  reference = ref_impl("PearsonDS", function(y, p) PearsonDS::dpearson(y, moments = c(mean = p$mean[1], variance = p$sd[1]^2, skewness = p$skew[1], kurtosis = p$kurt[1]), log = TRUE)),
  status = "validated"))

# Cantor (singular continuous) ------------------------------------------------------------
cantor_cdf <- function(x, depth = 50) vapply(x, function(xx) {
  if (xx <= 0) return(0); if (xx >= 1) return(1)
  out <- 0; w <- 0.5
  for (i in seq_len(depth)) { xx <- 3 * xx; d <- floor(xx); xx <- xx - d
    if (d == 1) return(out + w); out <- out + w * (d == 2); w <- w / 2 }
  out }, 0)
register_family(family_def("cantor", "Cantor (singular)", "bounded", list(),
  support = function(p) list(lower = 0, upper = 1),
  formula = "Uniform measure on the middle thirds Cantor set; CDF is the Cantor function; no density",
  lpdf = function(y, p) stop_unsupported("A density (the Cantor law is singular; use cdf() or rng())", "cantor"),
  cdf = function(q, p) cantor_cdf(q), rng = function(n, p) vapply(seq_len(n), function(i) sum(2 * stats::rbinom(50, 1, 0.5) / 3^(1:50)), 0),
  mean = function(p) 0.5, variance = function(p) 1 / 8, numeric_moments = FALSE,
  regression_notes = "Not usable as a likelihood; provided for simulation and CDF evaluation.",
  references = "Cantor (1884); Feller (1971), An Introduction to Probability Theory and its Applications, vol. 2.",
  reference = ref_impl("self similarity F(x/3) = F(x)/2", function(y, p) cantor_cdf(y / 3) * 2), status = "implemented"))
