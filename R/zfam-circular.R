# Circular and directional distributions. Angles are on (-pi, pi].

circ_support <- function(p) list(lower = -pi, upper = pi)
cfam <- function(id, display, pars, ..., support = circ_support) register_family(family_def(id, display, "circular", pars,
  support = support, numeric_moments = FALSE, ...))
circ_mu <- function() P("mu", "circular", role = "mean direction", offset_allowed = TRUE)
conc <- function(name = "kappa") P(name, "positive", role = "concentration")
wrap_angle <- function(x) x - 2 * pi * floor((x + pi) / (2 * pi))

cfam("von_mises", "von Mises", list(circ_mu(), conc()), aliases = c("circular_normal", "VM"),
  formula = "f(y) = exp(kappa cos(y - mu)) / (2 pi I0(kappa))",
  lpdf = function(y, p) p$kappa * (cos(y - p$mu) - 1) - log(2 * pi * besselI(p$kappa, 0, expon.scaled = TRUE)),
  rng = function(n, p) { r <- recycle_pars(c(list(.n = seq_len(n)), p)); wrap_angle(vapply(seq_len(n), function(i) {
    k <- r$kappa[i]; if (k < 1e-8) return(stats::runif(1, -pi, pi))
    a <- 1 + sqrt(1 + 4 * k^2); b <- (a - sqrt(2 * a)) / (2 * k); rr <- (1 + b^2) / (2 * b)
    repeat { u <- stats::runif(3); z <- cos(pi * u[1]); f <- (1 + rr * z) / (rr + z); cc <- k * (rr - f)
      if (cc * (2 - cc) - u[2] > 0 || log(cc / u[2]) + 1 - cc >= 0) break }
    r$mu[i] + sign(u[3] - 0.5) * acos(f) }, 0)) },
  references = "Mardia and Jupp (2000), Directional Statistics. Wiley.",
  reference = ref_impl("circular", function(y, p) log(mapply(function(x, m, k) circular::dvonmises(circular::circular(x), circular::circular(m), k), y, p$mu, p$kappa))),
  stan = stan_native("von_mises_lpdf({y} | {mu}, {kappa})", rng = "von_mises_rng({mu}, {kappa})"),
  regression_notes = "Use the tan_half link for the mean direction (Fisher and Lee 1992).",
  status = "validated")

cfam("axial_von_mises", "Axial von Mises", list(P("mu", "circular", role = "axis"), conc()),
  formula = "f(y) = cosh(kappa cos(y - mu)) / (2 pi I0(kappa)) on the circle (axial data)",
  lpdf = function(y, p) { c1 <- p$kappa * cos(y - p$mu); abs(c1) + log1p(exp(-2 * abs(c1))) - log(2) - log(2 * pi) - log(besselI(p$kappa, 0, expon.scaled = TRUE)) - p$kappa },
  references = "arXiv 2504.04681 (2025), Multimodal distributions for circular axial data.",
  reference = ref_impl("mixture of two von Mises", function(y, p) log(0.5 * exp(p$kappa * cos(y - p$mu)) / (2 * pi * besselI(p$kappa, 0)) +
    0.5 * exp(-p$kappa * cos(y - p$mu)) / (2 * pi * besselI(p$kappa, 0)))), status = "validated")

cfam("cardioid", "Cardioid", list(circ_mu(), P("rho", "unit", role = "concentration", interpretation = "rho in (0, 1/2), scaled by 1/2")),
  formula = "f(y) = (1 + 2 r cos(y - mu)) / (2 pi), r = rho / 2",
  lpdf = function(y, p) log1p(p$rho * cos(y - p$mu)) - log(2 * pi),
  references = "Jeffreys (1961); Mardia and Jupp (2000).",
  reference = ref_impl("circular", function(y, p) log(mapply(function(x, m, r) circular::dcardioid(circular::circular(x), circular::circular(m), r / 2), y, p$mu, p$rho))),
  status = "validated")

cfam("wrapped_cauchy", "Wrapped Cauchy", list(circ_mu(), P("rho", "unit", role = "concentration")),
  formula = "f(y) = (1 - rho^2) / (2 pi (1 + rho^2 - 2 rho cos(y - mu)))",
  lpdf = function(y, p) log1p(-p$rho^2) - log(2 * pi) - log1p(p$rho^2 - 2 * p$rho * cos(y - p$mu)),
  rng = function(n, p) wrap_angle(stats::rcauchy(n, p$mu, -log(p$rho))),
  references = "Levy (1939); Mardia and Jupp (2000).",
  reference = ref_impl("circular", function(y, p) log(mapply(function(x, m, r) circular::dwrappedcauchy(circular::circular(x), circular::circular(m), r), y, p$mu, p$rho))),
  stan = sfun("wrapped_cauchy", c("mu", "rho"), "    return log1m(square(rho)) - log(2 * pi()) - log1p(square(rho) - 2 * rho * cos(y - mu));"),
  status = "validated")

cfam("wrapped_normal", "Wrapped normal", list(circ_mu(), scale_par("sigma")),
  formula = "f(y) = sum_k phi((y - mu + 2 pi k)/sigma) / sigma",
  lpdf = function(y, p) log(pervec(y, p, function(x, pp) { K <- ceiling(6 * pp$sigma / (2 * pi)) + 2; k <- -K:K
    sum(stats::dnorm(x - pp$mu + 2 * pi * k, 0, pp$sigma)) })),
  rng = function(n, p) wrap_angle(stats::rnorm(n, p$mu, p$sigma)),
  references = "Mardia and Jupp (2000).",
  reference = ref_impl("circular", function(y, p) log(mapply(function(x, m, s) circular::dwrappednormal(circular::circular(x), circular::circular(m), sd = s), y, p$mu, p$sigma))),
  status = "validated")

#' Wrapped distributions
#'
#' Wraps a law on the real line (or the positive half line) onto the circle
#' of period `2 pi`: `f_W(y) = sum_k f(y + 2 pi k)`. The sum is truncated
#' once terms fall below `1e-14` of the running total (at most 2000 terms on
#' each side). Gives the wrapped normal, Cauchy, Laplace, asymmetric
#' Laplace, exponential, Levy, Lindley, xgamma and any other wrapped law.
#'
#' @param base Base family on the real or positive line.
#' @return A circular `regcomp_family`.
#' @references Mardia and Jupp (2000), Directional Statistics, sec. 3.5.
#' @export
#' @examples
#' we <- wrap("exponential")
#' d <- distribution(we, lambda = 0.8)
#' prob(d, c(-1, 0.5, 2))
wrap <- function(base) {
  base <- get_family(base)
  if (base$discrete || base$variate_type != "scalar") stop_unsupported("Wrapping a discrete or vector valued law", base$id)
  lpdf <- function(y, p) log(pervec(y, p, function(x, pp) {
    xx <- x + 2 * pi * (x < 0) * (base$category == "positive")
    # Sum windings in growing blocks until the base mass outside the covered
    # range (which bounds 2 pi times the remaining terms for decaying
    # densities) is negligible.
    dens <- function(t) exp(eval_lpdf(base, t, lapply(pp, rep_len, length(t))))
    tot <- dens(xx); K0 <- 0; K <- 64
    repeat {
      k <- (K0 + 1):K
      tot <- tot + sum(dens(xx + 2 * pi * k)) + sum(dens(xx - 2 * pi * k))
      hi <- xx + 2 * pi * K; lo <- xx - 2 * pi * K
      rest <- (1 - family_cdf(base, hi, pp)) + family_cdf(base, lo, pp)
      if (!is.finite(rest) || rest < 1e-12 * max(tot, 1e-300) * 2 * pi) break
      K0 <- K; K <- 4 * K
      if (K > 1e7) stop_numerical("The wrapped density series did not converge.")
    }
    tot }))
  f <- op_family(base, sprintf("wrap(%s)", base$id), sprintf("Wrapped %s", base$display_name), base$pars, lpdf,
                 rng = function(n, p) wrap_angle(family_rng(base, n, p)), support = circ_support,
                 operator = list(type = "wrap", base = base), formula = "f_W(y) = sum_k f(y + 2 pi k)")
  f$category <- "circular"; f$numeric_moments <- FALSE
  f
}

#' Sine-skewed circular distributions
#'
#' `f(y) = f0(y) (1 + lambda sin(y - mu))` for a symmetric circular base
#' law `f0` with location `mu` (Abe and Pewsey 2011), `lambda` in (-1, 1).
#'
#' @param base Symmetric circular family with a parameter named `mu`.
#' @return A `regcomp_family`.
#' @references Abe and Pewsey (2011), Sine-skewed circular distributions. Statistical Papers, 52, 683-707.
#' @export
#' @examples
#' ssvm <- sine_skew("von_mises")
sine_skew <- function(base) {
  base <- get_family(base)
  if (!"mu" %in% names(base$pars)) stop_input("The base law needs a location parameter 'mu'.")
  ln <- new_par_name(base, "lambda")
  pars <- c(base$pars, stats::setNames(list(P(ln, "symmetric_unit", role = "skewness")), ln))
  f <- op_family(base, sprintf("sine_skew(%s)", base$id), sprintf("Sine-skewed %s", base$display_name), pars,
    function(y, p) eval_lpdf(base, y, p[names(base$pars)]) + log1p(p[[ln]] * sin(y - p$mu)),
    operator = list(type = "sine_skew", base = base), formula = "f0(y) (1 + lambda sin(y - mu))")
  f$numeric_moments <- FALSE
  f
}

cfam("henyey_greenstein", "Henyey-Greenstein (phase function on the circle of scattering angles)", list(P("g", "symmetric_unit", role = "asymmetry")),
  support = function(p) list(lower = 0, upper = pi),
  formula = "Scattering angle density f(theta) = (1 - g^2) sin(theta) / (2 (1 + g^2 - 2 g cos theta)^(3/2)), theta in (0, pi)",
  lpdf = function(y, p) ifelse(y < 0 | y > pi, -Inf, log1p(-p$g^2) + log(sin(y)) - log(2) - 1.5 * log1p(p$g^2 - 2 * p$g * cos(y))),
  cdf = function(q, p) { x <- pmin(pmax(q, 0), pi); g <- p$g
    ifelse(abs(g) < 1e-10, (1 - cos(x)) / 2, (1 - g^2) / (2 * g) * (1 / (1 - g) - 1 / sqrt(1 + g^2 - 2 * g * cos(x)))) },
  references = "Henyey and Greenstein (1941), Diffuse radiation in the galaxy. Astrophysical J., 93, 70-83.",
  reference = ref_impl("definition", function(y, p) log((1 - p$g^2) * sin(y) / (2 * (1 + p$g^2 - 2 * p$g * cos(y))^1.5))),
  status = "validated")

pn_lpdf <- function(y, p) { u <- cbind(cos(y), sin(y)); m <- cbind(p$m1, p$m2)
  a <- rowSums(u * m); -0.5 * (p$m1^2 + p$m2^2) - log(2 * pi) + log1p(a * stats::pnorm(a) / stats::dnorm(a)) }
cfam("projected_normal", "Projected normal (isotropic)", list(P("m1", "real", link = "identity"), P("m2", "real", link = "identity")),
  formula = "Direction of a bivariate N((m1, m2), I) vector",
  lpdf = pn_lpdf, rng = function(n, p) atan2(stats::rnorm(n, p$m2), stats::rnorm(n, p$m1)),
  regression_notes = "Regress m1 and m2 on covariates (Presnell, Morrison and Littell 1998).",
  references = "Presnell, Morrison and Littell (1998), Projected multivariate linear models for directional data. JASA, 93, 1068-1077.",
  reference = ref_impl("integral", function(y, p) log(mapply(function(x, a, b) stats::integrate(function(r) r * stats::dnorm(r * cos(x), a) * stats::dnorm(r * sin(x), b), 0, Inf)$value, y, p$m1, p$m2))),
  status = "validated")

cfam("nnts", "Nonnegative trigonometric sums (Fernandez-Duran)", list(P("c", "real", role = "coefficients", regression = FALSE)),
  variate_type = "scalar_vector_parameter",
  formula = "f(y) = |sum_{k=0}^M c_k exp(i k y)|^2 / (2 pi), sum |c_k|^2 = 1 (real coefficients here)",
  lpdf = function(y, p) { cc <- p$c / sqrt(sum(p$c^2)); k <- seq_along(cc) - 1
    vapply(y, function(x) log(Mod(sum(cc * exp(1i * k * x)))^2 / (2 * pi)), 0) },
  references = "Fernandez-Duran (2004), Circular distributions based on nonnegative trigonometric sums. Biometrics, 60, 499-503.",
  reference = ref_impl("definition", function(y, p) { cc <- p$c / sqrt(sum(p$c^2)); vapply(y, function(x)
    log((sum(cc * cos((seq_along(cc) - 1) * x))^2 + sum(cc * sin((seq_along(cc) - 1) * x))^2) / (2 * pi)), 0) }),
  status = "validated")

cfam("dirac_comb", "Dirac comb on the circle (discrete uniform on m equally spaced angles)", list(P("m", "positive_integer", data = TRUE),
  P("phase", "circular", role = "offset")),
  formula = "P(Y = phase + 2 pi j / m) = 1/m", discrete = TRUE,
  lpdf = function(y, p) { j <- (y - p$phase) * p$m / (2 * pi); ifelse(abs(j - round(j)) < 1e-9, -log(p$m), -Inf) },
  rng = function(n, p) wrap_angle(p$phase + 2 * pi * sample.int(p$m[1], n, TRUE) / p$m[1]),
  extra = list(grid = 0), numeric_cdf = FALSE,
  references = "Wikipedia, Dirac comb; Mardia and Jupp (2000), lattice distributions.",
  reference = ref_impl("definition", function(y, p) rep(-log(p$m), length(y))), status = "validated")

# Directional laws on spheres --------------------------------------------------------
sph <- function(id, display, pars, ...) register_family(family_def(id, display, "directional", pars, variate_type = "unit_vector",
  numeric_cdf = FALSE, numeric_moments = FALSE, ...))
log_vmf_c <- function(k, d) (d / 2 - 1) * log(k) - (d / 2) * log(2 * pi) - (log(besselI(k, d / 2 - 1, expon.scaled = TRUE)) + k)
sph("von_mises_fisher", "von Mises-Fisher", list(P("mu", "any", role = "mean direction (unit vector)", regression = FALSE), conc()),
  aliases = "vmf", formula = "f(x) = C_d(kappa) exp(kappa mu'x), x on the unit sphere S^(d-1)",
  lpdf = function(y, p) { y <- as.matrix(y); m <- p$mu / sqrt(sum(p$mu^2)); log_vmf_c(p$kappa[1], ncol(y)) + p$kappa[1] * drop(y %*% m) },
  references = "Fisher (1953), Dispersion on a sphere. Proc. Roy. Soc. A, 217, 295-305.",
  reference = ref_impl("Directional", function(y, p) Directional::dvmf(as.matrix(y), p$mu / sqrt(sum(p$mu^2)), p$kappa[1], logden = TRUE)),
  status = "validated")

kent_logc <- function(k, b) { j <- if (b == 0) 0 else 0:200
  log(2 * pi) + log(sum(exp(lgamma(j + 0.5) - lgamma(j + 1) + ifelse(j == 0, 0, 2 * j * log(b)) + (-2 * j - 0.5) * log(k / 2) +
    log(besselI(k, 2 * j + 0.5, expon.scaled = TRUE)) + k))) }
sph("kent", "Kent (FB5)", list(P("gamma1", "any", role = "mean direction"), P("gamma2", "any", role = "major axis"), P("gamma3", "any", role = "minor axis"),
  conc(), P("beta", "nonneg", role = "ovalness")), aliases = "fisher_bingham_5",
  formula = "f(x) = exp(kappa g1'x + beta ((g2'x)^2 - (g3'x)^2)) / c(kappa, beta) on S^2, 2 beta < kappa for unimodality",
  lpdf = function(y, p) { y <- as.matrix(y); k <- p$kappa[1]; b <- p$beta[1]
    k * drop(y %*% p$gamma1) + b * (drop(y %*% p$gamma2)^2 - drop(y %*% p$gamma3)^2) - kent_logc(k, b) },
  numerical_notes = "Normalizing constant by Kent's (1982) series in Bessel functions.",
  references = "Kent (1982), The Fisher-Bingham distribution on the sphere. JRSS B, 44, 71-80.",
  reference = ref_impl("Directional", function(y, p) Directional::dkent(as.matrix(y), G = cbind(p$gamma1, p$gamma2, p$gamma3), param = c(p$kappa[1], p$beta[1]), logden = TRUE)),
  status = "validated")

bingham_logc <- function(lam) {
  # Normalizing constant by quadrature in spherical coordinates (d = 2, 3).
  d <- length(lam)
  if (d == 2) return(log(stats::integrate(function(t) exp(lam[1] * cos(t)^2 + lam[2] * sin(t)^2), 0, 2 * pi)$value))
  if (d != 3) stop_unsupported("Bingham normalizing constant for dimension above 3", "bingham")
  f <- function(th) vapply(th, function(t) stats::integrate(function(ph) exp(lam[1] * (sin(t) * cos(ph))^2 + lam[2] * (sin(t) * sin(ph))^2 +
    lam[3] * cos(t)^2) * sin(t), 0, 2 * pi)$value, 0)
  log(stats::integrate(f, 0, pi, rel.tol = 1e-10)$value)
}
sph("bingham", "Bingham (axial)", list(P("A", "any", role = "symmetric parameter matrix", regression = FALSE)),
  formula = "f(x) = exp(x' A x) / c(A) on S^(d-1), antipodally symmetric",
  lpdf = function(y, p) { y <- as.matrix(y); A <- as.matrix(p$A); e <- eigen(A, symmetric = TRUE)$values
    rowSums((y %*% A) * y) - bingham_logc(e) },
  numerical_notes = "Normalizing constant by quadrature for d = 2 and d = 3.",
  references = "Bingham (1974), An antipodally symmetric distribution on the sphere. Ann. Stat., 2, 1201-1225.",
  reference = ref_impl("Directional", function(y, p) { A <- as.matrix(p$A); e <- eigen(A, symmetric = TRUE)$values
    rowSums((as.matrix(y) %*% A) * as.matrix(y)) - log(Directional::fb.saddle(rep(0, 3), -e)[3]) }, note = "Saddlepoint approximation, so agreement is to about 1e-2."),
  status = "implemented")

sph("esag", "Elliptically symmetric angular Gaussian (ESAG, S^2)", list(P("mu", "any", role = "mean vector"), P("gamma1", "real", link = "identity"), P("gamma2", "real", link = "identity")),
  formula = "Projected N(mu, V) with V^-1 = I + g1 (x1 x1' - x2 x2') + g2 (x1 x2' + x2 x1') + (sqrt(1 + g1^2 + g2^2) - 1)(x1 x1' + x2 x2'), Paine et al. (2018)",
  lpdf = function(y, p) { y <- as.matrix(y); mu <- p$mu; g1 <- p$gamma1[1]; g2 <- p$gamma2[1]
    m0 <- sqrt(mu[2]^2 + mu[3]^2); x1 <- c(-m0^2, mu[1] * mu[2], mu[1] * mu[3]) / (m0 * sqrt(sum(mu^2))); x2 <- c(0, -mu[3], mu[2]) / m0
    Vi <- diag(3) + g1 * (tcrossprod(x1) - tcrossprod(x2)) + g2 * (tcrossprod(x1, x2) + tcrossprod(x2, x1)) +
      (sqrt(1 + g1^2 + g2^2) - 1) * (tcrossprod(x1) + tcrossprod(x2))
    C <- rowSums((y %*% Vi) * y); a <- drop(y %*% mu) / sqrt(C)
    M2 <- (1 + a^2) * stats::pnorm(a) + a * stats::dnorm(a)
    -log(2 * pi) - 1.5 * log(C) + 0.5 * (a^2 - sum(mu^2)) + log(M2) },
  references = c("Paine, Preston, Tsagris and Wood (2018), An elliptically symmetric angular Gaussian distribution. Stat. Comput., 28, 689-697.",
                 "arXiv 2212.05634 (2022), Elliptically symmetric distributions for directional data of arbitrary dimension."),
  reference = ref_impl("Directional", function(y, p) log(Directional::desag(as.matrix(y), p$mu, c(p$gamma1[1], p$gamma2[1])))),
  status = "implemented")
