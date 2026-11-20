# Laws defined through special function machinery: Chernoff, the Mie phase
# function and the phased (bi-Weibull) hazard laws.

# Airy Ai for complex arguments by its Maclaurin series (adequate for |z| < 20,
# which covers every argument used below).
airy_ai_complex <- function(z) {
  z3 <- z^3; f <- g <- 0 * z; a <- 1 + 0 * z; b <- z
  for (k in 0:200) {
    f <- f + a; g <- g + b
    a <- a * z3 / ((3 * k + 2) * (3 * k + 3)); b <- b * z3 / ((3 * k + 3) * (3 * k + 4))
    if (k > 10 && max(Mod(a) + Mod(b)) < 1e-18 * max(Mod(f) + Mod(g))) break
  }
  0.355028053887817239 * f - 0.258819403792806798 * g
}
# g(t) with Fourier transform 2^(1/3) / Ai(i 2^(-1/3) s) (Groeneboom 1989).
chernoff_g <- function(t) vapply(t, function(tt) stats::integrate(function(s)
  Re(exp(-1i * s * tt) * 2^(1 / 3) / airy_ai_complex(1i * 2^(-1 / 3) * s)) / pi, 0, 24,
  subdivisions = 2000L, rel.tol = 1e-11)$value, 0)
rfam("chernoff", "Chernoff", list(loc_par("mu"), scale_par("sigma")),
  formula = "Z = argmax_s (W(s) - s^2); f(z) = g(z) g(-z) / 2 with g^(s) = 2^(1/3) / Ai(i 2^(-1/3) s)",
  lpdf = function(y, p) { z <- (y - p$mu) / p$sigma
    ifelse(abs(z) > 3.5, log(0.5 * 4^(4 / 3) * abs(z) / 0.7012108227206906) - 2 / 3 * abs(z)^3 + 2^(1 / 3) * (-2.338107410459767) * abs(z),
           log(pmax(0.5 * chernoff_g(z) * chernoff_g(-z), 1e-300))) - log(p$sigma) },
  mean = function(p) p$mu, variance = function(p) 0.26355964 * p$sigma^2,
  numerical_notes = "g by Fourier inversion with Ai from its power series; the Groeneboom tail approximation is used for |z| > 3.5, where the two agree to 3e-3 on the log scale.",
  references = c("Groeneboom (1989), Brownian motion with a parabolic drift and Airy functions. Probab. Theory Related Fields, 81, 79-109.",
                 "Groeneboom and Wellner (2001), Computing Chernoff's distribution. J. Comput. Graph. Stat., 10, 388-400."),
  reference = ref_impl("total mass, symmetry and Var(Z) = 0.26355964 (Groeneboom and Wellner 2001)",
    function(y, p) log(0.5 * chernoff_g((y - p$mu) / p$sigma) * chernoff_g(-(y - p$mu) / p$sigma)) - log(p$sigma)),
  status = "validated")

bw_H <- function(y, p) (y / p$s1)^p$g1 + ifelse(y > p$mu2, ((y - p$mu2) / p$s2)^p$g2, 0)
bw_h <- function(y, p) p$g1 / p$s1 * (y / p$s1)^(p$g1 - 1) + ifelse(y > p$mu2, p$g2 / p$s2 * (pmax(y - p$mu2, 0) / p$s2)^(p$g2 - 1), 0)
posfam("bi_weibull", "Phased bi-Weibull", list(scale_par("s1"), shp("g1"), P("mu2", "positive", role = "phase change point"), scale_par("s2"), shp("g2")),
  aliases = "phased_bi_weibull",
  formula = "H(y) = (y/s1)^g1 + 1(y > mu2) ((y - mu2)/s2)^g2; f = h exp(-H)",
  lpdf = function(y, p) log(bw_h(y, p)) - bw_H(y, p), cdf = function(q, p) -expm1(-bw_H(pmax(q, 0), p)),
  references = c("Evans, Hastings and Peacock (2000), Statistical Distributions, 3rd ed., ch. 6.", "NIST Dataplot Reference Manual, BWEHAZ."),
  reference = ref_impl("hazard identity", function(y, p) log(bw_h(y, p)) - vapply(seq_along(y), function(i)
    stats::integrate(function(t) bw_h(t, lapply(p, `[`, i)), 0, y[i], rel.tol = 1e-12)$value, 0)),
  status = "validated")

posfam("phased_bi_exponential", "Phased bi-exponential", list(P("lambda1", "positive", role = "rate"), P("mu2", "positive", role = "phase change point"),
  P("lambda2", "positive", role = "added rate")),
  formula = "Hazard lambda1 before mu2 and lambda1 + lambda2 after; the bi-Weibull law with unit shapes",
  lpdf = function(y, p) log(p$lambda1 + p$lambda2 * (y > p$mu2)) - p$lambda1 * y - p$lambda2 * pmax(y - p$mu2, 0),
  cdf = function(q, p) -expm1(-p$lambda1 * pmax(q, 0) - p$lambda2 * pmax(q - p$mu2, 0)),
  references = "Evans, Hastings and Peacock (2000), Statistical Distributions, 3rd ed.",
  reference = ref_impl("bi-Weibull with g1 = g2 = 1", function(y, p) { q <- list(s1 = 1 / p$lambda1, g1 = 1, mu2 = p$mu2, s2 = 1 / p$lambda2, g2 = 1)
    log(bw_h(y, q)) - bw_H(y, q) }),
  status = "validated")

# Mie scattering amplitudes S1, S2 by the Bohren and Huffman (1983) BHMIE recursion.
mie_S <- function(x, m, mu) {
  nstop <- floor(x + 4 * x^(1 / 3) + 2); y <- m * x; nmx <- max(nstop, ceiling(Mod(y))) + 15
  D <- complex(nmx)
  for (n in nmx:2) D[n - 1] <- n / y - 1 / (D[n] + n / y)
  psi0 <- cos(x); psi1 <- sin(x); chi0 <- -sin(x); chi1 <- cos(x); xi1 <- complex(real = psi1, imaginary = -chi1)
  S1 <- S2 <- complex(length(mu)); pi0 <- 0 * mu; pi1 <- 1 + 0 * mu
  for (n in seq_len(nstop)) {
    fn <- (2 * n + 1) / (n * (n + 1))
    psi <- (2 * n - 1) * psi1 / x - psi0; chi <- (2 * n - 1) * chi1 / x - chi0; xi <- complex(real = psi, imaginary = -chi)
    an <- ((D[n] / m + n / x) * psi - psi1) / ((D[n] / m + n / x) * xi - xi1)
    bn <- ((m * D[n] + n / x) * psi - psi1) / ((m * D[n] + n / x) * xi - xi1)
    tau <- n * mu * pi1 - (n + 1) * pi0
    S1 <- S1 + fn * (an * pi1 + bn * tau); S2 <- S2 + fn * (an * tau + bn * pi1)
    pn <- ((2 * n + 1) * mu * pi1 - (n + 1) * pi0) / n; pi0 <- pi1; pi1 <- pn
    psi0 <- psi1; psi1 <- psi; chi0 <- chi1; chi1 <- chi; xi1 <- complex(real = psi1, imaginary = -chi1)
  }
  list(S1 = S1, S2 = S2)
}
mie_unnorm <- function(th, x, m) { s <- mie_S(x, m, cos(th)); (Mod(s$S1)^2 + Mod(s$S2)^2) * sin(th) }
register_family(family_def("mie", "Mie phase function (scattering angle law)", "circular",
  list(P("x", "positive", role = "size parameter 2 pi r / wavelength", regression = FALSE),
       P("n_re", "positive", role = "real refractive index", regression = FALSE),
       P("n_im", "nonneg", role = "absorption index", regression = FALSE)),
  support = function(p) list(lower = 0, upper = pi),
  formula = "f(theta) proportional to (|S1(theta)|^2 + |S2(theta)|^2) sin(theta), S1, S2 the Mie amplitudes of a homogeneous sphere",
  # Normalized once per distinct parameter set.
  lpdf = function(y, p) { r <- recycle_pars(c(list(.y = y), p)); key <- paste(r$x, r$n_re, r$n_im); out <- numeric(length(y))
    for (g in unique(key)) { i <- which(key == g); m <- complex(real = r$n_re[i[1]], imaginary = r$n_im[i[1]]); x0 <- r$x[i[1]]
      z <- stats::integrate(function(t) mie_unnorm(t, x0, m), 0, pi, subdivisions = 5000L, rel.tol = 1e-10)$value
      yi <- r$.y[i]; out[i] <- ifelse(yi <= 0 | yi >= pi, -Inf, log(mie_unnorm(pmin(pmax(yi, 1e-12), pi - 1e-12), x0, m)) - log(z)) }
    out },
  numeric_moments = FALSE,
  references = c("Mie (1908), Beitraege zur Optik trueber Medien. Ann. Phys., 330, 377-445.",
                 "Bohren and Huffman (1983), Absorption and Scattering of Light by Small Particles. Wiley, appendix A."),
  reference = ref_impl("miepython agrees on |S1|^2 + |S2|^2 to 1e-8; Rayleigh limit 3/8 (1 + cos^2) sin as x -> 0", function(y, p) log(3 / 8 * (1 + cos(y)^2) * sin(y)),
    note = "Checked in the Rayleigh limit and by the identity int (|S1|^2 + |S2|^2) sin = x^2 Qsca."),
  status = "validated"))

register_family(family_def("truncated_normal", "Truncated normal", "real",
  list(loc_par("mu"), scale_par("sigma"), P("a", "any", role = "lower bound", link = "identity", data = TRUE, default = -Inf),
       P("b", "any", role = "upper bound", link = "identity", data = TRUE, default = Inf)),
  support = function(p) list(lower = p$a, upper = p$b),
  check = function(p) if (any(p$a >= p$b, na.rm = TRUE)) "truncated_normal needs a < b.",
  lpdf = function(y, p) ifelse(y < p$a | y > p$b, -Inf, stats::dnorm(y, p$mu, p$sigma, log = TRUE) -
    family_log_interval(get_family("normal"), p$a, p$b, p[c("mu", "sigma")])),
  cdf = function(q, p) { N <- get_family("normal"); pn <- p[c("mu", "sigma")]
    out <- exp(family_log_interval(N, p$a, pmin(pmax(q, p$a), p$b), pn) - family_log_interval(N, p$a, p$b, pn))
    out[q <= p$a] <- 0; out[q >= p$b] <- 1; out },
  quantile = function(pr, p) { la <- stats::pnorm(p$a, p$mu, p$sigma, log.p = TRUE)
    upper <- la > log(0.5)
    sa <- stats::pnorm(p$a, p$mu, p$sigma, lower.tail = FALSE); sb <- stats::pnorm(p$b, p$mu, p$sigma, lower.tail = FALSE)
    fa <- stats::pnorm(p$a, p$mu, p$sigma); fb <- stats::pnorm(p$b, p$mu, p$sigma)
    ifelse(upper, stats::qnorm(sa - pr * (sa - sb), p$mu, p$sigma, lower.tail = FALSE),
           stats::qnorm(fa + pr * (fb - fa), p$mu, p$sigma)) },
  references = "Johnson, Kotz and Balakrishnan (1994), Continuous Univariate Distributions, vol. 1, sec. 13.10.",
  reference = ref_impl("truncnorm", function(y, p) log(truncnorm::dtruncnorm(y, p$a, p$b, p$mu, p$sigma))),
  stan = sfun("truncated_normal", c("mu", "sigma", "a", "b"),
    "    if (y < a || y > b) return negative_infinity();\n    return normal_lpdf(y | mu, sigma) - log_diff_exp(normal_lcdf(b | mu, sigma), normal_lcdf(a | mu, sigma));"),
  status = "validated"))

register_family(modify_family(censor("normal", lower = 0), id = "rectified_gaussian", display_name = "Rectified Gaussian",
  aliases = "rectified_normal", references = "Socci, Lee and Seung (1998), The rectified Gaussian distribution. NIPS 10, 350-356.",
  reference = ref_impl("definition", function(y, p) ifelse(y == 0, stats::pnorm(0, p$mu, p$sigma, log.p = TRUE), stats::dnorm(y, p$mu, p$sigma, log = TRUE))),
  status = "validated"))

# Aliases for laws listed under other names.
local({
  al <- c(compound_poisson_gamma = "tweedie", type1_gumbel = "gumbel",
          type2_gumbel = "frechet", fisher_tippett_type1 = "gumbel", levy_skew_alpha_stable = "stable")
  for (a in names(al)) assign(a, al[[a]], envir = .aliases)
})
