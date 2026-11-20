# Discrete distributions with finite support.

fin <- "finite_discrete"
trials_par <- function(name = "size") P(name, "nonneg_integer", role = "trials", data = TRUE,
                                         interpretation = "Number of trials (data)")

register_family(family_def(
  "bernoulli", "Bernoulli", fin, discrete = TRUE, aliases = "BI1",
  pars = list(P("prob", "unit", role = "probability", is_mean = TRUE)),
  formula = "p(y) = p^y (1 - p)^(1 - y), y in {0, 1}",
  mean_formula = "p", variance_formula = "p (1 - p)",
  lpdf = function(y, p) ifelse(y == 1, log(p$prob), log1p(-p$prob)),
  cdf = function(q, p) ifelse(q < 0, 0, ifelse(q < 1, 1 - p$prob, 1)),
  quantile = function(pr, p) as.numeric(pr > 1 - p$prob),
  rng = function(n, p) stats::rbinom(n, 1, p$prob),
  mean = function(p) p$prob, variance = function(p) p$prob * (1 - p$prob),
  support = function(p) list(lower = 0, upper = 1),
  references = JKK, reference = ref_impl("extraDistr", function(y, p) extraDistr::dbern(y, p$prob, log = TRUE)),
  stan = stan_native("bernoulli_lpmf({y} | {prob})", rng = "bernoulli_rng({prob})"),
  bugs = list(native = "dbern({prob})", loglik = "{y} * log({prob}) + (1 - {y}) * log(1 - {prob})"),
  start = function(y) list(prob = min(max(mean(y), 0.01), 0.99)), status = "validated"))

register_family(family_def(
  "rademacher", "Rademacher (P(Y = 1) = p, P(Y = -1) = 1 - p; p = 1/2 is the classical law)", fin, discrete = TRUE,
  pars = list(P("prob", "unit", role = "probability", default = 0.5)),
  formula = "P(Y = 1) = p, P(Y = -1) = 1 - p", mean_formula = "2p - 1", variance_formula = "4 p (1 - p)",
  lpdf = function(y, p) ifelse(y == 1, log(p$prob), log1p(-p$prob)),
  response_check = function(y, p) y %in% c(-1, 1),
  cdf = function(q, p) ifelse(q < -1, 0, ifelse(q < 1, 1 - p$prob, 1)),
  quantile = function(pr, p) ifelse(pr > 1 - p$prob, 1, -1),
  mean = function(p) 2 * p$prob - 1, variance = function(p) 4 * p$prob * (1 - p$prob),
  support = function(p) list(lower = -1, upper = 1),
  references = "Montgomery-Smith (1990), The distribution of Rademacher sums. Proc. AMS, 109, 517-522.",
  reference = ref_impl("closed form", function(y, p) log(ifelse(y == 1, p$prob, 1 - p$prob))),
  stan = sfun("rademacher", "prob", "    return y == 1 ? log(prob) : log1m(prob);", discrete = TRUE),
  bugs = list(loglik = "equals({y}, 1) * log({prob}) + equals({y}, -1) * log(1 - {prob})"),
  status = "validated"))

register_family(family_def(
  "binomial", "Binomial", fin, discrete = TRUE, aliases = c("BI", "binom"),
  pars = list(trials_par(), P("prob", "unit", role = "probability", interpretation = "Success probability")),
  formula = "p(y) = choose(n, y) p^y (1 - p)^(n - y)",
  mean_formula = "n p", variance_formula = "n p (1 - p)",
  lpdf = function(y, p) stats::dbinom(y, p$size, p$prob, log = TRUE),
  cdf = function(q, p) stats::pbinom(q, p$size, p$prob),
  quantile = function(pr, p) stats::qbinom(pr, p$size, p$prob),
  rng = function(n, p) stats::rbinom(n, p$size, p$prob),
  mean = function(p) p$size * p$prob, variance = function(p) p$size * p$prob * (1 - p$prob),
  support = function(p) list(lower = 0, upper = p$size),
  regression_notes = "The logit of prob is the usual logistic regression predictor; size is data.",
  references = JKK, reference = ref_impl("stats", function(y, p) stats::dbinom(y, p$size, p$prob, log = TRUE),
                                         function(q, p) stats::pbinom(q, p$size, p$prob)),
  stan = stan_native("binomial_lpmf({y} | {size}, {prob})", "binomial_lcdf({y} | {size}, {prob})",
                     "binomial_lccdf({y} | {size}, {prob})", "binomial_rng({size}, {prob})"),
  bugs = list(native = "dbin({prob}, {size})",
              loglik = "logfact({size}) - logfact({y}) - logfact({size} - {y}) + {y} * log({prob}) + ({size} - {y}) * log(1 - {prob})"),
  tmb = "binomial", start = function(y) list(prob = 0.5), status = "validated"))

register_family(family_def(
  "beta_binomial", "Beta binomial (mean and precision)", fin, discrete = TRUE, aliases = c("BB", "polya_eggenberger"),
  pars = list(trials_par(), P("mu", "unit", role = "location", interpretation = "Mean success probability"),
              P("phi", "positive", role = "precision", interpretation = "alpha + beta")),
  formula = "p(y) = choose(n, y) B(y + mu phi, n - y + (1 - mu) phi) / B(mu phi, (1 - mu) phi)",
  mean_formula = "n mu", variance_formula = "n mu (1 - mu) (phi + n) / (phi + 1)",
  regression_notes = "Also the Polya-Eggenberger (Polya urn) distribution with alpha = mu phi, beta = (1 - mu) phi.",
  lpdf = function(y, p) { a <- p$mu * p$phi; b <- (1 - p$mu) * p$phi
    lchoose(p$size, y) + lbeta(y + a, p$size - y + b) - lbeta(a, b) },
  rng = function(n, p) stats::rbinom(n, p$size, stats::rbeta(n, p$mu * p$phi, (1 - p$mu) * p$phi)),
  mean = function(p) p$size * p$mu,
  variance = function(p) p$size * p$mu * (1 - p$mu) * (p$phi + p$size) / (p$phi + 1),
  support = function(p) list(lower = 0, upper = p$size),
  references = c(JKK, "Skellam (1948), JRSS B, 10, 257-261."),
  reference = ref_impl("extraDistr", function(y, p) extraDistr::dbbinom(y, p$size, p$mu * p$phi, (1 - p$mu) * p$phi, log = TRUE)),
  stan = stan_native("beta_binomial_lpmf({y} | {size}, {mu} * {phi}, (1 - {mu}) * {phi})",
                     "beta_binomial_lcdf({y} | {size}, {mu} * {phi}, (1 - {mu}) * {phi})",
                     "beta_binomial_lccdf({y} | {size}, {mu} * {phi}, (1 - {mu}) * {phi})",
                     "beta_binomial_rng({size}, {mu} * {phi}, (1 - {mu}) * {phi})"),
  bugs = list(loglik = "logfact({size}) - logfact({y}) - logfact({size} - {y}) + loggam({y} + {mu} * {phi}) + loggam({size} - {y} + (1 - {mu}) * {phi}) - loggam({size} + {phi}) - loggam({mu} * {phi}) - loggam((1 - {mu}) * {phi}) + loggam({phi})"),
  start = function(y) list(mu = 0.5, phi = 5), status = "validated"))

register_family(family_def(
  "degenerate", "Degenerate (point mass, Dirac delta)", fin, discrete = TRUE, aliases = c("dirac", "dirac_delta", "point_mass"),
  pars = list(P("x0", "real", role = "location", is_mean = TRUE, regression = FALSE)),
  formula = "P(Y = x0) = 1", mean_formula = "x0", variance_formula = "0",
  lpdf = function(y, p) ifelse(y == p$x0, 0, -Inf),
  cdf = function(q, p) as.numeric(q >= p$x0), quantile = function(pr, p) p$x0 + 0 * pr,
  rng = function(n, p) rep_len(p$x0, n), mean = function(p) p$x0, variance = function(p) 0 * p$x0,
  support = function(p) list(lower = p$x0, upper = p$x0), numeric_moments = FALSE,
  references = "Feller (1971), An Introduction to Probability Theory, vol. 2.",
  reference = ref_impl("closed form", function(y, p) ifelse(y == p$x0, 0, -Inf)),
  extra = list(grid = 0), status = "validated"))

register_family(family_def(
  "discrete_uniform", "Discrete uniform on {a, ..., b}", fin, discrete = TRUE, aliases = "dunif",
  pars = list(P("a", "integer", role = "lower", regression = FALSE), P("b", "integer", role = "upper", regression = FALSE)),
  check = function(p) if (any(p$b < p$a, na.rm = TRUE)) "b must be at least a.",
  formula = "p(y) = 1 / (b - a + 1)", mean_formula = "(a + b) / 2", variance_formula = "((b - a + 1)^2 - 1) / 12",
  lpdf = function(y, p) -log(p$b - p$a + 1),
  cdf = function(q, p) pmin(pmax((floor(q) - p$a + 1) / (p$b - p$a + 1), 0), 1),
  quantile = function(pr, p) p$a + pmax(ceiling(pr * (p$b - p$a + 1)) - 1, 0),
  mean = function(p) (p$a + p$b) / 2, variance = function(p) ((p$b - p$a + 1)^2 - 1) / 12,
  support = function(p) list(lower = p$a, upper = p$b),
  references = JKK, reference = ref_impl("extraDistr", function(y, p) extraDistr::ddunif(y, p$a, p$b, log = TRUE)),
  stan = sfun("discrete_uniform", c("a", "b"), "    return -log(b - a + 1);", discrete = TRUE),
  bugs = list(loglik = "-log({b} - {a} + 1)"), status = "validated"))

register_family(family_def(
  "hypergeometric", "Hypergeometric", fin, discrete = TRUE,
  pars = list(P("m", "nonneg_integer", role = "successes", data = TRUE, interpretation = "Success states in the population"),
              P("n", "nonneg_integer", role = "failures", data = TRUE, interpretation = "Failure states in the population"),
              P("k", "nonneg_integer", role = "draws", data = TRUE, interpretation = "Number of draws")),
  formula = "p(y) = choose(m, y) choose(n, k - y) / choose(m + n, k)",
  mean_formula = "k m / (m + n)", variance_formula = "k m n (m + n - k) / ((m + n)^2 (m + n - 1))",
  lpdf = function(y, p) stats::dhyper(y, p$m, p$n, p$k, log = TRUE),
  cdf = function(q, p) stats::phyper(q, p$m, p$n, p$k), quantile = function(pr, p) stats::qhyper(pr, p$m, p$n, p$k),
  rng = function(n, p) stats::rhyper(n, p$m, p$n, p$k),
  mean = function(p) p$k * p$m / (p$m + p$n),
  variance = function(p) { N <- p$m + p$n; p$k * p$m * p$n * (N - p$k) / (N^2 * (N - 1)) },
  support = function(p) list(lower = pmax(0, p$k - p$n), upper = pmin(p$k, p$m)),
  regression_notes = "All arguments are data; use fisher_noncentral_hypergeometric or wallenius_noncentral_hypergeometric to regress an odds ratio.",
  references = JKK, reference = ref_impl("stats", function(y, p) stats::dhyper(y, p$m, p$n, p$k, log = TRUE)),
  stan = stan_native("hypergeometric_lpmf({y} | {k}, {m}, {n})", rng = "hypergeometric_rng({k}, {m}, {n})"),
  status = "validated"))

register_family(family_def(
  "negative_hypergeometric", "Negative hypergeometric (failures before the r-th success, no replacement)", fin, discrete = TRUE,
  pars = list(P("N", "positive_integer", role = "population", data = TRUE),
              P("K", "nonneg_integer", role = "successes", data = TRUE),
              P("r", "positive_integer", role = "target", data = TRUE)),
  formula = "p(y) = choose(y + r - 1, y) choose(N - r - y, K - r) / choose(N, K), y failures before r successes drawn from N items with K successes",
  mean_formula = "r (N - K) / (K + 1)",
  lpdf = function(y, p) lchoose(y + p$r - 1, y) + lchoose(p$N - p$r - y, p$K - p$r) - lchoose(p$N, p$K),
  mean = function(p) p$r * (p$N - p$K) / (p$K + 1),
  variance = function(p) p$r * (p$N + 1) * (p$N - p$K) / ((p$K + 1) * (p$K + 2)) * (1 - p$r / (p$K + 1)),
  support = function(p) list(lower = 0, upper = p$N - p$K),
  check = function(p) if (any(p$K > p$N | p$r > p$K, na.rm = TRUE)) "Need r <= K <= N.",
  references = JKK, reference = ref_impl("extraDistr", function(y, p) extraDistr::dnhyper(y + p$r, p$N - p$K, p$K, p$r, log = TRUE), note = "extraDistr counts total draws, y + r."),
  stan = sfun("negative_hypergeometric", c("N", "K", "r"),
              "    return lchoose(y + r - 1, y) + lchoose(N - r - y, K - r) - lchoose(N, K);", discrete = TRUE),
  status = "validated"))

register_family(family_def(
  "poisson_binomial", "Poisson binomial", fin, discrete = TRUE, variate_type = "scalar_vector_parameter",
  pars = list(P("probs", "unit", role = "probability", regression = FALSE, interpretation = "Vector of success probabilities")),
  formula = "Y = sum of independent Bernoulli(p_i); pmf by the discrete Fourier transform or recursion",
  mean_formula = "sum p_i", variance_formula = "sum p_i (1 - p_i)",
  lpdf = function(y, p) {
    pr <- p$probs; f <- c(1, numeric(length(pr)))
    for (q in pr) f <- c(f[1] * (1 - q), f[-1] * (1 - q) + f[-length(f)] * q)
    ifelse(y >= 0 & y <= length(pr) & y == round(y), log(f[pmin(pmax(y, 0), length(pr)) + 1]), -Inf)
  },
  rng = function(n, p) colSums(matrix(stats::runif(n * length(p$probs)) < p$probs, ncol = n)),
  references = c("Hong (2013), On computing the distribution function for the Poisson binomial distribution. CSDA, 59, 41-51."),
  reference = ref_impl("PoissonBinomial", function(y, p) log(PoissonBinomial::dpbinom(y, p$probs))),
  status = "validated"))

fnc_logZ <- function(m1, m2, n, omega) {
  lo <- max(0, n - m2); hi <- min(n, m1); x <- lo:hi
  l <- lchoose(m1, x) + lchoose(m2, n - x) + x * log(omega)
  list(x = x, l = l, lz = log_sum_exp(l))
}
register_family(family_def(
  "fisher_noncentral_hypergeometric", "Fisher's noncentral hypergeometric", fin, discrete = TRUE,
  pars = list(P("m1", "nonneg_integer", role = "type 1 items", data = TRUE), P("m2", "nonneg_integer", role = "type 2 items", data = TRUE),
              P("n", "nonneg_integer", role = "draws", data = TRUE),
              P("omega", "positive", role = "odds ratio", interpretation = "Odds ratio of drawing type 1 against type 2")),
  formula = "p(y) = choose(m1, y) choose(m2, n - y) omega^y / sum_x choose(m1, x) choose(m2, n - x) omega^x",
  lpdf = function(y, p) vapply(seq_along(y), function(i) {
    z <- fnc_logZ(p$m1[i], p$m2[i], p$n[i], p$omega[i])
    lchoose(p$m1[i], y[i]) + lchoose(p$m2[i], p$n[i] - y[i]) + y[i] * log(p$omega[i]) - z$lz }, 0),
  support = function(p) list(lower = pmax(0, p$n - p$m2), upper = pmin(p$n, p$m1)),
  regression_notes = "log(omega) is the conditional log odds ratio; the conditional likelihood of a 2x2 table.",
  references = "Fog (2008), Calculation methods for Wallenius' noncentral hypergeometric distribution. Comm. Stat. Sim. Comp., 37, 258-273.",
  reference = ref_impl("BiasedUrn", function(y, p) log(mapply(BiasedUrn::dFNCHypergeo, y, p$m1, p$m2, p$n, p$omega))),
  stan = sfun("fisher_nchg", c("m1", "m2", "n", "omega"), paste0(
    "    int lo = max(0, n - m2);\n    int hi = min(n, m1);\n",
    "    vector[hi - lo + 1] l;\n    for (x in lo:hi) l[x - lo + 1] = lchoose(m1, x) + lchoose(m2, n - x) + x * log(omega);\n",
    "    return lchoose(m1, y) + lchoose(m2, n - y) + y * log(omega) - log_sum_exp(l);"), discrete = TRUE, ints = c("m1", "m2", "n")),
  status = "validated"))

wnc_pmf1 <- function(y, m1, m2, n, omega) {
  if (y < max(0, n - m2) || y > min(n, m1)) return(-Inf)
  D <- omega * (m1 - y) + (m2 - (n - y))
  f <- function(t) (1 - t^(omega / D))^y * (1 - t^(1 / D))^(n - y)
  v <- stats::integrate(f, 0, 1, rel.tol = 1e-12, subdivisions = 500L)$value
  lchoose(m1, y) + lchoose(m2, n - y) + log(v)
}
register_family(family_def(
  "wallenius_noncentral_hypergeometric", "Wallenius' noncentral hypergeometric", fin, discrete = TRUE,
  pars = list(P("m1", "nonneg_integer", role = "type 1 items", data = TRUE), P("m2", "nonneg_integer", role = "type 2 items", data = TRUE),
              P("n", "nonneg_integer", role = "draws", data = TRUE), P("omega", "positive", role = "weight ratio")),
  formula = "p(y) = choose(m1, y) choose(m2, n - y) int_0^1 (1 - t^(omega/D))^y (1 - t^(1/D))^(n - y) dt, D = omega (m1 - y) + (m2 - n + y)",
  lpdf = function(y, p) vapply(seq_along(y), function(i) wnc_pmf1(y[i], p$m1[i], p$m2[i], p$n[i], p$omega[i]), 0),
  support = function(p) list(lower = pmax(0, p$n - p$m2), upper = pmin(p$n, p$m1)),
  numerical_notes = "The pmf uses adaptive quadrature of Wallenius' integral representation.",
  references = c("Wallenius (1963), Biased Sampling: the Noncentral Hypergeometric Probability Distribution. PhD thesis, Stanford.",
                 "Fog (2008), Comm. Stat. Sim. Comp., 37, 258-273."),
  reference = ref_impl("BiasedUrn", function(y, p) log(mapply(BiasedUrn::dWNCHypergeo, y, p$m1, p$m2, p$n, p$omega, precision = 1e-12))),
  status = "validated"))

register_family(family_def(
  "benford", "Benford's law (first significant digit in base b)", fin, discrete = TRUE, aliases = "benford_law",
  pars = list(P("base", "positive_integer", role = "base", regression = FALSE, default = 10)),
  formula = "P(D = d) = log_b(1 + 1/d), d = 1, ..., b - 1",
  lpdf = function(y, p) log(log1p(1 / y) / log(p$base)),
  cdf = function(q, p) pmin(pmax(log(1 + pmax(floor(q), 0)) / log(p$base), 0), 1),
  quantile = function(pr, p) pmax(ceiling(p$base^pr - 1 - 1e-12), 1),
  support = function(p) list(lower = 1, upper = p$base - 1),
  references = c("Benford (1938), Proc. Am. Philos. Soc., 78, 551-572.", "Hill (1995), Statistical Science, 10, 354-363."),
  reference = ref_impl("VGAM", function(y, p) VGAM::dbenf(y, log = TRUE), note = "VGAM covers base 10."),
  stan = sfun("benford", "base", "    return log(log1p(1.0 / y) / log(base));", discrete = TRUE),
  status = "validated"))

soliton_ideal <- function(N) c(1 / N, 1 / ((2:N) * (2:N - 1)))
soliton_robust <- function(N, c0, delta) {
  R <- c0 * log(N / delta) * sqrt(N); kk <- floor(N / R)
  tau <- numeric(N); i <- seq_len(N)
  if (kk >= 2) tau[seq_len(kk - 1)] <- R / (seq_len(kk - 1) * N)
  if (kk >= 1 && kk <= N) tau[kk] <- R * log(R / delta) / N
  m <- soliton_ideal(N) + tau; m / sum(m)
}
register_family(family_def(
  "soliton_ideal", "Ideal soliton", fin, discrete = TRUE,
  pars = list(P("N", "positive_integer", role = "size", regression = FALSE)),
  formula = "p(1) = 1/N, p(k) = 1 / (k (k - 1)), k = 2, ..., N",
  lpdf = function(y, p) log(ifelse(y == 1, 1 / p$N, 1 / (y * (y - 1)))),
  support = function(p) list(lower = 1, upper = p$N),
  references = "Luby (2002), LT codes. Proc. 43rd IEEE FOCS, 271-280.",
  reference = ref_impl("closed form", function(y, p) log(mapply(function(k, N) soliton_ideal(N)[k], y, p$N))),
  status = "validated"))

register_family(family_def(
  "soliton_robust", "Robust soliton", fin, discrete = TRUE,
  pars = list(P("N", "positive_integer", role = "size", regression = FALSE), P("c", "positive", role = "tuning", regression = FALSE),
              P("delta", "unit", role = "failure probability", regression = FALSE)),
  formula = "mu(k) proportional to rho(k) + tau(k): ideal soliton plus a spike at N/R, R = c log(N/delta) sqrt(N)",
  lpdf = function(y, p) log(mapply(function(k, N, c0, d) soliton_robust(N, c0, d)[k], y, p$N, p$c, p$delta)),
  support = function(p) list(lower = 1, upper = p$N),
  references = "Luby (2002), LT codes. Proc. 43rd IEEE FOCS, 271-280.",
  reference = ref_impl("definition", function(y, p) log(mapply(function(k, N, c0, d) { R <- c0 * log(N / d) * sqrt(N)
    kk <- floor(N / R); t <- vapply(1:N, function(i) if (i < kk) R / (i * N) else if (i == kk) R * log(R / d) / N else 0, 0)
    r <- soliton_ideal(N); ((r + t) / sum(r + t))[k] }, y, p$N, p$c, p$delta))),
  status = "validated"))

gen_harmonic <- function(N, s, q = 0) vapply(seq_along(N), function(i) sum((seq_len(N[i]) + q[i])^-s[i]), 0)
register_family(family_def(
  "zipf_finite", "Zipf (finite N)", fin, discrete = TRUE, aliases = "zipf_law",
  pars = list(P("N", "positive_integer", role = "size", data = TRUE), P("s", "positive", role = "exponent")),
  formula = "p(k) = k^(-s) / H(N, s), k = 1, ..., N",
  lpdf = function(y, p) { n <- length(y); N <- rep_len(p$N, n); s <- rep_len(p$s, n)
    -s * log(y) - log(gen_harmonic(N, s, rep(0, n))) },
  support = function(p) list(lower = 1, upper = p$N),
  references = "Zipf (1949), Human Behavior and the Principle of Least Effort.",
  reference = ref_impl("VGAM", function(y, p) VGAM::dzipf(y, p$N, p$s, log = TRUE)),
  stan = sfun("zipf_finite", c("N", "s"), "    real h = 0;\n    for (k in 1:N) h += pow(k, -s);\n    return -s * log(y) - log(h);", discrete = TRUE, ints = "N"),
  start = function(y) list(s = 1), status = "validated"))

register_family(family_def(
  "zipf_mandelbrot", "Zipf-Mandelbrot", fin, discrete = TRUE,
  pars = list(P("N", "positive_integer", role = "size", data = TRUE), P("s", "positive", role = "exponent"),
              P("q", "nonneg", role = "shift")),
  formula = "p(k) = (k + q)^(-s) / H(N, q, s), k = 1, ..., N",
  lpdf = function(y, p) { n <- length(y); -p$s * log(y + p$q) - log(gen_harmonic(rep_len(p$N, n), rep_len(p$s, n), rep_len(p$q, n))) },
  support = function(p) list(lower = 1, upper = p$N),
  references = "Mandelbrot (1953), An informational theory of the statistical structure of language.",
  reference = ref_impl("definition", function(y, p) log(mapply(function(k, N, s, q) ((1:N + q)^-s / sum((1:N + q)^-s))[k], y, p$N, p$s, p$q))),
  stan = sfun("zipf_mandelbrot", c("N", "s", "q"), "    real h = 0;\n    for (k in 1:N) h += pow(k + q, -s);\n    return -s * log(y + q) - log(h);", discrete = TRUE, ints = "N"),
  start = function(y) list(s = 1, q = 1), status = "validated"))

register_family(family_def(
  "categorical", "Categorical", "categorical", discrete = TRUE, variate_type = "categorical",
  pars = list(P("probs", "unit", role = "probabilities", regression = FALSE,
                interpretation = "Vector (or matrix, one row per observation) of category probabilities summing to one")),
  formula = "P(Y = k) = p_k, k = 1, ..., K",
  lpdf = function(y, p) { P0 <- if (is.matrix(p$probs)) p$probs else matrix(p$probs, nrow = length(y), ncol = length(p$probs), byrow = TRUE)
    ok <- y >= 1 & y <= ncol(P0) & y == round(y); out <- rep(-Inf, length(y))
    out[ok] <- log(P0[cbind(which(ok), as.integer(y[ok]))]); out },
  rng = function(n, p) sample.int(length(p$probs), n, replace = TRUE, prob = p$probs),
  check = function(p) if (!is.matrix(p$probs) && abs(sum(p$probs) - 1) > 1e-8) "Category probabilities must sum to one.",
  regression_notes = "Regression on category probabilities uses the multinomial logit family 'multinomial_logit' (see the categorical and ordinal families).",
  references = "Agresti (2013), Categorical Data Analysis, 3rd ed. Wiley.",
  reference = ref_impl("extraDistr", function(y, p) extraDistr::dcat(y, p$probs, log = TRUE)),
  status = "validated"))
