# Generated families: F(y) = T(G(y)) for a baseline CDF G and a CDF T on
# [0, 1]. Covers the beta-G, Kumaraswamy-G, Marshall-Olkin-G, transmuted-G,
# power series and many other "new distribution" constructions.

gen_family <- function(base, tag, display, newpars, T, logt, Tinv = NULL, formula = "", refs = character()) {
  base <- get_family(base)
  reject_atoms(base, "A CDF generator")
  if (base$variate_type != "scalar") stop_unsupported("A generator on a vector valued law", base$id)
  wanted <- vapply(newpars, `[[`, "", "name")
  actual <- character(0)
  for (k in seq_along(newpars)) {
    nm <- new_par_name(list(pars = c(base$pars, stats::setNames(vector("list", length(actual)), actual))), wanted[k])
    newpars[[k]]$name <- nm; actual <- c(actual, nm)
  }
  bn <- names(base$pars)
  gp <- function(p) stats::setNames(p[actual], wanted)
  bp <- function(p) p[bn]
  G <- function(y, p) pmin(pmax(family_cdf(base, y, bp(p)), 0), 1)
  lpdf <- if (base$discrete) {
    function(y, p) log(pmax(T(G(y, p), gp(p)) - T(G(y - 1, p), gp(p)), 0))
  } else {
    function(y, p) {
      lb <- eval_lpdf(base, y, bp(p))
      u <- pmin(pmax(G(y, p), 1e-300), 1 - .Machine$double.eps / 2)
      ifelse(lb == -Inf, -Inf, logt(u, gp(p)) + lb)
    }
  }
  # Every generator maps the endpoints 0 and 1 to themselves; set them
  # exactly, as the formulas give 0 log 0 or Inf / Inf there.
  ends <- function(f) function(u, g) { v <- f(u, g); v[!is.na(u) & u <= 0] <- 0; v[!is.na(u) & u >= 1] <- 1; v }
  cdf <- function(q, p) ends(T)(G(q, p), gp(p))
  quantile <- if (!is.null(Tinv)) function(pr, p) family_quantile(base, pmin(pmax(ends(Tinv)(pr, gp(p)), 0), 1), bp(p)) else NULL
  rng <- if (!is.null(Tinv)) function(n, p) family_quantile(base, ends(Tinv)(stats::runif(n), gp(p)), bp(p)) else NULL
  f <- op_family(base, sprintf("%s(%s)", tag, base$id), sprintf("%s %s", display, base$display_name),
                 c(base$pars, stats::setNames(newpars, actual)), lpdf, cdf = cdf, quantile = quantile, rng = rng,
                 operator = list(type = "generator", generator = tag, base = base),
                 formula = paste0("F(y) = T(G(y)); ", formula))
  f$references <- c(refs, base$references)
  f
}

#' Generated distribution families
#'
#' Build a new family from a baseline family `G` through a transformation of
#' its CDF, `F(y) = T(G(y))`, where `T` is a CDF on `[0, 1]`. The density is
#' `f(y) = t(G(y)) g(y)` for a continuous baseline and
#' `P(Y = y) = T(G(y)) - T(G(y - 1))` for a discrete one. New shape
#' parameters are appended to the baseline parameters (renamed with a
#' trailing underscore on a clash). All of them can take regression
#' predictors. Quantiles use the baseline quantile function where `T` has a
#' closed form inverse.
#'
#' * `exponentiated()`: `T(u) = u^a` (Lehmann alternative, exponentiated-G).
#' * `kumaraswamy_g()`: `T(u) = 1 - (1 - u^a)^b` (Cordeiro and de Castro 2011).
#' * `beta_g()`: `T(u) = I_u(a, b)` (Eugene, Lee and Famoye 2002).
#' * `mcdonald_g()`: `T(u) = I_{u^c}(a, b)` (generalized beta-G, GB1-G).
#' * `exp_generalized_g()`: `T(u) = (1 - (1 - u)^a)^b` (Cordeiro, Ortega and da Cunha 2013).
#' * `marshall_olkin()`: `T(u) = u / (alpha + (1 - alpha) u)` (Marshall and Olkin 1997).
#' * `transmuted()`: `T(u) = (1 + lambda) u - lambda u^2`, `|lambda| <= 1` (Shaw and Buckley 2007).
#' * `record_transmuted()`: `T(u) = u + p (1 - u) log(1 - u)` (Balakrishnan and He 2021).
#' * `odd_log_logistic()`: `T(u) = u^a / (u^a + (1 - u)^a)` (Gleaton and Lynch 2006).
#' * `odd_burr()`: `T(u) = 1 - (1 + (u / (1 - u))^a)^(-b)` (Alizadeh et al. 2017).
#' * `weibull_g()`: `T(u) = 1 - exp(-a (u / (1 - u))^b)` (Bourguignon, Silva and Cordeiro 2014).
#' * `gamma_g()`: Zografos and Balakrishnan (2009) `T(u) = P(a, -log(1 - u))`, or
#'   Ristic and Balakrishnan (2012) `T(u) = 1 - P(a, -log u)`.
#' * `topp_leone_g()`: `T(u) = (u (2 - u))^a`.
#' * `lindley_g()`: `T(u) = u^theta (1 + theta - theta log u) / (1 + theta)`.
#' * `harris_g()`: survival `(theta S^k / (1 - (1 - theta) S^k))^(1/k)` (Aly and Benkherouf 2011).
#' * `sine_g()`: `T(u) = sin(pi u / 2)` (Kumar, Singh and Singh 2015).
#' * `dus()`: `T(u) = (e^u - 1) / (e - 1)` (Kumar, Singh and Singh 2015).
#' * `km()`: `T(u) = e (1 - e^-u) / (e - 1)` (Kavya and Manoharan 2021).
#' * `power_series_max()` and `power_series_min()`: the maximum or minimum of a
#'   zero truncated power series number of baseline draws, with series
#'   `"geometric"`, `"poisson"`, `"logarithmic"` or `"binomial"`;
#'   `T(u) = C(theta u) / C(theta)` and `T(u) = 1 - C(theta (1 - u)) / C(theta)`.
#'
#' Generators nest, so for example `beta_g(kumaraswamy_g("weibull"))` gives
#' the beta Kumaraswamy Weibull law.
#'
#' @param base Baseline family (identifier or family object).
#' @param type For `gamma_g()`: `"zografos"` or `"ristic"`.
#' @param series Power series type.
#' @param size Binomial series size.
#' @return A `regcomp_family`.
#' @references
#' Eugene, Lee and Famoye (2002), Beta-normal distribution and its applications. Comm. Stat. Theory Methods, 31, 497-512.
#' Cordeiro and de Castro (2011), A new family of generalized distributions. J. Stat. Comput. Simul., 81, 883-898.
#' Marshall and Olkin (1997), A new method for adding a parameter to a family of distributions. Biometrika, 84, 641-652.
#' Shaw and Buckley (2007), The alchemy of probability distributions. arXiv 0901.0434.
#' Tahir and Cordeiro (2016), Compounding of distributions: a survey and new generalized classes. J. Stat. Distrib. Appl., 3, 13.
#' @export
#' @examples
#' bw <- beta_g("weibull")
#' d <- distribution(bw, shape = 1.5, scale = 2, a = 2, b = 0.5)
#' cdf(d, 1)
#' kwn <- kumaraswamy_g("normal")
#' mow <- marshall_olkin("weibull")
#' ewp <- power_series_max(exponentiated("weibull"), "poisson")
exponentiated <- function(base) gen_family(base, "exponentiated", "Exponentiated",
  list(P("a", "positive", role = "power")),
  T = function(u, g) u^g$a, logt = function(u, g) log(g$a) + (g$a - 1) * log(u),
  Tinv = function(v, g) v^(1 / g$a), formula = "T(u) = u^a",
  refs = "Mudholkar and Srivastava (1993); Gupta, Gupta and Gupta (1998), Comm. Stat. Theory Methods, 27, 887-904.")

#' @rdname exponentiated
#' @export
kumaraswamy_g <- function(base) gen_family(base, "kumaraswamy_g", "Kumaraswamy",
  list(P("a", "positive"), P("b", "positive")),
  T = function(u, g) 1 - (1 - u^g$a)^g$b,
  logt = function(u, g) log(g$a) + log(g$b) + (g$a - 1) * log(u) + (g$b - 1) * log1p(-u^g$a),
  Tinv = function(v, g) (1 - (1 - v)^(1 / g$b))^(1 / g$a), formula = "T(u) = 1 - (1 - u^a)^b",
  refs = "Cordeiro and de Castro (2011), J. Stat. Comput. Simul., 81, 883-898.")

#' @rdname exponentiated
#' @export
beta_g <- function(base) gen_family(base, "beta_g", "Beta",
  list(P("a", "positive"), P("b", "positive")),
  T = function(u, g) stats::pbeta(u, g$a, g$b), logt = function(u, g) stats::dbeta(u, g$a, g$b, log = TRUE),
  Tinv = function(v, g) stats::qbeta(v, g$a, g$b), formula = "T(u) = I_u(a, b)",
  refs = "Eugene, Lee and Famoye (2002), Comm. Stat. Theory Methods, 31, 497-512.")

#' @rdname exponentiated
#' @export
mcdonald_g <- function(base) gen_family(base, "mcdonald_g", "McDonald (GB1)",
  list(P("a", "positive"), P("b", "positive"), P("c", "positive")),
  T = function(u, g) stats::pbeta(u^g$c, g$a, g$b),
  logt = function(u, g) log(g$c) + (g$a * g$c - 1) * log(u) + (g$b - 1) * log1p(-u^g$c) - lbeta(g$a, g$b),
  Tinv = function(v, g) stats::qbeta(v, g$a, g$b)^(1 / g$c), formula = "T(u) = I_{u^c}(a, b)",
  refs = "Alexander, Cordeiro, Ortega and Sarabia (2012), Generalized beta-generated distributions. CSDA, 56, 1880-1897.")

#' @rdname exponentiated
#' @export
exp_generalized_g <- function(base) gen_family(base, "exp_generalized_g", "Exponentiated generalized",
  list(P("a", "positive"), P("b", "positive")),
  T = function(u, g) (1 - (1 - u)^g$a)^g$b,
  logt = function(u, g) log(g$a) + log(g$b) + (g$a - 1) * log1p(-u) + (g$b - 1) * log1p(-(1 - u)^g$a),
  Tinv = function(v, g) 1 - (1 - v^(1 / g$b))^(1 / g$a), formula = "T(u) = (1 - (1 - u)^a)^b",
  refs = "Cordeiro, Ortega and da Cunha (2013), The exponentiated generalized class of distributions. J. Data Science, 11, 1-27.")

#' @rdname exponentiated
#' @export
marshall_olkin <- function(base) gen_family(base, "marshall_olkin", "Marshall-Olkin",
  list(P("alpha", "positive", role = "tilt")),
  T = function(u, g) u / (g$alpha + (1 - g$alpha) * u),
  logt = function(u, g) log(g$alpha) - 2 * log(g$alpha + (1 - g$alpha) * u),
  Tinv = function(v, g) g$alpha * v / (1 - (1 - g$alpha) * v), formula = "T(u) = u / (alpha + (1 - alpha) u)",
  refs = "Marshall and Olkin (1997), Biometrika, 84, 641-652.")

#' @rdname exponentiated
#' @export
transmuted <- function(base) gen_family(base, "transmuted", "Transmuted",
  list(P("lambda", "symmetric_unit", role = "transmutation")),
  T = function(u, g) (1 + g$lambda) * u - g$lambda * u^2,
  logt = function(u, g) log(1 + g$lambda - 2 * g$lambda * u),
  Tinv = function(v, g) ifelse(abs(g$lambda) < 1e-10, v, ((1 + g$lambda) - sqrt((1 + g$lambda)^2 - 4 * g$lambda * v)) / (2 * g$lambda)),
  formula = "T(u) = (1 + lambda) u - lambda u^2",
  refs = "Shaw and Buckley (2007), arXiv 0901.0434.")

#' @rdname exponentiated
#' @export
record_transmuted <- function(base) gen_family(base, "record_transmuted", "Record-based transmuted",
  list(P("prob", "unit", role = "transmutation")),
  T = function(u, g) u + g$prob * (1 - u) * log1p(-u),
  logt = function(u, g) log(1 - g$prob - g$prob * log1p(-u)), formula = "T(u) = u + p (1 - u) log(1 - u)",
  refs = "Balakrishnan and He (2021), A record-based transmuted family of distributions. In Advances in Statistics, Springer.")

#' @rdname exponentiated
#' @export
odd_log_logistic <- function(base) gen_family(base, "odd_log_logistic", "Odd log-logistic",
  list(P("a", "positive")),
  T = function(u, g) u^g$a / (u^g$a + (1 - u)^g$a),
  logt = function(u, g) log(g$a) + (g$a - 1) * (log(u) + log1p(-u)) - 2 * log(u^g$a + (1 - u)^g$a),
  Tinv = function(v, g) { r <- (v / (1 - v))^(1 / g$a); r / (1 + r) }, formula = "T(u) = u^a / (u^a + (1 - u)^a)",
  refs = "Gleaton and Lynch (2006), Properties of generalized log-logistic families of lifetime distributions. J. Probability and Statistical Science, 4, 51-64.")

#' @rdname exponentiated
#' @export
odd_burr <- function(base) gen_family(base, "odd_burr", "Odd Burr XII",
  list(P("a", "positive"), P("b", "positive")),
  T = function(u, g) 1 - (1 + (u / (1 - u))^g$a)^(-g$b),
  logt = function(u, g) log(g$a) + log(g$b) + (g$a - 1) * log(u) - (g$a + 1) * log1p(-u) - (g$b + 1) * log1p((u / (1 - u))^g$a),
  Tinv = function(v, g) { r <- ((1 - v)^(-1 / g$b) - 1)^(1 / g$a); r / (1 + r) },
  formula = "T(u) = 1 - (1 + (u/(1 - u))^a)^(-b)",
  refs = "Alizadeh, Cordeiro, Nascimento, Lima and Ortega (2017), Odd-Burr generalized family of distributions with some applications. J. Stat. Comput. Simul., 87, 367-389.")

#' @rdname exponentiated
#' @export
weibull_g <- function(base) gen_family(base, "weibull_g", "Weibull",
  list(P("a", "positive"), P("b", "positive")),
  T = function(u, g) 1 - exp(-g$a * (u / (1 - u))^g$b),
  logt = function(u, g) log(g$a) + log(g$b) + (g$b - 1) * log(u) - (g$b + 1) * log1p(-u) - g$a * (u / (1 - u))^g$b,
  Tinv = function(v, g) { r <- (-log1p(-v) / g$a)^(1 / g$b); r / (1 + r) },
  formula = "T(u) = 1 - exp(-a (u/(1 - u))^b)",
  refs = "Bourguignon, Silva and Cordeiro (2014), The Weibull-G family of probability distributions. J. Data Science, 12, 53-68.")

#' @rdname exponentiated
#' @export
gamma_g <- function(base, type = c("zografos", "ristic")) {
  type <- match.arg(type)
  if (type == "zografos") gen_family(base, "gamma_g", "Gamma (Zografos-Balakrishnan)", list(P("a", "positive")),
    T = function(u, g) stats::pgamma(-log1p(-u), g$a), logt = function(u, g) stats::dgamma(-log1p(-u), g$a, log = TRUE) - log1p(-u),
    Tinv = function(v, g) -expm1(-stats::qgamma(v, g$a)), formula = "T(u) = P(a, -log(1 - u))",
    refs = "Zografos and Balakrishnan (2009), On families of beta- and generalized gamma-generated distributions. Statistical Methodology, 6, 344-362.")
  else gen_family(base, "gamma_g_ristic", "Gamma (Ristic-Balakrishnan)", list(P("a", "positive")),
    T = function(u, g) stats::pgamma(-log(u), g$a, lower.tail = FALSE), logt = function(u, g) stats::dgamma(-log(u), g$a, log = TRUE) - log(u),
    Tinv = function(v, g) exp(-stats::qgamma(v, g$a, lower.tail = FALSE)), formula = "T(u) = 1 - P(a, -log u)",
    refs = "Ristic and Balakrishnan (2012), The gamma-exponentiated exponential distribution. J. Stat. Comput. Simul., 82, 1191-1206.")
}

#' @rdname exponentiated
#' @export
topp_leone_g <- function(base) gen_family(base, "topp_leone_g", "Topp-Leone",
  list(P("a", "positive")),
  T = function(u, g) (u * (2 - u))^g$a, logt = function(u, g) log(2 * g$a) + log1p(-u) + (g$a - 1) * log(u * (2 - u)),
  Tinv = function(v, g) 1 - sqrt(1 - v^(1 / g$a)), formula = "T(u) = (u (2 - u))^a",
  refs = "Al-Shomrani, Arif, Shawky, Hanif and Shahbaz (2016), Topp-Leone family of distributions. Pakistan J. Statistics, 32, 1-20.")

#' @rdname exponentiated
#' @export
lindley_g <- function(base) gen_family(base, "lindley_g", "Lindley",
  list(P("theta", "positive")),
  T = function(u, g) u^g$theta * (1 + g$theta - g$theta * log(u)) / (1 + g$theta),
  logt = function(u, g) 2 * log(g$theta) - log1p(g$theta) + log1p(-log(u)) + (g$theta - 1) * log(u),
  formula = "T(u) = u^theta (1 + theta - theta log u) / (1 + theta)",
  refs = "Cakmakyapan and Ozel (2017), The Lindley family of distributions. Hacettepe J. Math. Stat., 46, 1113-1137.")

#' @rdname exponentiated
#' @export
harris_g <- function(base) gen_family(base, "harris_g", "Harris extended",
  list(P("theta", "positive"), P("k", "positive")),
  T = function(u, g) { S <- 1 - u; 1 - (g$theta * S^g$k / (1 - (1 - g$theta) * S^g$k))^(1 / g$k) },
  logt = function(u, g) { w <- (1 - u)^g$k; log(g$theta) / g$k - (1 / g$k + 1) * log1p(-(1 - g$theta) * w) },
  Tinv = function(v, g) { s <- (1 - v)^g$k; 1 - (s / (g$theta + (1 - g$theta) * s))^(1 / g$k) },
  formula = "survival (theta S^k / (1 - (1 - theta) S^k))^(1/k)",
  refs = "Aly and Benkherouf (2011), A new family of distributions based on probability generating functions. Sankhya B, 73, 70-80.")

#' @rdname exponentiated
#' @export
sine_g <- function(base) gen_family(base, "sine_g", "Sine", list(),
  T = function(u, g) sin(pi * u / 2), logt = function(u, g) log(pi / 2) + log(cos(pi * u / 2)),
  Tinv = function(v, g) 2 * asin(v) / pi, formula = "T(u) = sin(pi u / 2)",
  refs = "Kumar, Singh and Singh (2015), A method of proposing new distribution and its application to bladder cancer patients data. J. Stat. Appl. Prob. Lett., 2, 235-245.")

#' @rdname exponentiated
#' @export
dus <- function(base) gen_family(base, "dus", "DUS transformed", list(),
  T = function(u, g) expm1(u) / (exp(1) - 1), logt = function(u, g) u - log(exp(1) - 1),
  Tinv = function(v, g) log1p(v * (exp(1) - 1)), formula = "T(u) = (e^u - 1)/(e - 1)",
  refs = "Kumar, Singh and Singh (2015), A new method for generating distributions with an application to exponential distribution. Int. J. Stat. Appl. Prob.")

#' @rdname exponentiated
#' @export
km <- function(base) gen_family(base, "km", "Kavya-Manoharan transformed", list(),
  T = function(u, g) exp(1) / (exp(1) - 1) * (1 - exp(-u)), logt = function(u, g) 1 - u - log(exp(1) - 1),
  Tinv = function(v, g) -log1p(-v * (exp(1) - 1) / exp(1)), formula = "T(u) = e (1 - e^-u)/(e - 1)",
  refs = "Kavya and Manoharan (2021), Some parsimonious models for lifetimes and applications. J. Stat. Comput. Simul., 91, 3693-3708.")

ps_series <- function(series, size = 2) switch(series,
  geometric = list(C = function(t) t / (1 - t), dC = function(t) 1 / (1 - t)^2, dom = "unit"),
  poisson = list(C = function(t) expm1(t), dC = function(t) exp(t), dom = "positive"),
  logarithmic = list(C = function(t) -log1p(-t), dC = function(t) 1 / (1 - t), dom = "unit"),
  binomial = list(C = function(t) (1 + t)^size - 1, dC = function(t) size * (1 + t)^(size - 1), dom = "positive"),
  stop_input("series must be geometric, poisson, logarithmic or binomial."))

#' @rdname exponentiated
#' @export
power_series_max <- function(base, series = c("geometric", "poisson", "logarithmic", "binomial"), size = 2) {
  series <- match.arg(series); s <- ps_series(series, size)
  gen_family(base, paste0(series, "_max"), sprintf("%s maximum (power series) of", series),
    list(P("theta", s$dom, role = "series")),
    T = function(u, g) s$C(g$theta * u) / s$C(g$theta),
    logt = function(u, g) log(g$theta) + log(s$dC(g$theta * u)) - log(s$C(g$theta)),
    formula = "T(u) = C(theta u) / C(theta): maximum of N baseline draws, N zero truncated power series",
    refs = "Tahir and Cordeiro (2016), J. Stat. Distrib. Appl., 3, 13; Morais and Barreto-Souza (2011), A compound class of Weibull and power series distributions. CSDA, 55, 1410-1425.")
}

#' @rdname exponentiated
#' @export
power_series_min <- function(base, series = c("geometric", "poisson", "logarithmic", "binomial"), size = 2) {
  series <- match.arg(series); s <- ps_series(series, size)
  gen_family(base, paste0(series, "_min"), sprintf("%s minimum (power series) of", series),
    list(P("theta", s$dom, role = "series")),
    T = function(u, g) 1 - s$C(g$theta * (1 - u)) / s$C(g$theta),
    logt = function(u, g) log(g$theta) + log(s$dC(g$theta * (1 - u))) - log(s$C(g$theta)),
    formula = "T(u) = 1 - C(theta (1 - u)) / C(theta): minimum of N baseline draws",
    refs = c("Adamidis and Loukas (1998), A lifetime distribution with decreasing failure rate. Stat. Prob. Letters, 39, 35-42.",
             "Chahkandi and Ganjali (2009), On some lifetime distributions with decreasing failure rate. CSDA, 53, 4433-4440."))
}

#' Weighted and contaminated distributions
#'
#' `weighted()` gives the weighted law `f_w(y) = w(y) f(y) / E[w(Y)]`
#' with `w(y) = y^k` (length biased for `k = 1`, size biased in general).
#' The normalizing moment is computed by quadrature (or summation for a
#' discrete base).
#'
#' `contaminated()` gives `delta f(y; theta) + (1 - delta) f(y; theta*)`,
#' where `theta*` multiplies the parameter `par` by an inflation factor
#' `eta > 1` (Punzo and coauthors' contaminated models). Parameters `delta`
#' (logit link) and `eta` (log of `eta - 1`) are added.
#'
#' @param base Base family.
#' @param k Weight exponent.
#' @param par Name of the parameter inflated in the contaminating component.
#' @return A `regcomp_family`.
#' @references
#' Patil and Rao (1978), Weighted distributions and size-biased sampling. Biometrics, 34, 179-189.
#' Punzo, Mazza and Maruotti (2018), Fitting insurance and economic data with outliers: a flexible approach based on finite mixtures of contaminated gamma distributions. J. Appl. Stat., 45, 2563-2584.
#' @export
#' @examples
#' lb <- weighted("birnbaum_saunders")
#' cg <- contaminated("gamma", par = "sigma")
weighted <- function(base, k = 1) {
  base <- get_family(base)
  if (k == 0) return(base)
  if (k < 0) stop_input("k must be nonnegative.")
  logm <- function(p) {
    n <- max(lengths(p))
    vapply(seq_len(n), function(i) {
      pi_ <- lapply(p, function(v) v[(i - 1) %% length(v) + 1])
      m <- if (base$discrete) {
        s <- base$support(pi_)
        hi <- if (is.finite(s$upper)) s$upper else family_quantile(base, 1 - 1e-16, pi_) * 2 + 100
        if (hi > 1e7) stop_numerical("The weighted normalizer needs more than 1e7 terms.")
        ks <- max(s$lower, 0):hi; sum(ks^k * exp(eval_lpdf(base, ks, lapply(pi_, rep, length(ks)))))
      } else {
        s <- base$support(pi_)
        stats::integrate(function(y) y^k * exp(eval_lpdf(base, y, lapply(pi_, rep, length(y)))),
                         max(s$lower, 0), s$upper, rel.tol = 1e-10, subdivisions = 500L)$value
      }
      log(m) }, 0)
  }
  op_family(base, sprintf("weighted(%s, %s)", base$id, k), sprintf("Weighted (y^%s) %s", k, base$display_name),
            base$pars, function(y, p) k * log(y) + eval_lpdf(base, y, p) - logm(p),
            operator = list(type = "weighted", base = base, k = k),
            formula = "f_w(y) = y^k f(y) / E[Y^k]")
}

#' @rdname weighted
#' @export
contaminated <- function(base, par) {
  base <- get_family(base)
  if (!par %in% names(base$pars)) stop_input(sprintf("'%s' is not a parameter of %s.", par, base$id))
  dn <- new_par_name(base, "delta"); en <- new_par_name(base, "eta")
  pars <- c(base$pars, stats::setNames(list(P(dn, "unit", role = "share of good points"),
                                            P(en, "greater_than_one", role = "inflation")), c(dn, en)))
  bn <- names(base$pars)
  lpdf <- function(y, p) {
    p2 <- p[bn]; p2[[par]] <- p2[[par]] * p[[en]]
    log_add_exp(log(p[[dn]]) + eval_lpdf(base, y, p[bn]), log1p(-p[[dn]]) + eval_lpdf(base, y, p2))
  }
  cdf <- if (!is.null(base$cdf)) function(q, p) { p2 <- p[bn]; p2[[par]] <- p2[[par]] * p[[en]]
    p[[dn]] * family_cdf(base, q, p[bn]) + (1 - p[[dn]]) * family_cdf(base, q, p2) } else NULL
  op_family(base, sprintf("contaminated(%s, %s)", base$id, par), sprintf("Contaminated %s", base$display_name),
            pars, lpdf, cdf = cdf, operator = list(type = "contaminated", base = base, par = par),
            support = function(p) { p2 <- p[bn]; p2[[par]] <- p2[[par]] * p[[en]]
              s1 <- base$support(p[bn]); s2 <- base$support(p2)
              list(lower = pmin(s1$lower, s2$lower), upper = pmax(s1$upper, s2$upper)) },
            formula = sprintf("delta f(y) + (1 - delta) f(y; %s * eta)", par))
}

#' Quantile parameterized families
#'
#' Replace a parameter of a family by the `tau` quantile of the law, so a
#' regression on that parameter is a parametric quantile regression. For
#' given values of the remaining parameters, the original parameter is
#' found by root finding on its link scale so that `Q(tau) = q`.
#'
#' @param base Base family.
#' @param par Parameter to replace (default: the first one).
#' @param tau Quantile level in (0, 1).
#' @return A `regcomp_family` whose first parameter `q` is the quantile.
#' @references Noufaily and Jones (2013), Parametric quantile regression based on the generalized gamma distribution. JRSS C, 62, 723-740.
#' @export
#' @examples
#' qw <- quantile_param("weibull", par = "scale", tau = 0.5)
quantile_param <- function(base, par = names(get_family(base)$pars)[1], tau = 0.5) {
  base <- get_family(base)
  if (base$discrete) stop_unsupported("Quantile parameterization of a discrete law", base$id)
  ps <- base$pars[[par]]; lk <- link(ps$link)
  others <- setdiff(names(base$pars), par)
  solve_par <- function(p) {
    n <- max(lengths(p))
    vapply(seq_len(n), function(i) {
      pi_ <- lapply(p, function(v) v[(i - 1) %% length(v) + 1])
      f <- function(e) { q <- pi_[others]; q[[par]] <- lk$linkinv(e); family_quantile(base, tau, q[names(base$pars)]) - pi_$q }
      lo <- -2; hi <- 2; s0 <- sign(f(0))
      for (j in 1:60) { if (sign(f(lo)) != sign(f(hi))) break; lo <- lo * 1.6; hi <- hi * 1.6 }
      lk$linkinv(stats::uniroot(f, c(lo, hi), tol = 1e-12)$root) }, 0)
  }
  full <- function(p) { q <- p[others]; q[[par]] <- solve_par(p); q[names(base$pars)] }
  qpar <- P("q", if (base$category == "bounded") "unit" else if (base$category == "positive") "positive" else "real",
            role = "quantile", offset_allowed = TRUE)
  pars <- c(list(q = qpar), base$pars[others])
  op_family(base, sprintf("quantile_param(%s, %s, %s)", base$id, par, tau),
            sprintf("%s (tau = %s quantile parameterization)", base$display_name, tau), pars,
            function(y, p) eval_lpdf(base, y, full(p)),
            cdf = function(q, p) family_cdf(base, q, full(p)),
            quantile = function(pr, p) family_quantile(base, pr, full(p)),
            operator = list(type = "quantile_param", base = base, par = par, tau = tau),
            formula = sprintf("%s chosen so that Q(%s) = q", par, tau))
}
