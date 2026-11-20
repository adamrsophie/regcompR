# Count laws from the count regression literature (Hilbe 2014, Modeling Count
# Data, ch. 5 to 9): three parameter models, truncated variants and the
# multivariate hypergeometric law.

cfam3 <- function(id, display, pars, ...) register_family(family_def(id, display, "count", pars, discrete = TRUE,
  support = function(p) list(lower = 0, upper = Inf), ...))

gw_lpmf <- function(y, a, k, r) lgamma(a + r) + lgamma(k + r) - lgamma(r) - lgamma(a + k + r) +
  lgamma(a + y) - lgamma(a) + lgamma(k + y) - lgamma(k) - (lgamma(a + k + r + y) - lgamma(a + k + r)) - lgamma(y + 1)
cfam3("generalized_waring", "Generalized Waring (GWRM mean parameterization)", list(log_mean_par("mu"), shp("k"), P("rho", "greater_than_one", role = "shape")),
  aliases = c("gwrm", "GW"),
  formula = "p(y) = Gamma(a + rho) Gamma(k + rho) / (Gamma(rho) Gamma(a + k + rho)) (a)_y (k)_y / ((a + k + rho)_y y!), a = mu (rho - 1) / k",
  lpdf = function(y, p) gw_lpmf(y, p$mu * (p$rho - 1) / p$k, p$k, p$rho),
  mean = function(p) p$mu,
  variance = function(p) { a <- p$mu * (p$rho - 1) / p$k; ifelse(p$rho > 2, p$mu * (a + p$rho - 1) * (p$k + p$rho - 1) / ((p$rho - 1) * (p$rho - 2)), Inf) },
  regression_notes = "Log link on mu with constant k and rho is the generalized Waring regression of Rodriguez-Avi et al. (2009); the variance splits into randomness, liability and proneness.",
  references = c("Xekalaki (1983), The univariate generalized Waring distribution in relation to accident theory. Biometrics, 39, 887-895.",
                 "Rodriguez-Avi, Conde-Sanchez, Saez-Castillo, Olmo-Jimenez and Martinez-Rodriguez (2009), A generalized Waring regression model for count data. Comput. Stat. Data Anal., 53, 3717-3725."),
  reference = ref_impl("hypergeometric series ratio", function(y, p) { a <- p$mu * (p$rho - 1) / p$k
    vapply(seq_along(y), function(i) { j <- seq_len(y[i]) - 1
      gw_lpmf(0, a[i], p$k[i], p$rho[i]) + sum(log((a[i] + j) * (p$k[i] + j) / ((a[i] + p$k[i] + p$rho[i] + j) * (j + 1)))) }, 0) }),
  stan = sfun("gwaring", c("mu", "k", "rho"),
    "    real a = mu * (rho - 1) / k;\n    return lgamma(a + rho) + lgamma(k + rho) - lgamma(rho) - lgamma(a + k + rho) + lgamma(a + y) - lgamma(a) + lgamma(k + y) - lgamma(k) - lgamma(a + k + rho + y) + lgamma(a + k + rho) - lgamma(y + 1);",
    discrete = TRUE),
  status = "validated")

gnb_lpmf <- function(y, m, b, th) log(m) - log(m + b * y) + lgamma(m + b * y + 1) - lgamma(y + 1) - lgamma(m + b * y - y + 1) +
  y * log(th) + (m + b * y - y) * log1p(-th)
cfam3("generalized_negative_binomial", "Generalized negative binomial (Jain-Consul; Famoye NBF)", list(log_mean_par("mu"), P("m", "positive", role = "shape"),
  P("beta", "positive", role = "generalization", interpretation = "beta >= 1; beta = 1 is the negative binomial, beta > 1 adds overdispersion")),
  check = function(p) if (any(p$beta < 1, na.rm = TRUE)) "Need beta >= 1.",
  aliases = c("jain_consul_negative_binomial", "famoye_negative_binomial", "NBF_famoye"),
  formula = "p(y) = m / (m + beta y) C(m + beta y, y) theta^y (1 - theta)^(m + beta y - y), theta = mu / (m + beta mu)",
  lpdf = function(y, p) gnb_lpmf(y, p$m, p$beta, p$mu / (p$m + p$beta * p$mu)),
  mean = function(p) p$mu,
  variance = function(p) { th <- p$mu / (p$m + p$beta * p$mu); p$m * th * (1 - th) / (1 - p$beta * th)^3 },
  references = c("Jain and Consul (1971), A generalized negative binomial distribution. SIAM J. Appl. Math., 21, 501-513.",
                 "Famoye (1995), Generalized binomial regression model. Biometrical J., 37, 581-594.",
                 "Harris, Hilbe and Hardin (2014), Modeling count data with generalized distributions. Stata J., 14, 562-579."),
  reference = ref_impl("beta = 1 is NB(size m, prob 1 - theta); mass and mean by summation", function(y, p) gnb_lpmf(y, p$m, p$beta, p$mu / (p$m + p$beta * p$mu))),
  status = "validated")

nbl_mix <- function(y, r, th) vapply(seq_along(y), function(i) log(stats::integrate(function(l)
  exp(lgamma(r[i] + y[i]) - lgamma(r[i]) - lgamma(y[i] + 1) - r[i] * l + (if (y[i] == 0) 0 else y[i] * log(-expm1(-pmax(l, 1e-300)))) +
    2 * log(th[i]) - log1p(th[i]) + log1p(l) - th[i] * l), 0, Inf, rel.tol = 1e-12)$value), 0)
cfam3("negative_binomial_lindley", "Negative binomial-Lindley", list(shp("r"), P("theta", "greater_than_one", role = "Lindley shape")),
  aliases = c("nb_lindley", "NBL"),
  formula = "Y | lambda ~ NB(r, p = exp(-lambda)), lambda ~ Lindley(theta); p(y) = Gamma(r + y)/(Gamma(r) y!) theta^2/(theta + 1) sum_j C(y, j) (-1)^j (theta + r + j + 1)/(theta + r + j)^2",
  lpdf = function(y, p) nbl_mix(y, p$r, p$theta),
  mean = function(p) p$r * (p$theta^3 / ((p$theta + 1) * (p$theta - 1)^2) - 1),
  regression_notes = "Lord and Geedipally (2011) link the mean through r with theta constant; use the log link on r.",
  references = c("Zamani and Ismail (2010), Negative binomial-Lindley distribution and its application. J. Math. Stat., 6, 4-9.",
                 "Lord and Geedipally (2011), The negative binomial-Lindley distribution as a tool for analyzing crash data characterized by a large amount of zeros. Accid. Anal. Prev., 43, 1738-1742."),
  reference = ref_impl("alternating closed form", function(y, p) vapply(seq_along(y), function(i) { j <- 0:y[i]; th <- p$theta[i]; r <- p$r[i]
    lgamma(r + y[i]) - lgamma(r) - lgamma(y[i] + 1) + log(th^2 / (th + 1) * sum(choose(y[i], j) * (-1)^j * (th + r + j + 1) / (th + r + j)^2)) }, 0),
    note = "The alternating sum is used for small y only."),
  status = "validated")

eppm_probs <- function(a, b, cc, nmax) { n <- 0:nmax; lam <- pmin(a * (b + n)^cc, 745)
  Q <- matrix(0, nmax + 1, nmax + 1); diag(Q) <- -lam; Q[cbind(1:nmax, 2:(nmax + 1))] <- lam[1:nmax]; Q[nmax + 1, nmax + 1] <- 0
  pmax(as.numeric(as.matrix(Matrix::expm(Matrix::Matrix(Q)))[1, ]), 0) }
cfam3("eppm_faddy", "Extended Poisson process model (Faddy and Smith)", list(P("a", "positive", role = "rate"), P("b", "positive", role = "shift"),
  P("c", "real", role = "rate power", link = "identity", interpretation = "c < 1; c = 0 is Poisson, c = 1 the negative binomial limit")),
  aliases = c("eppm", "faddy"),
  formula = "State at time 1 of a pure birth process from 0 with rates lambda_n = a (b + n)^c",
  # One matrix exponential per distinct parameter set, covering all its counts.
  lpdf = function(y, p) { r <- recycle_pars(c(list(.y = y), p)); key <- paste(r$a, r$b, r$c); out <- numeric(length(r$.y))
    for (g in unique(key)) { i <- which(key == g); a <- r$a[i[1]]; b <- r$b[i[1]]; cc <- r$c[i[1]]; x <- r$.y[i]
      if (abs(cc) < 1e-12) { out[i] <- stats::dpois(x, a, log = TRUE); next }
      xm <- max(x); nmax <- max(50, xm + 40, ceiling(5 * a * (b + xm)^max(cc, 0) + 50))
      out[i] <- log(eppm_probs(a, b, cc, nmax)[x + 1]) }
    out },
  numerical_notes = "Probabilities from the matrix exponential of the generator truncated at a state far beyond the data (the last state absorbs).",
  references = c("Faddy (1997), Extended Poisson process modelling and analysis of count data. Biometrical J., 39, 431-440.",
                 "Faddy and Smith (2011), Analysis of count data with covariate dependence in both mean and variance. J. Appl. Stat., 38, 2683-2694."),
  reference = ref_impl("CountsEPPM", function(y, p) vapply(seq_along(y), function(i)
    log(CountsEPPM::Faddyprob.general(list(p$a[i], p$b[i], p$c[i]), max(y[i] + 60, 120))[y[i] + 1]), 0)),
  status = "validated")

register_family(family_def("multivariate_hypergeometric", "Multivariate hypergeometric", "multivariate",
  list(P("K", "positive_integer", role = "population counts per category", regression = FALSE),
       P("size", "positive_integer", data = TRUE, default = NA, interpretation = "Number of draws; NA conditions on each observed total")),
  discrete = TRUE, variate_type = "vector",
  numeric_cdf = FALSE, numeric_moments = FALSE,
  check = function(p) if (!all(is.na(p$size)) && any(p$size > sum(p$K), na.rm = TRUE)) "The number of draws exceeds the population size sum(K).",
  formula = "P(y) = prod C(K_i, y_i) / C(N, n), N = sum K, n = sum y (or the given size)",
  lpdf = function(y, p) { y <- if (is.matrix(y)) y else matrix(y, nrow = 1)
    out <- apply(y, 1, function(r) if (any(r > p$K | r < 0)) -Inf else sum(lchoose(p$K, r)) - lchoose(sum(p$K), sum(r)))
    if (!all(is.na(p$size))) out[rowSums(y) != rep_len(p$size, nrow(y))] <- -Inf
    out },
  rng = function(n, p) { if (anyNA(p$size)) stop_input("Simulating a multivariate hypergeometric law needs the number of draws: give size.")
    K <- p$K; t(vapply(seq_len(n), function(i) tabulate(sample(rep(seq_along(K), K), p$size[1]), length(K)), numeric(length(K)))) },
  mean = function(p) p$size[1] * p$K / sum(p$K),
  references = "Johnson, Kotz and Balakrishnan (1997), Discrete Multivariate Distributions, ch. 39.",
  reference = ref_impl("extraDistr", function(y, p) { y <- if (is.matrix(y)) y else matrix(y, nrow = 1); extraDistr::dmvhyper(y, p$K, rowSums(y)[1], log = TRUE) }),
  status = "validated"))

# Truncated variants used for hurdle and zero-free data (Hilbe 2014, sec. 7.1).
local({
  zt <- list(zero_truncated_pig = c("poisson_inverse_gaussian", "Zero-truncated Poisson inverse Gaussian", "ZTPIG"),
             zero_truncated_nbp = c("negative_binomial_family", "Zero-truncated NB-P", "ZTNBP"),
             zero_truncated_poisson_lognormal = c("poisson_lognormal", "Zero-truncated Poisson lognormal", "ZTPLN"),
             zero_truncated_nb1 = c("negative_binomial_nb1", "Zero-truncated NB1", "ZTNB1"),
             zero_truncated_generalized_poisson = c("generalized_poisson_mean", "Zero-truncated generalized Poisson", "ZTGP"))
  for (id in names(zt)) { base <- get_family(zt[[id]][1])
    register_family(modify_family(truncate_dist(base, lower = 1), id = id, display_name = zt[[id]][2], aliases = zt[[id]][3],
      reference = ref_impl("renormalized base", local({ b <- base; function(y, p) eval_lpdf(b, y, p) - log1p(-exp(eval_lpdf(b, rep(0, length(y)), p))) })),
      status = "validated")) }
  assign("nbp", "negative_binomial_family", envir = .aliases); assign("negative_binomial_p", "negative_binomial_family", envir = .aliases)
})
