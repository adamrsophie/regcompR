# Ordinal and nominal response laws (McCullagh 1980; McCullagh and Nelder
# 1989, ch. 5; Agresti 2010). Categories are coded 1, ..., K. The linear
# predictor `mu` is the location; `cuts` holds the K - 1 increasing
# thresholds (or intercepts) and `scale` gives location-scale models.

ord_pars <- function(cuts = "ordered") list(P("mu", "real", role = "location", link = "identity", offset_allowed = TRUE),
  P("cuts", "real", role = if (cuts == "ordered") "thresholds (K - 1, increasing)" else "intercepts (K - 1)", regression = FALSE, estimate = cuts),
  P("scale", "positive", role = "scale", default = 1))
# The location has no intercept (the thresholds absorb it); the scale is 1
# unless it gets a formula, which then has no intercept either.
ord_extra <- list(no_intercept = c("mu", "scale"), fixed_unless_formula = "scale")
ord_obs <- function(y, p) { n <- length(y); K <- length(p$cuts) + 1
  list(n = n, K = K, mu = rep_len(p$mu, n), s = if (is.null(p$scale)) rep(1, n) else rep_len(p$scale, n)) }
ord_check <- function(p) if (any(diff(p$cuts) <= 0)) "Thresholds must be strictly increasing."
ord_reg <- function(id, display, lpmf, formula, refs, reference, check = ord_check, cuts = "ordered") register_family(family_def(id, display, "ordinal", ord_pars(cuts),
  discrete = TRUE, variate_type = "categorical", check = check, extra = ord_extra, formula = formula, lpdf = lpmf,
  support = function(p) list(lower = 1, upper = length(p$cuts) + 1), numeric_moments = FALSE,
  references = refs, reference = reference, status = "validated"))

cum_lpmf <- function(F) function(y, p) { o <- ord_obs(y, p)
  up <- ifelse(y >= o$K, 1, F((p$cuts[pmin(y, o$K - 1)] - o$mu) / o$s))
  lo <- ifelse(y <= 1, 0, F((p$cuts[pmax(y - 1, 1)] - o$mu) / o$s))
  ifelse(y >= 1 & y <= o$K & y == round(y), log(pmax(up - lo, 0)), -Inf) }

#' Cumulative link ordinal families
#'
#' `P(Y <= k) = F((theta_k - mu) / scale)` with `F` the inverse of the
#' chosen link. Registered as `cumulative_logit` (proportional odds),
#' `cumulative_probit`, `cumulative_cloglog` (proportional hazards),
#' `cumulative_loglog` and `cumulative_cauchit`. A regression on `scale`
#' gives the location-scale model of McCullagh (1980).
#' @name ordinal-families
#' @references McCullagh (1980), Regression models for ordinal data. JRSS B, 42, 109-142.
NULL
for (lk in c("logit", "probit", "cloglog", "loglog", "cauchit")) local({
  F <- link(lk)$linkinv
  ord_reg(paste0("cumulative_", lk), sprintf("Cumulative %s (ordinal)", lk), cum_lpmf(F),
    formula = sprintf("P(Y <= k) = %s^-1((theta_k - mu) / scale)", lk),
    refs = c("McCullagh (1980), Regression models for ordinal data. JRSS B, 42, 109-142.", "McCullagh and Nelder (1989), Generalized Linear Models, 2nd ed., ch. 5."),
    reference = ref_impl("definition (MASS::polr convention)",
      local({ L <- lk; function(y, p) { cdf <- c(0, switch(L, logit = stats::plogis, probit = stats::pnorm, cauchit = stats::pcauchy,
        cloglog = function(x) -expm1(-exp(x)), loglog = function(x) exp(-exp(-x)))((p$cuts - p$mu[1]) / p$scale[1]), 1)
        log(diff(cdf))[y] } })))
})
assign("proportional_odds", "cumulative_logit", envir = .aliases)
assign("ordered_probit", "cumulative_probit", envir = .aliases)
assign("proportional_hazards_ordinal", "cumulative_cloglog", envir = .aliases)

# Continuation ratio (sequential): P(Y = k | Y >= k) = F(theta_k - mu).
cr_lpmf <- function(F) function(y, p) { o <- ord_obs(y, p)
  vapply(seq_len(o$n), function(i) { k <- y[i]; if (k < 1 || k > o$K || k != round(k)) return(-Inf)
    h <- F((p$cuts - o$mu[i]) / o$s[i]); sum(log1p(-h[seq_len(k - 1)])) + if (k < o$K) log(h[k]) else 0 }, 0) }
for (lk in c("logit", "probit", "cloglog")) local({
  F <- link(lk)$linkinv
  ord_reg(paste0("continuation_ratio_", lk), sprintf("Continuation ratio %s (sequential ordinal)", lk), cr_lpmf(F), check = NULL, cuts = "free",
    formula = sprintf("P(Y = k | Y >= k) = %s^-1(theta_k - mu)", lk),
    refs = c("Fienberg (1980), The Analysis of Cross-Classified Categorical Data.", "Tutz (1991), Sequential models in categorical regression. Comput. Stat. Data Anal., 11, 275-295."),
    reference = ref_impl("product of conditional probabilities", local({ G <- F; function(y, p) { h <- G((p$cuts - p$mu[1]) / p$scale[1])
      pr <- c(h, 1) * cumprod(c(1, 1 - h)); log(pr)[y] } })))
})
assign("sequential_ordinal", "continuation_ratio_logit", envir = .aliases)

# Adjacent category logit: log(P(k + 1) / P(k)) = (mu - theta_k) / scale.
ac_lpmf <- function(y, p) { o <- ord_obs(y, p)
  vapply(seq_len(o$n), function(i) { k <- y[i]; if (k < 1 || k > o$K || k != round(k)) return(-Inf)
    eta <- c(0, cumsum((o$mu[i] - p$cuts) / o$s[i])); eta[k] - log_sum_exp(eta) }, 0) }
ord_reg("adjacent_category_logit", "Adjacent category logit (ordinal)", ac_lpmf,
  formula = "log(P(Y = k + 1) / P(Y = k)) = (mu - theta_k) / scale",
  refs = "Agresti (2010), Analysis of Ordinal Categorical Data, 2nd ed., ch. 4.", check = NULL, cuts = "free",
  reference = ref_impl("VGAM acat", function(y, p) { eta <- c(0, cumsum((p$mu[1] - p$cuts) / p$scale[1])); (eta - log(sum(exp(eta))))[y] }))

# Stereotype model (Anderson 1984): P(Y = k) proportional to exp(alpha_k + phi_k mu).
register_family(family_def("stereotype", "Stereotype (ordinal)", "ordinal",
  list(P("mu", "real", role = "location", link = "identity"), P("alpha", "real", role = "category intercepts (K, first 0)", regression = FALSE),
       P("phi", "real", role = "category scores (K, 0 = phi_1 <= ... <= phi_K = 1)", regression = FALSE)),
  discrete = TRUE, variate_type = "categorical", numeric_moments = FALSE, support = function(p) list(lower = 1, upper = length(p$alpha)),
  formula = "P(Y = k) = exp(alpha_k + phi_k mu) / sum_j exp(alpha_j + phi_j mu)",
  lpdf = function(y, p) { mu <- rep_len(p$mu, length(y)); vapply(seq_along(y), function(i) { e <- p$alpha + p$phi * mu[i]; e[y[i]] - log_sum_exp(e) }, 0) },
  references = "Anderson (1984), Regression and ordered categorical variables. JRSS B, 46, 1-30.", extra = list(no_regression = TRUE),
  reference = ref_impl("softmax", function(y, p) { e <- exp(p$alpha + p$phi * p$mu[1]); log(e / sum(e))[y] }), status = "validated"))

# Multinomial (baseline category) logit: one linear predictor per non baseline category.
register_family(family_def("multinomial_logit", "Multinomial logit (nominal)", "categorical",
  list(P("eta", "real", role = "log odds versus category 1 (K - 1 values, or a matrix with one row per observation)", regression = FALSE)),
  discrete = TRUE, variate_type = "categorical", numeric_moments = FALSE,
  formula = "P(Y = k) = exp(eta_k) / (1 + sum_j exp(eta_j)), eta_1 = 0",
  lpdf = function(y, p) { E <- if (is.matrix(p$eta)) cbind(0, p$eta) else matrix(c(0, p$eta), length(y), length(p$eta) + 1, byrow = TRUE)
    vapply(seq_along(y), function(i) E[i, y[i]] - log_sum_exp(E[i, ]), 0) },
  regression_notes = "Each eta_k takes its own linear predictor (baseline category logit; McCullagh and Nelder 1989, sec. 5.2.4). In regcomp() the family expands to parameters eta2, ..., etaK that share the main formula unless given their own.",
  extra = list(expand = function(K) mlogit_family(K)),
  references = "McCullagh and Nelder (1989), Generalized Linear Models, 2nd ed., ch. 5.",
  reference = ref_impl("VGAM multinomial", function(y, p) { e <- exp(c(0, p$eta)); log(e / sum(e))[y] }), status = "validated"))

# Dale (1986) / Palmgren bivariate binary model with a global odds ratio.
plackett_p11 <- function(p1, p2, psi) ifelse(abs(psi - 1) < 1e-10, p1 * p2, {
  a <- 1 + (p1 + p2) * (psi - 1); (a - sqrt(a^2 - 4 * psi * (psi - 1) * p1 * p2)) / (2 * (psi - 1)) })
register_family(family_def("bivariate_binary_odds_ratio", "Bivariate binary with odds ratio (Dale; Palmgren)", "multivariate",
  list(P("p1", "unit", role = "P(Y1 = 1)", link = "logit"), P("p2", "unit", role = "P(Y2 = 1)", link = "logit"),
       P("psi", "positive", role = "odds ratio", link = "log")),
  aliases = c("bivariate_logistic", "dale_bivariate_binary", "palmgren"),
  discrete = TRUE, variate_type = "vector", numeric_cdf = FALSE, numeric_moments = FALSE,
  formula = "Marginal logits for Y1 and Y2; log odds ratio log(p11 p00 / (p10 p01)) = log(psi); p11 from the Plackett solution",
  lpdf = function(y, p) { y <- if (is.matrix(y)) y else matrix(y, nrow = 1); n <- nrow(y)
    p1 <- rep_len(p$p1, n); p2 <- rep_len(p$p2, n); ps <- rep_len(p$psi, n); p11 <- plackett_p11(p1, p2, ps)
    log(ifelse(y[, 1] == 1, ifelse(y[, 2] == 1, p11, p1 - p11), ifelse(y[, 2] == 1, p2 - p11, 1 - p1 - p2 + p11))) },
  references = c("Dale (1986), Global cross-ratio models for bivariate, discrete, ordered responses. Biometrics, 42, 909-917.",
                 "Palmgren (1989), Regression models for bivariate binary responses. UW Biostatistics Working Paper 101.",
                 "McCullagh and Nelder (1989), Generalized Linear Models, 2nd ed., sec. 6.5."),
  reference = ref_impl("odds ratio equation solved numerically", function(y, p) { y <- if (is.matrix(y)) y else matrix(y, nrow = 1)
    vapply(seq_len(nrow(y)), function(i) { a <- p$p1[i]; b <- p$p2[i]
      q <- stats::uniroot(function(x) log(x * (1 - a - b + x)) - log((a - x) * (b - x)) - log(p$psi[i]), c(max(0, a + b - 1) + 1e-12, min(a, b) - 1e-12), tol = 1e-14)$root
      log(c(1 - a - b + q, b - q, a - q, q))[1 + y[i, 2] + 2 * y[i, 1]] }, 0) }),
  status = "validated"))

# Multinomial logit with K categories as K - 1 regression parameters.
mlogit_family <- function(K) {
  nm <- paste0("eta", 2:K)
  family_def(paste0("multinomial_logit"), sprintf("Multinomial logit (%d categories)", K), "categorical",
    lapply(nm, function(k) P(k, "real", role = sprintf("log odds of category %s versus 1", sub("eta", "", k)), link = "identity")),
    discrete = TRUE, variate_type = "categorical", numeric_moments = FALSE,
    lpdf = function(y, p) { E <- cbind(0, do.call(cbind, p[nm]))
      m <- apply(E, 1, max); E[cbind(seq_along(y), y)] - m - log(rowSums(exp(E - m))) },
    extra = list(share_formula = TRUE), status = "validated")
}

# Heckman (1979) sample selection model as a law for the pair (s, y): s is
# the selection indicator and y the outcome, missing when s = 0.
register_family(family_def("heckman_selection", "Heckman sample selection (bivariate normal)", "multivariate",
  list(P("mu", "real", role = "location", link = "identity", offset_allowed = TRUE, interpretation = "Outcome mean"),
       P("gamma", "real", role = "selection index", link = "identity", interpretation = "Probit index of selection"),
       P("sigma", "positive", role = "scale"), P("rho", "symmetric_unit", role = "correlation")),
  aliases = c("heckman", "sample_selection", "tobit2"),
  discrete = FALSE, variate_type = "vector", numeric_cdf = FALSE, numeric_moments = FALSE,
  formula = "P(s = 0) = Phi(-gamma); f(s = 1, y) = phi((y - mu)/sigma)/sigma Phi((gamma + rho (y - mu)/sigma) / sqrt(1 - rho^2))",
  lpdf = function(y, p) { y <- if (is.matrix(y)) y else matrix(y, nrow = 1); n <- nrow(y)
    mu <- rep_len(p$mu, n); g <- rep_len(p$gamma, n); s <- rep_len(p$sigma, n); r <- rep_len(p$rho, n)
    z <- (y[, 2] - mu) / s
    ifelse(y[, 1] == 0, stats::pnorm(-g, log.p = TRUE),
           stats::dnorm(z, log = TRUE) - log(s) + stats::pnorm((g + r * z) / sqrt(1 - r^2), log.p = TRUE)) },
  regression_notes = "Response cbind(s, y) with y missing (any value) when s = 0; the main formula models the outcome mean and gamma ~ z the selection equation (Heckman 1979; Type 2 Tobit). Fitted by full maximum likelihood.",
  references = c("Heckman (1979), Sample selection bias as a specification error. Econometrica, 47, 153-161.",
                 "Amemiya (1985), Advanced Econometrics, ch. 10."),
  reference = ref_impl("bivariate normal integral", function(y, p) { y <- if (is.matrix(y)) y else matrix(y, nrow = 1)
    vapply(seq_len(nrow(y)), function(i) if (y[i, 1] == 0) log(stats::pnorm(-p$gamma[i])) else
      log(stats::integrate(function(v) mvtnorm::dmvnorm(cbind(v, y[i, 2]), c(p$gamma[i], p$mu[i]),
        matrix(c(1, p$rho[i] * p$sigma[i], p$rho[i] * p$sigma[i], p$sigma[i]^2), 2)), 0, Inf, rel.tol = 1e-12)$value), 0) }),
  status = "validated"))
