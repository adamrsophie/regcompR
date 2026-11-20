# Further operators: mixture cure, splicing, convolution and copulas.

rename_pars <- function(fam, suffix) {
  out <- lapply(fam$pars, function(p) { p$name <- paste0(p$name, suffix); p$symbol <- p$name; p$is_mean <- FALSE; p })
  stats::setNames(out, paste0(names(fam$pars), suffix))
}
suffix_pars <- function(p, fam, suffix) stats::setNames(lapply(names(fam$pars), function(k) p[[paste0(k, suffix)]]), names(fam$pars))
# Constraint check of a multi component family: each component's own check
# on its renamed parameters.
component_checks <- function(fams, suffixes) {
  if (all(vapply(fams, function(f) is.null(f$check), TRUE))) return(NULL)
  function(p) {
    for (k in seq_along(fams)) if (!is.null(fams[[k]]$check)) {
      m <- fams[[k]]$check(suffix_pars(p, fams[[k]], suffixes[k]))
      if (is.character(m) && length(m) && nzchar(m[1])) return(m[1])
    }
    NULL
  }
}

#' Mixture cure models
#'
#' A fraction `cure` of subjects never has the event: `S(t) = cure + (1 -
#' cure) S0(t)` and `f(t) = (1 - cure) f0(t)` (Boag 1949; Berkson and Gage
#' 1952; Farewell 1982). The cure fraction is a parameter in (0, 1) with a
#' logit link, so `cure ~ z` gives the incidence model and the base
#' parameters give the latency model. Use with right censored responses
#' (`surv()`); the law is improper (its distribution function tends to
#' `1 - cure`).
#'
#' @param base A family for the event time of uncured subjects.
#' @return A `regcomp_family`.
#' @export
#' @examples
#' cw <- cure("weibull")
#' set.seed(1)
#' n <- 400; cured <- runif(n) < 0.3
#' t <- ifelse(cured, Inf, rweibull(n, 1.5, 2)); cens <- runif(n, 0, 6)
#' d <- data.frame(time = pmin(t, cens), event = as.numeric(t <= cens))
#' fit <- regcomp(surv(time, event) ~ 1, family = cw, data = d)
#' plogis(coef(fit)["cure:(Intercept)"])
cure <- function(base) {
  base <- get_family(base)
  if (base$discrete) stop_input("cure() needs a continuous event time law.")
  pars <- c(base$pars, list(cure = P("cure", "unit", role = "cure fraction", link = "logit")))
  bp <- function(p) p[names(base$pars)]
  lpdf <- function(y, p) log1p(-p$cure) + eval_lpdf(base, y, bp(p))
  cdf <- function(q, p) (1 - p$cure) * family_cdf(base, q, bp(p))
  lccdf <- function(q, p) log(p$cure + (1 - p$cure) * (1 - family_cdf(base, q, bp(p))))
  op_family(base, sprintf("cure(%s)", base$id), sprintf("Mixture cure %s", base$display_name), pars, lpdf,
    cdf = cdf, lccdf = lccdf, numeric_moments = FALSE,
    # The cured fraction is mass at infinity: the support extends to Inf, so
    # F(t) + S(t) = 1 also above a finite upper bound of the uncured law.
    support = function(p) list(lower = base$support(bp(p))$lower, upper = Inf),
    quantile = function(pr, p) { r <- pr / (1 - p$cure); out <- rep(Inf, length(pr)); ok <- r < 1
      if (any(ok)) out[ok] <- family_quantile(base, r[ok], lapply(bp(p), function(v) if (length(v) == length(pr)) v[ok] else v)); out },
    rng = function(n, p) ifelse(stats::runif(n) < p$cure, Inf, family_rng(base, n, bp(p))),
    operator = list(type = "cure", base = base),
    formula = "S(t) = cure + (1 - cure) S0(t); f(t) = (1 - cure) f0(t)",
    extra = list(primary = names(base$pars)[vapply(base$pars, function(p) p$role %in% c("scale", "location", "rate") || p$is_mean, TRUE)][1]))
}

#' Spliced (composite) distributions
#'
#' Joins a body law below a threshold `u` and a tail law above it:
#' `f(y) = w f_b(y) / F_b(u)` for `y <= u` and `(1 - w) f_t(y) / S_t(u)` for
#' `y > u` (Cooray and Ananda 2005; Scollnik 2007). The threshold is a
#' parameter (`threshold`, identity link unless the body is positive) and
#' `w` the body weight (logit link). Parameters of the body get suffix
#' `.b` and of the tail `.t`.
#'
#' @param body,tail Families.
#' @return A `regcomp_family`.
#' @export
#' @examples
#' sp <- splice("lognormal", "lomax")
#' names(sp$pars)
splice <- function(body, tail) {
  body <- get_family(body); tail <- get_family(tail)
  if (body$discrete || tail$discrete) stop_input("splice() needs continuous laws.")
  # Point masses of the body at or below the threshold and of the tail above
  # it are kept, rescaled like the densities.
  pars <- c(rename_pars(body, ".b"), rename_pars(tail, ".t"),
            list(threshold = P("threshold", "real", role = "threshold", link = "identity"),
                 w = P("w", "unit", role = "body weight", link = "logit")))
  pb <- function(p) suffix_pars(p, body, ".b"); pt <- function(p) suffix_pars(p, tail, ".t")
  lpdf <- function(y, p) { u <- p$threshold
    lb <- log(p$w) + eval_lpdf(body, y, pb(p)) - log(family_cdf(body, u, pb(p)))
    lt <- log1p(-p$w) + eval_lpdf(tail, y, pt(p)) - log1p(-family_cdf(tail, u, pt(p)))
    ifelse(y <= u, lb, lt) }
  cdf <- function(q, p) { u <- p$threshold
    Fb <- family_cdf(body, pmin(q, u), pb(p)) / family_cdf(body, u, pb(p))
    Ftu <- family_cdf(tail, u, pt(p))
    ifelse(q <= u, p$w * Fb, p$w + (1 - p$w) * (family_cdf(tail, pmax(q, u), pt(p)) - Ftu) / (1 - Ftu)) }
  support <- function(p) list(lower = body$support(pb(p))$lower, upper = tail$support(pt(p))$upper)
  atomic <- has_atoms(body) || has_atoms(tail)
  atom <- if (atomic) function(y, p) { u <- p$threshold
    ifelse(y <= u, p$w * family_atom(body, y, pb(p)) / family_cdf(body, u, pb(p)),
           (1 - p$w) * family_atom(tail, y, pt(p)) / (1 - family_cdf(tail, u, pt(p)))) }
  known <- (!has_atoms(body) || !is.null(body$extra$atom_points)) && (!has_atoms(tail) || !is.null(tail$extra$atom_points))
  f <- family_def(sprintf("splice(%s, %s)", body$id, tail$id), sprintf("Spliced %s body and %s tail", body$display_name, tail$display_name),
    body$category, pars, lpdf = lpdf, cdf = cdf, support = support, status = "implemented", atom = atom,
    extra = if (atomic && known) list(atom_points = sort(unique(c(body$extra$atom_points, tail$extra$atom_points)))) else list(),
    check = component_checks(list(body, tail), c(".b", ".t")),
    formula = "w f_b(y)/F_b(u) for y <= u; (1 - w) f_t(y)/S_t(u) for y > u", references = c("Cooray and Ananda (2005), Modeling actuarial data with a composite lognormal-Pareto model. Scand. Actuar. J., 2005, 321-334.",
      "Scollnik (2007), On composite lognormal-Pareto models. Scand. Actuar. J., 2007, 20-33."))
  f$operator <- list(type = "splice", body = body, tail = tail)
  f
}

#' Convolution of two independent laws
#'
#' The law of `X1 + X2` for independent `X1` and `X2`: the density is
#' `int f1(y - x) f2(x) dx` by adaptive quadrature for continuous laws and
#' a finite sum for counts. Parameters get suffixes `.1` and `.2`. Examples
#' include the exponentially modified laws, the normal plus Laplace law of
#' Reed and Jorgensen and Poisson plus negative binomial counts.
#'
#' @param a,b Families, both continuous (possibly with point masses at
#'   known places, which give point masses of the sum at the sums of their
#'   locations) or both discrete.
#' @return A `regcomp_family`.
#' @export
#' @examples
#' cv <- convolve_dist("normal", "exponential")
#' d <- distribution(cv, mu.1 = 0, sigma.1 = 1, lambda.2 = 2)
#' prob(d, 1)
convolve_dist <- function(a, b) {
  a <- get_family(a); b <- get_family(b)
  if (a$discrete != b$discrete) stop_input("convolve_dist() needs two continuous or two discrete laws.")
  # Laws with point masses at known places c_j (masses A_j) and d_k (B_k):
  # the sum has the continuous part int f1c(y - x) f2c(x) dx +
  # sum_j A_j f2c(y - c_j) + sum_k B_k f1c(y - d_k) and point masses
  # A_j B_k at c_j + d_k.
  if ((has_atoms(a) && is.null(a$extra$atom_points)) || (has_atoms(b) && is.null(b$extra$atom_points)))
    stop_unsupported("Convolution of a law whose point mass locations are unknown")
  pts <- a$extra$atom_points; ptsb <- b$extra$atom_points
  both <- length(pts) && length(ptsb)
  spts <- if (both) sort(unique(as.vector(outer(pts, ptsb, `+`)))) else numeric()
  if (a$discrete && (lattice_step(a) != 1 || lattice_step(b) != 1)) stop_unsupported("Convolution of discrete laws off the integer lattice")
  sa_ <- ".1"; sb_ <- ".2"
  pars <- c(rename_pars(a, ".1"), rename_pars(b, ".2"))
  pa <- function(p) suffix_pars(p, a, sa_); pb <- function(p) suffix_pars(p, b, sb_)
  one <- function(y, p) {
    if (a$discrete) {
      sb <- b$support(pb(p)); sa <- a$support(pa(p))
      lo <- max(sb$lower, y - sa$upper); hi <- min(sb$upper, y - sa$lower)
      if (!is.finite(lo) || !is.finite(hi)) stop_unsupported("Convolution of two laws unbounded on the same side")
      if (hi < lo) return(-Inf)
      x <- seq(lo, hi)
      return(log_sum_exp(eval_lpdf(a, y - x, lapply(pa(p), rep, length(x))) + eval_lpdf(b, x, lapply(pb(p), rep, length(x)))))
    }
    if (both && any(abs(spts - y) <= 1e-12)) return(log(atom1(y, p)))
    sb <- b$support(pb(p))
    # Continuous parts: the density with the point masses removed.
    fc <- function(f, pts_, x, pp) { l <- eval_lpdf(f, x, lapply(pp, rep, length(x))); if (length(pts_)) l[x %in% pts_] <- -Inf; exp(l) }
    dens <- function(x) fc(a, pts, y - x, pa(p)) * fc(b, ptsb, x, pb(p))
    sa <- a$support(pa(p))
    lo <- max(sb$lower, y - sa$upper); hi <- min(sb$upper, y - sa$lower)
    v <- if (hi > lo) integrate_safe(dens, lo, hi) else 0
    if (length(pts)) v <- v + sum(family_atom(a, pts, lapply(pa(p), rep, length(pts))) * fc(b, ptsb, y - pts, pb(p)))
    if (length(ptsb)) v <- v + sum(family_atom(b, ptsb, lapply(pb(p), rep, length(ptsb))) * fc(a, pts, y - ptsb, pa(p)))
    log(v)
  }
  atom1 <- function(y, p) { if (!both) return(0)
    A <- family_atom(a, pts, lapply(pa(p), rep, length(pts))); B <- family_atom(b, ptsb, lapply(pb(p), rep, length(ptsb)))
    sum(outer(A, B)[abs(outer(pts, ptsb, `+`) - y) <= 1e-12]) }
  # P(X1 + X2 <= q) = E F2(q - X1) over the continuous part and atoms of X1.
  cdf <- if (both) function(q, p) pervec(q, p, function(x, pp) {
    sa <- a$support(pa(pp))
    fa <- function(t) { l <- eval_lpdf(a, t, lapply(pa(pp), rep, length(t))); l[t %in% pts] <- -Inf; exp(l) }
    g <- function(t) fa(t) * family_cdf(b, x - t, lapply(pb(pp), rep, length(t)))
    integrate_safe(g, sa$lower, sa$upper) + sum(family_atom(a, pts, lapply(pa(pp), rep, length(pts))) *
      family_cdf(b, x - pts, lapply(pb(pp), rep, length(pts)))) })
  lpdf <- function(y, p) pervec(y, p, function(x, pp) one(x, pp))
  support <- function(p) { s1 <- a$support(pa(p)); s2 <- b$support(pb(p)); list(lower = s1$lower + s2$lower, upper = s1$upper + s2$upper) }
  f <- family_def(sprintf("convolve(%s, %s)", a$id, b$id), sprintf("Convolution of %s and %s", a$display_name, b$display_name),
    a$category, pars, lpdf = lpdf, support = support, discrete = a$discrete, status = "implemented", cdf = cdf,
    atom = if (both) function(y, p) pervec(y, p, function(x, pp) atom1(x, pp)),
    extra = if (both) list(atom_points = spts) else list(),
    rng = function(n, p) family_rng(a, n, pa(p)) + family_rng(b, n, pb(p)),
    check = component_checks(list(a, b), c(sa_, sb_)),
    mean = if (!is.null(a$mean) && !is.null(b$mean)) function(p) a$mean(pa(p)) + b$mean(pb(p)),
    formula = "f(y) = int f1(y - x) f2(x) dx", references = "Feller (1971), An Introduction to Probability Theory and Its Applications, vol. 2, sec. V.4.")
  f$operator <- list(type = "convolve", a = a, b = b)
  f
}

#' Bivariate copula models
#'
#' Joint law of a pair `(y1, y2)` with margins `m1`, `m2` and a copula:
#' `"gaussian"` (correlation `rho`), `"clayton"` (`theta > 0`), `"gumbel"`
#' (`theta > 1`, parameterized as `1 + exp(eta)`) or `"frank"` (`theta`
#' real, not 0). For continuous margins the log density is
#' `log c(F1(y1), F2(y2)) + log f1(y1) + log f2(y2)`; for discrete margins
#' the rectangle probability `C(F1, F2) - C(F1-, F2) - C(F1, F2-) +
#' C(F1-, F2-)` is used (Trivedi and Zimmer 2007). Margins may differ in
#' type, and continuous margins may have point masses: at a jump of one
#' margin and a continuous value of the other the density is
#' `(h(F1, F2) - h(F1-, F2)) f2(y2)` with `h = dC/dv`. Margin parameters get
#' suffixes `.1` and `.2`; every parameter can take a formula, giving
#' copula regression. The response is `cbind(y1, y2)`.
#'
#' @param m1,m2 Scalar margin families (continuous, discrete, or continuous
#'   with point masses).
#' @param copula Copula name.
#' @return A `regcomp_family`.
#' @export
#' @examples
#' cp <- copula_dist("normal", "gamma", "clayton")
#' names(cp$pars)
copula_dist <- function(m1, m2, copula = c("gaussian", "clayton", "gumbel", "frank")) {
  copula <- match.arg(copula)
  m1 <- get_family(m1); m2 <- get_family(m2)
  if (m1$variate_type != "scalar" || m2$variate_type != "scalar") stop_input("copula_dist() needs scalar margins.")
  dep <- switch(copula, gaussian = P("rho", "symmetric_unit", role = "copula correlation"),
    clayton = P("theta", "positive", role = "copula dependence"),
    gumbel = P("theta", "greater_than_one", role = "copula dependence"),
    frank = P("theta", "real", role = "copula dependence", link = "identity"))
  pars <- c(rename_pars(m1, ".1"), rename_pars(m2, ".2"), stats::setNames(list(dep), dep$name))
  p1 <- function(p) suffix_pars(p, m1, ".1"); p2 <- function(p) suffix_pars(p, m2, ".2")
  Cfun <- function(u, v, th) { u <- pmin(pmax(u, 1e-15), 1 - 1e-15); v <- pmin(pmax(v, 1e-15), 1 - 1e-15)
    switch(copula,
      gaussian = vapply(seq_along(u), function(i) { if (i == 1 && !requireNamespace("mvtnorm", quietly = TRUE))
        stop_unsupported("A Gaussian copula with discrete margins without the 'mvtnorm' package; it"); mvtnorm::pmvnorm(upper = stats::qnorm(c(u[i], v[i])), corr = matrix(c(1, th[i], th[i], 1), 2))[1] }, 0),
      clayton = pmax(u^-th + v^-th - 1, 0)^(-1 / th),
      gumbel = exp(-((-log(u))^th + (-log(v))^th)^(1 / th)),
      frank = ifelse(th == 0, u * v, -log1p(expm1(-th * u) * expm1(-th * v) / expm1(-th)) / th)) }
  logc <- function(u, v, th) switch(copula,
    gaussian = { a <- stats::qnorm(u); b <- stats::qnorm(v); -0.5 * log1p(-th^2) - (th^2 * (a^2 + b^2) - 2 * th * a * b) / (2 * (1 - th^2)) },
    clayton = log1p(th) - (1 + th) * (log(u) + log(v)) - (2 + 1 / th) * log(u^-th + v^-th - 1),
    gumbel = { x <- -log(u); y <- -log(v); A <- (x^th + y^th)^(1 / th)
      -A - log(u) - log(v) + (th - 1) * (log(x) + log(y)) + (1 / th - 2) * log(x^th + y^th) + log(A + th - 1) },
    # theta and 1 - exp(-theta) share their sign, so take absolute values;
    # theta = 0 is independence.
    frank = ifelse(th == 0, 0, log(abs(th)) + log(abs(expm1(-th))) - th * (u + v) - 2 * log(abs(-expm1(-th) - expm1(-th * u) * expm1(-th * v)))))
  # h(u, v) = dC(u, v)/dv; the four copulas are exchangeable, so
  # dC/du at (u, v) is h(v, u).
  hfun <- function(u, v, th) { v <- pmin(pmax(v, 1e-15), 1 - 1e-15); uc <- pmin(pmax(u, 1e-15), 1 - 1e-15)
    out <- switch(copula,
      gaussian = stats::pnorm((stats::qnorm(uc) - th * stats::qnorm(v)) / sqrt(1 - th^2)),
      clayton = v^(-th - 1) * pmax(uc^-th + v^-th - 1, 0)^(-1 / th - 1),
      gumbel = { x <- -log(uc); z <- -log(v); s <- x^th + z^th; exp(-s^(1 / th)) * s^(1 / th - 1) * z^(th - 1) / v },
      frank = ifelse(th == 0, uc, exp(-th * v) * expm1(-th * uc) / (expm1(-th) + expm1(-th * uc) * expm1(-th * v))))
    out[u <= 0] <- 0; out[u >= 1] <- 1; out }
  # Each coordinate is either at a jump of its margin (a discrete value or
  # a point mass), where the copula is differenced across the jump, or
  # continuous, where it is differentiated: rectangle probabilities, h
  # function differences times a density, or the copula density.
  lpdf <- function(y, p) { y <- if (is.matrix(y)) y else matrix(y, nrow = 1); n <- nrow(y)
    pp <- lapply(p, rep_len, n); a <- p1(pp); b <- p2(pp); th <- pp[[dep$name]]
    l1 <- eval_lpdf(m1, y[, 1], a); l2 <- eval_lpdf(m2, y[, 2], b)
    u1 <- family_cdf(m1, y[, 1], a); v1 <- family_cdf(m2, y[, 2], b)
    j1 <- family_atom(m1, y[, 1], a) > 0; j2 <- family_atom(m2, y[, 2], b) > 0
    u0 <- ifelse(j1, family_cdf_left(m1, y[, 1], a), u1); v0 <- ifelse(j2, family_cdf_left(m2, y[, 2], b), v1)
    out <- rep(-Inf, n); ok <- is.finite(l1) & is.finite(l2)
    cc <- ok & !j1 & !j2; jj <- ok & j1 & j2; jc <- ok & j1 & !j2; cj <- ok & !j1 & j2
    cl <- function(x) pmin(pmax(x, 1e-15), 1 - 1e-15)
    if (any(cc)) out[cc] <- logc(cl(u1[cc]), cl(v1[cc]), th[cc]) + l1[cc] + l2[cc]
    if (any(jj)) out[jj] <- log(pmax(Cfun(u1[jj], v1[jj], th[jj]) - Cfun(u0[jj], v1[jj], th[jj]) -
                                     Cfun(u1[jj], v0[jj], th[jj]) + Cfun(u0[jj], v0[jj], th[jj]), 1e-300))
    if (any(jc)) out[jc] <- log(pmax(hfun(u1[jc], v1[jc], th[jc]) - hfun(u0[jc], v1[jc], th[jc]), 1e-300)) + l2[jc]
    if (any(cj)) out[cj] <- log(pmax(hfun(v1[cj], u1[cj], th[cj]) - hfun(v0[cj], u1[cj], th[cj]), 1e-300)) + l1[cj]
    out }
  f <- family_def(sprintf("copula_%s(%s, %s)", copula, m1$id, m2$id), sprintf("%s copula of %s and %s", copula, m1$display_name, m2$display_name),
    "multivariate", pars, lpdf = lpdf, discrete = m1$discrete && m2$discrete, variate_type = "vector", numeric_cdf = FALSE, numeric_moments = FALSE,
    check = component_checks(list(m1, m2), c(".1", ".2")),
    status = "implemented", formula = "log c(F1(y1), F2(y2)) + log f1(y1) + log f2(y2); rectangle probabilities for discrete margins",
    references = c("Nelsen (2006), An Introduction to Copulas, 2nd ed.", "Trivedi and Zimmer (2007), Copula modeling: an introduction for practitioners. Found. Trends Econom., 1, 1-111."),
    extra = list(primary = paste0(primary_par(m1), ".1")))
  f$operator <- list(type = "copula", copula = copula, margins = list(m1, m2))
  f
}
