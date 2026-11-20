# Operator algebra: families built from other families.

new_par_name <- function(base, nm) {
  while (nm %in% names(base$pars)) nm <- paste0(nm, "_")
  nm
}

sub_pars <- function(p, names) p[intersect(names(p), names)]

op_family <- function(base, id, display, pars, lpdf, cdf = NULL, lccdf = NULL,
                      quantile = NULL, rng = NULL, mean = NULL, variance = NULL,
                      support = NULL, discrete = base$discrete, stan = NULL, bugs = NULL,
                      operator, formula = "", extra = list(), numeric_moments = TRUE, check = NULL,
                      atom = NULL) {
  # Base constraints carry over when the base parameters keep their names.
  if (is.null(check) && !is.null(base$check) && all(names(base$pars) %in% names(pars))) {
    bn <- names(base$pars); check <- function(p) base$check(p[bn])
  }
  # A discrete result on the same support keeps the base lattice.
  if (discrete && base$discrete && is.null(extra$grid) && !is.null(base$extra$grid)) extra$grid <- base$extra$grid
  f <- family_def(
    id, display, base$category, pars, lpdf = lpdf, cdf = cdf, lccdf = lccdf,
    quantile = quantile, rng = rng, mean = mean, variance = variance,
    support = if (is.null(support)) base$support else support,
    discrete = discrete, variate_type = base$variate_type, formula = formula,
    references = base$references, stan = stan, bugs = bugs,
    status = if (base$status %in% c("validated", "production")) "implemented" else base$status,
    numeric_moments = numeric_moments, extra = extra, check = check, atom = atom)
  f$operator <- operator
  f
}

#' Inflated distributions
#'
#' Adds point masses to a base law: with weights `pi_j` at points `c_j`,
#' `P(Y = c_j) = pi_j + (1 - sum(pi)) P0(Y = c_j)` and elsewhere the base law
#' is scaled by `1 - sum(pi)`. For a continuous base law the result is a
#' mixed distribution whose log likelihood uses the point mass at `c_j` and
#' the scaled density elsewhere. One inflation parameter per point is added
#' (named `zi`, `oi` or `pi1`, `pi2`, ...), each in (0, 1) with a logit link.
#' The parameters are the point mass weights themselves; with several points
#' their sum must stay below one. Points must be distinct, and for a discrete
#' base law they must lie on its integer support.
#'
#' Zero inflation can only add zeros; compare [hurdle()].
#'
#' @param base Base family (identifier or object).
#' @param at Inflation points.
#' @param names Optional parameter names.
#' @return A `regcomp_family`.
#' @export
#' @examples
#' zip <- inflate("poisson", at = 0)
#' d <- distribution(zip, lambda = 2, zi = 0.3)
#' prob(d, 0:3)
#' mean(d)
#' zoib <- inflate("beta", at = c(0, 1))
#' zoib
inflate <- function(base, at = 0, names = NULL) {
  base <- get_family(base)
  if (base$variate_type != "scalar") stop_unsupported("Inflation of a vector valued law", base$id)
  k <- length(at)
  if (anyDuplicated(at)) stop_input("Inflation points must be distinct.")
  g <- lattice_step(base)
  if (base$discrete && (g == 0 || any(abs(at / g - round(at / g)) > 1e-9)))
    stop_input("Inflation points of a discrete law must lie on its lattice (integers for counts).")
  atomic <- base$discrete || has_atoms(base)
  if (is.null(names)) names <- if (k == 1 && at == 0) "zi" else if (k == 1 && at == 1) "oi" else paste0("pi", seq_len(k))
  names <- vapply(names, function(n) new_par_name(base, n), "")
  bnames <- names(base$pars)
  wts <- function(p) {
    if (k == 1) return(matrix(p[[names]], ncol = 1))
    W <- sapply(names, function(n) p[[n]])
    matrix(W, ncol = k)
  }
  lpdf <- function(y, p) {
    W <- wts(p); bp <- sub_pars(p, bnames)
    rest <- 1 - rowSums(W)
    lb <- eval_lpdf(base, y, bp)
    out <- log(rest) + lb
    for (j in seq_len(k)) {
      hit <- !is.na(y) & y == at[j]
      if (any(hit)) out[hit] <- if (atomic)
        log(W[hit, j] + rest[hit] * family_atom(base, y[hit], lapply(bp, function(v) rep_len(v, length(y))[hit])))
        else log(W[hit, j])
    }
    out
  }
  atom <- if (!base$discrete) function(y, p) {
    W <- wts(p); bp <- sub_pars(p, bnames)
    out <- (1 - rowSums(W)) * family_atom(base, y, bp)
    for (j in seq_len(k)) out <- out + W[, j] * (y == at[j])
    out
  }
  cdf <- function(q, p) {
    W <- wts(p); bp <- sub_pars(p, bnames)
    out <- (1 - rowSums(W)) * family_cdf(base, q, bp)
    for (j in seq_len(k)) out <- out + W[, j] * (q >= at[j])
    out
  }
  support <- function(p) {
    s <- base$support(sub_pars(p, bnames))
    list(lower = pmin(s$lower, min(at)), upper = pmax(s$upper, max(at)))
  }
  rng <- function(n, p) {
    W <- wts(p); bp <- sub_pars(p, bnames)
    y <- family_rng(base, n, bp)
    u <- stats::runif(n)
    cw <- t(apply(cbind(0, W), 1, cumsum))
    for (j in seq_len(k)) {
      hit <- u >= cw[, j] & u < cw[, j + 1]
      y[hit] <- at[j]
    }
    y
  }
  mean <- function(p) {
    W <- wts(p); m <- family_moment(base, "mean", sub_pars(p, bnames), quiet = TRUE)
    (1 - rowSums(W)) * m + drop(W %*% at)
  }
  variance <- function(p) {
    W <- wts(p); bp <- sub_pars(p, bnames)
    m <- family_moment(base, "mean", bp, quiet = TRUE)
    v <- family_moment(base, "variance", bp, quiet = TRUE)
    r <- 1 - rowSums(W)
    ey2 <- r * (v + m^2) + drop(W %*% at^2)
    ey <- r * m + drop(W %*% at)
    ey2 - ey^2
  }
  newpars <- c(base$pars, stats::setNames(lapply(names, function(n)
    P(n, "unit", role = "inflation", interpretation = "Probability of the added point mass")), names))
  stan <- NULL
  if (!is.null(base$stan) && k == 1 && !has_atoms(base)) {
    fn <- stan_fn_name(paste0("inflate_", base$id), paste(base$id, format(at)))
    args <- stan_args(base$pars, bnames, stan_ident(bnames))
    pm <- stats::setNames(as.list(stan_ident(bnames)), bnames)
    yt <- if (base$discrete) "int" else "real"
    blp <- stan_fill(base$stan$lpdf, "y", pm)
    body <- if (base$discrete) sprintf(
      "    if (y == %s) return log_sum_exp(log(pi), log1m(pi) + %s);\n    return log1m(pi) + %s;",
      format(at), blp, blp) else sprintf(
      "    if (y == %s) return log(pi);\n    return log1m(pi) + %s;", format(at), blp)
    lc <- if (!is.null(base$stan$lcdf)) {
      bl <- stan_fill(base$stan$lcdf, "y", pm)
      sprintf("
  real %1$s_lcdf(%2$s y, %3$s, real pi) {
    if (y >= %4$s) return log_sum_exp(log(pi), log1m(pi) + %5$s);
    return log1m(pi) + %5$s;
  }
  real %1$s_lccdf(%2$s y, %3$s, real pi) {
    return log1m_exp(%1$s_lcdf(y | %6$s, pi));
  }", fn, yt, args, format(at), bl, paste(stan_ident(bnames), collapse = ", "))
    } else ""
    stan <- list(native = FALSE, needs = base$stan$needs, constraint = base$stan$constraint,
      functions = paste0(base$stan$functions, sprintf("
  real %s_%s(%s y, %s, real pi) {
%s
  }", fn, if (base$discrete) "lpmf" else "lpdf", yt, args, body), lc),
      lpdf = sprintf("%s_%s({y} | %s, {%s})", fn, if (base$discrete) "lpmf" else "lpdf",
                     paste(sprintf("{%s}", bnames), collapse = ", "), names),
      lcdf = if (nzchar(lc)) sprintf("%s_lcdf({y} | %s, {%s})", fn, paste(sprintf("{%s}", bnames), collapse = ", "), names),
      lccdf = if (nzchar(lc)) sprintf("%s_lccdf({y} | %s, {%s})", fn, paste(sprintf("{%s}", bnames), collapse = ", "), names))
  }
  bugs <- NULL
  if (!is.null(base$bugs$loglik) && k == 1 && !has_atoms(base)) {
    bl <- base$bugs$loglik
    bugs <- list(loglik = if (base$discrete)
      sprintf("log(equals({y}, %s) * {%s} + (1 - {%s}) * exp(%s))", format(at), names, names, bl)
      else sprintf("equals({y}, %s) * log({%s}) + (1 - equals({y}, %s)) * (log(1 - {%s}) + %s)",
                   format(at), names, format(at), names,
                   gsub("{y}", sprintf("({y} + equals({y}, %s))", format(at)), bl, fixed = TRUE)))
  }
  f <- op_family(base, sprintf("inflate(%s, at = %s)", base$id, paste(format(at), collapse = ",")),
            sprintf("%s inflated at %s", base$display_name, paste(format(at), collapse = ", ")),
            newpars, lpdf, cdf = cdf, rng = rng, mean = mean, variance = variance,
            support = support, discrete = base$discrete, stan = stan, bugs = bugs, atom = atom,
            extra = if (!base$discrete && (!has_atoms(base) || !is.null(base$extra$atom_points)))
              list(atom_points = sort(unique(c(at, base$extra$atom_points)))) else list(),
            operator = list(type = "inflate", base = base, at = at, pars = names),
            check = function(p) { m <- if (is.null(base$check)) NULL else base$check(sub_pars(p, bnames))
              if (!is.null(m)) return(m)
              if (k > 1 && any(Reduce(`+`, p[names]) >= 1, na.rm = TRUE)) "Inflation weights must sum to less than 1." },
            formula = sprintf("P(Y = c) = pi + (1 - pi) P0(Y = c); elsewhere (1 - pi) f0(y)"))
  # Starts: observed proportions at the points (sum below 0.9), base law
  # starts from the remaining observations.
  f$start <- function(y) {
    pr <- pmin(pmax(vapply(at, function(a) base::mean(y == a), 0), 0.01), 0.9 / k)
    yr <- y[!y %in% at]; if (!length(yr)) yr <- y
    st <- if (!is.null(base$start)) tryCatch(base$start(yr), error = function(e) list()) else list()
    for (b in bnames) if (is.null(st[[b]])) st[[b]] <- generic_start(base$pars[[b]], yr)
    c(st, stats::setNames(as.list(pr), names))
  }
  f
}

#' Hurdle distributions
#'
#' Replaces the base probability at `at` by a separate hurdle probability:
#' `P(Y = at) = pi` and, for `y != at`, `P(Y = y) = (1 - pi) f0(y) / (1 -
#' P0(Y = at))`. Unlike inflation, a hurdle can produce fewer as well as
#' more values at the hurdle than the base law. For a continuous base law
#' the conditioning term is one, giving a two part (zero adjusted) model such
#' as the zero adjusted gamma.
#'
#' @param base Base family.
#' @param at Hurdle point.
#' @param name Name of the hurdle probability parameter.
#' @return A `regcomp_family`.
#' @export
#' @examples
#' hnb <- hurdle("negative_binomial_nb2", at = 0)
#' d <- distribution(hnb, mu = 3, phi = 1, hu = 0.1)
#' prob(d, 0:3)
#' zaga <- hurdle("gamma", at = 0)
hurdle <- function(base, at = 0, name = "hu") {
  base <- get_family(base)
  if (base$variate_type != "scalar") stop_unsupported("Hurdles for a vector valued law", base$id)
  g <- lattice_step(base)
  if (base$discrete && (g == 0 || abs(at / g - round(at / g)) > 1e-9))
    stop_input("The hurdle point of a discrete law must lie on its lattice.")
  name <- new_par_name(base, name)
  bnames <- names(base$pars)
  # Log base probability of the hurdle point (an atom of a mixed base too).
  lp_at <- function(bp) { n <- max(lengths(bp))
    if (base$discrete) eval_lpdf(base, rep(at, n), bp) else log(family_atom(base, rep(at, n), bp)) }
  lpdf <- function(y, p) {
    bp <- sub_pars(p, bnames); pi <- p[[name]]
    out <- log1p(-pi) + eval_lpdf(base, y, bp) - log1mexp(-lp_at(bp))
    hit <- !is.na(y) & y == at
    out[hit] <- log(pi[hit])
    out
  }
  cdf <- function(q, p) {
    bp <- sub_pars(p, bnames); pi <- p[[name]]
    lpa <- lp_at(bp); l1 <- log1mexp(-lpa)
    # Conditional law off the hurdle point; above it use the survival
    # function so a base law concentrated near the point keeps precision.
    above <- q >= at
    cond <- ifelse(above, 1 - exp(family_lccdf(base, q, bp) - l1), exp(log(family_cdf(base, q, bp)) - l1))
    pi * above + (1 - pi) * pmin(pmax(cond, 0), 1)
  }
  atom <- if (!base$discrete) function(y, p) {
    bp <- sub_pars(p, bnames); pi <- p[[name]]
    ifelse(y == at, pi, (1 - pi) * family_atom(base, y, bp) / -expm1(lp_at(bp)))
  }
  support <- function(p) {
    s <- base$support(sub_pars(p, bnames))
    list(lower = pmin(s$lower, at), upper = pmax(s$upper, at))
  }
  rng <- function(n, p) {
    bp <- sub_pars(p, bnames); pi <- p[[name]]
    y <- rep(at, n)
    todo <- which(stats::runif(n) >= pi); tries <- 0
    # Rows whose base law has almost no mass off the hurdle point: draw from
    # the conditional law directly by summing masses along the lattice.
    lpa <- rep_len(lp_at(bp), n)
    direct <- todo[base$discrete & is.finite(base$support(lapply(bp, `[`, 1))$lower) & lpa[todo] > log(0.99)]
    for (i in direct) {
      bi <- lapply(bp, `[`, i); l1 <- log(-expm1(lpa[i])); u <- stats::runif(1)
      k <- base$support(bi)$lower; acc <- 0
      for (s in 1:100000) {
        if (k != at) { acc <- acc + exp(eval_lpdf(base, k, bi) - l1); if (acc >= u) break }
        k <- k + g
      }
      y[i] <- k
    }
    todo <- setdiff(todo, direct)
    while (length(todo)) {
      tries <- tries + 1
      if (tries > 1000) stop_numerical("Hurdle sampling: the base law has almost no mass away from the hurdle point.")
      d <- family_rng(base, length(todo), lapply(bp, `[`, todo))
      ok <- d != at
      y[todo[ok]] <- d[ok]
      todo <- todo[!ok]
    }
    y
  }
  mean <- function(p) {
    bp <- sub_pars(p, bnames); pi <- p[[name]]
    m <- family_moment(base, "mean", bp, quiet = TRUE); la <- lp_at(bp); pa <- exp(la)
    pi * at + (1 - pi) * (m - pa * at) / -expm1(la)
  }
  variance <- function(p) {
    bp <- sub_pars(p, bnames); pi <- p[[name]]
    m <- family_moment(base, "mean", bp, quiet = TRUE)
    v <- family_moment(base, "variance", bp, quiet = TRUE); la <- lp_at(bp); pa <- exp(la); q1 <- -expm1(la)
    m1 <- (m - pa * at) / q1; m2 <- (v + m^2 - pa * at^2) / q1
    ey <- pi * at + (1 - pi) * m1; ey2 <- pi * at^2 + (1 - pi) * m2
    ey2 - ey^2
  }
  hcheck <- function(p) { bp <- sub_pars(p, bnames)
    m <- if (is.null(base$check)) NULL else base$check(bp)
    if (!is.null(m)) return(m)
    if (any(lp_at(bp) >= 0, na.rm = TRUE)) "The base law has no mass away from the hurdle point." }
  newpars <- c(base$pars, stats::setNames(list(P(name, "unit", role = "hurdle",
    interpretation = "Probability of the hurdle value")), name))
  stan <- NULL
  if (!is.null(base$stan) && !has_atoms(base)) {
    fn <- stan_fn_name(paste0("hurdle_", base$id), paste(base$id, format(at)))
    args <- stan_args(base$pars, bnames, stan_ident(bnames))
    pm <- stats::setNames(as.list(stan_ident(bnames)), bnames)
    yt <- if (base$discrete) "int" else "real"
    blp <- stan_fill(base$stan$lpdf, "y", pm)
    b0 <- if (base$discrete) sprintf("log1m_exp(%s)", stan_fill(base$stan$lpdf, format(at), pm)) else "0"
    stan <- list(native = FALSE, needs = base$stan$needs, constraint = base$stan$constraint,
      functions = paste0(base$stan$functions, sprintf("
  real %s_%s(%s y, %s, real pi) {
    if (y == %s) return log(pi);
    return log1m(pi) + %s - %s;
  }", fn, if (base$discrete) "lpmf" else "lpdf", yt, args, format(at), blp, b0)),
      lpdf = sprintf("%s_%s({y} | %s, {%s})", fn, if (base$discrete) "lpmf" else "lpdf",
                     paste(sprintf("{%s}", bnames), collapse = ", "), name))
  }
  bugs <- NULL
  if (!is.null(base$bugs$loglik) && base$discrete) {
    bl <- base$bugs$loglik
    b0 <- gsub("{y}", format(at), bl, fixed = TRUE)
    bugs <- list(loglik = sprintf(
      "equals({y}, %s) * log({%s}) + (1 - equals({y}, %s)) * (log(1 - {%s}) + %s - log(1 - exp(%s)))",
      format(at), name, format(at), name, gsub("{y}", sprintf("({y} + equals({y}, %s))", format(at)), bl, fixed = TRUE), b0))
  }
  op_family(base, sprintf("hurdle(%s, at = %s)", base$id, format(at)),
            sprintf("%s hurdle at %s", base$display_name, format(at)),
            newpars, lpdf, cdf = cdf, rng = rng, mean = mean, variance = variance,
            support = support, stan = stan, bugs = bugs, atom = atom,
            extra = if (!base$discrete && (!has_atoms(base) || !is.null(base$extra$atom_points)))
              list(atom_points = sort(unique(c(at, base$extra$atom_points)))) else list(),
            operator = list(type = "hurdle", base = base, at = at, pars = name), check = hcheck,
            formula = "P(Y = c) = pi; P(Y = y) = (1 - pi) f0(y) / (1 - P0(c)) otherwise")
}

#' Truncated distributions
#'
#' Restricts a base law to `[lower, upper]` (discrete; both end points kept)
#' or `(lower, upper]` (continuous) and renormalizes:
#' `f_T(y) = f(y) / (F(upper) - F(lower-))`, computed with log scale
#' differences. For observation specific truncation use [cens()].
#'
#' @param base Base family.
#' @param lower,upper Truncation bounds.
#' @return A `regcomp_family`.
#' @export
#' @examples
#' ztp <- truncate_dist("poisson", lower = 1)
#' d <- distribution(ztp, lambda = 1.5)
#' prob(d, 0:3); mean(d)
#' tn <- truncate_dist("normal", lower = 0)
truncate_dist <- function(base, lower = -Inf, upper = Inf) {
  base <- get_family(base)
  atomic <- has_atoms(base)
  if (lower == -Inf && upper == Inf) return(base)
  if (base$variate_type != "scalar") stop_unsupported("Truncation of a vector valued law", base$id)
  if (lower >= upper) stop_input("lower must be smaller than upper.")
  bn <- names(base$pars)
  if (base$discrete) { lower <- ceiling(lower); upper <- floor(upper) }
  lo_q <- if (base$discrete) lower - lattice_step(base) else lower
  # Log normalizer and conditional probabilities in log space, so a region
  # far in a tail keeps its mass.
  # A law with point masses keeps the atoms at both bounds: the region is
  # [lower, upper], with mass F(upper) - F(lower-).
  Fleft <- function(n, pp) if (is.finite(lower)) family_cdf_left(base, rep(lower, n), pp) else rep(0, n)
  lnorm <- function(p) { n <- max(lengths(p)); pp <- recycle_pars(p, n)
    if (atomic) return(log(family_cdf(base, rep(upper, n), pp) - Fleft(n, pp)))
    family_log_interval(base, rep(lo_q, n), rep(upper, n), pp) }
  lpdf <- function(y, p) eval_lpdf(base, y, p) - lnorm(p)
  cdf <- function(q, p) {
    n <- length(q); pp <- recycle_pars(p, n)
    qq <- pmin(pmax(q, lower), upper)
    out <- if (atomic) (family_cdf(base, qq, pp) - Fleft(n, pp)) / exp(lnorm(pp))
           else exp(family_log_interval(base, rep(lo_q, n), qq, pp) - lnorm(pp))
    out[q < lower] <- 0; out[q >= upper] <- 1
    pmin(pmax(out, 0), 1)
  }
  quantile <- function(pr, p) {
    n <- length(pr); pp <- recycle_pars(p, n)
    Fl <- if (atomic) Fleft(n, pp) else if (is.finite(lo_q)) family_cdf(base, rep(lo_q, n), pp) else rep(0, n)
    Fu <- if (is.finite(upper)) family_cdf(base, rep(upper, n), pp) else rep(1, n)
    out <- family_quantile(base, Fl + pr * (Fu - Fl), pp)
    # Where the region is too far in a tail for that, solve on the log scale.
    bad <- which(!(Fu - Fl > 1e-10) | !is.finite(out) | out < lower | out > upper)
    for (i in bad) {
      pi <- lapply(pp, `[`, i)
      g <- function(q) family_log_interval(base, lo_q, q, pi) - lnorm(pi) - log(pr[i])
      if (pr[i] <= 0) { out[i] <- lower; next }
      if (pr[i] >= 1) { out[i] <- upper; next }
      hi <- if (is.finite(upper)) upper else lower + 1
      while (!is.finite(upper) && g(hi) < 0) hi <- lower + 2 * (hi - lower)
      lo <- if (is.finite(lower)) lower else hi - 1
      while (!is.finite(lower) && g(lo) > 0) lo <- hi - 2 * (hi - lo)
      if (base$discrete) {
        lo <- floor(lo); hi <- ceiling(hi)
        while (hi - lo > 1) { mid <- floor((lo + hi) / 2); if (g(mid) >= 0) hi <- mid else lo <- mid }
        out[i] <- if (g(lo) >= 0) lo else hi
      } else out[i] <- stats::uniroot(g, c(lo, hi), tol = 1e-12)$root
    }
    out
  }
  support <- function(p) {
    s <- base$support(p)
    list(lower = pmax(s$lower, lower), upper = pmin(s$upper, upper))
  }
  stan <- NULL
  # The Stan normalizer F(upper) - F(lower) would drop a point mass at lower.
  if (!is.null(base$stan) && !is.null(base$stan$lcdf) && !atomic) {
    fn <- stan_fn_name(paste0("trunc_", base$id), paste(base$id, format(lower), format(upper)))
    pm <- stats::setNames(as.list(stan_ident(bn)), bn)
    args <- stan_args(base$pars, bn, stan_ident(bn))
    yt <- if (base$discrete) "int" else "real"
    nrm <- if (is.finite(lower) && is.finite(upper))
      sprintf("log_diff_exp(%s, %s)", stan_fill(base$stan$lcdf, format(upper), pm), stan_fill(base$stan$lcdf, format(lo_q), pm))
      else if (is.finite(lower)) stan_fill(base$stan$lccdf, format(lo_q), pm)
      else stan_fill(base$stan$lcdf, format(upper), pm)
    guard <- paste(c(if (is.finite(lower)) sprintf("y < %s", format(lower)),
                     if (is.finite(upper)) sprintf("y > %s", format(upper))), collapse = " || ")
    stan <- list(native = FALSE, needs = base$stan$needs, constraint = base$stan$constraint,
      functions = paste0(base$stan$functions, sprintf("
  real %s_%s(%s y, %s) {
    if (%s) return negative_infinity();
    return %s - %s;
  }", fn, if (base$discrete) "lpmf" else "lpdf", yt, args, guard, stan_fill(base$stan$lpdf, "y", pm), nrm)),
      lpdf = sprintf("%s_%s({y} | %s)", fn, if (base$discrete) "lpmf" else "lpdf",
                     paste(sprintf("{%s}", bn), collapse = ", ")))
  }
  op_family(base, sprintf("truncate(%s, %s, %s)", base$id, format(lower), format(upper)),
            sprintf("%s truncated to [%s, %s]", base$display_name, format(lower), format(upper)),
            base$pars, lpdf, cdf = cdf, quantile = quantile,
            rng = function(n, p) quantile(stats::runif(n), p),
            support = support, stan = stan,
            atom = if (atomic) function(y, p) ifelse(y >= lower & y <= upper, family_atom(base, y, p) / exp(lnorm(p)), 0),
            extra = if (atomic && !is.null(base$extra$atom_points)) {
              pts <- base$extra$atom_points; list(atom_points = pts[pts >= lower & pts <= upper]) } else list(),
            operator = list(type = "truncate", base = base, lower = lower, upper = upper),
            formula = "f(y) / (F(upper) - F(lower-)) on the retained region")
}

#' Finite mixtures
#'
#' Mixture of `K` component families with weights `w_k`:
#' `f(y) = sum_k w_k f_k(y)`. Component parameters are renamed with suffixes
#' `.1`, `.2`, ... (for example `mu.1`). Weights are parameters `w2, ...,
#' wK` on the multinomial logit scale relative to the first component
#' (`w1` is implied), each with an identity link. Likelihoods are
#' marginalized over the component labels; label switching is possible and
#' is reported by `identifiability_check()`.
#'
#' @param ... Component families.
#' @return A `regcomp_family`.
#' @export
#' @examples
#' mx <- mixture("normal", "normal")
#' d <- distribution(mx, mu.1 = -2, sigma.1 = 1, mu.2 = 2, sigma.2 = 0.5, w2 = 0)
#' plot(d)
mixture <- function(...) {
  comps <- lapply(list(...), get_family)
  K <- length(comps)
  if (K < 2) stop_input("A mixture needs at least two components.")
  disc <- vapply(comps, `[[`, TRUE, "discrete")
  if (length(unique(disc)) > 1) stop_input("Mixture components must all be discrete or all continuous.")
  grids <- unique(vapply(comps, lattice_step, 0))
  if (disc[1] && length(grids) > 1) stop_unsupported("A mixture of discrete laws on different lattices")
  atomic <- !disc[1] && any(vapply(comps, has_atoms, TRUE))
  pars <- list(); cmap <- vector("list", K)
  for (k in seq_len(K)) {
    for (p in comps[[k]]$pars) {
      nm <- paste0(p$name, ".", k)
      q <- p; q$name <- nm; q$symbol <- nm; q$is_mean <- FALSE
      pars[[nm]] <- q
      cmap[[k]] <- c(cmap[[k]], stats::setNames(nm, p$name))
    }
  }
  for (k in 2:K) pars[[paste0("w", k)]] <- P(paste0("w", k), "real", role = "mixing",
    interpretation = sprintf("Log odds of component %d against component 1", k))
  weights <- function(p) {
    n <- max(lengths(p))
    E <- cbind(0, matrix(vapply(2:K, function(k) rep_len(p[[paste0("w", k)]], n), numeric(n)), nrow = n))
    E <- matrix(E, nrow = n)
    E <- E - apply(E, 1, max)
    W <- exp(E); W / rowSums(W)
  }
  comp_pars <- function(p, k) {
    m <- cmap[[k]]
    stats::setNames(lapply(m, function(nm) p[[nm]]), names(m))
  }
  atom <- if (atomic) function(y, p) {
    W <- weights(p)
    rowSums(matrix(vapply(seq_len(K), function(k) W[, k] * family_atom(comps[[k]], y, comp_pars(p, k)), numeric(length(y))), ncol = K))
  }
  lpdf <- function(y, p) {
    W <- weights(p)
    L <- sapply(seq_len(K), function(k) log(W[, k]) + eval_lpdf(comps[[k]], y, comp_pars(p, k)))
    L <- matrix(L, ncol = K)
    m <- apply(L, 1, max)
    out <- m + log(rowSums(exp(L - m)))
    out[!is.finite(m)] <- m[!is.finite(m)]
    # At an atom of any component the law has a point mass there; the
    # densities of the other components do not add to it.
    if (atomic) { a <- atom(y, p); out[!is.na(a) & a > 0] <- log(a[!is.na(a) & a > 0]) }
    out
  }
  cdf <- function(q, p) {
    W <- weights(p)
    rowSums(matrix(sapply(seq_len(K), function(k) W[, k] * family_cdf(comps[[k]], q, comp_pars(p, k))), ncol = K))
  }
  rng <- function(n, p) {
    W <- weights(p)
    z <- apply(W, 1, function(w) sample.int(K, 1, prob = w))
    y <- numeric(n)
    for (k in seq_len(K)) {
      idx <- which(z == k)
      if (length(idx)) y[idx] <- family_rng(comps[[k]], length(idx), lapply(comp_pars(p, k), `[`, idx))
    }
    y
  }
  mean <- function(p) {
    W <- weights(p)
    rowSums(matrix(sapply(seq_len(K), function(k) W[, k] * family_moment(comps[[k]], "mean", comp_pars(p, k), quiet = TRUE)), ncol = K))
  }
  variance <- function(p) {
    W <- weights(p)
    m <- sapply(seq_len(K), function(k) family_moment(comps[[k]], "mean", comp_pars(p, k), quiet = TRUE))
    v <- sapply(seq_len(K), function(k) family_moment(comps[[k]], "variance", comp_pars(p, k), quiet = TRUE))
    m <- matrix(m, ncol = K); v <- matrix(v, ncol = K)
    mu <- rowSums(W * m)
    rowSums(W * (v + m^2)) - mu^2
  }
  support <- function(p) {
    s <- lapply(seq_len(K), function(k) comps[[k]]$support(comp_pars(p, k)))
    list(lower = do.call(pmin, lapply(s, `[[`, "lower")), upper = do.call(pmax, lapply(s, `[[`, "upper")))
  }
  stan <- NULL
  if (all(vapply(comps, function(f) !is.null(f$stan), TRUE)) && !atomic) {
    fn <- stan_fn_name(paste0("mix_", paste(vapply(comps, `[[`, "", "id"), collapse = "_")),
                       paste(vapply(comps, `[[`, "", "id"), collapse = "\r"))
    yt <- if (disc[1]) "int" else "real"
    args <- stan_args(pars, names(pars), stan_ident(names(pars)))
    terms <- vapply(seq_len(K), function(k) {
      pm <- stats::setNames(as.list(stan_ident(cmap[[k]])), names(cmap[[k]]))
      sprintf("log(w[%d]) + %s", k, stan_fill(comps[[k]]$stan$lpdf, "y", pm))
    }, "")
    funs <- unique(vapply(comps, function(f) f$stan$functions, ""))
    # Component constraints, with the component parameters renamed.
    cons <- unlist(lapply(seq_len(K), function(k) if (!is.null(comps[[k]]$stan$constraint))
      stan_fill(comps[[k]]$stan$constraint, "y", stats::setNames(as.list(sprintf("{%s}", cmap[[k]])), names(cmap[[k]])))))
    stan <- list(native = FALSE, needs = unique(unlist(lapply(comps, function(f) f$stan$needs))),
      constraint = if (length(cons)) paste(sprintf("(%s)", cons), collapse = " && "),
      functions = paste0(paste(funs, collapse = "\n"), sprintf("
  real %s_%s(%s y, %s) {
    vector[%d] w = softmax([0, %s]');
    vector[%d] lp;
%s
    return log_sum_exp(lp);
  }", fn, if (disc[1]) "lpmf" else "lpdf", yt, args, K,
        paste(stan_ident(paste0("w", 2:K)), collapse = ", "), K,
        paste(sprintf("    lp[%d] = %s;", seq_len(K), terms), collapse = "\n"))),
      lpdf = sprintf("%s_%s({y} | %s)", fn, if (disc[1]) "lpmf" else "lpdf",
                     paste(sprintf("{%s}", names(pars)), collapse = ", ")))
  }
  f <- family_def(
    paste0("mixture(", paste(vapply(comps, `[[`, "", "id"), collapse = ", "), ")"),
    paste("Mixture of", paste(vapply(comps, `[[`, "", "display_name"), collapse = " and ")),
    comps[[1]]$category, pars, lpdf = lpdf, cdf = cdf, rng = rng, mean = mean,
    variance = variance, support = support, discrete = disc[1], stan = stan, atom = atom,
    check = function(p) { for (k in seq_len(K)) if (!is.null(comps[[k]]$check)) {
      m <- comps[[k]]$check(comp_pars(p, k)); if (is.character(m) && length(m) && nzchar(m[1])) return(sprintf("Component %d: %s", k, m[1])) }
      NULL },
    extra = if (disc[1]) list(grid = grids[1]) else if (atomic && all(vapply(comps, function(f) !has_atoms(f) || !is.null(f$extra$atom_points), TRUE)))
      list(atom_points = sort(unique(unlist(lapply(comps, function(f) f$extra$atom_points))))) else list(),
    status = "implemented", formula = "f(y) = sum_k w_k f_k(y)",
    identifiability_notes = "Component labels are exchangeable when components share a family (label switching).")
  f$operator <- list(type = "mixture", components = comps)
  f
}

#' Monotone transformations
#'
#' The law of `Y = g(X)` for a monotone, differentiable transformation `g`
#' of a continuous base law: `log f_Y(y) = log f_X(g^{-1}(y)) + log |d
#' g^{-1}(y) / dy|`. Named transformations: `"log"` (`Y = log X`), `"exp"`
#' (`Y = exp X`), `"logit"` (`Y = logit X`), `"inv_logit"`, `"reciprocal"`
#' (`Y = 1 / X` on positive support) and `"sqrt"`; or give `forward`,
#' `inverse` and `log_jacobian` functions yourself.
#'
#' @param base Base family.
#' @param transform Name of a built in transformation, or `NULL` when the
#'   functions are given.
#' @param forward,inverse Functions `g` and `g^{-1}`.
#' @param log_jacobian Function of `y` returning `log |d g^{-1}(y) / dy|`.
#' @param increasing Whether `g` is increasing.
#' @param support Optional function of the parameters returning the new
#'   support bounds.
#' @return A `regcomp_family`.
#' @export
#' @examples
#' ln <- transform_dist("normal", "exp")   # a lognormal law
#' d <- distribution(ln, mu = 0, sigma = 1)
#' prob(d, 1); dlnorm(1)
transform_dist <- function(base, transform = c("log", "exp", "logit", "inv_logit", "reciprocal", "sqrt"),
                           forward = NULL, inverse = NULL, log_jacobian = NULL,
                           increasing = TRUE, support = NULL) {
  base <- get_family(base)
  if (base$discrete) stop_unsupported("Monotone transformation of a discrete law (use discretize())", base$id)
  atomic <- has_atoms(base)
  if (atomic && is.null(base$extra$atom_points))
    stop_unsupported("Monotone transformation of a law whose point mass locations are unknown", base$id)
  label <- "custom"
  if (is.null(forward)) {
    transform <- match.arg(transform)
    label <- transform
    tr <- switch(transform,
      log = list(f = log, i = exp, j = function(y) y, inc = TRUE),
      exp = list(f = exp, i = log, j = function(y) -log(y), inc = TRUE),
      logit = list(f = stats::qlogis, i = stats::plogis, j = function(y) stats::dlogis(y, log = TRUE), inc = TRUE),
      inv_logit = list(f = stats::plogis, i = stats::qlogis, j = function(y) -log(y) - log1p(-y), inc = TRUE),
      reciprocal = list(f = function(x) 1 / x, i = function(y) 1 / y, j = function(y) -2 * log(y), inc = FALSE),
      sqrt = list(f = sqrt, i = function(y) y^2, j = function(y) log(2 * y), inc = TRUE))
    forward <- tr$f; inverse <- tr$i; log_jacobian <- tr$j; increasing <- tr$inc
  }
  sup <- if (!is.null(support)) support else function(p) {
    s <- base$support(p)
    a <- suppressWarnings(forward(s$lower)); b <- suppressWarnings(forward(s$upper))
    if (anyNA(c(a, b))) stop_domain("The transformation is not defined on the whole support of the base law.")
    list(lower = pmin(a, b), upper = pmax(a, b))
  }
  # Point masses move to forward(c) with their mass and no Jacobian; the
  # base atom is looked up at its exact location (inverse() may round).
  pts <- base$extra$atom_points
  at_pts <- function(y) { if (!atomic) return(rep(NA_real_, length(y)))
    fp <- suppressWarnings(forward(pts))
    vapply(y, function(v) { k <- which(abs(fp - v) <= 1e-10 * pmax(1, abs(v)) | fp == v); if (length(k)) pts[k[1]] else NA_real_ }, 0) }
  atom <- if (atomic) function(y, p) { c0 <- at_pts(y); out <- rep(0, length(y)); k <- !is.na(c0)
    if (any(k)) out[k] <- family_atom(base, c0[k], lapply(p, function(v) rep_len(v, length(y))[k])); out }
  lpdf <- function(y, p) {
    out <- eval_lpdf(base, inverse(y), p) + log_jacobian(y)
    if (atomic) { c0 <- at_pts(y); k <- !is.na(c0)
      if (any(k)) out[k] <- eval_lpdf(base, c0[k], lapply(p, function(v) rep_len(v, length(y))[k])) }
    out
  }
  cdf <- function(q, p) {
    x <- inverse(q)
    if (increasing) family_cdf(base, x, p) else 1 - family_cdf_left(base, x, p)
  }
  quantile <- function(pr, p) forward(family_quantile(base, if (increasing) pr else 1 - pr, p))
  rng <- function(n, p) forward(family_rng(base, n, p))
  op_family(base, sprintf("transform(%s, %s)", base$id, label),
            sprintf("%s transformed (%s)", base$display_name, label), base$pars,
            lpdf, cdf = cdf, quantile = quantile, rng = rng, support = sup, atom = atom,
            extra = if (atomic) list(atom_points = sort(suppressWarnings(forward(pts)))) else list(),
            operator = list(type = "transform", base = base, transform = label),
            formula = "log f_Y(y) = log f_X(g^-1(y)) + log |d g^-1(y)/dy|")
}

#' Discretized and rounded distributions
#'
#' `discretize()` turns a continuous law into a discrete one on the integers
#' by interval probabilities `P(Y = k) = F(k + 1 - offset) - F(k - offset)`
#' (the default `offset = 0` maps `[k, k + 1)` to `k`; `offset = 0.5` rounds
#' to the nearest integer). `round_dist()` gives the law of a value rounded
#' to `digits` decimals, `P(Y = y) = F(y + h/2) - F(y - h/2)` with
#' `h = 10^-digits`, on the grid of multiples of `h`.
#'
#' @param base Continuous base family.
#' @param offset Shift of the interval boundaries.
#' @param digits Number of decimals.
#' @return A `regcomp_family`.
#' @export
#' @examples
#' dg <- discretize("gamma")
#' d <- distribution(dg, mu = 3, sigma = 0.5)
#' sum(prob(d, 0:200))
#' r1 <- round_dist("normal", digits = 1)
discretize <- function(base, offset = 0) {
  base <- get_family(base)
  if (base$discrete) stop_input("The base law is already discrete.")
  # floor(X + offset) = k exactly when k - offset <= X < k + 1 - offset, so
  # P(Y = k) = F((k + 1 - offset)-) - F((k - offset)-): left limits keep an
  # atom of a mixed base in its own bin. Without atoms the log interval
  # probability keeps tail masses.
  lpdf <- function(y, p) {
    if (!has_atoms(base)) return(family_log_interval(base, y - offset, y + 1 - offset, recycle_pars(p, length(y))))
    log(pmax(family_cdf_left(base, y + 1 - offset, p) - family_cdf_left(base, y - offset, p), 0))
  }
  cdf <- function(q, p) family_cdf_left(base, floor(q) + 1 - offset, p)
  support <- function(p) {
    s <- base$support(p)
    list(lower = floor(s$lower + offset), upper = ifelse(is.finite(s$upper), floor(s$upper + offset), Inf))
  }
  rng <- function(n, p) floor(family_rng(base, n, p) + offset)
  op_family(base, sprintf("discretize(%s)", base$id), sprintf("Discretized %s", base$display_name),
            base$pars, lpdf, cdf = cdf, rng = rng, support = support, discrete = TRUE,
            operator = list(type = "discretize", base = base, offset = offset),
            formula = "P(Y = k) = F(k + 1 - offset) - F(k - offset)")
}

#' @rdname discretize
#' @export
round_dist <- function(base, digits = 0) {
  base <- get_family(base)
  if (base$discrete) stop_input("The base law is already discrete.")
  # Bins are (y - h/2, y + h/2]; point masses of the base fall in their bin.
  h <- 10^-digits
  f <- family_def(
    sprintf("round(%s, %d)", base$id, as.integer(digits)), sprintf("Rounded %s", base$display_name),
    base$category, base$pars, discrete = TRUE,
    lpdf = function(y, p) family_log_interval(base, y - h / 2, y + h / 2, recycle_pars(p, length(y))),
    cdf = function(q, p) family_cdf(base, floor(q / h + 1e-9) * h + h / 2, p),
    quantile = function(pr, p) ceiling((family_quantile(base, pr, p) - h / 2) / h - 1e-9) * h,
    rng = function(n, p) ceiling(family_rng(base, n, p) / h - 0.5) * h,
    support = function(p) { s <- base$support(p); g <- function(x) ifelse(is.finite(x), ceiling((x - h / 2) / h - 1e-9) * h, x)
      list(lower = g(s$lower), upper = g(s$upper)) },
    status = "implemented", numeric_moments = FALSE,
    numeric_cdf = FALSE,
    formula = "P(Y = y) = F(y + h/2) - F(y - h/2), y on the grid of multiples of h")
  f$extra$grid <- h
  f$operator <- list(type = "round", base = base, digits = digits)
  f
}
