# Operator algebra: families built from other families.

new_par_name <- function(base, nm) {
  while (nm %in% names(base$pars)) nm <- paste0(nm, "_")
  nm
}

sub_pars <- function(p, names) p[intersect(names(p), names)]

op_family <- function(base, id, display, pars, lpdf, cdf = NULL, lccdf = NULL,
                      quantile = NULL, rng = NULL, mean = NULL, variance = NULL,
                      support = NULL, discrete = base$discrete, stan = NULL, bugs = NULL,
                      operator, formula = "", extra = list(), numeric_moments = TRUE) {
  f <- family_def(
    id, display, base$category, pars, lpdf = lpdf, cdf = cdf, lccdf = lccdf,
    quantile = quantile, rng = rng, mean = mean, variance = variance,
    support = if (is.null(support)) base$support else support,
    discrete = discrete, variate_type = base$variate_type, formula = formula,
    references = base$references, stan = stan, bugs = bugs,
    status = if (base$status %in% c("validated", "production")) "implemented" else base$status,
    numeric_moments = numeric_moments, extra = extra)
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
#' (named `zi`, `oi` or `pi1`, `pi2`, ...), each in (0, 1) with a logit link;
#' with several points the parameters are mapped to weights through a
#' multinomial logit (stick breaking is not used), so
#' `pi_j = exp(eta_j) / (1 + sum(exp(eta)))`.
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
  if (is.null(names)) names <- if (k == 1 && at == 0) "zi" else if (k == 1 && at == 1) "oi" else paste0("pi", seq_len(k))
  names <- vapply(names, function(n) new_par_name(base, n), "")
  bnames <- names(base$pars)
  wts <- function(p) {
    if (k == 1) return(matrix(p[[names]], ncol = 1))
    W <- sapply(names, function(n) p[[n]])
    W <- matrix(W, ncol = k)
    s <- rowSums(W)
    if (any(s >= 1, na.rm = TRUE)) W[s >= 1, ] <- NA
    W
  }
  lpdf <- function(y, p) {
    W <- wts(p); bp <- sub_pars(p, bnames)
    rest <- 1 - rowSums(W)
    lb <- eval_lpdf(base, y, bp)
    out <- log(rest) + lb
    for (j in seq_len(k)) {
      hit <- !is.na(y) & y == at[j]
      if (any(hit)) {
        if (base$discrete) out[hit] <- log(W[hit, j] + rest[hit] * exp(lb[hit]))
        else out[hit] <- log(W[hit, j])
      }
    }
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
  if (!is.null(base$stan) && k == 1) {
    fn <- sprintf("rc_inflate_%s_%s", stan_ident(base$id), gsub("[^0-9a-z]", "", tolower(format(at))))
    args <- paste(sprintf("real %s", bnames), collapse = ", ")
    pm <- stats::setNames(lapply(bnames, identity), bnames)
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
  }", fn, yt, args, format(at), bl, paste(bnames, collapse = ", "))
    } else ""
    stan <- list(native = FALSE, needs = base$stan$needs,
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
  if (!is.null(base$bugs$loglik) && k == 1) {
    bl <- base$bugs$loglik
    bugs <- list(loglik = if (base$discrete)
      sprintf("log(equals({y}, %s) * {%s} + (1 - {%s}) * exp(%s))", format(at), names, names, bl)
      else sprintf("equals({y}, %s) * log({%s}) + (1 - equals({y}, %s)) * (log(1 - {%s}) + %s)",
                   format(at), names, format(at), names, gsub("\\{y\\}", "({y} + equals({y}, %s))", bl)))
    if (!base$discrete) bugs$loglik <- sprintf(bugs$loglik, format(at))
  }
  op_family(base, sprintf("inflate(%s, at = %s)", base$id, paste(format(at), collapse = ",")),
            sprintf("%s inflated at %s", base$display_name, paste(format(at), collapse = ", ")),
            newpars, lpdf, cdf = cdf, rng = rng, mean = mean, variance = variance,
            support = support, discrete = base$discrete, stan = stan, bugs = bugs,
            operator = list(type = "inflate", base = base, at = at, pars = names),
            formula = sprintf("P(Y = c) = pi + (1 - pi) P0(Y = c); elsewhere (1 - pi) f0(y)"))
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
  name <- new_par_name(base, name)
  bnames <- names(base$pars)
  lp_at <- function(bp) if (base$discrete) eval_lpdf(base, rep(at, max(lengths(bp))), bp) else rep(-Inf, max(lengths(bp)))
  lpdf <- function(y, p) {
    bp <- sub_pars(p, bnames); pi <- p[[name]]
    out <- log1p(-pi) + eval_lpdf(base, y, bp) - log1mexp(-lp_at(bp))
    hit <- !is.na(y) & y == at
    out[hit] <- log(pi[hit])
    out
  }
  cdf <- function(q, p) {
    bp <- sub_pars(p, bnames); pi <- p[[name]]
    F0 <- family_cdf(base, q, bp); pa <- exp(lp_at(bp))
    cond <- (F0 - pa * (q >= at)) / (1 - pa)
    pi * (q >= at) + (1 - pi) * cond
  }
  support <- function(p) {
    s <- base$support(sub_pars(p, bnames))
    list(lower = pmin(s$lower, at), upper = pmax(s$upper, at))
  }
  rng <- function(n, p) {
    bp <- sub_pars(p, bnames); pi <- p[[name]]
    y <- rep(at, n)
    todo <- which(stats::runif(n) >= pi)
    while (length(todo)) {
      d <- family_rng(base, length(todo), lapply(bp, `[`, todo))
      ok <- d != at
      y[todo[ok]] <- d[ok]
      todo <- todo[!ok]
    }
    y
  }
  mean <- function(p) {
    bp <- sub_pars(p, bnames); pi <- p[[name]]
    m <- family_moment(base, "mean", bp, quiet = TRUE); pa <- exp(lp_at(bp))
    pi * at + (1 - pi) * (m - pa * at) / (1 - pa)
  }
  variance <- function(p) {
    bp <- sub_pars(p, bnames); pi <- p[[name]]
    m <- family_moment(base, "mean", bp, quiet = TRUE)
    v <- family_moment(base, "variance", bp, quiet = TRUE); pa <- exp(lp_at(bp))
    m1 <- (m - pa * at) / (1 - pa); m2 <- (v + m^2 - pa * at^2) / (1 - pa)
    ey <- pi * at + (1 - pi) * m1; ey2 <- pi * at^2 + (1 - pi) * m2
    ey2 - ey^2
  }
  newpars <- c(base$pars, stats::setNames(list(P(name, "unit", role = "hurdle",
    interpretation = "Probability of the hurdle value")), name))
  stan <- NULL
  if (!is.null(base$stan)) {
    fn <- sprintf("rc_hurdle_%s", stan_ident(base$id))
    args <- paste(sprintf("real %s", bnames), collapse = ", ")
    pm <- stats::setNames(lapply(bnames, identity), bnames)
    yt <- if (base$discrete) "int" else "real"
    blp <- stan_fill(base$stan$lpdf, "y", pm)
    b0 <- if (base$discrete) sprintf("log1m_exp(%s)", stan_fill(base$stan$lpdf, format(at), pm)) else "0"
    stan <- list(native = FALSE, needs = base$stan$needs,
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
            support = support, stan = stan, bugs = bugs,
            operator = list(type = "hurdle", base = base, at = at, pars = name),
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
  if (base$variate_type != "scalar") stop_unsupported("Truncation of a vector valued law", base$id)
  if (lower >= upper) stop_input("lower must be smaller than upper.")
  bn <- names(base$pars)
  lo_q <- if (base$discrete) lower - 1 else lower
  lnorm <- function(p) {
    Fu <- if (is.finite(upper)) family_cdf(base, rep(upper, max(lengths(p))), p) else 1
    Fl <- if (is.finite(lo_q)) family_cdf(base, rep(lo_q, max(lengths(p))), p) else 0
    log(pmax(Fu - Fl, 0))
  }
  lpdf <- function(y, p) eval_lpdf(base, y, p) - lnorm(p)
  cdf <- function(q, p) {
    Fl <- if (is.finite(lo_q)) family_cdf(base, rep(lo_q, length(q)), p) else 0
    Fq <- family_cdf(base, pmin(pmax(q, lower), upper), p)
    out <- (Fq - Fl) / exp(lnorm(p))
    out[q < lower] <- 0; out[q >= upper] <- 1
    pmin(pmax(out, 0), 1)
  }
  quantile <- function(pr, p) {
    Fl <- if (is.finite(lo_q)) family_cdf(base, rep(lo_q, length(pr)), p) else 0
    Fu <- if (is.finite(upper)) family_cdf(base, rep(upper, length(pr)), p) else 1
    family_quantile(base, Fl + pr * (Fu - Fl), p)
  }
  support <- function(p) {
    s <- base$support(p)
    list(lower = pmax(s$lower, lower), upper = pmin(s$upper, upper))
  }
  stan <- NULL
  if (!is.null(base$stan) && !is.null(base$stan$lcdf)) {
    fn <- sprintf("rc_trunc_%s", stan_ident(base$id))
    pm <- stats::setNames(lapply(bn, identity), bn)
    args <- paste(sprintf("real %s", bn), collapse = ", ")
    yt <- if (base$discrete) "int" else "real"
    nrm <- if (is.finite(lower) && is.finite(upper))
      sprintf("log_diff_exp(%s, %s)", stan_fill(base$stan$lcdf, format(upper), pm), stan_fill(base$stan$lcdf, format(lo_q), pm))
      else if (is.finite(lower)) stan_fill(base$stan$lccdf, format(lo_q), pm)
      else stan_fill(base$stan$lcdf, format(upper), pm)
    stan <- list(native = FALSE, needs = base$stan$needs,
      functions = paste0(base$stan$functions, sprintf("
  real %s_%s(%s y, %s) {
    return %s - %s;
  }", fn, if (base$discrete) "lpmf" else "lpdf", yt, args, stan_fill(base$stan$lpdf, "y", pm), nrm)),
      lpdf = sprintf("%s_%s({y} | %s)", fn, if (base$discrete) "lpmf" else "lpdf",
                     paste(sprintf("{%s}", bn), collapse = ", ")))
  }
  op_family(base, sprintf("truncate(%s, %s, %s)", base$id, format(lower), format(upper)),
            sprintf("%s truncated to [%s, %s]", base$display_name, format(lower), format(upper)),
            base$pars, lpdf, cdf = cdf, quantile = quantile,
            rng = function(n, p) quantile(stats::runif(n), p),
            support = support, stan = stan,
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
    E <- cbind(0, sapply(2:K, function(k) rep_len(p[[paste0("w", k)]], n)))
    E <- matrix(E, nrow = n)
    E <- E - apply(E, 1, max)
    W <- exp(E); W / rowSums(W)
  }
  comp_pars <- function(p, k) {
    m <- cmap[[k]]
    stats::setNames(lapply(m, function(nm) p[[nm]]), names(m))
  }
  lpdf <- function(y, p) {
    W <- weights(p)
    L <- sapply(seq_len(K), function(k) log(W[, k]) + eval_lpdf(comps[[k]], y, comp_pars(p, k)))
    L <- matrix(L, ncol = K)
    m <- apply(L, 1, max)
    out <- m + log(rowSums(exp(L - m)))
    out[!is.finite(m)] <- m[!is.finite(m)]
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
  if (all(vapply(comps, function(f) !is.null(f$stan), TRUE))) {
    fn <- paste0("rc_mix_", paste(vapply(comps, function(f) stan_ident(f$id), ""), collapse = "_"))
    yt <- if (disc[1]) "int" else "real"
    args <- paste(c(sprintf("real %s", stan_ident(names(pars)))), collapse = ", ")
    terms <- vapply(seq_len(K), function(k) {
      pm <- stats::setNames(as.list(stan_ident(cmap[[k]])), names(cmap[[k]]))
      sprintf("log(w[%d]) + %s", k, stan_fill(comps[[k]]$stan$lpdf, "y", pm))
    }, "")
    funs <- unique(vapply(comps, function(f) f$stan$functions, ""))
    stan <- list(native = FALSE, needs = unique(unlist(lapply(comps, function(f) f$stan$needs))),
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
    variance = variance, support = support, discrete = disc[1], stan = stan,
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
    b <- sort(c(suppressWarnings(forward(s$lower)), suppressWarnings(forward(s$upper))))
    list(lower = b[1], upper = b[2])
  }
  lpdf <- function(y, p) eval_lpdf(base, inverse(y), p) + log_jacobian(y)
  cdf <- function(q, p) {
    v <- family_cdf(base, inverse(q), p)
    if (increasing) v else 1 - v
  }
  quantile <- function(pr, p) forward(family_quantile(base, if (increasing) pr else 1 - pr, p))
  rng <- function(n, p) forward(family_rng(base, n, p))
  op_family(base, sprintf("transform(%s, %s)", base$id, label),
            sprintf("%s transformed (%s)", base$display_name, label), base$pars,
            lpdf, cdf = cdf, quantile = quantile, rng = rng, support = sup,
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
  lpdf <- function(y, p) {
    log(pmax(family_cdf(base, y + 1 - offset, p) - family_cdf(base, y - offset, p), 0))
  }
  cdf <- function(q, p) family_cdf(base, floor(q) + 1 - offset, p)
  support <- function(p) {
    s <- base$support(p)
    list(lower = floor(s$lower + offset), upper = if (is.finite(s$upper)) floor(s$upper + offset) else Inf)
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
  h <- 10^-digits
  f <- family_def(
    sprintf("round(%s, %d)", base$id, as.integer(digits)), sprintf("Rounded %s", base$display_name),
    base$category, base$pars, discrete = TRUE,
    lpdf = function(y, p) log(pmax(family_cdf(base, y + h / 2, p) - family_cdf(base, y - h / 2, p), 0)),
    cdf = function(q, p) family_cdf(base, floor(q / h + 1e-9) * h + h / 2, p),
    quantile = function(pr, p) ceiling((family_quantile(base, pr, p) - h / 2) / h - 1e-9) * h,
    rng = function(n, p) round(family_rng(base, n, p) / h) * h,
    support = base$support, status = "implemented", numeric_moments = FALSE,
    numeric_cdf = FALSE,
    formula = "P(Y = y) = F(y + h/2) - F(y - h/2), y on the grid of multiples of h")
  f$extra$grid <- h
  f$operator <- list(type = "round", base = base, digits = digits)
  f
}
