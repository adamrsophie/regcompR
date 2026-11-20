# Distribution objects: a family plus parameter values, with the canonical
# probability interface.

#' Create a distribution object
#'
#' Combines a registered family (or one built by an operator such as
#' [inflate()]) with parameter values. Parameters are vectorized: a
#' distribution object with vectors of length `n` represents `n`
#' distributions, one per position, as in regression predictions.
#'
#' @param family Identifier, alias or `regcomp_family` object.
#' @param ... Named parameter values.
#' @param parameterization Optional parameterization name, see
#'   [get_family()].
#' @return An object of class `regcomp_dist`.
#' @seealso [log_prob()], [cdf()], [quantile.regcomp_dist()], [rng()],
#'   [capabilities()], [plot.regcomp_dist()].
#' @export
#' @examples
#' d <- distribution("negative_binomial", parameterization = "nb2",
#'                   mu = 4, phi = 2)
#' d
#' prob(d, 0:5)
#' cdf(d, 3)
#' quantile(d, c(0.1, 0.5, 0.9))
#' mean(d); variance(d)
#' distribution("conway_maxwell_poisson_mean", mu = 3, nu = 1.5)
distribution <- function(family, ..., parameterization = NULL) {
  fam <- get_family(family, parameterization)
  vals <- list(...)
  if (length(vals) && (is.null(names(vals)) || any(!nzchar(names(vals)))))
    stop_input("Distribution parameters must be named.")
  unknown <- setdiff(names(vals), names(fam$pars))
  if (length(unknown))
    stop_input(sprintf("Unknown parameter%s for '%s': %s. Parameters are: %s.",
                       if (length(unknown) > 1) "s" else "", fam$id,
                       paste(unknown, collapse = ", "),
                       paste(names(fam$pars), collapse = ", ")))
  for (p in fam$pars) {
    if (is.null(vals[[p$name]])) {
      if (is.null(p$default))
        stop_input(sprintf("Parameter '%s' of '%s' is missing.", p$name, fam$id))
      vals[[p$name]] <- p$default
    }
  }
  vals <- vals[names(fam$pars)]
  if (fam$variate_type == "scalar") vals <- recycle_pars(vals)
  validate_parameters(fam, vals)
  structure(list(family = fam, pars = vals), class = "regcomp_dist")
}

#' Validate parameters and responses
#'
#' `validate_parameters()` signals a `regcomp_domain_error` when a value lies
#' outside its domain or violates a family constraint. `validate_response()`
#' returns a logical vector that is `TRUE` where a response value lies in the
#' support.
#'
#' @param x A family, a distribution or an identifier.
#' @param pars Named list of parameter values.
#' @param y Response values.
#' @return `validate_parameters()` returns `pars` invisibly;
#'   `validate_response()` returns a logical vector.
#' @export
#' @examples
#' validate_parameters("poisson", list(lambda = 2))
#' validate_response(distribution("poisson", lambda = 2), c(-1, 0, 1.5, 3))
validate_parameters <- function(x, pars) {
  fam <- if (inherits(x, "regcomp_dist")) x$family else get_family(x)
  for (p in fam$pars) {
    v <- pars[[p$name]]
    if (is.null(v)) next
    ok <- check_par_domain(p, v)
    if (!all(ok | is.na(v)))
      stop_domain(sprintf("Parameter '%s' of '%s' must be in its %s domain; got %s.",
                          p$name, fam$id, p$domain,
                          format(v[!ok & !is.na(v)][1], digits = 6)))
  }
  if (!is.null(fam$check)) {
    msg <- fam$check(pars)
    if (is.character(msg) && length(msg) && nzchar(msg[1])) stop_domain(msg[1])
  }
  if (fam$variate_type == "scalar") {
    s <- tryCatch(fam$support(pars), error = function(e) NULL)
    # Discrete laws and operators (censoring) may legitimately collapse to a point.
    if (!is.null(s) && any(s$lower > s$upper | (s$lower == s$upper & !fam$discrete & is.null(fam$operator)), na.rm = TRUE))
      stop_domain(sprintf("The parameters of '%s' give an empty support (lower bound not below the upper bound).", fam$id))
  }
  invisible(pars)
}

#' @rdname validate_parameters
#' @export
validate_response <- function(x, y) {
  d <- as_dist(x)
  in_support(d$family, y, d$pars)
}

as_dist <- function(x) {
  if (inherits(x, "regcomp_dist")) return(x)
  stop_input("Expected a distribution object created by distribution().")
}

# Shape checks for vector and matrix responses: integer counts, simplex
# rows, unit vectors, symmetric positive definite matrices.
mv_support <- function(fam, y) {
  vt <- fam$variate_type
  if (vt == "matrix" && fam$id %in% c("wishart", "inverse_wishart", "lkj", "matrix_beta")) {
    mats <- if (is.list(y)) y else if (length(dim(y)) == 3) lapply(seq_len(dim(y)[3]), function(i) y[, , i]) else list(y)
    pd <- function(m) all(eigen(m, symmetric = TRUE, only.values = TRUE)$values > 0)
    return(vapply(mats, function(m) { m <- unname(as.matrix(m))
      all(is.finite(m)) && isSymmetric(m, tol = 1e-8) && pd(m) &&
        (fam$id != "lkj" || all(abs(diag(m) - 1) < 1e-8)) &&
        (fam$id != "matrix_beta" || pd(diag(nrow(m)) - m)) }, TRUE))
  }
  if (!vt %in% c("vector", "unit_vector")) return(TRUE)
  Y <- as_rows(y); ok <- rep(TRUE, nrow(Y))
  if (fam$discrete) ok <- ok & !apply(!is.na(Y) & (Y < 0 | abs(Y - round(Y)) > 1e-9), 1, any)
  if (vt == "unit_vector") ok <- ok & abs(sqrt(rowSums(Y^2)) - 1) < 1e-6
  if (fam$id %in% c("dirichlet", "continuous_categorical")) ok <- ok & rowSums(Y <= 0) == 0 & abs(rowSums(Y) - 1) < 1e-8
  if (grepl("binary", fam$id)) ok <- ok & rowSums(Y != 0 & Y != 1) == 0
  ok
}

in_support <- function(fam, y, pars) {
  if (fam$variate_type != "scalar") {
    ok <- mv_support(fam, y)
    if (!is.null(fam$response_check)) ok <- ok & fam$response_check(y, pars)
    return(ok)
  }
  r <- recycle_pars(c(list(.y = y), pars))
  y <- r$.y
  s <- fam$support(r[-1])
  lo <- rep_len(s$lower, length(y)); hi <- rep_len(s$upper, length(y))
  ok <- !is.na(y)
  ok[ok] <- y[ok] >= lo[ok] & y[ok] <= hi[ok]
  if (isTRUE(s$open_lower)) ok[ok] <- y[ok] > lo[ok]
  if (isTRUE(s$open_upper)) ok[ok] <- y[ok] < hi[ok]
  if (fam$discrete) {
    g <- if (is.null(fam$extra$grid)) 1 else fam$extra$grid
    ok[ok] <- is.finite(y[ok])
    if (g > 0) ok[ok] <- abs(y[ok] / g - round(y[ok] / g)) < 1e-9
  }
  if (!is.null(fam$response_check)) ok[ok] <- fam$response_check(y, r[-1])[ok]
  ok
}

# Core evaluation shared by every method: recycle, mask out of support, call
# the family function on valid positions.
eval_lpdf <- function(fam, y, pars) {
  if (fam$variate_type != "scalar") {
    # A plain vector is one observation of a vector valued law.
    if (fam$variate_type %in% c("unit_vector", "vector", "vector_mixed") && is.null(dim(y))) y <- matrix(y, nrow = 1)
    out <- fam$lpdf(y, pars)
    ok <- in_support(fam, y, pars)
    if (length(ok) == 1L || length(ok) == length(out)) out[!rep_len(ok, length(out))] <- -Inf
    return(out)
  }
  r <- recycle_pars(c(list(.y = y), pars))
  y <- r$.y; p <- r[-1]
  out <- rep(-Inf, length(y))
  out[is.na(y)] <- NA
  ok <- in_support(fam, y, p)
  pok <- rep(TRUE, length(y))
  for (sp in fam$pars) if (!sp$data || TRUE) pok <- pok & check_par_domain(sp, p[[sp$name]])
  out[!pok] <- NaN
  sel <- ok & pok
  if (any(sel)) out[sel] <- fam$lpdf(y[sel], lapply(p, `[`, sel))
  # A continuous density that is NaN at a finite support bound (0 log 0 and
  # similar forms) takes its limit from inside the support.
  nb <- which(sel & is.nan(out))
  if (!fam$discrete && length(nb)) {
    pb <- lapply(p, `[`, nb)
    s <- tryCatch(fam$support(pb), error = function(e) NULL)
    if (!is.null(s)) {
      lo <- rep_len(s$lower, length(nb)); hi <- rep_len(s$upper, length(nb)); yb <- y[nb]
      eps <- 1e-10 * pmax(1, abs(yb))
      inner <- ifelse(is.finite(lo) & yb == lo, yb + eps, ifelse(is.finite(hi) & yb == hi, yb - eps, NA))
      v <- suppressWarnings(fam$lpdf(ifelse(is.na(inner), yb, inner), pb))
      # A limit of zero density shows as a large negative value.
      v[!is.na(inner) & v < -1e6] <- -Inf
      out[nb[!is.na(inner)]] <- v[!is.na(inner)]
    }
  }
  out
}

#' Probability interface
#'
#' Evaluate a distribution object. `log_prob()` returns the log density (or
#' log mass for discrete laws), `prob()` the density or mass, `cdf()` the
#' distribution function, `survival()` the survival function, and `rng()`
#' random draws. Values outside the support give log density `-Inf`.
#' Methods without a closed form use a declared numerical fallback
#' (quadrature, summation or inversion with tolerance checks); when none is
#' justified a `regcomp_unsupported` error is raised.
#'
#' @param x A distribution object.
#' @param y,q Response values or quantiles.
#' @param lower_tail Logical; `FALSE` gives upper tail probabilities.
#' @param log_p Logical; return log probabilities.
#' @param n Number of draws. With vectorized parameters of length `m > 1`,
#'   `n` must equal `m` or be omitted.
#' @param ... Unused.
#' @return Numeric vector.
#' @name probability-interface
#' @examples
#' d <- distribution("gamma", mu = 2, sigma = 0.5)
#' log_prob(d, 1:3)
#' prob(d, 1)
#' cdf(d, 2)
#' survival(d, 2, log_p = TRUE)
#' set.seed(1); rng(d, 5)
NULL

#' @rdname probability-interface
#' @export
log_prob <- function(x, y, ...) UseMethod("log_prob")

#' @export
log_prob.regcomp_dist <- function(x, y, ...) eval_lpdf(x$family, y, x$pars)

#' @rdname probability-interface
#' @export
prob <- function(x, y, ...) UseMethod("prob")

#' @export
prob.regcomp_dist <- function(x, y, ...) exp(log_prob(x, y))

#' @rdname probability-interface
#' @export
cdf <- function(x, q, ...) UseMethod("cdf")

#' @rdname probability-interface
#' @export
cdf.regcomp_dist <- function(x, q, lower_tail = TRUE, log_p = FALSE, ...) {
  fam <- x$family
  if (fam$variate_type == "scalar") { r <- recycle_pars(c(list(.q = q), x$pars)); q <- r$.q; p <- r[-1] }
  else if (scalar_vecpar(fam)) p <- x$pars
  else stop_unsupported("A univariate CDF", fam$id)
  if (!log_p && lower_tail) return(family_cdf(fam, q, p))
  out <- if (lower_tail) family_lcdf(fam, q, p) else family_lccdf(fam, q, p)
  if (log_p) out else exp(out)
}

# Category probabilities of an ordinal or nominal law (n x K matrix).
cat_probs <- function(fam, p, n) {
  K <- fam$support(p)$upper
  matrix(vapply(seq_len(K), function(k) exp(eval_lpdf(fam, rep(k, n), p)), numeric(n)), n, K)
}

# Log distribution and survival functions, exact in the tails when the
# family supplies them, and the log probability of an interval
# P(a < Y <= b) chosen from the tail with usable precision.
family_lcdf <- function(fam, q, p) {
  if (!is.null(fam$lcdf)) return(fam$lcdf(q, p))
  F <- family_cdf(fam, q, p)
  tail_refine(fam, q, p, log(F), F < 1e-12, upper = FALSE)
}
family_lccdf <- function(fam, q, p) {
  if (!is.null(fam$lccdf)) return(fam$lccdf(q, p))
  F <- family_cdf(fam, q, p)
  tail_refine(fam, q, p, log1p(-F), F > 1 - 1e-12, upper = TRUE)
}

# Where a tail probability from the CDF has lost its precision, recompute it
# in log space from the density: an integral scaled by the density at q
# (continuous laws) or a log sum of masses (discrete laws).
tail_refine <- function(fam, q, p, out, bad, upper) {
  bad <- which(bad & is.finite(q))
  if (!length(bad) || fam$variate_type != "scalar" || scalar_vecpar(fam)) return(out)
  n <- length(q); pr <- lapply(p, function(v) if (length(v) == n && n > 1) v else rep_len(v, n))
  for (i in bad) out[i] <- tryCatch({
    pi_ <- lapply(pr, `[`, i)
    lf <- function(t) eval_lpdf(fam, t, lapply(pi_, rep_len, length(t)))
    if (fam$discrete) {
      # Masses beyond q (upper) or up to q (lower), in blocks until negligible.
      acc <- -Inf; start <- if (upper) floor(q[i]) + 1 else floor(q[i])
      for (b in 0:199) {
        t <- if (upper) start + b * 1000 + 0:999 else start - b * 1000 - 0:999
        l <- lf(t); m <- suppressWarnings(max(l))
        acc <- log_sum_exp(c(acc, l))
        if (!is.finite(m) && b > 0 || is.finite(m) && m < acc - 40) break
      }
      acc
    } else {
      l0 <- lf(q[i]); if (!is.finite(l0)) stop("no density")
      I <- stats::integrate(function(t) exp(lf(t) - l0), if (upper) q[i] else -Inf, if (upper) Inf else q[i],
                            rel.tol = 1e-8)$value
      if (I > 0) l0 + log(I) else out[i]
    }
  }, error = function(e) out[i])
  out
}
# Measure contract. A discrete law has mass only on its support points: a
# lattice of multiples of `lattice_step()` (0 means an irregular finite set
# described by the density). A continuous law has a density with respect to
# Lebesgue measure plus, when `fam$atom` is set, point masses given by
# `fam$atom(y, p)`; at an atom `lpdf` returns the log mass. `family_atom()`
# gives the point mass at y and `family_cdf_left()` the left limit F(y-).
lattice_step <- function(fam) if (is.null(fam$extra$grid)) 1 else fam$extra$grid
has_atoms <- function(f) !f$discrete && !is.null(f$atom)
family_atom <- function(fam, y, p) {
  if (fam$discrete) { v <- exp(eval_lpdf(fam, y, p)); v[is.na(v)] <- 0; return(v) }
  if (is.null(fam$atom)) return(rep(0, length(y)))
  r <- recycle_pars(c(list(.y = y), p))
  v <- rep(0, length(r$.y)); ok <- !is.na(r$.y)
  if (any(ok)) v[ok] <- fam$atom(r$.y[ok], lapply(r[-1], `[`, ok))
  v
}
family_cdf_left <- function(fam, q, p) {
  F <- family_cdf(fam, q, p)
  if (!fam$discrete && is.null(fam$atom)) return(F)
  pmax(F - family_atom(fam, q, p), 0)
}

family_log_interval <- function(fam, a, b, p) {
  n <- length(a); out <- rep(-Inf, n)
  pos <- !is.na(a) & !is.na(b) & b > a
  out[is.na(a) | is.na(b)] <- NA
  if (!any(pos)) return(out)
  i <- which(pos); pp <- lapply(p, function(v) if (length(v) == n) v[i] else v)
  ai <- a[i]; bi <- b[i]
  lFb <- ifelse(is.infinite(bi) & bi > 0, 0, family_lcdf(fam, ifelse(is.finite(bi), bi, 0), pp))
  lFa <- ifelse(is.infinite(ai) & ai < 0, -Inf, family_lcdf(fam, ifelse(is.finite(ai), ai, 0), pp))
  lSa <- ifelse(is.infinite(ai) & ai < 0, 0, family_lccdf(fam, ifelse(is.finite(ai), ai, 0), pp))
  lSb <- ifelse(is.infinite(bi) & bi > 0, -Inf, family_lccdf(fam, ifelse(is.finite(bi), bi, 0), pp))
  upper <- lFa > log(0.5)
  v <- ifelse(upper, log_diff_exp(lSa, lSb), log_diff_exp(lFb, lFa))
  v[is.nan(v)] <- -Inf
  out[i] <- v
  out
}

family_cdf <- function(fam, q, p) {
  if (fam$variate_type == "categorical") {
    P <- cat_probs(fam, p, length(q))
    k <- pmin(pmax(floor(q), 0), ncol(P))
    return(rowSums(P * (col(P) <= k)))
  }
  vec <- scalar_vecpar(fam)
  n <- length(q)
  # Exterior points get 0 or 1 from the support, so interior formulas are
  # never evaluated outside their domain.
  s <- tryCatch(fam$support(p), error = function(e) NULL)
  lo <- if (is.null(s) || vec && length(s$lower) != 1) rep(-Inf, n) else rep_len(s$lower, n)
  hi <- if (is.null(s) || vec && length(s$upper) != 1) rep(Inf, n) else rep_len(s$upper, n)
  out <- rep(NA_real_, n)
  ok <- !is.na(q) & !is.na(lo) & !is.na(hi)
  below <- ok & (q < lo | q == -Inf); above <- ok & !below & is.finite(hi) & q >= hi
  out[below] <- 0; out[above] <- 1
  inn <- !is.na(q) & !below & !above
  if (!any(inn)) return(out)
  pp <- if (vec) p else lapply(p, function(v) if (length(v) == n && n > 1) v[inn] else v)
  qq <- q[inn]
  out[inn] <- if (!is.null(fam$cdf)) fam$cdf(qq, pp) else {
    if (!isTRUE(fam$numeric_cdf)) stop_unsupported("The CDF", fam$id)
    l <- lo[inn]
    if (fam$discrete) {
      if (any(!is.finite(l))) stop_unsupported("A numerical CDF on an unbounded integer support", fam$id)
      numeric_cdf_disc(function(y, pv) eval_lpdf(fam, y, pv), qq, pp, l, whole = vec)
    } else numeric_cdf_cont(function(y, pv) eval_lpdf(fam, y, pv), qq, pp, l, whole = vec)
  }
  out[inn & q == Inf & is.nan(out)] <- 1
  out
}
#' @rdname probability-interface
#' @export
survival <- function(x, q, ...) UseMethod("survival")

#' @rdname probability-interface
#' @export
survival.regcomp_dist <- function(x, q, log_p = FALSE, ...)
  cdf(x, q, lower_tail = FALSE, log_p = log_p)

#' Quantiles of a distribution object
#'
#' For discrete laws returns the smallest support point `y` with
#' `F(y) >= p`. Without a closed form the CDF is inverted numerically with
#' bracketing; failures raise `regcomp_numerical_failure`.
#'
#' @param x A distribution object.
#' @param probs Probabilities.
#' @param lower_tail,log_p As in [cdf()].
#' @param ... Unused.
#' @return Numeric vector.
#' @export
#' @examples
#' quantile(distribution("weibull", shape = 2, scale = 1), c(0.25, 0.5))
quantile.regcomp_dist <- function(x, probs = c(0.025, 0.5, 0.975),
                                  lower_tail = TRUE, log_p = FALSE, ...) {
  fam <- x$family
  if (fam$variate_type != "scalar" && !scalar_vecpar(fam)) stop_unsupported("A univariate quantile function", fam$id)
  lp <- if (log_p) probs else log(probs)
  pr <- exp(lp)
  if (!lower_tail) pr <- -expm1(lp)
  if (scalar_vecpar(fam)) return(family_quantile(fam, pr, x$pars))
  r <- recycle_pars(c(list(.p = pr, .lp = lp), x$pars))
  out <- family_quantile(fam, r$.p, r[-(1:2)])
  # Tail probabilities that do not survive the conversion to an ordinary
  # probability are inverted on the log scale.
  tiny <- which(is.finite(r$.lp) & r$.lp < log(1e-12))
  if (!fam$discrete) for (i in tiny) {
    pi <- lapply(r[-(1:2)], `[`, i)
    g <- function(q) (if (lower_tail) family_lcdf(fam, q, pi) else family_lccdf(fam, q, pi)) - r$.lp[i]
    s0 <- fam$support(pi); start <- family_quantile(fam, if (lower_tail) 1e-10 else 1 - 1e-10, pi)
    step <- 1; a <- start; b <- start
    if (lower_tail) { while (is.finite(g(a)) && g(a) > 0 && a > s0$lower) { a <- max(s0$lower, a - step); step <- 2 * step } }
    else { while (is.finite(g(b)) && g(b) > 0 && b < s0$upper) { b <- min(s0$upper, b + step); step <- 2 * step } }
    out[i] <- tryCatch(stats::uniroot(g, sort(c(a, b)), tol = 1e-12)$root, error = function(e) out[i])
  }
  out
}

family_quantile <- function(fam, p, pars, grid = TRUE) {
  if (!is.null(fam$quantile)) return(fam$quantile(p, pars))
  if (is.null(fam$cdf) && !isTRUE(fam$numeric_cdf)) stop_unsupported("The quantile function", fam$id)
  s <- fam$support(pars)
  n <- length(p)
  # Many probabilities with one parameter set (random draws): invert on a
  # grid once, then refine every value at once inside its grid cell.
  if (grid && !fam$discrete && !scalar_vecpar(fam) && n > 40 && length(pars) &&
      all(vapply(pars, function(v) length(unique(v)) == 1, TRUE)))
    return(quantile_grid(fam, p, lapply(pars, `[`, 1)))
  start <- tryCatch({
    m <- family_moment(fam, "mean", pars, quiet = TRUE)
    ifelse(is.finite(m), m, 0)
  }, error = function(e) rep(0, n))
  start <- rep_len(start, n)
  lo <- rep_len(s$lower, n); hi <- rep_len(s$upper, n)
  start <- pmin(pmax(start, ifelse(is.finite(lo), lo + 1e-8, start)),
                ifelse(is.finite(hi), hi - 1e-8, start))
  if (scalar_vecpar(fam))
    return(invert_cdf(function(q, pp) family_cdf(fam, q, pars), p, list(), fam$discrete, lo, hi, start))
  invert_cdf(function(q, pp) family_cdf(fam, q, pp), p, pars, fam$discrete,
             lo, hi, start)
}

#' @rdname probability-interface
#' @export
rng <- function(x, n, ...) UseMethod("rng")

#' @rdname probability-interface
#' @export
rng.regcomp_dist <- function(x, n = NULL, ...) {
  if (...length()) stop_input("Unused arguments to rng(); give sample sizes and similar values as distribution parameters.")
  fam <- x$family
  m <- if (fam$variate_type == "scalar") max(vapply(x$pars, length, 1L)) else 1L
  if (is.null(n)) n <- m
  if (m > 1 && n != m)
    stop_input("With vectorized parameters, n must equal the parameter length.")
  p <- if (fam$variate_type == "scalar") lapply(x$pars, rep_len, n) else x$pars
  family_rng(fam, n, p)
}

# Scalar variates whose parameters include a vector or matrix (metalog
# coefficients, phase-type generators): never recycle the parameters.
scalar_vecpar <- function(fam) grepl("^scalar_.*_parameter$", fam$variate_type)

family_rng <- function(fam, n, pars) {
  if (!is.null(fam$rng)) return(fam$rng(n, pars))
  if (fam$variate_type == "categorical") {
    P <- cat_probs(fam, pars, n)
    return(apply(P, 1, function(pr) sample.int(length(pr), 1L, prob = pr)))
  }
  if (fam$variate_type != "scalar" && !scalar_vecpar(fam)) stop_unsupported("Random generation", fam$id)
  family_quantile(fam, stats::runif(n), pars)
}

#' Moments of a distribution object
#'
#' `mean()`, `variance()`, `skewness()` and `kurtosis()` (excess kurtosis)
#' return `NA` with a `regcomp_undefined_moment` warning where the moment does
#' not exist for the parameter values. Without closed forms, moments are
#' computed by quadrature or summation when the family allows it.
#'
#' @param x A distribution object.
#' @param ... Unused.
#' @return Numeric vector.
#' @name moments
#' @examples
#' d <- distribution("student_t", mu = 0, sigma = 1, nu = c(1.5, 3, 10))
#' suppressWarnings(variance(d))
NULL

#' @rdname moments
#' @export
mean.regcomp_dist <- function(x, ...) family_moment(x$family, "mean", x$pars)

#' @rdname moments
#' @export
variance <- function(x, ...) UseMethod("variance")

#' @rdname moments
#' @export
variance.regcomp_dist <- function(x, ...) family_moment(x$family, "variance", x$pars)

#' @rdname moments
#' @export
skewness <- function(x, ...) UseMethod("skewness")

#' @rdname moments
#' @export
skewness.regcomp_dist <- function(x, ...) family_moment(x$family, "skewness", x$pars)

#' @rdname moments
#' @export
kurtosis <- function(x, ...) UseMethod("kurtosis")

#' @rdname moments
#' @export
kurtosis.regcomp_dist <- function(x, ...) family_moment(x$family, "kurtosis", x$pars)

family_moment <- function(fam, what, pars, quiet = FALSE) {
  f <- fam[[what]]
  if (fam$variate_type == "scalar") pars <- recycle_pars(pars)
  if (is.null(f)) {
    if (!isTRUE(fam$numeric_moments) || fam$variate_type != "scalar")
      stop_unsupported(sprintf("The %s", what), fam$id)
    out <- numeric_family_moment(fam, what, pars)
  } else out <- f(pars)
  bad <- is.na(out) & !Reduce(`|`, lapply(pars, function(v) is.na(v)), FALSE)
  if (!quiet && any(bad)) warn_undefined_moment(what, sum(bad))
  out
}

numeric_family_moment <- function(fam, what, pars) {
  lp <- function(y, pp) eval_lpdf(fam, y, pp)
  s <- fam$support(pars)
  n <- max(vapply(pars, length, 1L))
  lo <- rep_len(s$lower, n); hi <- rep_len(s$upper, n)
  mom <- function(fun) {
    if (fam$discrete) {
      if (any(!is.finite(lo))) return(rep(NA_real_, n))
      vapply(seq_len(n), function(i)
        numeric_moment_disc(lp, lapply(pars, `[`, i), fun, lo[i]), 0)
    } else numeric_moment_cont(lp, pars, fun, lo, hi)
  }
  m1 <- mom(identity)
  if (what == "mean") return(m1)
  c2 <- vapply(seq_len(n), function(i) {
    pp <- lapply(pars, `[`, i)
    if (fam$discrete) numeric_moment_disc(lp, pp, function(y) (y - m1[i])^2, lo[i])
    else numeric_moment_cont(lp, pp, function(y) (y - m1[i])^2, lo[i], hi[i])
  }, 0)
  if (what == "variance") return(c2)
  k <- if (what == "skewness") 3 else 4
  ck <- vapply(seq_len(n), function(i) {
    pp <- lapply(pars, `[`, i)
    if (fam$discrete) numeric_moment_disc(lp, pp, function(y) (y - m1[i])^k, lo[i])
    else numeric_moment_cont(lp, pp, function(y) (y - m1[i])^k, lo[i], hi[i])
  }, 0)
  if (what == "skewness") ck / c2^1.5 else ck / c2^2 - 3
}

#' Support of a distribution object
#'
#' @param x A distribution object.
#' @param ... Unused.
#' @return A list with `lower` and `upper` bounds (vectors when the support
#'   depends on parameters) and `discrete`.
#' @export
#' @examples
#' support(distribution("generalized_pareto", mu = 0, sigma = 1, xi = -0.5))
support <- function(x, ...) UseMethod("support")

#' @export
support.regcomp_dist <- function(x, ...) {
  s <- x$family$support(x$pars)
  c(s, list(discrete = x$family$discrete))
}

#' Score and Hessian of the log likelihood
#'
#' Derivatives of `log_prob(x, y)` with respect to the natural parameters,
#' per observation. Analytic expressions are used when a family provides
#' them; otherwise central finite differences with a step scaled to each
#' parameter.
#'
#' @param x A distribution object.
#' @param y Response values.
#' @param ... Unused.
#' @return `score()`: a matrix with one column per parameter. `hessian()`: an
#'   array `n x k x k`.
#' @export
#' @examples
#' d <- distribution("normal", mu = 0, sigma = 1)
#' score(d, c(-1, 0, 2))
score <- function(x, y, ...) UseMethod("score")

#' @export
score.regcomp_dist <- function(x, y, ...) {
  fam <- x$family
  free <- regression_pars(fam)
  r <- recycle_pars(c(list(.y = y), x$pars))
  y <- r$.y; p <- r[-1]
  out <- matrix(NA_real_, length(y), length(free), dimnames = list(NULL, free))
  for (k in free) {
    h <- 1e-5 * ifelse(p[[k]] == 0, 1e-3, abs(p[[k]]))
    up <- p; dn <- p
    up[[k]] <- p[[k]] + h; dn[[k]] <- p[[k]] - h
    out[, k] <- (eval_lpdf(fam, y, up) - eval_lpdf(fam, y, dn)) / (2 * h)
  }
  out
}

#' @rdname score
#' @export
hessian <- function(x, y, ...) UseMethod("hessian")

#' @export
hessian.regcomp_dist <- function(x, y, ...) {
  fam <- x$family
  free <- regression_pars(fam)
  r <- recycle_pars(c(list(.y = y), x$pars))
  y <- r$.y; p <- r[-1]
  K <- length(free)
  out <- array(NA_real_, c(length(y), K, K), list(NULL, free, free))
  f <- function(pp) eval_lpdf(fam, y, pp)
  h <- lapply(p[free], function(v) 1e-4 * ifelse(v == 0, 1e-2, abs(v)))
  for (a in seq_len(K)) for (b in a:K) {
    ka <- free[a]; kb <- free[b]
    shift <- function(sa, sb) {
      pp <- p
      pp[[ka]] <- pp[[ka]] + sa * h[[ka]]
      pp[[kb]] <- pp[[kb]] + sb * h[[kb]]
      f(pp)
    }
    v <- if (a == b) {
      (shift(1, 0) - 2 * f(p) + shift(-1, 0)) / h[[ka]]^2
    } else {
      (shift(1, 1) - shift(1, -1) - shift(-1, 1) + shift(-1, -1)) /
        (4 * h[[ka]] * h[[kb]])
    }
    out[, a, b] <- v; out[, b, a] <- v
  }
  out
}

#' @export
print.regcomp_dist <- function(x, ...) {
  fam <- x$family
  n <- if (fam$variate_type == "scalar") max(vapply(x$pars, length, 1L)) else 1L
  cat(sprintf("<regcomp_dist> %s%s\n", fam$display_name,
              if (n > 1) sprintf(" [%d distributions]", n) else ""))
  for (nm in names(x$pars)) {
    v <- x$pars[[nm]]
    s <- if (is.null(dim(v))) paste(format(utils::head(v, 6), digits = 4), collapse = ", ")
         else sprintf("%s array", paste(dim(v), collapse = " x "))
    if (is.null(dim(v)) && length(v) > 6) s <- paste0(s, ", ...")
    cat(sprintf("  %s = %s\n", nm, s))
  }
  invisible(x)
}

#' @export
length.regcomp_dist <- function(x) {
  switch(x$family$variate_type, scalar = max(vapply(x$pars, length, 1L)),
         categorical = length(x$pars[[1]]), 1L)
}

#' @export
format.regcomp_dist <- function(x, ...) sprintf("%s(%s)", x$family$id,
  paste(names(x$pars), "=", vapply(x$pars, function(v)
    paste(format(utils::head(as.numeric(v), 3), digits = 3), collapse = ","), ""),
    collapse = ", "))

quantile_grid <- function(fam, p, p1) {
  g <- c(1e-8, 1e-5, 1e-3, seq(0.01, 0.99, by = 0.02), 1 - 1e-3, 1 - 1e-5, 1 - 1e-8)
  qg <- family_quantile(fam, g, p1, grid = FALSE)
  out <- rep(NA_real_, length(p))
  F <- function(q) family_cdf(fam, q, lapply(p1, rep, length(q)))
  inn <- which(!is.na(p) & p >= g[1] & p <= g[length(g)])
  if (length(inn)) {
    j <- findInterval(p[inn], g, rightmost.closed = TRUE)
    a <- qg[j]; b <- qg[j + 1]; fa <- g[j] - p[inn]; fb <- g[j + 1] - p[inn]
    # Illinois (modified regula falsi), vectorized over all values.
    x <- a; side <- 0L
    for (it in 1:60) {
      x <- ifelse(fb != fa, (a * fb - b * fa) / (fb - fa), (a + b) / 2)
      fx <- F(x) - p[inn]
      if (max(abs(fx)) < 1e-13 || max(abs(b - a) / pmax(1, abs(x))) < 1e-14) break
      lft <- fx * fa > 0
      a[lft] <- x[lft]; fa[lft] <- fx[lft]; fb[lft & side == 1L] <- fb[lft & side == 1L] / 2
      b[!lft] <- x[!lft]; fb[!lft] <- fx[!lft]; fa[!lft & side == -1L] <- fa[!lft & side == -1L] / 2
      side <- ifelse(lft, 1L, -1L)
    }
    out[inn] <- x
  }
  rest <- setdiff(seq_along(p), inn)
  if (length(rest)) {
    pr <- p[rest]; ok <- !is.na(pr)
    if (any(ok)) out[rest[ok]] <- vapply(pr[ok], function(u) family_quantile(fam, u, p1), 0)
  }
  out
}
