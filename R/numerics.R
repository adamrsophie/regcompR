# Numerical utilities shared by every family. All functions are vectorized.

#' Stable numerical building blocks
#'
#' Log scale helpers used throughout the package. `log1mexp(x)` computes
#' `log(1 - exp(-x))` for `x > 0` (Maechler 2012); `log1pexp(x)` computes
#' `log(1 + exp(x))`; `log_sum_exp()` and `log_diff_exp()` combine log scale
#' values without overflow.
#'
#' @param x,a,b Numeric vectors.
#' @return Numeric vector.
#' @name numerics
#' @examples
#' log1mexp(1e-20)
#' log1pexp(800)
#' log_sum_exp(c(1000, 1000))
#' log_diff_exp(log(3), log(2))
NULL

#' @rdname numerics
#' @export
log1mexp <- function(x) {
  out <- x
  small <- !is.na(x) & x <= log(2)
  out[small] <- log(-expm1(-x[small]))
  out[!small] <- log1p(-exp(-x[!small]))
  out
}

#' @rdname numerics
#' @export
log1pexp <- function(x) {
  out <- x
  a <- !is.na(x) & x <= -37
  b <- !is.na(x) & x > -37 & x <= 18
  c <- !is.na(x) & x > 18 & x <= 33.3
  out[a] <- exp(x[a])
  out[b] <- log1p(exp(x[b]))
  out[c] <- x[c] + exp(-x[c])
  out
}

#' @rdname numerics
#' @export
log_sum_exp <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(-Inf)
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}

# Elementwise log(exp(a) + exp(b)).
log_add_exp <- function(a, b) {
  m <- pmax(a, b)
  out <- m + log1p(exp(-abs(a - b)))
  out[is.infinite(m) & m < 0] <- -Inf
  out[is.infinite(m) & m > 0] <- Inf
  out
}

#' @rdname numerics
#' @export
log_diff_exp <- function(a, b) {
  out <- a + log1mexp(a - b)
  out[!is.na(a) & !is.na(b) & a == b] <- -Inf
  out[!is.na(a) & !is.na(b) & a < b] <- NaN
  out[is.infinite(b) & b < 0 & !is.na(a)] <- a[is.infinite(b) & b < 0 & !is.na(a)]
  out
}

# x * log(y) with 0 * log(0) = 0.
xlogy <- function(x, y) {
  out <- x * log(y)
  out[!is.na(x) & x == 0 & !is.na(y)] <- 0
  out
}

lchoose_safe <- function(n, k) lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)

# Recycle a list of parameter vectors (and optionally y) to a common length.
recycle_pars <- function(pars, n = NULL) {
  lens <- vapply(pars, length, 1L)
  if (any(lens == 0L)) return(lapply(pars, function(p) p[0]))
  L <- max(c(lens, n))
  lapply(pars, function(p) if (length(p) == L) p else rep_len(p, L))
}

# Log of the modified Bessel function of the second kind, real order.
log_besselK <- function(x, nu) {
  r <- recycle_pars(list(x = x, nu = nu))
  out <- suppressWarnings(log(besselK(r$x, r$nu, expon.scaled = TRUE)) - r$x)
  bad <- which(!is.finite(out) & is.finite(r$x) & r$x > 0)
  for (i in bad) {
    # Upward recurrence K_{v+1} = K_{v-1} + (2v/x) K_v in log space from the
    # fractional order, stable for K.
    v <- abs(r$nu[i]); x <- r$x[i]; v0 <- v - floor(v)
    lk0 <- log(besselK(x, v0, expon.scaled = TRUE)) - x
    lk1 <- log(besselK(x, v0 + 1, expon.scaled = TRUE)) - x
    if (floor(v) == 0) { out[i] <- lk0; next }
    cur <- v0 + 1
    while (cur < v - 1e-9) {
      lk2 <- log_add_exp(lk0, log(2 * cur / x) + lk1)
      lk0 <- lk1; lk1 <- lk2; cur <- cur + 1
    }
    out[i] <- lk1
  }
  out
}

# Log of the modified Bessel function of the first kind, real order.
log_besselI <- function(x, nu) {
  r <- recycle_pars(list(x = x, nu = nu))
  log(besselI(r$x, r$nu, expon.scaled = TRUE)) + r$x
}

# Sum a series of log terms term(k), k = k0, k0 + 1, ..., until the terms are
# negligible relative to the running sum and decreasing. Errors when the
# maximum number of terms is reached without convergence.
log_series <- function(term, k0 = 0, rel_tol = 1e-12, max_terms = 1e5,
                       chunk = 256L, what = "series") {
  total <- -Inf
  k <- k0
  prev_max <- -Inf
  repeat {
    ks <- k:(k + chunk - 1)
    lt <- term(ks)
    lt[is.na(lt)] <- -Inf
    total <- log_sum_exp(c(total, lt))
    last <- lt[length(lt)]
    if (is.finite(total) && last < total + log(rel_tol) &&
        last <= max(lt) && lt[length(lt)] <= lt[length(lt) - 1L]) break
    if (!is.finite(total) && all(lt == -Inf) && k - k0 > 4 * chunk) break
    k <- k + chunk
    if (k - k0 > max_terms)
      stop_numerical(sprintf("The %s did not converge within %d terms.",
                             what, as.integer(max_terms)))
  }
  total
}

# Numerical integration of a vectorized function on (lower, upper) with an
# error on failure.
integrate_safe <- function(f, lower, upper, rel_tol = 1e-10, ...) {
  res <- tryCatch(stats::integrate(f, lower, upper, rel.tol = rel_tol,
                                   subdivisions = 1000L, stop.on.error = FALSE, ...),
                  error = function(e) NULL)
  if (is.null(res) || !res$message %in% "OK")
    stop_numerical("Numerical integration failed to reach the requested tolerance.")
  res$value
}

# Generic numerical CDF of a continuous family by integrating its density from
# the lower end of the support.
numeric_cdf_cont <- function(lpdf, q, pars, lower) {
  vapply(seq_along(q), function(i) {
    if (is.na(q[i])) return(NA_real_)
    if (q[i] <= lower[i]) return(0)
    pi_ <- lapply(pars, `[`, i)
    f <- function(t) exp(lpdf(t, lapply(pi_, rep_len, length(t))))
    min(1, max(0, integrate_safe(f, lower[i], q[i])))
  }, 0)
}

# Generic CDF of a discrete family on the non-negative integers by summing
# the mass function.
numeric_cdf_disc <- function(lpmf, q, pars, lower = 0) {
  q <- floor(q)
  vapply(seq_along(q), function(i) {
    if (is.na(q[i])) return(NA_real_)
    lo <- if (length(lower) > 1) lower[i] else lower
    if (q[i] < lo) return(0)
    if (!is.finite(q[i])) return(1)
    pi_ <- lapply(pars, `[`, i)
    ys <- lo:q[i]
    min(1, sum(exp(lpmf(ys, lapply(pi_, rep_len, length(ys))))))
  }, 0)
}

# Quantile by inversion of a CDF. For discrete families returns the smallest
# integer y with F(y) >= p. Brackets the root and fails with a typed error
# when it cannot.
invert_cdf <- function(cdf, p, pars, discrete, lower = -Inf, upper = Inf,
                       start = 0, tol = 1e-10) {
  n <- max(length(p), vapply(pars, length, 1L))
  p <- rep_len(p, n)
  lower <- rep_len(lower, n); upper <- rep_len(upper, n); start <- rep_len(start, n)
  out <- numeric(n)
  for (i in seq_len(n)) {
    pi_ <- lapply(pars, function(v) rep_len(v, n)[i])
    if (is.na(p[i]) || any(is.na(unlist(pi_)))) { out[i] <- NA; next }
    if (p[i] < 0 || p[i] > 1) { out[i] <- NaN; next }
    if (p[i] == 0) { out[i] <- lower[i]; next }
    if (p[i] == 1) { out[i] <- upper[i]; next }
    F <- function(x) cdf(x, pi_) - p[i]
    if (discrete) {
      lo <- if (is.finite(lower[i])) lower[i] else floor(start[i]) - 1
      hi <- max(lo + 1, ceiling(start[i]))
      step <- 1
      while (F(hi) < 0) {
        lo <- hi; hi <- hi + step; step <- step * 2
        if (hi > upper[i]) { hi <- upper[i]; break }
        if (step > 2^60) stop_numerical("Quantile search did not bracket the root.")
      }
      if (!is.finite(lower[i])) {
        step <- 1
        while (F(lo) >= 0) { hi <- lo; lo <- lo - step; step <- step * 2
          if (step > 2^60) stop_numerical("Quantile search did not bracket the root.") }
      }
      while (hi - lo > 1) {
        mid <- floor((lo + hi) / 2)
        if (F(mid) >= 0) hi <- mid else lo <- mid
      }
      out[i] <- if (is.finite(lower[i]) && F(lower[i]) >= 0) lower[i] else hi
    } else {
      lo <- if (is.finite(lower[i])) lower[i] else start[i] - 1
      hi <- if (is.finite(upper[i])) upper[i] else start[i] + 1
      w <- 1
      while (!is.finite(lower[i]) && F(lo) > 0) {
        w <- w * 2; lo <- start[i] - w
        if (w > 1e300) stop_numerical("Quantile search did not bracket the root.")
      }
      w <- 1
      while (!is.finite(upper[i]) && F(hi) < 0) {
        w <- w * 2; hi <- start[i] + w
        if (w > 1e300) stop_numerical("Quantile search did not bracket the root.")
      }
      r <- tryCatch(stats::uniroot(F, c(lo, hi), tol = tol * max(1, abs(start[i])),
                                   maxiter = 2000L),
                    error = function(e) NULL)
      if (is.null(r)) stop_numerical("Quantile inversion did not converge.")
      out[i] <- r$root
    }
  }
  out
}

# Inverse transform sampling from a quantile function.
rng_by_inversion <- function(n, quantile, pars) {
  pars <- lapply(pars, rep_len, n)
  quantile(stats::runif(n), pars)
}

# Moment of a discrete family by summation of the mass function up to a
# tail threshold.
numeric_moment_disc <- function(lpmf, pars, fun, lower = 0, max_terms = 1e6) {
  n <- max(vapply(pars, length, 1L))
  pars <- lapply(pars, rep_len, n)
  vapply(seq_len(n), function(i) {
    pi_ <- lapply(pars, `[`, i)
    total <- 0; k <- lower; mass <- 0
    repeat {
      ys <- k:(k + 999)
      lp <- lpmf(ys, lapply(pi_, rep_len, length(ys)))
      pr <- exp(lp)
      total <- total + sum(fun(ys) * pr)
      mass <- mass + sum(pr)
      if (mass > 1 - 1e-14 && pr[length(pr)] < 1e-300) break
      if (mass > 1 - 1e-12 && all(diff(pr[900:1000]) <= 0) &&
          abs(fun(ys[1000]) * pr[1000]) < 1e-16 * max(1, abs(total))) break
      k <- k + 1000
      if (k - lower > max_terms) return(NA_real_)
    }
    total
  }, 0)
}

# Moment of a continuous family by quadrature over its support.
numeric_moment_cont <- function(lpdf, pars, fun, lower, upper) {
  n <- max(vapply(pars, length, 1L))
  pars <- lapply(pars, rep_len, n)
  lower <- rep_len(lower, n); upper <- rep_len(upper, n)
  vapply(seq_len(n), function(i) {
    pi_ <- lapply(pars, `[`, i)
    f <- function(t) fun(t) * exp(lpdf(t, lapply(pi_, rep_len, length(t))))
    tryCatch(integrate_safe(f, lower[i], upper[i], rel_tol = 1e-8),
             regcomp_numerical_failure = function(e) NA_real_)
  }, 0)
}
