# Helpers used while defining families. Loaded first (file name aaa-).

P <- function(...) par_spec(...)

# Stan argument list for a wrapper function: integer data arguments (such
# as binomial trials) stay int so native lpmfs accept them.
stan_args <- function(pars, nm, ids = nm) {
  int <- vapply(nm, function(k) isTRUE(pars[[k]]$data) &&
    isTRUE(pars[[k]]$domain %in% c("positive_integer", "nonneg_integer", "integer")), TRUE)
  paste(ifelse(int, "int ", "real "), ids, sep = "", collapse = ", ")
}

# Stan template helper: native function call with {par} placeholders.
stan_native <- function(lpdf, lcdf = NULL, lccdf = NULL, rng = NULL,
                        functions = "") {
  list(native = TRUE, lpdf = lpdf, lcdf = lcdf, lccdf = lccdf, rng = rng,
       functions = functions)
}

# Stan template helper for package generated functions.
stan_generated <- function(functions, lpdf, lcdf = NULL, lccdf = NULL,
                           rng = NULL, needs = character()) {
  list(native = FALSE, lpdf = lpdf, lcdf = lcdf, lccdf = lccdf, rng = rng,
       functions = functions, needs = needs)
}

# Shared Stan helper functions, included on demand by name.
stan_helpers <- list(
  log_besselK = "
  // log K_nu(x) for real order by the integral representation
  // K_nu(x) = int_0^inf exp(-x cosh t) cosh(nu t) dt, trapezoid rule.
  real rc_log_besselK(real nu, real x) {
    real h = 0.02;
    real tmax = 1;
    real lmax;
    real acc;
    while (x * cosh(tmax) - abs(nu) * tmax < x + 50) tmax += 1;
    acc = log(0.5) - x;
    lmax = acc;
    real t = h;
    while (t <= tmax) {
      real lt = -x * cosh(t) + log_sum_exp(nu * t, -nu * t) - log(2);
      acc = log_sum_exp(acc, lt);
      t += h;
    }
    return acc + log(h);
  }",
  zeta = "
  // Riemann zeta function for s > 1 by Euler Maclaurin summation.
  real rc_zeta(real s) {
    int N = 20;
    real acc = 0;
    vector[6] b2k = [1.0 / 6, -1.0 / 30, 1.0 / 42, -1.0 / 30, 5.0 / 66, -691.0 / 2730]';
    real term;
    for (k in 1:(N - 1)) acc += pow(k, -s);
    acc += pow(N, 1 - s) / (s - 1) + 0.5 * pow(N, -s);
    term = s * pow(N, -s - 1);
    for (j in 1:6) {
      acc += b2k[j] / tgamma(2 * j + 1) * term;
      term = term * (s + 2 * j) * (s + 2 * j - 1) / (N * N);
    }
    return acc;
  }"
)

# Riemann zeta for s > 1 (vectorized) by Euler Maclaurin summation.
zeta_fn <- function(s) {
  N <- 20
  b2k <- c(1 / 6, -1 / 30, 1 / 42, -1 / 30, 5 / 66, -691 / 2730)
  vapply(s, function(si) {
    if (is.na(si) || si <= 1) return(NA_real_)
    acc <- sum((1:(N - 1))^(-si)) + N^(1 - si) / (si - 1) + 0.5 * N^(-si)
    term <- si * N^(-si - 1)
    for (j in 1:6) {
      acc <- acc + b2k[j] / factorial(2 * j) * term
      term <- term * (si + 2 * j) * (si + 2 * j - 1) / N^2
    }
    acc
  }, 0)
}

# Row wise log normalizer for count laws defined by an unnormalized log term
# lterm(k, pars) on k = 0, 1, .... Grows the truncation point until the
# terms beyond the mode are below rel_tol of the total for every row.
series_lognorm <- function(lterm, pars, k0 = 0, rel_tol = 1e-15,
                           max_k = 2e6, cap = NULL) {
  n <- max(vapply(pars, length, 1L))
  pars <- lapply(pars, rep_len, n)
  out <- rep(NA_real_, n)
  todo <- seq_len(n)
  K <- 64
  while (length(todo)) {
    ks <- k0:(k0 + K - 1)
    pp <- lapply(pars, function(v) rep(v[todo], each = length(ks)))
    lt <- matrix(lterm(rep(ks, length(todo)), pp), nrow = length(ks))
    lt[is.na(lt)] <- -Inf
    m <- apply(lt, 2, max)
    tot <- m + log(colSums(exp(lt - rep(m, each = length(ks)))))
    tot[!is.finite(m)] <- -Inf
    last <- lt[length(ks), ]
    prev <- lt[length(ks) - 1L, ]
    done <- (last < tot + log(rel_tol) & last <= prev) | !is.finite(m) |
      (last == -Inf & prev == -Inf)
    if (!is.null(cap)) done <- done | (k0 + K - 1 >= cap(lapply(pars, `[`, todo)))
    out[todo[done]] <- tot[done]
    todo <- todo[!done]
    K <- K * 4
    if (K > max_k && length(todo))
      stop_numerical("A normalizing series did not converge; parameters may be extreme.")
  }
  out
}

# Build lpmf, cdf, mean and variance for a count law from its unnormalized
# log term.
normalized_count <- function(lterm, k0 = 0) {
  logZ <- function(p) series_lognorm(lterm, p, k0)
  lpmf <- function(y, p) lterm(y, p) - logZ(p)
  cdf <- function(q, p) {
    n <- length(q)
    lz <- logZ(p)
    vapply(seq_len(n), function(i) {
      if (is.na(q[i])) return(NA_real_)
      if (q[i] < k0) return(0)
      if (!is.finite(q[i])) return(1)
      ks <- k0:floor(q[i])
      pi_ <- lapply(p, function(v) rep(v[i], length(ks)))
      min(1, sum(exp(lterm(ks, pi_) - lz[i])))
    }, 0)
  }
  moment <- function(p, fun) {
    n <- max(vapply(p, length, 1L))
    p <- lapply(p, rep_len, n)
    lz <- logZ(p)
    vapply(seq_len(n), function(i) {
      K <- 256
      repeat {
        ks <- k0:(k0 + K - 1)
        pi_ <- lapply(p, function(v) rep(v[i], length(ks)))
        pr <- exp(lterm(ks, pi_) - lz[i])
        if (sum(pr) > 1 - 1e-13 || K > 4e6) return(sum(fun(ks, i) * pr))
        K <- K * 4
      }
    }, 0)
  }
  mean <- function(p) moment(p, function(k, i) k)
  variance <- function(p) {
    m <- mean(p)
    moment(p, function(k, i) (k - m[i])^2)
  }
  list(lpmf = lpmf, cdf = cdf, mean = mean, variance = variance, logZ = logZ)
}

# Stan code for a series normalized count law. `lterm_body` is Stan code
# for the unnormalized log term with argument names k and the parameters.
stan_series_count <- function(name, args, lterm_body, k0 = 0, cap = NULL) {
  sig <- paste(sprintf("real %s", args), collapse = ", ")
  call <- paste(args, collapse = ", ")
  cap_line <- if (is.null(cap)) "" else
    sprintf("    if (k > %s) break;\n", cap)
  sprintf("
  real %1$s_lterm(int k, %2$s) {
%3$s
  }
  real %1$s_lognorm(%2$s) {
    real lz = negative_infinity();
    real prev = negative_infinity();
    int k = %5$d;
    while (1) {
      real lt;
%6$s      lt = %1$s_lterm(k, %4$s);
      lz = log_sum_exp(lz, lt);
      if (k > %5$d + 10 && lt < lz - 37 && lt <= prev) break;
      if (k > 100000) reject(\"series did not converge\");
      prev = lt;
      k += 1;
    }
    return lz;
  }
  real %1$s_lpmf(int y, %2$s) {
    return %1$s_lterm(y, %4$s) - %1$s_lognorm(%4$s);
  }", name, sig, lterm_body, call, k0, cap_line)
}

gamlss_ref <- function(fun, map, discrete = TRUE, extra = list()) {
  list(package = "gamlss.dist", fun = fun, map = map, extra = extra)
}

# Reference implementation descriptor: a package and functions computing the
# log density and (optionally) CDF on regcompR's parameterization.
ref_impl <- function(package, lpdf, cdf = NULL, note = "") {
  list(package = package, lpdf = lpdf, cdf = cdf, note = note)
}

# Generated Stan functions for a scalar family. Bodies are Stan statements
# using `y` and the parameter names; they must return a real.
sfun <- function(id, args, lpdf, lcdf = NULL, lccdf = NULL, discrete = FALSE,
                 rng = NULL, needs = character(), ints = character()) {
  sid <- paste0("rc_", gsub("[^A-Za-z0-9]", "_", id))
  yt <- if (discrete) "int" else "real"
  suf <- if (discrete) "lpmf" else "lpdf"
  sig <- paste(ifelse(args %in% ints, "int", "real"), args, collapse = ", ")
  call <- paste(sprintf("{%s}", args), collapse = ", ")
  f <- sprintf("\n  real %s_%s(%s y, %s) {\n%s\n  }", sid, suf, yt, sig, lpdf)
  if (!is.null(lcdf)) f <- paste0(f, sprintf("\n  real %s_lcdf(%s y, %s) {\n%s\n  }", sid, yt, sig, lcdf))
  if (!is.null(lcdf) && is.null(lccdf))
    lccdf <- sprintf("    return log1m_exp(%s_lcdf(y | %s));", sid, paste(args, collapse = ", "))
  if (!is.null(lccdf)) f <- paste0(f, sprintf("\n  real %s_lccdf(%s y, %s) {\n%s\n  }", sid, yt, sig, lccdf))
  if (!is.null(rng)) f <- paste0(f, sprintf("\n  %s %s_rng(%s) {\n%s\n  }", if (discrete) "int" else "real", sid, sig, rng))
  list(native = FALSE, functions = f, needs = needs,
       lpdf = sprintf("%s_%s({y} | %s)", sid, suf, call),
       lcdf = if (!is.null(lcdf)) sprintf("%s_lcdf({y} | %s)", sid, call),
       lccdf = if (!is.null(lccdf)) sprintf("%s_lccdf({y} | %s)", sid, call),
       rng = if (!is.null(rng)) sprintf("%s_rng(%s)", sid, call))
}

loc_par <- function(name = "mu", mean = FALSE) P(name, "real", role = "location", is_mean = mean,
                                                  offset_allowed = TRUE)
scale_par <- function(name = "sigma") P(name, "positive", role = "scale")

JKK <- "Johnson, Kemp and Kotz (2005), Univariate Discrete Distributions, 3rd ed. Wiley."

# Copy a family object with some fields replaced.
modify_family <- function(f, ...) {
  a <- list(...)
  for (n in names(a)) f[n] <- list(a[[n]])
  f
}

GS <- "Rigby, Stasinopoulos, Heller and De Bastiani (2019), Distributions for Modeling Location, Scale, and Shape: Using GAMLSS in R. Chapman and Hall/CRC."
# Reference wrappers for gamlss.dist (Suggests): resolved only when called.
gd <- function(f, ...) function(y, p) { fn <- getExportedValue("gamlss.dist", f)
  do.call(fn, c(list(y), p[intersect(names(p), names(formals(fn)))], list(log = TRUE))) }
gp <- function(f) function(q, p) { fn <- getExportedValue("gamlss.dist", f)
  do.call(fn, c(list(q), p[intersect(names(p), names(formals(fn)))])) }

mean_par <- function(name = "mu") P(name, "positive", role = "location", is_mean = TRUE, offset_allowed = TRUE)
shp <- function(name, ...) P(name, "positive", role = "shape", ...)
JKB <- "Johnson, Kotz and Balakrishnan (1994, 1995), Continuous Univariate Distributions, vols 1 and 2, 2nd ed. Wiley."
KLUG <- "Klugman, Panjer and Willmot (2012), Loss Models: From Data to Decisions, 4th ed. Wiley."


# Density by Fourier inversion of a characteristic function phi(t) for a law
# with real values: f(x) = (1/pi) int_0^inf Re(exp(-i t x) phi(t)) dt.
cf_density <- function(x, cf, upper = Inf) vapply(x, function(xx)
  max(stats::integrate(function(t) Re(exp(-1i * t * xx) * cf(t)) / pi, 0, upper,
                       subdivisions = 2000L, rel.tol = 1e-9, stop.on.error = FALSE)$value, 0), 0)
pervec <- function(y, p, fun) { r <- recycle_pars(c(list(.y = y), p))
  vapply(seq_along(r$.y), function(i) fun(r$.y[i], lapply(r[-1], `[`, i)), 0) }


`%||%` <- function(a, b) if (is.null(a)) b else a
