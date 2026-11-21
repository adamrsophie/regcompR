#' Random intercepts with non-normal distributions
#'
#' Fits a model with one random intercept per level of `group` in the
#' first parameter's linear predictor, whose distribution is not normal:
#' a log-gamma effect (`exp(b)` gamma with mean 1 and variance `v`, the
#' gamma frailty of a Poisson or survival model), a Student t effect with
#' scale `s` and estimated degrees of freedom, a logistic or a Laplace
#' effect with scale `s`. The marginal likelihood is computed by adaptive
#' Gauss-Hermite quadrature over `b = F^-1(Phi(z))` with standard normal
#' `z`, so any family of the package can be used.
#'
#' @param formula Model formula for the first parameter, without random
#'   effects; other parameters (`sigma ~ x`) go in `...`.
#' @param group A one sided formula naming the grouping factor, `~ g`.
#' @param family Family (any regcompR family).
#' @param data A data frame.
#' @param re `"gamma"`, `"t"`, `"logistic"` or `"laplace"`.
#' @param nq Number of quadrature nodes.
#' @param ... Formulas for the other parameters.
#' @return An object of class `regcomp_nongauss` with coefficients
#'   (including `re:log_v` or `re:log_scale`, and `re:log_df_minus_2` for t
#'   effects), their covariance, the log likelihood, AIC and the predicted
#'   random effects (posterior means).
#' @export
#' @examples
#' set.seed(1)
#' d <- data.frame(g = factor(rep(1:30, each = 5)), x = rnorm(150))
#' d$y <- rpois(150, exp(0.3 + 0.4 * d$x) * rgamma(30, 2, 2)[d$g])
#' re_nongauss(y ~ x, ~ g, family = "poisson", data = d, re = "gamma")
re_nongauss <- function(formula, group, family, data, re = c("gamma", "t", "logistic", "laplace"), nq = 20, ...) {
  re <- match.arg(re)
  cl <- match.call()
  f0 <- regcomp(formula, ..., family = family, data = data)
  spec <- f0$spec
  if (spec$Q > 0) stop_unsupported("Other random effects in re_nongauss() (use one grouping factor)")
  k1 <- spec$rpars[1]
  gv <- factor(eval(group[[2]], f0$model_data, environment(group)))
  gi <- as.integer(gv); G <- nlevels(gv)
  rq <- switch(re,
    gamma = function(u, th) { v <- exp(th[1]); log(stats::qgamma(u, 1 / v, 1 / v)) },
    t = function(u, th) exp(th[1]) * stats::qt(u, 2 + exp(th[2])),
    logistic = function(u, th) exp(th[1]) * stats::qlogis(u),
    laplace = function(u, th) exp(th[1]) * ifelse(u < 0.5, log(2 * u), -log(2 * (1 - u))))
  thn <- c(if (re == "gamma") "re:log_v" else "re:log_scale", if (re == "t") "re:log_df_minus_2")
  gh <- gauss_hermite(nq)
  nb <- length(f0$coefficients)
  # Group sums of the log likelihood with effect b(z_g) for each group.
  hz <- function(beta, th, z) {
    u <- stats::pnorm(z); u <- pmin(pmax(u, 1e-15), 1 - 1e-15)
    eta <- compute_eta(spec, beta); eta[[k1]] <- eta[[k1]] + rq(u, th)[gi]
    ll <- spec$weights * loglik_eta(spec, eta); ll[is.na(ll)] <- -Inf
    drop(rowsum(ll, gi, reorder = TRUE)) + stats::dnorm(z, log = TRUE)
  }
  marg <- function(p, keep = FALSE) {
    beta <- p[seq_len(nb)]; th <- p[-seq_len(nb)]
    z <- rep(0, G); dz <- 1e-4
    for (it in 1:30) {
      h0 <- hz(beta, th, z); hp <- hz(beta, th, z + dz); hm <- hz(beta, th, z - dz)
      g1 <- (hp - hm) / (2 * dz); g2 <- (hp - 2 * h0 + hm) / dz^2
      st <- ifelse(is.finite(g2) & g2 < 0, -g1 / g2, sign(g1) * 0.5)
      st[!is.finite(st)] <- 0; st <- pmax(pmin(st, 1), -1); z <- z + st
      if (max(abs(st)) < 1e-8) break
    }
    h0 <- hz(beta, th, z); hp <- hz(beta, th, z + dz); hm <- hz(beta, th, z - dz)
    s <- 1 / sqrt(pmax(-(hp - 2 * h0 + hm) / dz^2, 1e-8))
    H <- vapply(seq_len(nq), function(q) hz(beta, th, z + sqrt(2) * s * gh$x[q]) + gh$x[q]^2 + log(gh$w[q]), numeric(G))
    H <- matrix(H, G)
    mx <- apply(H, 1, max); lg <- mx + log(rowSums(exp(H - mx))) + log(sqrt(2) * s)
    if (keep) {
      W <- exp(H - mx); W <- W / rowSums(W)
      zn <- outer(z, rep(1, nq)) + sqrt(2) * outer(s, gh$x)
      b <- rowSums(W * matrix(rq(pmin(pmax(stats::pnorm(zn), 1e-15), 1 - 1e-15), th), G))
      return(list(ll = sum(lg), b = stats::setNames(b, levels(gv))))
    }
    s <- sum(lg); if (is.finite(s)) s else -1e300
  }
  p0 <- c(unname(f0$coefficients), if (re == "gamma") log(0.3) else log(0.5), if (re == "t") log(5))
  o <- stats::nlminb(p0, function(p) -marg(p), control = list(eval.max = 2000, iter.max = 1000, rel.tol = 1e-10))
  nm <- c(names(f0$coefficients), thn)
  H <- num_hessian(function(p) marg(p), o$par)
  V <- tryCatch(solve(-H), error = function(e) matrix(NA_real_, length(o$par), length(o$par)))
  dimnames(V) <- list(nm, nm)
  fin <- marg(o$par, keep = TRUE)
  structure(list(coefficients = stats::setNames(o$par, nm), vcov = V, logLik = fin$ll, df = length(o$par),
                 nobs = spec$n, ranef = fin$b, re = re, family = spec$fam$id, convergence = o$convergence, call = cl),
            class = "regcomp_nongauss")
}

#' @export
print.regcomp_nongauss <- function(x, digits = 4, ...) {
  cat(sprintf("<regcomp_nongauss> %s model with a %s random intercept, n = %d\n", x$family, x$re, x$nobs))
  se <- sqrt(diag(x$vcov))
  print(round(data.frame(estimate = x$coefficients, std_error = se, z = x$coefficients / se), digits))
  cat(sprintf("log likelihood %.3f, AIC %.3f\n", x$logLik, -2 * x$logLik + 2 * x$df))
  invisible(x)
}
#' @export
coef.regcomp_nongauss <- function(object, ...) object$coefficients
#' @export
vcov.regcomp_nongauss <- function(object, ...) object$vcov
#' @export
logLik.regcomp_nongauss <- function(object, ...) structure(object$logLik, df = object$df, nobs = object$nobs, class = "logLik")
#' @export
nobs.regcomp_nongauss <- function(object, ...) object$nobs
