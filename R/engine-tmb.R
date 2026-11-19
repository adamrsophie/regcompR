# The TMB engine: the same specification fitted with automatic
# differentiation and the Laplace approximation of 'TMB', through 'RTMB'.

# Log densities written with operations that RTMB can tape. Each takes the
# response and a list of natural parameters (advectors) and returns the
# vector of log densities.
tmb_lpdf <- list(
  normal = function(y, p) RTMB::dnorm(y, p$mu, p$sigma, log = TRUE),
  poisson = function(y, p) RTMB::dpois(y, p[[1]], log = TRUE),
  binomial = function(y, p) RTMB::dbinom(y, p$size, p$prob, log = TRUE),
  bernoulli = function(y, p) RTMB::dbinom(y, 1, p[[1]], log = TRUE),
  nb2 = function(y, p) RTMB::dnbinom2(y, p$mu, p$mu + p$mu^2 / p$phi, log = TRUE),
  nb1 = function(y, p) RTMB::dnbinom2(y, p$mu, p$mu * (1 + p$phi), log = TRUE),
  student_t = function(y, p) RTMB::dt((y - p$mu) / p$sigma, p$nu, log = TRUE) - log(p$sigma),
  gamma = function(y, p) RTMB::dgamma(y, 1 / p$sigma^2, scale = p$mu * p$sigma^2, log = TRUE),
  lognormal = function(y, p) RTMB::dnorm(log(y), p$mu, p$sigma, log = TRUE) - log(y),
  weibull = function(y, p) RTMB::dweibull(y, p$shape, p$scale, log = TRUE),
  beta = function(y, p) RTMB::dbeta(y, p$mu * p$phi, (1 - p$mu) * p$phi, log = TRUE),
  exponential = function(y, p) RTMB::dexp(y, p[[1]], log = TRUE)
)

# Families checked against the registry density at random points before
# the engine accepts them (the parameter names and meanings must agree).
tmb_family_key <- function(fam) {
  key <- fam$tmb %||% switch(fam$id, bernoulli = "bernoulli", normal = "normal", gamma = "gamma",
    lognormal = "lognormal", weibull = "weibull", beta = "beta", exponential = "exponential", NULL)
  if (is.null(key) || is.null(tmb_lpdf[[key]])) return(NULL)
  key
}

tmb_linkinv <- function(name) switch(name,
  identity = function(x) x, log = function(x) exp(x),
  logit = function(x) 1 / (1 + exp(-x)), probit = function(x) RTMB::pnorm(x),
  cloglog = function(x) 1 - exp(-exp(x)), inverse = function(x) 1 / x,
  sqrt = function(x) x^2, softplus = function(x) log(1 + exp(x)),
  fisherz = function(x) tanh(x), stop_unsupported(sprintf("The '%s' link in the TMB engine", name)))

fit_tmb <- function(spec, control) {
  if (!requireNamespace("RTMB", quietly = TRUE)) stop_unsupported("The TMB engine needs the 'RTMB' package; it")
  fam <- spec$fam
  key <- tmb_family_key(fam)
  if (is.null(key)) stop_unsupported("The TMB engine", fam$id)
  if (!is.null(spec$cmat)) stop_unsupported("Censored or truncated responses in the TMB engine (use engine = \"ml\")")
  if (length(spec$aux)) stop_unsupported("Estimated thresholds in the TMB engine")
  # Check the taped density against the registry density.
  pr <- eta_to_pars(spec, compute_eta(spec, start_values(spec, control)))
  i <- which(is.finite(spec$y))[seq_len(min(5, sum(is.finite(spec$y))))]
  sub <- lapply(pr, `[`, i)
  dif <- tryCatch(max(abs(tmb_lpdf[[key]](spec$y[i], sub) - eval_lpdf(fam, spec$y[i], sub))), error = function(e) Inf)
  if (!is.finite(dif) || dif > 1e-8)
    stop_unsupported("The TMB engine for this parameterization", fam$id)
  np <- length(spec$beta_names); Q <- spec$Q
  lk <- lapply(spec$links, function(l) tmb_linkinv(l$name))
  y <- spec$y; w <- spec$weights
  nll <- function(par) {
    beta <- par$beta; u <- if (Q > 0) par$u else NULL
    p <- spec$data_pars
    for (k in spec$rpars) {
      pk <- spec$preds[[k]]
      eta <- drop(pk$X %*% beta[spec$beta_idx[[k]]]) + pk$offset
      if (Q > 0 && !is.null(spec$Z[[k]])) eta <- eta + drop(spec$Z[[k]] %*% u)
      p[[k]] <- lk[[k]](eta)
    }
    for (k in names(spec$fixed)) p[[k]] <- rep(spec$fixed[[k]], spec$n)
    out <- -sum(w * tmb_lpdf[[key]](y, p))
    if (Q > 0) for (j in seq_along(spec$blocks)) {
      b <- spec$blocks[[j]]
      th <- par$theta[spec$theta_idx[[j]]]
      uj <- u[spec$u_idx[[j]]]
      if (b$struct %in% c("diag", "homdiag")) {
        sdv <- if (b$struct == "diag") exp(th) else rep(exp(th[1]), b$q)
        out <- out - sum(RTMB::dnorm(uj, 0, rep(sdv, b$L), log = TRUE))
      } else {
        X <- RTMB::AD(matrix(0, b$L, b$q))
        X[cbind(rep(seq_len(b$L), each = b$q), rep(seq_len(b$q), b$L))] <- uj
        out <- out - sum(RTMB::dmvnorm(X, Sigma = block_cov(b, th), log = TRUE))
      }
    }
    out
  }
  beta0 <- start_values(spec, control)
  pl <- list(beta = beta0)
  if (Q > 0) { pl$u <- rep(0, Q); pl$theta <- unlist(lapply(spec$blocks, block_theta_start)) }
  obj <- RTMB::MakeADFun(nll, pl, random = if (Q > 0) "u", silent = !control$trace)
  opt <- stats::nlminb(obj$par, obj$fn, obj$gr, control = list(iter.max = control$maxit,
    eval.max = 2 * control$maxit, rel.tol = control$rel_tol))
  rep <- if (control$hessian) RTMB::sdreport(obj) else NULL
  est <- opt$par
  beta <- est[names(est) == "beta"]
  theta <- est[names(est) == "theta"]
  u <- if (Q > 0) obj$env$last.par.best[obj$env$random] else numeric(0)
  cov <- if (!is.null(rep)) rep$cov.fixed else NULL
  fit <- new_fit(spec, unname(beta), unname(theta), unname(u), if (!is.null(cov)) unname(cov), -opt$objective,
                 opt, if (Q > 0) "laplace" else "none", engine = "tmb")
  fit$tmb <- obj
  fit
}
