# Model fitting entry point and the maximum likelihood engine.

#' Fit a distributional regression model
#'
#' Every regression parameter of the response distribution gets its own
#' linear predictor on its link scale. The first formula gives the response
#' and the predictor of the first regression parameter (usually the
#' location); further formulas such as `sigma ~ x` give predictors for other
#' parameters, which otherwise are constant. Predictors may contain
#' offsets, factors and interactions, random effects in `lme4` syntax
#' (`(1 | g)`, `(1 + x | g)`, uncorrelated `(1 + x || g)`), and smooth terms
#' `s()`, `te()`, `ti()` and `t2()` from 'mgcv', which are fitted in their
#' mixed model representation.
#'
#' @section Engines:
#' * `"ml"`: maximum likelihood in R. Without random effects or smooths the
#'   log likelihood is maximized by `nlminb()` with an analytic chain rule
#'   gradient. With random effects or smooths the marginal likelihood is
#'   approximated by the Laplace method (joint over all random effects, so
#'   crossed and nested factors and random effects in several parameters
#'   are allowed); the fitted object records `approximation = "laplace"`.
#'   Penalties (see [penalty()]) give penalized likelihood estimates.
#' * `"tmb"`: the same model through 'TMB' automatic differentiation for the
#'   families that declare a TMB implementation (see [capabilities()]).
#' * `"bayes"`: Bayesian inference. A 'Stan' program is generated from the
#'   same specification (see [stan_code()]) and run by 'bayescompR'
#'   (default, no compiler needed), 'cmdstanr' or 'rstan'. Priors are set
#'   with [prior()].
#'
#' @param formula Two sided formula `response ~ predictors` for the first
#'   regression parameter. The response may be wrapped in [cens()] or
#'   [surv()].
#' @param ... Further formulas for other parameters, such as `sigma ~ x`,
#'   and data arguments of the family (for example `trials = n` for the
#'   binomial), evaluated in `data`.
#' @param family A family identifier, family object or a family built by an
#'   operator such as [inflate()].
#' @param data A data frame.
#' @param engine `"ml"`, `"tmb"` or `"bayes"`.
#' @param weights Optional prior weights multiplying each log likelihood
#'   contribution.
#' @param links Optional named list or character vector of link names per
#'   parameter.
#' @param fixed Optional named list of parameters held fixed at given
#'   values (not estimated).
#' @param penalty Optional penalty from [penalty()] (engine `"ml"`).
#' @param priors Optional priors from [prior()] (engine `"bayes"`).
#' @param control List of control settings, see [regcomp_control()].
#' @param parameterization Optional parameterization name passed to
#'   [get_family()].
#' @param subset Optional logical or index vector selecting rows of `data`.
#' @return An object of class `regcomp_fit`.
#' @export
#' @examples
#' set.seed(1)
#' d <- data.frame(x = rnorm(300), z = rnorm(300))
#' d$y <- rnbinom(300, mu = exp(1 + 0.5 * d$x), size = exp(0.5 + 0.4 * d$z))
#' fit <- regcomp(y ~ x, phi ~ z, family = "nb2", data = d)
#' summary(fit)
#' head(predict(fit, type = "mean"))
#' predict(fit, newdata = data.frame(x = 0, z = 0), type = "quantile",
#'         probs = c(0.1, 0.9))
regcomp <- function(formula, ..., family, data, engine = c("ml", "tmb", "bayes"),
                    weights = NULL, links = NULL, fixed = NULL, penalty = NULL,
                    priors = NULL, control = regcomp_control(),
                    parameterization = NULL, subset = NULL) {
  engine <- match.arg(engine)
  cl <- match.call()
  fam <- get_family(family, parameterization)
  if (isTRUE(fam$extra$no_regression))
    stop_unsupported("Regression", fam$id)
  if (missing(data)) data <- environment(formula)
  data <- as.data.frame(data)
  if (!is.null(subset)) data <- data[subset, , drop = FALSE]
  dots <- eval(substitute(alist(...)))
  spec <- build_spec(formula, dots, fam, data, weights = substitute(weights),
                     links = links, fixed = fixed, env = parent.frame())
  spec$call <- cl
  fit <- switch(engine,
    ml = fit_ml(spec, control, penalty),
    tmb = fit_tmb(spec, control),
    bayes = fit_bayes(spec, priors, control))
  fit$call <- cl
  fit$engine <- engine
  fit
}

#' Control settings for model fitting
#'
#' @param maxit Maximum outer iterations of `nlminb()`.
#' @param rel_tol Relative tolerance of the outer optimization.
#' @param inner_maxit Maximum Newton iterations of the Laplace inner problem.
#' @param inner_tol Convergence tolerance of the inner problem.
#' @param restarts Number of random restarts tried when the first
#'   optimization fails.
#' @param hessian Compute the covariance of the estimates.
#' @param trace Print progress.
#' @param chains,iter_warmup,iter_sampling,seed,runner,cores Bayesian
#'   settings: chains, iterations, seed, the Stan runner (`"bayescompR"`,
#'   `"cmdstanr"` or `"rstan"`) and parallel chains.
#' @param adapt_delta Target acceptance for NUTS.
#' @param start Optional named list of starting values for the coefficients
#'   of each parameter (vectors on the link scale).
#' @return A list.
#' @export
#' @examples
#' regcomp_control(maxit = 500)
regcomp_control <- function(maxit = 1000, rel_tol = 1e-10, inner_maxit = 50,
                            inner_tol = 1e-8, restarts = 2, hessian = TRUE,
                            trace = FALSE, chains = 4, iter_warmup = 1000,
                            iter_sampling = 1000, seed = 1,
                            runner = c("bayescompR", "cmdstanr", "rstan"),
                            cores = 1, adapt_delta = 0.8, start = NULL) {
  list(maxit = maxit, rel_tol = rel_tol, inner_maxit = inner_maxit,
       inner_tol = inner_tol, restarts = restarts, hessian = hessian,
       trace = trace, chains = chains, iter_warmup = iter_warmup,
       iter_sampling = iter_sampling, seed = seed, runner = match.arg(runner),
       cores = cores, adapt_delta = adapt_delta, start = start)
}

# Assemble everything the engines need.
build_spec <- function(formula, dots, fam, data, weights = NULL, links = NULL,
                       fixed = NULL, env = parent.frame(), xlev = NULL,
                       prev = NULL) {
  is_formula <- vapply(dots, function(d) {
    v <- tryCatch(eval(d, env), error = function(e) NULL)
    inherits(v, "formula")
  }, TRUE)
  forms <- lapply(dots[is_formula], eval, env)
  if (length(forms)) names(forms) <- NULL
  dargs <- dots[!is_formula]
  fl <- parse_formulas(formula, forms, fam)
  for (k in names(fixed)) {
    if (!k %in% names(fam$pars)) stop_input(sprintf("Unknown fixed parameter '%s'.", k))
    fl[[k]] <- NULL
  }
  rpars <- names(fl)
  n <- nrow(data)
  y <- if (is.null(prev)) get_response(formula, data) else NULL
  cmat <- NULL
  if (inherits(y, "regcomp_cens")) { cmat <- unclass(y); y <- cmat[, "lower"] }
  # Links.
  lk <- list()
  for (k in rpars) {
    nm <- if (!is.null(links[[k]])) links[[k]] else fam$pars[[k]]$link
    if (!nm %in% fam$pars[[k]]$links)
      stop_input(sprintf("Link '%s' is not allowed for parameter '%s' (%s).", nm, k,
                         paste(fam$pars[[k]]$links, collapse = ", ")))
    lk[[k]] <- link(nm)
  }
  preds <- list()
  for (k in rpars) {
    preds[[k]] <- build_predictor(expand_dot(rhs_of(fl[[k]]), data, all.vars(formula[[2]])), data, k,
      xlev = if (!is.null(prev)) prev$preds[[k]]$xlev,
      smooth_specs = if (!is.null(prev)) prev$preds[[k]]$smooths,
      re_levels = if (!is.null(prev)) stats::setNames(
        lapply(prev$preds[[k]]$blocks, `[[`, "levels"),
        vapply(prev$preds[[k]]$blocks, `[[`, "", "label")))
    if (any(preds[[k]]$offset != 0) && !fam$pars[[k]]$offset_allowed &&
        is.null(prev))
      stop_input(sprintf("An offset is not allowed for parameter '%s' of '%s'.", k, fam$id))
  }
  # Data arguments (binomial trials and similar).
  dp <- list()
  for (k in data_pars(fam)) {
    if (!is.null(dargs[[k]])) dp[[k]] <- rep_len(eval(dargs[[k]], data, env), n)
    else if (!is.null(fixed[[k]])) dp[[k]] <- rep_len(fixed[[k]], n)
    else if (!is.null(fam$pars[[k]]$default)) dp[[k]] <- rep_len(fam$pars[[k]]$default, n)
    else stop_input(sprintf("The '%s' family needs the data argument '%s'.", fam$id, k))
  }
  unknown <- setdiff(names(dargs), c(data_pars(fam), ""))
  if (length(unknown)) stop_input(sprintf("Unknown argument%s: %s.",
    if (length(unknown) > 1) "s" else "", paste(unknown, collapse = ", ")))
  fixed <- fixed[setdiff(names(fixed), data_pars(fam))]
  w <- if (is.null(weights)) rep(1, n) else {
    wv <- eval(weights, data, env)
    if (is.null(wv)) rep(1, n) else as.numeric(wv)
  }
  # Coefficient layout.
  beta_idx <- list(); off <- 0
  for (k in rpars) {
    p <- ncol(preds[[k]]$X)
    beta_idx[[k]] <- off + seq_len(p); off <- off + p
  }
  beta_names <- unlist(lapply(rpars, function(k)
    paste0(k, ":", colnames(preds[[k]]$X))))
  # Random effect layout.
  blocks <- unlist(lapply(rpars, function(k) preds[[k]]$blocks), recursive = FALSE)
  Z <- list(); theta_idx <- list(); u_idx <- list(); toff <- 0; uoff <- 0
  Q <- sum(vapply(blocks, function(b) b$q * b$L, 1))
  for (j in seq_along(blocks)) {
    b <- blocks[[j]]
    nt <- b$q + if (b$corr) b$q * (b$q - 1) / 2 else 0
    theta_idx[[j]] <- toff + seq_len(nt); toff <- toff + nt
    u_idx[[j]] <- uoff + seq_len(b$q * b$L); uoff <- uoff + b$q * b$L
  }
  for (k in rpars) {
    js <- which(vapply(blocks, `[[`, "", "par") == k)
    if (!length(js)) { Z[[k]] <- NULL; next }
    Zk <- Matrix::Matrix(0, n, Q, sparse = TRUE)
    for (j in js) Zk[, u_idx[[j]]] <- block_Z(blocks[[j]], n)
    Z[[k]] <- Zk
  }
  structure(list(
    fam = fam, formulas = fl, rpars = rpars, links = lk, preds = preds,
    y = y, cmat = cmat, n = n, weights = w, data_pars = dp, fixed = fixed,
    beta_idx = beta_idx, beta_names = beta_names, blocks = blocks,
    theta_idx = theta_idx, u_idx = u_idx, Z = Z, Q = Q, ntheta = toff,
    data_args = dargs, response_name = paste(deparse(formula[[2]]), collapse = "")
  ), class = "regcomp_spec")
}

# Covariance matrix of one block from its unconstrained parameters.
block_cov <- function(b, th) {
  sd <- exp(th[seq_len(b$q)])
  if (!b$corr) return(diag(sd^2, b$q))
  L <- diag(b$q)
  L[lower.tri(L)] <- th[-seq_len(b$q)]
  L <- L / sqrt(rowSums(L^2))
  R <- tcrossprod(L)
  outer(sd, sd) * R
}

# Sparse precision of all random effects.
re_precision <- function(spec, theta) {
  mats <- lapply(seq_along(spec$blocks), function(j) {
    b <- spec$blocks[[j]]
    S <- block_cov(b, theta[spec$theta_idx[[j]]])
    Si <- solve(S)
    Matrix::kronecker(Matrix::Diagonal(b$L), Matrix::Matrix(Si))
  })
  Matrix::bdiag(mats)
}

re_logdet_cov <- function(spec, theta) {
  sum(vapply(seq_along(spec$blocks), function(j) {
    b <- spec$blocks[[j]]
    S <- block_cov(b, theta[spec$theta_idx[[j]]])
    b$L * as.numeric(determinant(S, logarithm = TRUE)$modulus)
  }, 0))
}

# Weighted log likelihood and its gradient with respect to beta.
ml_objective <- function(spec, beta, u = NULL) {
  eta <- compute_eta(spec, beta, u)
  ll <- loglik_eta(spec, eta)
  s <- sum(spec$weights * ll)
  if (is.na(s)) -Inf else s
}

ml_gradient <- function(spec, beta, u = NULL) {
  eta <- compute_eta(spec, beta, u)
  d <- eta_derivs(spec, eta, second = FALSE)
  out <- numeric(length(beta))
  for (a in seq_along(spec$rpars)) {
    k <- spec$rpars[a]
    out[spec$beta_idx[[k]]] <- drop(crossprod(spec$preds[[k]]$X, spec$weights * d$g[, a]))
  }
  out
}

# Inner Laplace problem: mode of the random effects for given beta, theta.
laplace_inner <- function(spec, beta, theta, u0, control) {
  Lam <- re_precision(spec, theta)
  u <- if (is.null(u0)) rep(0, spec$Q) else u0
  K <- length(spec$rpars)
  fval <- function(u) {
    eta <- compute_eta(spec, beta, u)
    sum(spec$weights * loglik_eta(spec, eta)) - 0.5 * sum(u * as.vector(Lam %*% u))
  }
  f <- fval(u)
  H <- NULL
  for (it in seq_len(control$inner_maxit)) {
    eta <- compute_eta(spec, beta, u)
    d <- eta_derivs(spec, eta, second = TRUE)
    grad <- -as.vector(Lam %*% u)
    H <- -Lam
    for (a in seq_len(K)) {
      ka <- spec$rpars[a]
      if (is.null(spec$Z[[ka]])) next
      grad <- grad + as.vector(Matrix::crossprod(spec$Z[[ka]], spec$weights * d$g[, a]))
      for (b in seq_len(K)) {
        kb <- spec$rpars[b]
        if (is.null(spec$Z[[kb]])) next
        w <- spec$weights * d$W[, a, b]
        H <- H + Matrix::crossprod(spec$Z[[ka]], Matrix::Diagonal(x = w) %*% spec$Z[[kb]])
      }
    }
    negH <- Matrix::forceSymmetric(-H)
    ch <- tryCatch(Matrix::Cholesky(negH, perm = TRUE, LDL = FALSE), error = function(e) NULL)
    lam <- 0
    while (is.null(ch)) {
      lam <- if (lam == 0) 1e-6 else lam * 10
      ch <- tryCatch(Matrix::Cholesky(negH + lam * Matrix::Diagonal(spec$Q), perm = TRUE),
                     error = function(e) NULL)
      if (lam > 1e8) break
    }
    step <- as.vector(Matrix::solve(ch, grad))
    t <- 1
    repeat {
      un <- u + t * step
      fn <- fval(un)
      if (is.finite(fn) && fn >= f - 1e-10) break
      t <- t / 2
      if (t < 1e-8) { un <- u; fn <- f; break }
    }
    conv <- abs(fn - f) < control$inner_tol * (1 + abs(f)) && max(abs(un - u)) < 1e-6
    u <- un; f <- fn
    if (conv) break
  }
  # Hessian at the mode.
  eta <- compute_eta(spec, beta, u)
  d <- eta_derivs(spec, eta, second = TRUE)
  H <- -Lam
  for (a in seq_len(K)) for (b in seq_len(K)) {
    ka <- spec$rpars[a]; kb <- spec$rpars[b]
    if (is.null(spec$Z[[ka]]) || is.null(spec$Z[[kb]])) next
    H <- H + Matrix::crossprod(spec$Z[[ka]],
      Matrix::Diagonal(x = spec$weights * d$W[, a, b]) %*% spec$Z[[kb]])
  }
  negH <- Matrix::forceSymmetric(-H)
  ld <- tryCatch(as.numeric(Matrix::determinant(negH, logarithm = TRUE)$modulus),
                 error = function(e) NA_real_)
  list(u = u, f = f, logdet = ld, negH = negH)
}

# Laplace approximation of the marginal log likelihood.
laplace_objective <- function(spec, beta, theta, state, control) {
  inner <- laplace_inner(spec, beta, theta, state$u, control)
  if (!is.finite(inner$f) || is.na(inner$logdet)) return(-Inf)
  state$u <- inner$u
  inner$f - 0.5 * re_logdet_cov(spec, theta) - 0.5 * inner$logdet
}

# Starting values: intercepts at link transformed family start values.
start_values <- function(spec, control) {
  fam <- spec$fam
  beta <- numeric(length(spec$beta_names))
  y <- spec$y
  ok <- is.finite(y)
  st <- if (!is.null(fam$start)) tryCatch(fam$start(y[ok]), error = function(e) list())
        else list()
  for (k in spec$rpars) {
    X <- spec$preds[[k]]$X
    v <- st[[k]]
    if (is.null(v)) v <- generic_start(fam$pars[[k]], y[ok])
    lv <- spec$links[[k]]$linkfun(v)
    if (!is.finite(lv)) lv <- 0
    lv <- lv - mean(spec$preds[[k]]$offset)
    ic <- which(colnames(X) == "(Intercept)")
    if (length(ic)) beta[spec$beta_idx[[k]][ic]] <- lv
    if (!is.null(control$start[[k]])) beta[spec$beta_idx[[k]]] <- control$start[[k]]
  }
  beta
}

generic_start <- function(par, y) {
  if (!is.null(par$default)) return(par$default)
  if (par$role %in% c("location") || par$is_mean) {
    m <- mean(y)
    if (par$domain == "positive") return(max(m, 1e-2))
    if (par$domain == "unit") return(min(max(m, 0.05), 0.95))
    return(m)
  }
  if (par$role %in% c("scale")) return(max(stats::sd(y), 1e-2))
  switch(par$domain, positive = 1, unit = 0.5, greater_than_one = 2,
         greater_than_two = 4, nonneg = 1, real = 0, circular = 0, 1)
}

fit_ml <- function(spec, control, penalty = NULL) {
  beta0 <- start_values(spec, control)
  has_re <- spec$Q > 0
  pen <- make_penalty(penalty, spec)
  if (!is.null(pen) && !pen$smooth && pen$type != "flac") {
    if (has_re) stop_unsupported(sprintf("The %s penalty with random effects or smooth terms (use ridge, firth, log_f or cauchy)", pen$type))
    return(fit_proximal(spec, beta0, pen, control))
  }
  if (!is.null(pen) && pen$type == "flac") {
    if (has_re) stop_unsupported("FLAC with random effects")
    fl <- flac_fit(spec, control)
    return(new_fit(spec, fl$beta, numeric(0), numeric(0), fl$cov, ml_objective(spec, fl$beta),
                   fl$opt, "none", pen))
  }
  pen_val <- function(beta) if (is.null(pen)) 0 else pen$value(beta)
  pen_grad <- function(beta) if (is.null(pen)) 0 else pen$grad(beta)
  if (!has_re) {
    obj <- function(b) {
      v <- ml_objective(spec, b) - pen_val(b)
      if (!is.finite(v)) 1e300 else -v
    }
    gr <- function(b) -(ml_gradient(spec, b) - pen_grad(b))
    opt <- robust_nlminb(beta0, obj, gr, control)
    beta <- opt$par
    if (!is.null(pen) && pen$type == "flic") beta <- flic_adjust(spec, beta, control)
    theta <- numeric(0); u <- numeric(0)
    cov <- NULL
    if (control$hessian && !is.null(pen) && pen$type %in% c("firth", "flic")) {
      cov <- safe_inverse(info_matrix(spec, beta))
    } else if (control$hessian) {
      H <- num_jacobian(function(b) -gr(b), beta)
      H <- (H + t(H)) / 2
      cov <- safe_inverse(-H)
    }
    ll <- ml_objective(spec, beta)
    approx <- "none"
  } else {
    np <- length(beta0)
    # Stage 1: fixed effects only for good starting values.
    obj0 <- function(b) { v <- ml_objective(spec, b) - pen_val(b); if (!is.finite(v)) 1e300 else -v }
    gr0 <- function(b) -(ml_gradient(spec, b) - pen_grad(b))
    b1 <- tryCatch(robust_nlminb(beta0, obj0, gr0, control)$par, error = function(e) beta0)
    theta0 <- unlist(lapply(spec$blocks, function(b)
      c(rep(log(0.5), b$q), if (b$corr) rep(0, b$q * (b$q - 1) / 2))))
    state <- new.env(); state$u <- rep(0, spec$Q)
    obj_re <- function(par) {
      v <- laplace_objective(spec, par[seq_len(np)], par[-seq_len(np)], state, control) -
        pen_val(par[seq_len(np)])
      if (!is.finite(v)) 1e300 else -v
    }
    opt <- robust_nlminb(c(b1, theta0), obj_re, NULL, control)
    beta <- opt$par[seq_len(np)]; theta <- opt$par[-seq_len(np)]
    inner <- laplace_inner(spec, beta, theta, state$u, control)
    u <- inner$u
    cov <- NULL
    if (control$hessian) {
      H <- num_hessian(function(p) -obj_re(p), opt$par)
      cov <- safe_inverse(-H)
    }
    ll <- -opt$objective + pen_val(beta)
    approx <- "laplace"
  }
  new_fit(spec, beta, theta, u, cov, ll, opt, approx, pen)
}

robust_nlminb <- function(start, obj, gr, control) {
  best <- NULL
  starts <- list(start)
  for (r in seq_len(control$restarts)) starts[[r + 1]] <- start + stats::rnorm(length(start), 0, 0.1 * r)
  for (s in starts) {
    o <- tryCatch(stats::nlminb(s, obj, gr, control = list(
      iter.max = control$maxit, eval.max = 2 * control$maxit, rel.tol = control$rel_tol,
      trace = as.integer(control$trace))), error = function(e) NULL)
    if (is.null(o)) next
    if (is.null(best) || o$objective < best$objective) best <- o
    if (o$convergence == 0 && o$objective < 1e299) break
  }
  if (is.null(best) || best$objective >= 1e299)
    stop_numerical("The optimizer could not find a finite log likelihood from the starting values.")
  best
}

# Forward difference Jacobian of a vector function.
num_jacobian <- function(f, x) {
  f0 <- f(x)
  J <- matrix(0, length(f0), length(x))
  for (j in seq_along(x)) {
    h <- 1e-5 * max(1, abs(x[j]))
    xp <- x; xp[j] <- xp[j] + h
    xm <- x; xm[j] <- xm[j] - h
    J[, j] <- (f(xp) - f(xm)) / (2 * h)
  }
  J
}

# Central difference Hessian of a scalar function.
num_hessian <- function(f, x) {
  d <- length(x)
  H <- matrix(0, d, d)
  h <- 1e-4 * pmax(1, abs(x))
  f0 <- f(x)
  for (i in seq_len(d)) {
    for (j in i:d) {
      if (i == j) {
        xp <- x; xp[i] <- xp[i] + h[i]; xm <- x; xm[i] <- xm[i] - h[i]
        H[i, i] <- (f(xp) - 2 * f0 + f(xm)) / h[i]^2
      } else {
        e <- function(si, sj) { z <- x; z[i] <- z[i] + si * h[i]; z[j] <- z[j] + sj * h[j]; f(z) }
        H[i, j] <- H[j, i] <- (e(1, 1) - e(1, -1) - e(-1, 1) + e(-1, -1)) / (4 * h[i] * h[j])
      }
    }
  }
  H
}

safe_inverse <- function(A) {
  out <- tryCatch(solve(A), error = function(e) NULL)
  if (is.null(out)) {
    ev <- eigen((A + t(A)) / 2, symmetric = TRUE)
    if (any(ev$values <= 0)) return(structure(matrix(NA_real_, nrow(A), ncol(A)),
                                              singular = TRUE))
    out <- ev$vectors %*% diag(1 / ev$values, length(ev$values)) %*% t(ev$vectors)
  }
  out
}

new_fit <- function(spec, beta, theta, u, cov, ll, opt, approx, pen = NULL,
                    engine = "ml") {
  names(beta) <- spec$beta_names
  tn <- theta_names(spec)
  names(theta) <- tn
  if (!is.null(cov)) {
    nn <- c(spec$beta_names, tn)[seq_len(nrow(cov))]
    dimnames(cov) <- list(nn, nn)
  }
  edf <- length(beta) + length(theta)
  fit <- structure(list(
    spec = spec, family = spec$fam, coefficients = beta, theta = theta, u = u,
    vcov = cov, logLik = ll, df = edf, nobs = sum(spec$weights > 0 & !is.na(spec$y)),
    optimizer = list(convergence = opt$convergence, message = opt$message,
                     iterations = opt$iterations, objective = opt$objective),
    approximation = approx, penalty = pen, engine = engine,
    inference = "frequentist"
  ), class = "regcomp_fit")
  fit$diagnostics <- fit_diagnostics(fit)
  fit
}

theta_names <- function(spec) {
  unlist(lapply(seq_along(spec$blocks), function(j) {
    b <- spec$blocks[[j]]
    nm <- paste0(b$par, ":", b$label, ":log_sd(", b$terms, ")")
    if (b$corr) {
      idx <- which(lower.tri(diag(b$q)), arr.ind = TRUE)
      nm <- c(nm, paste0(b$par, ":", b$label, ":chol(", b$terms[idx[, 1]], ",",
                         b$terms[idx[, 2]], ")"))
    }
    nm
  }))
}

fit_diagnostics <- function(fit) {
  spec <- fit$spec
  out <- list(converged = fit$optimizer$convergence == 0)
  if (fit$approximation == "none" && length(fit$coefficients)) {
    g <- tryCatch(ml_gradient(spec, unname(fit$coefficients)), error = function(e) NA)
    out$max_abs_gradient <- max(abs(g))
  }
  if (!is.null(fit$vcov)) {
    ev <- tryCatch(eigen(fit$vcov, symmetric = TRUE, only.values = TRUE)$values,
                   error = function(e) NA)
    out$hessian_positive_definite <- isTRUE(all(ev > 0))
    out$condition_number <- if (isTRUE(all(ev > 0))) max(ev) / min(ev) else NA
  }
  if (length(fit$theta)) {
    sds <- fit$theta[grepl("log_sd", names(fit$theta))]
    out$boundary_variance <- names(sds)[sds < log(1e-4)]
  }
  out
}
