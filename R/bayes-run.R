# Running generated programs and mapping draws back to the specification.

fit_bayes <- function(spec, priors, control) {
  runner <- control$runner
  st <- vapply(spec$blocks, function(b) b$struct, "")
  if (runner %in% c("jags", "openbugs") && (length(spec$aux) || any(!st %in% c("us", "diag")) ||
      any(vapply(spec$blocks, function(b) !is.null(b$W), TRUE))))
    stop_unsupported("Ordinal thresholds, structured or multiple membership random effects with the BUGS runners (use a Stan runner)")
  if (runner %in% c("jags", "openbugs")) {
    res <- run_bugs(spec, priors, control, runner)
  } else {
    code <- stan_program(spec, priors)
    data <- stan_data(spec)
    res <- run_stan(code, data, control, runner)
    res$code <- code
  }
  draws <- map_draws(spec, res$draws, res$language)
  fit <- structure(list(
    spec = spec, family = spec$fam, draws = draws$fixed, u_draws = draws$u,
    theta = { tn <- theta_names(spec); stats::setNames(if (all(tn %in% colnames(draws$fixed))) colMeans(draws$fixed[, tn, drop = FALSE]) else rep(NA_real_, length(tn)), tn) },
    coefficients = colMeans(draws$fixed[, spec$beta_names, drop = FALSE]),
    rhat = draws$rhat, ess = draws$ess, runner = runner, code = res$code,
    approximation = "none", inference = "bayesian", engine = "bayes",
    nobs = obs_count(spec), sampler = res$sampler,
    chain = rep(seq_len(dim(res$draws)[2]), each = dim(res$draws)[1]),
    df = length(spec$beta_names) + spec$ntheta
  ), class = "regcomp_fit")
  fit$diagnostics <- list(max_rhat = suppressWarnings(max(draws$rhat, na.rm = TRUE)),
                          min_ess_bulk = suppressWarnings(min(draws$ess, na.rm = TRUE)),
                          divergences = res$divergences)
  if (isTRUE(fit$diagnostics$max_rhat > 1.01))
    warning("Some R-hat values exceed 1.01; run longer chains before using the posterior.", call. = FALSE)
  fit
}

run_stan <- function(code, data, control, runner) {
  if (runner == "bayescompR") {
    if (!requireNamespace("bayescompR", quietly = TRUE))
      stop_unsupported("The 'bayescompR' runner needs the 'bayescompR' package; it")
    m <- bayescompR::bayescomp_model(code, data = data, seed = control$seed)
    f <- bayescompR::bayescomp_sample(m, chains = control$chains, iter_warmup = control$iter_warmup,
                                      iter_sampling = control$iter_sampling, seed = control$seed,
                                      refresh = 0, adapt_delta = control$adapt_delta,
                                      cores = control$cores, init = control$init %||% 2)
    arr <- bayescompR::bayescomp_draws(f, format = "array")
    div <- tryCatch(sum(bayescompR::bayescomp_draws(f, variables = "divergent__", format = "matrix")),
                    error = function(e) NA)
    return(list(draws = arr, language = "stan", sampler = f, divergences = div))
  }
  if (runner == "cmdstanr") {
    if (!requireNamespace("cmdstanr", quietly = TRUE))
      stop_unsupported("The 'cmdstanr' runner needs the 'cmdstanr' package; it")
    file <- cmdstanr::write_stan_file(code)
    m <- cmdstanr::cmdstan_model(file)
    f <- m$sample(data = data, chains = control$chains, iter_warmup = control$iter_warmup,
                  iter_sampling = control$iter_sampling, seed = control$seed, refresh = 0,
                  adapt_delta = control$adapt_delta, parallel_chains = control$cores,
                  init = control$init, show_messages = FALSE, show_exceptions = FALSE)
    arr <- unclass(f$draws())
    div <- tryCatch(sum(f$diagnostic_summary(quiet = TRUE)$num_divergent), error = function(e) NA)
    return(list(draws = arr, language = "stan", sampler = f, divergences = div))
  }
  if (!requireNamespace("rstan", quietly = TRUE))
    stop_unsupported("The 'rstan' runner needs the 'rstan' package; it")
  ini <- control$init %||% 2
  f <- rstan::stan(model_code = code, data = data, chains = control$chains,
                   warmup = control$iter_warmup, iter = control$iter_warmup + control$iter_sampling,
                   seed = control$seed, refresh = 0, cores = control$cores,
                   init = if (is.numeric(ini)) "random" else ini, init_r = if (is.numeric(ini)) ini else 2,
                   control = list(adapt_delta = control$adapt_delta))
  arr <- as.array(f)
  div <- tryCatch(sum(rstan::get_divergent_iterations(f)), error = function(e) NA)
  list(draws = arr, language = "stan", sampler = f, divergences = div)
}

run_bugs <- function(spec, priors, control, runner) {
  dialect <- if (runner == "jags") "jags" else "openbugs"
  code <- bugs_program(spec, priors, dialect)
  data <- bugs_data(spec)
  params <- c(paste0("b_", spec$rpars))
  for (j in seq_along(spec$blocks)) {
    b <- spec$blocks[[j]]
    params <- c(params, if (b$type == "smooth") c(paste0("sd_", j), paste0("s_", j))
                else if (b$corr) c(paste0("Sigma_", j), paste0("r_", j))
                else c(paste0("sd_", j), paste0("r_", j)))
  }
  file <- tempfile(fileext = ".bug")
  writeLines(code, file)
  iter <- control$iter_warmup + control$iter_sampling
  inits <- function() {
    out <- list()
    for (k in spec$rpars) out[[paste0("b_", k)]] <- rep(0, ncol(spec$preds[[k]]$X))
    out
  }
  # R2jags and R2OpenBUGS also draw from R's generator.
  if (!is.null(control$seed)) local_seed(control$seed)
  if (runner == "jags") {
    if (!requireNamespace("R2jags", quietly = TRUE))
      stop_unsupported("The 'jags' runner needs the 'R2jags' package and JAGS; it")
    # One JAGS random number generator seed per chain makes runs reproducible.
    jinits <- lapply(seq_len(control$chains), function(ch) c(inits(), list(.RNG.name = "base::Mersenne-Twister",
      .RNG.seed = as.integer(control$seed %||% sample.int(1e6, 1)) + ch)))
    f <- R2jags::jags(data = data, inits = jinits, parameters.to.save = params,
                      model.file = file, n.chains = control$chains, n.iter = iter,
                      n.burnin = control$iter_warmup, n.thin = 1, DIC = FALSE,
                      jags.seed = control$seed, progress.bar = "none", quiet = TRUE)
    arr <- f$BUGSoutput$sims.array
  } else {
    if (!requireNamespace("R2OpenBUGS", quietly = TRUE))
      stop_unsupported("The 'openbugs' runner needs the 'R2OpenBUGS' package and OpenBUGS; it")
    f <- R2OpenBUGS::bugs(data = data, inits = inits, parameters.to.save = params,
                          model.file = file, n.chains = control$chains, n.iter = iter,
                          n.burnin = control$iter_warmup, n.thin = 1, DIC = FALSE,
                          bugs.seed = control$seed)
    arr <- f$sims.array
  }
  list(draws = arr, language = "bugs", sampler = f, divergences = NA, code = code)
}

# Convert a draws array (iterations x chains x variables) to coefficient,
# variance parameter and random effect draws in the package's naming.
map_draws <- function(spec, arr, language) {
  if (language == "stan") spec <- stan_spec(spec)
  vn <- dimnames(arr)[[3]]
  get <- function(name) {
    i <- match(name, vn)
    if (is.na(i) && grepl("\\[1\\]$", name)) i <- match(sub("\\[1\\]$", "", name), vn)
    if (is.na(i)) stop_numerical(sprintf("Variable '%s' missing from the sampler output.", name))
    as.vector(arr[, , i])
  }
  S <- dim(arr)[1] * dim(arr)[2]
  fixed <- matrix(NA_real_, S, length(spec$beta_names) + spec$ntheta,
                  dimnames = list(NULL, c(spec$beta_names, theta_names(spec))))
  rh <- es <- stats::setNames(rep(NA_real_, ncol(fixed)), colnames(fixed))
  put <- function(col, name) {
    fixed[, col] <<- get(name)
    a <- matrix(get(name), dim(arr)[1])
    if (is.null(dim(a))) a <- matrix(a, ncol = 1)
    rh[col] <<- rhat_basic(a); es[col] <<- ess_basic(a)
  }
  ix <- function(base, ...) {
    idx <- c(...)
    if (language == "bugs" || length(idx) > 0) sprintf("%s[%s]", base, paste(idx, collapse = ","))
    else base
  }
  for (k in spec$rpars) for (j in seq_along(spec$beta_idx[[k]]))
    put(spec$beta_idx[[k]][j], ix(paste0("b_", k), j))
  for (k in names(spec$aux)) for (j in seq_along(spec$beta_idx[[k]]))
    put(spec$beta_idx[[k]][j], sprintf("%s[%d]", k, j))
  u <- if (spec$Q > 0) matrix(0, S, spec$Q) else NULL
  np <- length(spec$beta_names)
  for (j in seq_along(spec$blocks)) {
    b <- spec$blocks[[j]]
    ti <- np + spec$theta_idx[[j]]
    if (b$type == "smooth") {
      sdname <- paste0("sd_", j)
      fixed[, ti[1]] <- log(get(sdname))
      sdv <- get(sdname)
      for (l in seq_len(b$L)) u[, spec$u_idx[[j]][l]] <- if (language == "stan")
        sdv * get(sprintf("z_%d[%d]", j, l)) else get(sprintf("s_%d[%d]", j, l))
      next
    }
    if (!b$struct %in% c("us", "diag")) {
      hom <- b$struct != "cs"
      sdl <- if (hom) log(get(sprintf("sdh_%d", j))) else NULL
      if (b$struct == "cs") for (e in seq_len(b$q)) fixed[, ti[e]] <- log(get(sprintf("sd_%d[%d]", j, e)))
      else fixed[, ti[1]] <- sdl
      last <- ti[length(ti)]
      if (b$struct %in% c("cs", "homcs")) { a <- 1 / (b$q - 1); fixed[, last] <- stats::qlogis((get(sprintf("rho_%d", j)) + a) / (1 + a)) }
      else if (b$struct == "ar1") fixed[, last] <- atanh(get(sprintf("rho_%d", j)))
      else if (b$struct != "homdiag") fixed[, last] <- log(get(sprintf("range_%d", j)))
    } else if (language == "bugs" && b$corr) {
      Sig <- array(NA_real_, c(S, b$q, b$q))
      for (a in seq_len(b$q)) for (c2 in seq_len(b$q)) Sig[, a, c2] <- get(sprintf("Sigma_%d[%d,%d]", j, a, c2))
      for (e in seq_len(b$q)) fixed[, ti[e]] <- 0.5 * log(Sig[, e, e])
      idx <- which(lower.tri(diag(b$q)), arr.ind = TRUE)
      fr <- t(vapply(seq_len(S), function(s) {
        Lc <- t(chol(stats::cov2cor(Sig[s, , ])))
        Lc[idx] / diag(Lc)[idx[, 1]]
      }, numeric(nrow(idx))))
      if (nrow(idx) == 1) fr <- t(fr)
      fixed[, ti[-seq_len(b$q)]] <- fr
    } else {
      for (e in seq_len(b$q)) {
        nm <- if (language == "stan") sprintf("sd_%d[%d]", j, e) else sprintf("sd_%d[%d]", j, e)
        fixed[, ti[e]] <- log(get(nm))
      }
      if (b$corr) {
        idx <- which(lower.tri(diag(b$q)), arr.ind = TRUE)
        for (r in seq_len(nrow(idx))) {
          a <- idx[r, 1]; c2 <- idx[r, 2]
          fixed[, ti[b$q + r]] <- get(sprintf("Lc_%d[%d,%d]", j, a, c2)) / get(sprintf("Lc_%d[%d,%d]", j, a, a))
        }
      }
    }
    for (l in seq_len(b$L)) for (e in seq_len(b$q)) {
      nm <- if (language == "stan") sprintf("r_%d[%d,%d]", j, e, l) else sprintf("r_%d[%d,%d]", j, l, e)
      u[, spec$u_idx[[j]][(l - 1) * b$q + e]] <- get(nm)
    }
  }
  for (c2 in seq_len(ncol(fixed))) if (is.na(rh[c2])) {
    a <- matrix(fixed[, c2], dim(arr)[1])
    rh[c2] <- rhat_basic(a); es[c2] <- ess_basic(a)
  }
  list(fixed = fixed, u = u, rhat = rh, ess = es)
}

# Rank normalized split R-hat (Vehtari et al. 2021), iterations x chains.
rhat_basic <- function(a) {
  a <- as.matrix(a)
  if (any(!is.finite(a))) return(NA_real_)
  n <- nrow(a)
  if (n < 4) return(NA_real_)
  h <- floor(n / 2)
  sp <- cbind(a[seq_len(h), , drop = FALSE], a[n - h + seq_len(h), , drop = FALSE])
  z <- matrix(stats::qnorm((rank(sp) - 3 / 8) / (length(sp) + 1 / 4)), nrow(sp))
  rh <- function(x) {
    m <- ncol(x); nn <- nrow(x)
    B <- nn * stats::var(colMeans(x)); W <- mean(apply(x, 2, stats::var))
    if (W == 0) return(NA_real_)
    sqrt(((nn - 1) / nn * W + B / nn) / W)
  }
  zf <- matrix(stats::qnorm((rank(abs(sp - stats::median(sp))) - 3 / 8) / (length(sp) + 1 / 4)), nrow(sp))
  max(rh(z), rh(zf))
}

# Bulk effective sample size (Vehtari et al. 2021): split chains, rank
# normalize, Geyer initial positive and monotone sequence estimators.
ess_basic <- function(a) {
  a <- as.matrix(a)
  if (any(!is.finite(a))) return(NA_real_)
  n <- nrow(a)
  if (n < 4) return(NA_real_)
  h <- floor(n / 2)
  sp <- cbind(a[seq_len(h), , drop = FALSE], a[n - h + seq_len(h), , drop = FALSE])
  if (diff(range(sp)) == 0) return(NA_real_)
  z <- matrix(stats::qnorm((rank(sp) - 3 / 8) / (length(sp) + 1 / 4)), nrow(sp))
  ess_raw(z)
}
ess_raw <- function(x) {
  n <- nrow(x); m <- ncol(x)
  acov <- vapply(seq_len(m), function(j) stats::acf(x[, j], lag.max = n - 1, type = "covariance",
                                                  plot = FALSE, demean = TRUE)$acf[, 1, 1], numeric(n))
  acov <- matrix(acov, n)
  mean_var <- mean(acov[1, ]) * n / (n - 1)
  var_plus <- mean_var * (n - 1) / n + if (m > 1) stats::var(colMeans(x)) else 0
  rho <- numeric(n); rho_of <- function(t) 1 - (mean_var - mean(acov[t, ])) / var_plus
  t <- 0; even <- 1; rho[1] <- even; odd <- rho_of(2); rho[2] <- odd
  while (t < n - 5 && !is.nan(even + odd) && even + odd > 0) {
    t <- t + 2; even <- rho_of(t + 1); odd <- rho_of(t + 2)
    if (even + odd >= 0) { rho[t + 1] <- even; rho[t + 2] <- odd }
  }
  max_t <- t
  if (even > 0) rho[max_t + 1] <- even
  t <- 0
  while (t <= max_t - 4) {
    t <- t + 2
    if (rho[t + 1] + rho[t + 2] > rho[t - 1] + rho[t]) { rho[t + 1] <- (rho[t - 1] + rho[t]) / 2; rho[t + 2] <- rho[t + 1] }
  }
  S <- n * m
  tau <- max(-1 + 2 * sum(rho[seq_len(max_t)]) + rho[max_t + 1], 1 / log10(S))
  S / tau
}
