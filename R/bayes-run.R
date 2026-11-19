# Running generated programs and mapping draws back to the specification.

fit_bayes <- function(spec, priors, control) {
  runner <- control$runner
  if (length(spec$aux))
    stop_unsupported("Bayesian fitting of ordinal thresholds (use engine = \"ml\")")
  if (any(vapply(spec$blocks, function(b) !is.null(b$W), TRUE)))
    stop_unsupported("Bayesian fitting of multiple membership random effects (use engine = \"ml\")")
  st <- vapply(spec$blocks, function(b) b$struct, "")
  if (any(!st %in% c("us", "diag")))
    stop_unsupported(sprintf("Bayesian fitting of '%s' random effect structures (use engine = \"ml\")",
                             paste(unique(st[!st %in% c("us", "diag")]), collapse = "', '")))
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
    theta = stats::setNames(numeric(spec$ntheta), theta_names(spec)),
    coefficients = colMeans(draws$fixed[, spec$beta_names, drop = FALSE]),
    rhat = draws$rhat, ess = draws$ess, runner = runner, code = res$code,
    approximation = "none", inference = "bayesian", engine = "bayes",
    nobs = sum(spec$weights > 0 & !is.na(spec$y)), sampler = res$sampler,
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
                                      cores = control$cores)
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
                  show_messages = FALSE, show_exceptions = FALSE)
    arr <- unclass(f$draws())
    div <- tryCatch(sum(f$diagnostic_summary(quiet = TRUE)$num_divergent), error = function(e) NA)
    return(list(draws = arr, language = "stan", sampler = f, divergences = div))
  }
  if (!requireNamespace("rstan", quietly = TRUE))
    stop_unsupported("The 'rstan' runner needs the 'rstan' package; it")
  f <- rstan::stan(model_code = code, data = data, chains = control$chains,
                   warmup = control$iter_warmup, iter = control$iter_warmup + control$iter_sampling,
                   seed = control$seed, refresh = 0, cores = control$cores,
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
  if (runner == "jags") {
    if (!requireNamespace("R2jags", quietly = TRUE))
      stop_unsupported("The 'jags' runner needs the 'R2jags' package and JAGS; it")
    f <- R2jags::jags(data = data, inits = inits, parameters.to.save = params,
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
    if (language == "bugs" && b$corr) {
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

# Bulk effective sample size from autocorrelations (Geyer initial
# monotone sequence) of rank normalized draws.
ess_basic <- function(a) {
  a <- as.matrix(a)
  if (any(!is.finite(a))) return(NA_real_)
  n <- nrow(a); m <- ncol(a)
  if (n < 4) return(NA_real_)
  z <- matrix(stats::qnorm((rank(a) - 3 / 8) / (length(a) + 1 / 4)), n)
  ac <- sapply(seq_len(m), function(c2) {
    x <- z[, c2] - mean(z[, c2])
    v <- stats::acf(x, lag.max = n - 1, plot = FALSE, demean = FALSE)$acf[, 1, 1]
    v
  })
  rho <- rowMeans(as.matrix(ac))
  W <- mean(apply(z, 2, stats::var)); B <- n * stats::var(colMeans(z))
  vp <- (n - 1) / n * W + if (m > 1) B / n else 0
  rho_hat <- 1 - (W - rowMeans(as.matrix(ac)) * W) / vp
  rho_hat[1] <- 1
  s <- 0; k <- 1
  while (k + 1 <= length(rho_hat)) {
    pr <- rho_hat[k] + rho_hat[k + 1]
    if (pr < 0) break
    s <- s + pr; k <- k + 2
  }
  tau <- -1 + 2 * s
  min(n * m / max(tau, 1 / log10(n * m)), n * m * log10(n * m))
}
