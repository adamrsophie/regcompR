# Observation level log likelihood with censoring and truncation, and its
# derivatives with respect to the linear predictors.

# Log probability of the observation y given parameters, honoring cens
# matrices. `pars` is a named list of length n vectors (including data
# parameters).
obs_loglik <- function(fam, y, pars, cmat = NULL) {
  if (is.null(cmat)) return(eval_lpdf(fam, y, pars))
  n <- nrow(cmat)
  lo <- cmat[, "lower"]; up <- cmat[, "upper"]
  tl <- cmat[, "tlower"]; tu <- cmat[, "tupper"]
  out <- rep(NA_real_, n)
  exact <- !is.na(lo) & lo == up
  sub <- function(idx) lapply(pars, `[`, idx)
  if (any(exact)) out[exact] <- eval_lpdf(fam, lo[exact], sub(exact))
  ce <- !exact & !is.na(lo)
  # A censoring interval only counts where it meets the truncation region.
  lo[ce] <- pmax(lo[ce], tl[ce]); up[ce] <- pmin(up[ce], tu[ce])
  empty <- ce & (if (fam$discrete) lo > up else lo >= up)
  # Interval [lo, up] for discrete laws and (lo, up] for continuous ones,
  # as log P(a < Y <= b) with a = lo - 1 for discrete laws.
  lower_q <- function(q) if (fam$discrete) ceiling(q) - 1 else q
  if (any(ce)) out[ce] <- family_log_interval(fam, lower_q(lo[ce]), up[ce], sub(ce))
  tr <- (is.finite(tl) | is.finite(tu)) & !is.na(out)
  if (any(tr)) {
    out[tr] <- out[tr] - family_log_interval(fam, lower_q(tl[tr]), tu[tr], sub(tr))
    outside <- (exact & tr) & ((if (fam$discrete) lo < tl else lo <= tl) | lo > tu)
    out[outside] <- -Inf
  }
  out[empty] <- -Inf
  out
}

# Natural parameters from linear predictors.
eta_to_pars <- function(spec, eta) {
  pars <- spec$data_pars
  for (k in spec$rpars) pars[[k]] <- spec$links[[k]]$linkinv(eta[[k]])
  for (k in names(spec$fixed)) pars[[k]] <- if (k %in% names(spec$aux) || !is.null(spec$fam$pars[[k]]$estimate))
    spec$fixed[[k]] else rep_len(spec$fixed[[k]], spec$n)
  for (k in names(attr(eta, "aux"))) pars[[k]] <- attr(eta, "aux")[[k]]
  pars[names(spec$fam$pars)]
}

# Weighted per observation log likelihood from linear predictors.
loglik_eta <- function(spec, eta) {
  pars <- eta_to_pars(spec, eta)
  # Family wide constraints (beyond each parameter's own domain): a
  # violation means the law is not a probability distribution.
  if (!is.null(spec$fam$check)) {
    msg <- tryCatch(spec$fam$check(pars), error = function(e) conditionMessage(e))
    if (is.character(msg) && length(msg) && nzchar(msg[1])) return(rep(-Inf, spec$n))
  }
  obs_loglik(spec$fam, spec$y, pars, spec$cmat)
}

# First and second derivatives of per observation log likelihood with
# respect to each linear predictor, by central finite differences.
# Returns g (n x K) and W (n x K x K) together with ll.
# First and second derivatives of the log density in one parameter, for
# common families (uncensored rows); others use finite differences.
analytic_dpar <- list(
  normal = list(mu = function(y, p) list((y - p$mu) / p$sigma^2, rep(-1, length(y)) / p$sigma^2)),
  poisson = list(lambda = function(y, p) list(y / p$lambda - 1, -y / p$lambda^2)),
  bernoulli = list(prob = function(y, p) list(y / p$prob - (1 - y) / (1 - p$prob), -y / p$prob^2 - (1 - y) / (1 - p$prob)^2)),
  binomial = list(prob = function(y, p) list(y / p$prob - (p$size - y) / (1 - p$prob), -y / p$prob^2 - (p$size - y) / (1 - p$prob)^2)),
  negative_binomial_nb2 = list(mu = function(y, p) list(y / p$mu - (y + p$phi) / (p$mu + p$phi), -y / p$mu^2 + (y + p$phi) / (p$mu + p$phi)^2)))

# d theta / d eta and d2 theta / d eta2 of an inverse link.
link_d12 <- function(lk, e) {
  switch(lk$name,
    identity = list(rep(1, length(e)), rep(0, length(e))),
    log = { v <- exp(e); list(v, v) },
    logit = { v <- stats::plogis(e); d <- v * (1 - v); list(d, d * (1 - 2 * v)) },
    { h <- 1e-4 * pmax(1, abs(e)); i0 <- lk$linkinv(e); ip <- lk$linkinv(e + h); im <- lk$linkinv(e - h)
      list((ip - im) / (2 * h), (ip - 2 * i0 + im) / h^2) })
}

eta_derivs <- function(spec, eta, second = TRUE, which = spec$rpars) {
  rp <- which
  if (length(rp) == 1 && is.null(spec$cmat) && !is.null(f1 <- analytic_dpar[[spec$fam$id]][[rp]])) {
    ll0 <- loglik_eta(spec, eta)
    pars <- eta_to_pars(spec, eta); d <- f1(spec$y, pars); hd <- link_d12(spec$links[[rp]], eta[[rp]])
    g <- d[[1]] * hd[[1]]; w <- d[[2]] * hd[[1]]^2 + d[[1]] * hd[[2]]
    g[!is.finite(g)] <- 0; w[!is.finite(w)] <- 0
    return(list(ll = ll0, g = matrix(g, ncol = 1), W = if (second) array(w, c(spec$n, 1, 1))))
  }
  K <- length(rp)
  n <- spec$n
  ll0 <- loglik_eta(spec, eta)
  g <- matrix(0, n, K)
  W <- if (second) array(0, c(n, K, K)) else NULL
  h <- lapply(eta, function(e) 1e-4 * pmax(1, abs(e)))
  shifted <- function(k, s, k2 = NULL, s2 = 0) {
    e <- eta
    e[[k]] <- e[[k]] + s * h[[k]]
    if (!is.null(k2)) e[[k2]] <- e[[k2]] + s2 * h[[k2]]
    loglik_eta(spec, e)
  }
  plus <- vector("list", K); minus <- vector("list", K)
  for (a in seq_len(K)) {
    k <- rp[a]
    lp <- shifted(k, 1); lm <- shifted(k, -1)
    plus[[a]] <- lp; minus[[a]] <- lm
    ga <- (lp - lm) / (2 * h[[k]])
    # One sided differences where a shift leaves the support.
    bad <- !is.finite(ga)
    if (any(bad)) {
      fwd <- (lp - ll0) / h[[k]]; bwd <- (ll0 - lm) / h[[k]]
      ga[bad] <- ifelse(is.finite(fwd[bad]), fwd[bad], bwd[bad])
    }
    ga[!is.finite(ga)] <- 0
    g[, a] <- ga
    if (second) {
      wa <- (lp - 2 * ll0 + lm) / h[[k]]^2
      wa[!is.finite(wa)] <- 0
      W[, a, a] <- wa
    }
  }
  if (second && K > 1) for (a in 1:(K - 1)) for (b in (a + 1):K) {
    ka <- rp[a]; kb <- rp[b]
    v <- (shifted(ka, 1, kb, 1) - shifted(ka, 1, kb, -1) -
            shifted(ka, -1, kb, 1) + shifted(ka, -1, kb, -1)) / (4 * h[[ka]] * h[[kb]])
    v[!is.finite(v)] <- 0
    W[, a, b] <- v; W[, b, a] <- v
  }
  list(ll = ll0, g = g, W = W)
}

# Linear predictors from coefficient vector and random effects.
compute_eta <- function(spec, beta, u = NULL) {
  out <- list()
  for (k in spec$rpars) {
    pr <- spec$preds[[k]]
    e <- drop(pr$X %*% beta[spec$beta_idx[[k]]]) + pr$offset
    if (!is.null(u) && !is.null(spec$Z[[k]]))
      e <- e + if (!is.null(spec$Zd[[k]])) drop(spec$Zd[[k]] %*% u) else drop(as.matrix(spec$Z[[k]] %*% u))
    out[[k]] <- e
  }
  if (length(spec$aux)) attr(out, "aux") <- lapply(stats::setNames(nm = names(spec$aux)), function(k)
    aux_transform(spec$aux[[k]]$type, beta[spec$beta_idx[[k]]]))
  out
}
