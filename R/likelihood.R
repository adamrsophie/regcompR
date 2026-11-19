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
  lcdf <- function(q, idx) {
    if (!any(idx)) return(numeric())
    v <- family_cdf(fam, q, sub(idx))
    log(pmax(v, 0))
  }
  # P(Y <= q) and P(Y < q) for the interval bounds.
  below <- function(q, idx, strict) {
    qq <- if (fam$discrete && strict) ceiling(q) - 1 else q
    res <- rep(0, sum(idx))
    fin <- is.finite(qq[idx])
    res[!fin & qq[idx] > 0] <- 1
    if (any(fin)) {
      ii <- which(idx)[fin]
      res[fin] <- family_cdf(fam, qq[ii], sub(ii))
    }
    res
  }
  ce <- !exact & !is.na(lo)
  if (any(ce)) {
    # Interval [lo, up] for discrete; (lo, up] for continuous.
    Fu <- below(up, ce, FALSE)
    Fl <- below(lo, ce, fam$discrete)
    right <- is.infinite(up[ce]) & up[ce] > 0
    v <- log(pmax(Fu - Fl, 0))
    if (any(right) && !is.null(fam$lccdf)) {
      ii <- which(ce)[right]
      qq <- if (fam$discrete) ceiling(lo[ii]) - 1 else lo[ii]
      v[right] <- fam$lccdf(qq, sub(ii))
    }
    out[ce] <- v
  }
  tr <- (is.finite(tl) | is.finite(tu)) & !is.na(out)
  if (any(tr)) {
    Ftu <- below(tu, tr, FALSE)
    Ftl <- below(tl, tr, fam$discrete)
    out[tr] <- out[tr] - log(pmax(Ftu - Ftl, 0))
    outside <- (exact & tr) & ((if (fam$discrete) lo < tl else lo <= tl) | lo > tu)
    out[outside] <- -Inf
  }
  out
}

# Natural parameters from linear predictors.
eta_to_pars <- function(spec, eta) {
  pars <- spec$data_pars
  for (k in spec$rpars) pars[[k]] <- spec$links[[k]]$linkinv(eta[[k]])
  for (k in names(spec$fixed)) pars[[k]] <- rep_len(spec$fixed[[k]], spec$n)
  pars[names(spec$fam$pars)]
}

# Weighted per observation log likelihood from linear predictors.
loglik_eta <- function(spec, eta) {
  ll <- obs_loglik(spec$fam, spec$y, eta_to_pars(spec, eta), spec$cmat)
  ll
}

# First and second derivatives of per observation log likelihood with
# respect to each linear predictor, by central finite differences.
# Returns g (n x K) and W (n x K x K) together with ll.
eta_derivs <- function(spec, eta, second = TRUE) {
  K <- length(spec$rpars)
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
    k <- spec$rpars[a]
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
    ka <- spec$rpars[a]; kb <- spec$rpars[b]
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
    if (!is.null(u) && !is.null(spec$Z[[k]])) e <- e + drop(as.matrix(spec$Z[[k]] %*% u))
    out[[k]] <- e
  }
  out
}
