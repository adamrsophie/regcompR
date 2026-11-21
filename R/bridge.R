# Marginal likelihoods of Stan fits by bridge sampling.

# Rewrite a generated Stan program so that the log density keeps every
# normalizing constant: sampling statements become target += *_lpdf, and
# priors on bounded parameters are renormalized to their support
# (half normal and half t standard deviations, uniform correlations).
stan_normalize <- function(code) {
  L <- strsplit(code, "\n", fixed = TRUE)[[1]]
  ps <- which(L == "parameters {"); pe <- ps + which(L[(ps + 1):length(L)] == "}")[1]
  dec <- regmatches(L[ps:pe], regexec("^\\s*(real|vector\\[[^]]+\\])<([^>]*)>\\s+(\\w+);", L[ps:pe]))
  bnd <- list()
  for (m in dec) if (length(m)) {
    b <- trimws(strsplit(m[3], ",")[[1]]); get <- function(k) { x <- grep(paste0("^", k, "\\s*="), b, value = TRUE); if (length(x)) trimws(sub("^[a-z]+\\s*=", "", x)) }
    bnd[[m[4]]] <- list(lower = get("lower"), upper = get("upper"), n = if (m[2] == "real") "1" else sprintf("num_elements(%s)", m[4]))
  }
  ms <- which(L == "model {"); me <- ms + which(L[(ms + 1):length(L)] == "}")[1]
  seen <- character(); extra <- character()
  for (i in (ms + 1):(me - 1)) {
    m <- regmatches(L[i], regexec("^(\\s*)(.+?) ~ (\\w+)\\((.*)\\);\\s*$", L[i]))[[1]]
    if (!length(m)) next
    tg <- m[3]; d <- m[4]; a <- m[5]
    L[i] <- sprintf("%starget += %s_lpdf(%s%s);", m[2], d, tg, if (nzchar(a)) paste0(" | ", a) else "")
    if (!is.null(bnd[[tg]]) && d != "lkj_corr_cholesky") {
      seen <- c(seen, tg); b <- bnd[[tg]]; arg <- if (nzchar(a)) paste0(" | ", a) else ""
      if (!is.null(b$lower) && is.null(b$upper))
        extra <- c(extra, sprintf("  target += -%s * %s_lccdf(%s%s);", b$n, d, b$lower, arg))
      else if (!is.null(b$lower) && !is.null(b$upper))
        extra <- c(extra, sprintf("  target += -%s * log_diff_exp(%s_lcdf(%s%s), %s_lcdf(%s%s));", b$n, d, b$upper, arg, d, b$lower, arg))
    }
  }
  # Bounded parameters without a prior have an implicit uniform density.
  for (nm in setdiff(names(bnd), seen)) {
    b <- bnd[[nm]]
    if (!is.null(b$lower) && !is.null(b$upper)) extra <- c(extra, sprintf("  target += -%s * log((%s) - (%s));", b$n, b$upper, b$lower))
  }
  paste(c(L[1:(me - 1)], extra, L[me:length(L)]), collapse = "\n")
}

#' Marginal likelihood and Bayes factors by bridge sampling
#'
#' Estimates the log marginal likelihood of a Bayesian fit made with the
#' rstan or cmdstanr runner by bridge sampling (Meng and Wong 1996; Gronau
#' et al. 2017) with the bridgesampling package. The model is re-expressed
#' with all normalizing constants (sampling statements as `target +=`,
#' priors on standard deviations and correlations renormalized to their
#' support), so marginal likelihoods of different models are comparable.
#' `bayes_factor_regcomp()` gives the Bayes factor of two fits. Run long
#' chains (several thousand draws) for stable estimates.
#'
#' @param object,object1,object2 Bayesian fits (`engine = "bayes"`, Stan
#'   runner rstan or cmdstanr).
#' @param ... Passed to `bridgesampling::bridge_sampler()`.
#' @return `bridge_regcomp()`: a `bridge` object with the estimate
#'   `logml`; `bayes_factor_regcomp()`: a data frame with the Bayes factor
#'   of the first fit against the second and its logarithm.
#' @export
bridge_regcomp <- function(object, ...) {
  if (!identical(object$inference, "bayesian") || !object$runner %in% c("rstan", "cmdstanr"))
    stop_unsupported("Bridge sampling for fits other than Stan fits with the rstan or cmdstanr runner")
  for (p in c("bridgesampling", "rstan")) if (!requireNamespace(p, quietly = TRUE))
    stop_unsupported(sprintf("Bridge sampling without the %s package", p))
  samples <- if (object$runner == "rstan") object$sampler else rstan::read_stan_csv(object$sampler$output_files())
  m0 <- rstan::stan(model_code = stan_normalize(object$code), data = stan_data(object$spec), chains = 0, refresh = 0)
  bridgesampling::bridge_sampler(samples, stanfit_model = m0, silent = TRUE, ...)
}

#' @rdname bridge_regcomp
#' @export
bayes_factor_regcomp <- function(object1, object2, ...) {
  l1 <- bridge_regcomp(object1, ...)$logml; l2 <- bridge_regcomp(object2, ...)$logml
  data.frame(BF = exp(l1 - l2), log_BF = l1 - l2, logml1 = l1, logml2 = l2)
}
