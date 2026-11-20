# Proper scoring rules, cross validation, model comparison, user defined
# families and identifiability diagnostics.

#' Proper scoring rules for predictive distributions
#'
#' Scores each observation under the fitted predictive distribution
#' (Gneiting and Raftery 2007). Smaller is better for every rule.
#' * `"log"`: negative log predictive density or mass.
#' * `"crps"`: continuous ranked probability score
#'   `int (F(z) - 1(y <= z))^2 dz`, by quadrature for continuous laws and
#'   by summation for counts.
#' * `"dss"`: Dawid-Sebastiani score `(y - m)^2 / v + log v`.
#' * `"brier"` and `"rps"`: Brier and ranked probability scores for
#'   ordinal and nominal families (summed over the classes). For bernoulli
#'   and binomial families `"brier"` is the usual `(y / size - p)^2` with
#'   `p` the predicted success probability.
#'
#' For Bayesian fits every rule scores the posterior predictive
#' distribution, the mixture of the laws of the posterior draws (at most
#' 1000, evenly thinned): the log score uses the mean density, the CRPS the
#' mean CDF, the Dawid-Sebastiani score the total predictive mean and
#' variance, and the Brier and ranked probability scores the mean class
#' probabilities. The CRPS is infinite when the predictive law has no finite
#' mean; a quadrature that does not converge raises
#' `regcomp_numerical_failure`.
#'
#' @param object A fitted model.
#' @param newdata Optional data frame with the response; defaults to the
#'   fitting data.
#' @param rule Scoring rule.
#' @param re Random effects used for prediction, as in [predict()]:
#'   `"conditional"` uses the estimated effects of known groups,
#'   `"zero"` sets them to zero (new groups).
#' @return Numeric vector of scores, one per observation.
#' @export
#' @examples
#' set.seed(1)
#' d <- data.frame(x = rnorm(100)); d$y <- rnorm(100, 1 + d$x)
#' fit <- regcomp(y ~ x, family = "normal", data = d)
#' mean(scoring_rule(fit, rule = "crps"))
scoring_rule <- function(object, newdata = NULL, rule = c("log", "crps", "dss", "brier", "rps"),
                         re = c("conditional", "zero")) {
  rule <- match.arg(rule); re <- match.arg(re)
  spec <- new_spec(object, newdata, response = TRUE)
  # Parameter sets of the predictive law: the estimate, or posterior draws.
  PS <- if (object$inference == "bayesian") {
    S <- nrow(object$draws); idx <- if (S <= 1000) seq_len(S) else unique(round(seq(1, S, length.out = 1000)))
    lapply(idx, function(s) pars_from(spec, unname(object$draws[s, spec$beta_names]),
                                      if (length(object$u_draws)) object$u_draws[s, ] else NULL, re))
  } else list(pars_from(spec, unname(object$coefficients), object$u, re))
  pars <- PS[[1]]
  y <- spec$y; fam <- spec$fam
  if (rule == "log") {
    # Bayesian fits: the posterior predictive density, log mean over draws.
    if (object$inference == "bayesian") {
      ll <- log_lik(object, newdata, re = re)
      return(-apply(ll, 2, function(v) log_sum_exp(v) - log(length(v))))
    }
    return(-obs_loglik(fam, y, pars, spec$cmat))
  }
  if (!is.null(spec$cmat) && any(spec$cmat[, "lower"] != spec$cmat[, "upper"], na.rm = TRUE))
    stop_unsupported(sprintf("The '%s' score for censored observations", rule))
  if (rule %in% c("crps", "dss") && !is.null(spec$cmat) && any(is.finite(spec$cmat[, c("tlower", "tupper")])))
    stop_unsupported(sprintf("The '%s' score for truncated observations", rule))
  if (rule == "brier" && fam$id %in% c("bernoulli", "binomial")) {
    m <- if (fam$id == "binomial") rep_len(pars$size, length(y)) else 1
    p <- Reduce(`+`, lapply(PS, function(q) family_moment(fam, "mean", q) / m)) / length(PS)
    return((y / m - p)^2)
  }
  if (rule %in% c("brier", "rps")) {
    if (fam$variate_type != "categorical") stop_input("Brier and ranked probability scores need an ordinal or nominal family.")
    # All fitted classes, whatever occurs in the assessment set.
    P <- Reduce(`+`, lapply(PS, function(p) cat_probs(fam, p, length(y)))) / length(PS); K <- ncol(P)
    if (any(!y %in% seq_len(K))) stop_input("Responses must be fitted categories 1, ..., K.")
    O <- outer(y, seq_len(K), `==`) * 1
    if (rule == "brier") return(rowSums((P - O)^2))
    return(rowSums((t(apply(P, 1, cumsum)) - t(apply(O, 1, cumsum)))^2)[, drop = TRUE])
  }
  if (fam$variate_type != "scalar") stop_unsupported(sprintf("The '%s' score for this family", rule))
  if (rule == "dss") {
    # Total predictive moments: mean of means; mean variance plus the
    # variance of the means.
    M <- sapply(PS, function(p) mean(structure(list(family = fam, pars = p), class = "regcomp_dist")))
    V <- sapply(PS, function(p) variance(structure(list(family = fam, pars = p), class = "regcomp_dist")))
    M <- matrix(M, length(y)); V <- matrix(V, length(y))
    m <- rowMeans(M); v <- rowMeans(V) + rowMeans(M^2) - m^2
    return((y - m)^2 / v + log(v))
  }
  PP <- lapply(PS, recycle_pars)
  vapply(seq_along(y), function(i) {
    pis <- lapply(PP, function(pp) lapply(pp, `[`, i))
    # Mean CDF over the parameter sets, in one vectorized call.
    pst <- lapply(stats::setNames(names(pis[[1]]), names(pis[[1]])), function(k) vapply(pis, function(pi) pi[[k]], 0))
    Fz <- function(z) rowMeans(matrix(family_cdf(fam, rep(z, length(pis)), lapply(pst, rep, each = length(z))), length(z)))
    # The CRPS is finite exactly when the predictive mean is.
    mu <- vapply(pis, function(pi) tryCatch(family_moment(fam, "mean", pi, quiet = TRUE), error = function(e) NaN), 0)
    if (any(is.na(mu) & !is.nan(mu)) || any(is.infinite(mu))) return(Inf)
    lo <- min(vapply(pis, function(pi) family_quantile(fam, 1e-9, pi), 0))
    hi <- max(vapply(pis, function(pi) family_quantile(fam, 1 - 1e-9, pi), 0))
    if (fam$discrete) {
      # The CDF is a step function on the lattice: each step's squared
      # difference times the lattice spacing.
      g <- lattice_step(fam)
      if (g == 0) stop_unsupported("CRPS for a discrete law on an irregular support", fam$id)
      k <- seq(min(lo, y[i]), max(hi, y[i]), by = g)
      return(g * sum((Fz(k) - (y[i] <= k + g / 2))^2))
    }
    # Integrate over the whole support: finite cutoffs at extreme quantiles
    # make the quadrature miss the mass of heavy tails.
    sp <- lapply(pis, fam$support)
    a <- min(unlist(lapply(sp, `[[`, "lower")), y[i]); b <- max(unlist(lapply(sp, `[[`, "upper")), y[i])
    quad <- function(f, lo, up) tryCatch(stats::integrate(f, lo, up, rel.tol = 1e-8, subdivisions = 1000L, stop.on.error = FALSE),
      error = function(e) stop_numerical(sprintf("The CRPS integral of observation %d failed (%s).", i, conditionMessage(e))))
    i1 <- quad(function(z) Fz(z)^2, a, y[i])
    i2 <- quad(function(z) (1 - Fz(z))^2, y[i], b)
    # Roundoff messages on an accurate integral are harmless; a negative or
    # non finite value, or a large error estimate, is a failure.
    v <- i1$value + i2$value; err <- i1$abs.error + i2$abs.error
    if (!is.finite(v) || v < 0 || ((i1$message != "OK" || i2$message != "OK") && !(err <= 1e-4 * max(1, v))))
      stop_numerical(sprintf("The CRPS integral of observation %d did not converge (%s).", i,
                             if (i1$message != "OK") i1$message else i2$message))
    v
  }, 0)
}

new_spec_dots <- function(spec) {
  two_sided <- function(f, k) if (length(f) == 3) f else
    stats::as.formula(call("~", as.name(k), f[[2]]), env = environment(f))
  prim <- if (is.null(spec$formula)) primary_par(spec$fam, spec$rpars) else primary_par(spec$fam, regression_pars(spec$fam))
  oth <- setdiff(names(spec$formulas), prim)
  c(unname(Map(two_sided, spec$formulas[oth], oth)), spec$data_args)
}

#' Cross validation of a fitted model
#'
#' Refits the model on training folds and scores the held out
#' observations with a proper scoring rule. With `group` the folds leave out
#' whole groups (leave group out cross validation), and held out groups are
#' predicted with their random effects at zero.
#'
#' @param object A fitted model from [regcomp()].
#' @param K Number of folds.
#' @param rule Scoring rule passed to [scoring_rule()].
#' @param group Optional grouping variable name; folds then hold out whole
#'   groups.
#' @param folds Optional integer vector assigning each analysis row (the
#'   `nobs()` rows used in the fit) to a fold; a vector with one entry per
#'   row of the original data is also accepted and is restricted to the
#'   analysis rows.
#' @param seed Seed for the fold assignment.
#' @return A list with `scores` (per observation), `estimate` (mean score),
#'   `se` and `folds`.
#' @export
#' @examples
#' set.seed(1)
#' d <- data.frame(x = rnorm(120)); d$y <- rpois(120, exp(0.5 + 0.3 * d$x))
#' fit <- regcomp(y ~ x, family = "poisson", data = d)
#' cv_regcomp(fit, K = 5)$estimate
cv_regcomp <- function(object, K = 10, rule = "log", group = NULL, folds = NULL, seed = 1) {
  env <- environment(object$spec$formula %||% object$spec$formulas[[1]])
  # Folds split the stored analysis rows; external weights and data
  # arguments are columns of that data, so they are split with the rows.
  data <- object$model_data
  if (is.null(data)) stop_unsupported("Cross validation of a fit without stored analysis rows (refit with this version)")
  ctrl <- if (is.null(object$call$control)) regcomp_control() else eval(object$call$control, env)
  n <- nrow(data)
  if (!is.null(folds) && length(folds) != n) {
    if (length(folds) >= max(object$rows)) folds <- folds[object$rows]
    else stop_input(sprintf("folds must have %d entries (one per analysis row) or one per row of the data.", n))
  }
  if (is.null(folds)) {
    local_seed(seed)
    if (is.null(group)) folds <- sample(rep_len(seq_len(K), n))
    else { g <- as.character(data[[group]]); ug <- unique(g); fg <- sample(rep_len(seq_len(K), length(ug))); folds <- fg[match(g, ug)] }
  }
  sc <- rep(NA_real_, n)
  for (k in sort(unique(folds))) {
    test <- folds == k
    cl <- object$call; cl$data <- data[!test, , drop = FALSE]; cl$subset <- NULL
    cl$weights <- object$model_args$weights
    for (a in setdiff(names(object$model_args$dots), "")) if (!is.null(cl[[a]])) cl[[a]] <- object$model_args$dots[[a]]
    cl$control <- utils::modifyList(ctrl, list(hessian = FALSE))
    f <- eval(cl, env)
    sc[test] <- scoring_rule(f, data[test, , drop = FALSE], rule = rule, re = if (is.null(group)) "conditional" else "zero")
  }
  list(scores = sc, estimate = mean(sc), se = stats::sd(sc) / sqrt(n), folds = folds, rule = rule)
}

#' Compare fitted models
#'
#' Tabulates log likelihood, number of parameters, AIC and BIC, and
#' optionally a cross validated score, ordered by AIC.
#'
#' @param ... Fitted models, optionally named.
#' @param cv Optional scoring rule name; when given, `cv_regcomp()` with
#'   `K` folds is added.
#' @param K Number of folds for cross validation.
#' @return A data frame.
#' @export
#' @examples
#' set.seed(1)
#' d <- data.frame(x = rnorm(200)); d$y <- rnbinom(200, mu = exp(1 + 0.4 * d$x), size = 2)
#' compare_models(poisson = regcomp(y ~ x, family = "poisson", data = d),
#'                nb2 = regcomp(y ~ x, family = "nb2", data = d))
compare_models <- function(..., cv = NULL, K = 10) {
  fits <- list(...)
  if (is.null(names(fits))) names(fits) <- paste0("model", seq_along(fits))
  # Likelihood criteria compare models only on the same observations of the
  # same response; Bayesian fits have no maximized likelihood.
  if (any(vapply(fits, function(f) identical(f$inference, "bayesian"), TRUE)))
    stop_unsupported("AIC and BIC for Bayesian fits (compare them by held out scores or LOO)")
  y1 <- fits[[1]]$spec$y
  for (f in fits[-1]) if (!identical(dim(f$spec$y), dim(y1)) || NROW(f$spec$y) != NROW(y1) ||
                          !isTRUE(all.equal(unclass(f$spec$y), unclass(y1), check.attributes = FALSE)))
    stop_input("compare_models() needs fits to the same observations of the same response (same rows, same transformation).")
  out <- do.call(rbind, lapply(names(fits), function(nm) { f <- fits[[nm]]
    data.frame(model = nm, family = f$family$id, logLik = f$logLik, df = f$df,
               AIC = -2 * f$logLik + 2 * f$df, BIC = -2 * f$logLik + log(f$nobs) * f$df) }))
  if (!is.null(cv)) out[[paste0("cv_", cv)]] <- vapply(fits, function(f) cv_regcomp(f, K = K, rule = cv)$estimate, 0)
  out <- out[order(out$AIC), ]
  rownames(out) <- NULL
  out
}

#' User defined families
#'
#' Builds a family from a log density (or log mass) so that every engine
#' path that needs only the likelihood (maximum likelihood, random effects,
#' penalties, prediction by numerical inversion, scoring) works with it.
#' The normalization is checked numerically at `check` parameter values.
#'
#' @param id Identifier.
#' @param pars List of parameters from [par_spec()].
#' @param lpdf Function `(y, p)` returning log densities, vectorized over
#'   `y` and the parameter vectors in `p`.
#' @param discrete Whether the law is discrete.
#' @param support Optional function of the parameters returning
#'   `list(lower, upper)`.
#' @param cdf,quantile,rng,mean,variance Optional functions.
#' @param stan Optional Stan template from the registry helpers.
#' @param check Optional named list of parameter values at which total
#'   mass is verified.
#' @param ... Further arguments to the family constructor.
#' @return A `regcomp_family`.
#' @export
#' @examples
#' f <- custom_family("my_exp", list(par_spec("rate", "positive", role = "rate")),
#'   lpdf = function(y, p) dexp(y, p$rate, log = TRUE),
#'   support = function(p) list(lower = 0, upper = Inf), check = list(rate = 2))
#' set.seed(1)
#' d <- data.frame(y = rexp(100, 2))
#' coef(regcomp(y ~ 1, family = f, data = d))
custom_family <- function(id, pars, lpdf, discrete = FALSE, support = NULL, cdf = NULL,
                          quantile = NULL, rng = NULL, mean = NULL, variance = NULL,
                          stan = NULL, check = NULL, ...) {
  if (is.null(names(pars))) names(pars) <- vapply(pars, `[[`, "", "name")
  fam <- family_def(id, id, if (discrete) "count" else "real", pars, lpdf = lpdf, discrete = discrete,
                    support = support, cdf = cdf, quantile = quantile, rng = rng, mean = mean,
                    variance = variance, stan = stan, status = "experimental", ...)
  if (!is.null(check)) {
    d <- do.call(distribution, c(list(fam), check))
    sp <- if (!is.null(support)) support(check) else list(lower = -Inf, upper = Inf)
    tot <- if (discrete) sum(prob(d, seq(max(sp$lower, -1e4), min(sp$upper, 1e5)))) else
      stats::integrate(function(x) prob(d, x), sp$lower, sp$upper, rel.tol = 1e-8)$value
    if (abs(tot - 1) > 1e-4) stop_input(sprintf("The density of '%s' integrates to %.6f, not 1, at the check values.", id, tot))
  }
  fam
}

#' Identifiability diagnostics of a fit
#'
#' Reports the conditioning of the observed information, the parameter
#' combinations that are nearly flat (eigenvectors of the smallest
#' eigenvalues of the information), random effect variances at the
#' boundary, and for mixtures the possibility of label switching.
#'
#' @param object A frequentist fit.
#' @param tol Eigenvalue ratio below which a direction is reported.
#' @return A list with `condition_number`, `weak_directions` (a list of
#'   named loading vectors), `boundary_variances` and `notes`.
#' @export
identifiability_check <- function(object, tol = 1e-8) {
  if (object$inference != "frequentist" || is.null(object$vcov)) stop_unsupported("Identifiability diagnostics without a covariance matrix")
  info <- safe_inverse(object$vcov)
  ev <- eigen((info + t(info)) / 2, symmetric = TRUE)
  rel <- ev$values / max(ev$values)
  weak <- lapply(which(rel < tol | !is.finite(rel)), function(j) { v <- ev$vectors[, j]; names(v) <- rownames(object$vcov)
    sort(round(v[abs(v) > 0.1], 3), decreasing = TRUE) })
  notes <- character()
  if (!is.null(object$family$operator) && identical(object$family$operator$type, "mixture"))
    notes <- c(notes, "Mixture components are exchangeable: labels may switch between fits; order them by a parameter to compare.")
  if (any(!is.finite(ev$values)) || any(ev$values <= 0)) notes <- c(notes, "The information matrix is singular: some parameters are not identified.")
  list(condition_number = max(ev$values) / min(abs(ev$values)), weak_directions = weak,
       boundary_variances = object$diagnostics$boundary_variance, notes = notes)
}
