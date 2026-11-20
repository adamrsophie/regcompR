# Model fitting entry point and the maximum likelihood engine.

#' Fit a distributional regression model
#'
#' Every regression parameter of the response distribution gets its own
#' linear predictor on its link scale. The first formula gives the response
#' and the predictor of the main parameter: the parameter equal to the mean
#' if there is one, else the location, scale or rate parameter, else the
#' first regression parameter (for example `scale` of the Weibull law, so
#' that survival models are accelerated failure time models); further
#' formulas such as `sigma ~ x` or `shape ~ x` give predictors for other
#' parameters, which otherwise are constant. Predictors may contain
#' offsets, factors and interactions, random effects in `lme4` syntax
#' (`(1 | g)`, `(1 + x | g)`, uncorrelated `(1 + x || g)`, nested `(1 | a/b)`,
#' crossed `(1 | a) + (1 | b)`), structured random effects for longitudinal and
#' spatial data written as glmmTMB style wrappers, and smooth terms
#' `s()` and `t2()` from 'mgcv', which are fitted in their mixed model
#' representation (as in `mgcv::gamm()`). Tensor products `te()` and `ti()`
#' have no such representation; use `t2()` instead.
#'
#' @section Random effect structures:
#' A bar term may be wrapped in a covariance structure:
#' * `us(x | g)`: unstructured (the default for `|`); `diag(x | g)` and
#'   `homdiag(x | g)`: independent effects with separate or common variance
#'   (the default for `||`).
#' * `cs(0 + time | g)` and `homcs(0 + time | g)`: compound symmetry with
#'   heterogeneous or common variances.
#' * `ar1(0 + time | g)`: first order autoregressive correlation
#'   `rho^|i - j|` across the columns (time must be a factor with equally
#'   spaced, ordered levels), with a common variance.
#' * `ou(t | g)` or `exp(t | g)`: continuous time autoregressive
#'   (Ornstein-Uhlenbeck) correlation `exp(-|t_i - t_j| / range)` for
#'   unequally spaced numeric times; `gau(x + y | g)` and
#'   `matern32(x + y | g)` give Gaussian and Matern (smoothness 3/2)
#'   correlation in the coordinates on the left of the bar, for spatial
#'   effects. One random effect is created per distinct coordinate.
#'
#' Structured random effects are available in every parameter and for every
#' family, so a latent autoregressive process can enter a count, survival or
#' ordinal model. The Bayesian engine supports `us` and `diag` structures.
#'
#' @section Engines:
#' * `"ml"`: maximum likelihood in R. Without random effects or smooths the
#'   log likelihood is maximized by `nlminb()` with an analytic chain rule
#'   gradient. With random effects or smooths the marginal likelihood is
#'   approximated by the Laplace method (joint over all random effects, so
#'   crossed and nested factors and random effects in several parameters
#'   are allowed); the fitted object records `approximation = "laplace"`.
#'   Penalties (see [penalty()]) give penalized likelihood estimates.
#' * `"tmb"`: the same model, random effects and covariance structures
#'   through automatic differentiation and the Laplace approximation of
#'   'TMB' (via 'RTMB'), for the normal, Poisson, binomial, Bernoulli,
#'   NB1, NB2, Student t, gamma, lognormal, Weibull, beta and exponential
#'   families (see [family_capabilities()]). It gives the same estimates as
#'   `"ml"` and is much faster for large multilevel models; the taped
#'   density is checked against the registry density before fitting.
#' * `"bayes"`: Bayesian inference. A 'Stan' program is generated from the
#'   same specification (see [stan_code()]) and run by 'bayescompR'
#'   (default, no compiler needed), 'cmdstanr' or 'rstan'. Priors are set
#'   with [prior()].
#'
#' @param formula Two sided formula `response ~ predictors` for the first
#'   regression parameter. The response may be wrapped in [cens()] or
#'   [surv()]. A nonlinear predictor is written `response ~ nl(expr)`, as in
#'   `y ~ nl(a * exp(-b * x))`: the names of `expr` that have their own
#'   formula in `...` (`a ~ 1 + (1 | g)`, `b ~ z`) are nonlinear parameters,
#'   each with a linear predictor that may hold covariates and random
#'   effects; the other names are covariates. `expr` gives the first
#'   parameter on its link scale, as in `nls()`, `nlme()` and `nlmer()`.
#'   Nonlinear parameters need starting values (`control =
#'   regcomp_control(start = list(a = 10, b = 0.5))`, also the initial
#'   values of the Stan chains) and are fitted by `engine = "ml"` or by
#'   `engine = "bayes"` with a Stan runner.
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
  dots <- eval(substitute(alist(...)))
  if (missing(data)) data <- env_data(formula, dots, parent.frame())
  if (!is.null(penalty) && engine != "ml") stop_input("penalty is only used by engine = \"ml\".")
  if (!is.null(priors) && engine != "bayes") stop_input("priors are only used by engine = \"bayes\".")
  fam <- model_family(family, parameterization, formula, data)
  if (isTRUE(fam$extra$no_regression))
    stop_unsupported("Regression", fam$id)
  data <- as.data.frame(data)
  mr <- model_rows(formula, dots, data, substitute(weights), parent.frame(), subset = subset)
  # A shrinkage penalty identifies aliased or p > n designs (lasso
  # coefficients may then be non unique while predictions are unique).
  shrink <- engine == "ml" && inherits(penalty, "regcomp_penalty") && isTRUE(penalty$lambda > 0) &&
    penalty$type %in% c("ridge", "lasso", "elastic_net", "adaptive_lasso", "group_lasso", "scad", "mcp")
  spec <- build_spec(formula, mr$dots, fam, mr$data, weights = mr$weights,
                     links = links, fixed = fixed, env = parent.frame(), allow_aliased = shrink)
  if (!is.null(spec$nl) && !is.null(penalty)) stop_unsupported("A penalty with a nonlinear formula")
  if (!is.null(spec$nl) && engine == "tmb") stop_unsupported("A nonlinear formula with engine \"tmb\" (use engine = \"ml\" or \"bayes\")")
  spec$data_args <- dots[setdiff(names(spec$data_args), "")]
  spec$call <- cl
  fit <- switch(engine,
    ml = fit_ml(spec, control, penalty),
    tmb = fit_tmb(spec, control),
    bayes = fit_bayes(spec, priors, control))
  fit$call <- cl
  fit$engine <- engine
  # The analysis rows: their indices in the data, and the data with every
  # external vector stored as a column, so cross validation and clustered
  # covariances use exactly the fitted observations.
  fit$rows <- mr$rows
  fit$model_data <- mr$data
  fit$model_args <- list(dots = mr$dots, weights = mr$weights)
  fit
}

# Family of a model; families that depend on the number of response
# categories (multinomial logit) are expanded here.
model_family <- function(family, parameterization, formula, data) {
  fam <- get_family(family, parameterization)
  if (is.null(fam$extra$expand)) return(fam)
  y0 <- eval(formula[[2]], data, environment(formula))
  fam$extra$expand(if (is.matrix(y0)) ncol(y0) else if (is.factor(y0)) nlevels(y0) else max(y0, na.rm = TRUE))
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
#'   settings: chains, iterations, seed, the runner (`"bayescompR"`,
#'   `"cmdstanr"` or `"rstan"` for Stan, `"jags"` through 'R2jags' or
#'   `"openbugs"` through 'R2OpenBUGS') and parallel chains.
#' @param adapt_delta Target acceptance for NUTS.
#' @param init Initial values for the Stan runners: a radius for uniform
#'   initial values on the unconstrained scale (the default is 2), or a list
#'   with one named list per chain using the parameter names of the
#'   generated program (`b_mu`, `b_mu_1`, ...). Useful for multimodal
#'   posteriors such as mixtures, where a chain started in another mode
#'   stays there.
#' @param start Optional named list of starting values for the coefficients
#'   of each parameter (vectors on the link scale).
#' @param prox_maxit,prox_tol Iteration limit and relative step tolerance of
#'   the proximal gradient (FISTA) solver used for lasso type penalties
#'   (`maxit` and `rel_tol` apply to `nlminb()` only).
#' @return A list.
#' @export
#' @examples
#' regcomp_control(maxit = 500)
regcomp_control <- function(maxit = 1000, rel_tol = 1e-10, inner_maxit = 50,
                            inner_tol = 1e-8, restarts = 2, hessian = TRUE,
                            trace = FALSE, chains = 4, iter_warmup = 1000,
                            iter_sampling = 1000, seed = 1,
                            runner = c("bayescompR", "cmdstanr", "rstan", "jags", "openbugs"),
                            cores = 1, adapt_delta = 0.8, start = NULL,
                            prox_maxit = 5000, prox_tol = 1e-8, init = NULL) {
  list(maxit = maxit, rel_tol = rel_tol, inner_maxit = inner_maxit,
       inner_tol = inner_tol, restarts = restarts, hessian = hessian,
       trace = trace, chains = chains, iter_warmup = iter_warmup,
       iter_sampling = iter_sampling, seed = seed, runner = match.arg(runner),
       cores = cores, adapt_delta = adapt_delta, start = start,
       prox_maxit = prox_maxit, prox_tol = prox_tol, init = init)
}

# Rows used in a fit. External vectors (weights and data arguments given
# outside `data`) must have one value or one per row of `data`; they are
# stored as reserved columns before any row is dropped, so they stay
# aligned with the observations. Then `subset` is applied and rows with a
# missing or non finite value in anything the model evaluates (response,
# covariates after transformation and dot expansion, data arguments,
# weights) or with zero weight are omitted, as na.omit does for glm.
model_rows <- function(formula, dots, data, weights, env, subset = NULL) {
  n <- nrow(data)
  w <- if (is.null(weights)) NULL else eval(weights, data, env)
  if (!is.null(w)) {
    if (!is.numeric(w) || !length(w) %in% c(1L, n) || any(!is.finite(w) | w < 0))
      stop_input(sprintf("weights must be finite, nonnegative and have one value or %d values (one per row of data).", n))
    data[["(weights)"]] <- rep_len(w, n)
    # Weights given as a column expression of data keep it (so a dot in the
    # formula leaves that column out); anything else uses the stored column.
    wname <- if (length(all.vars(weights)) && !length(setdiff(all.vars(weights), names(data))) && length(w) == n) weights
             else as.name("(weights)")
  }
  is_formula <- vapply(dots, function(d) inherits(tryCatch(eval(d, env), error = function(e) NULL), "formula"), TRUE)
  for (k in names(dots)[!is_formula & nzchar(names(dots))]) {
    v <- tryCatch(eval(dots[[k]], data, env), error = function(e) NULL)
    if (is.null(v)) next
    if (!length(v) %in% c(1L, n))
      stop_input(sprintf("The data argument '%s' must have one value or %d values (one per row of data); got %d.", k, n, length(v)))
    if (!length(setdiff(all.vars(dots[[k]]), names(data)))) next
    if (length(v) == n && n > 1) { data[[paste0("(", k, ")")]] <- v; dots[[k]] <- as.name(paste0("(", k, ")")) }
  }
  rows <- seq_len(n)
  if (!is.null(subset)) {
    rows <- rows[subset]
    if (anyNA(rows)) stop_input("subset selects rows that are not in data.")
    data <- data[rows, , drop = FALSE]; n <- nrow(data)
  }
  forms <- c(list(formula), lapply(dots[is_formula], eval, env))
  bad <- rep(FALSE, n)
  for (f in forms) bad <- bad | nonfinite_rows(f, data, env)
  for (k in names(dots)[!is_formula & nzchar(names(dots))]) {
    v <- tryCatch(eval(dots[[k]], data, env), error = function(e) NULL)
    if (length(v) == n) bad <- bad | is.na(v)
  }
  keep <- !bad
  if (!is.null(w)) keep <- keep & data[["(weights)"]] > 0
  if (!any(keep)) stop_input("No rows remain after removing missing values and zero weights.")
  list(data = data[keep, , drop = FALSE], dots = dots, weights = if (!is.null(w)) wname, rows = rows[keep])
}

# Rows where a variable or term of a formula is missing or not finite after
# evaluation: the response, each fixed effect term (a dot expands to the
# other columns), the variables of random effect and smooth terms.
nonfinite_rows <- function(f, data, env) {
  n <- nrow(data); bad <- rep(FALSE, n)
  mark <- function(v) {
    if (is.null(v) || NROW(v) != n) return()
    b <- if (is.numeric(v)) !is.finite(v) else is.na(v)
    if (is.matrix(b)) b <- if (inherits(v, "regcomp_cens")) rowSums(is.na(unclass(v)[, c("lower", "upper"), drop = FALSE])) > 0 else rowSums(b) > 0
    bad <<- bad | b
  }
  ev <- function(e) tryCatch(suppressWarnings(eval(e, data, env)), error = function(e) NULL)
  if (length(f) == 3) {
    y <- ev(f[[2]])
    if (is.numeric(y) && is.matrix(y) && !inherits(y, "regcomp_cens")) mark(rowSums(is.na(y)) > 0) else mark(y)
  }
  rhs <- f[[length(f)]]
  walk <- function(e) {
    if (is.name(e)) {
      nm <- as.character(e)
      if (nm == ".") { for (cn in setdiff(names(data), all.vars(f[[2]]))) if (!startsWith(cn, "(")) mark(data[[cn]]) }
      else mark(ev(e))
      return()
    }
    if (!is.call(e)) return()
    fn <- as.character(e[[1]])[1]
    if (fn %in% c("+", "-", "*", ":", "/", "(", "^", "%in%")) { for (a in as.list(e)[-1]) walk(a); return() }
    if (fn %in% c("|", "||", "nl") || is_smooth_call(e) || is_struct_call(e)) {
      # Random effect and smooth terms: their variables must be present.
      for (v in all.vars(e)) if (v %in% names(data)) mark(data[[v]])
      return()
    }
    mark(ev(e))
  }
  walk(rhs)
  bad
}

# Data frame for a model fitted without `data`: the formula variables from
# its environment.
env_data <- function(formula, dots, env) {
  fe <- environment(formula) %||% env
  vars <- unique(c(all.vars(formula), unlist(lapply(dots, all.vars))))
  vals <- lapply(vars, function(v) tryCatch(eval(as.name(v), fe), error = function(e) NULL))
  keep <- vapply(vals, function(v) is.atomic(v) && length(v) > 1 || is.factor(v) || is.matrix(v), TRUE)
  as.data.frame(stats::setNames(vals[keep], vars[keep]), optional = TRUE, stringsAsFactors = FALSE)
}

# Assemble everything the engines need.
build_spec <- function(formula, dots, fam, data, weights = NULL, links = NULL,
                       fixed = NULL, env = parent.frame(), xlev = NULL,
                       prev = NULL, response = is.null(prev), allow_aliased = FALSE) {
  is_formula <- vapply(dots, function(d) {
    v <- tryCatch(eval(d, env), error = function(e) NULL)
    inherits(v, "formula")
  }, TRUE)
  forms <- lapply(dots[is_formula], eval, env)
  if (length(forms)) names(forms) <- NULL
  dargs <- dots[!is_formula]
  fl <- parse_formulas(formula, forms, fam)
  nl <- attr(fl, "nl")
  if (!is.null(nl)) {
    dv <- intersect(nl$pars, names(data))
    if (length(dv)) stop_input(sprintf("Nonlinear parameter names must differ from data columns: rename %s.", paste(dv, collapse = ", ")))
    vs <- setdiff(all.vars(nl$expr), nl$pars)
    nl$data <- stats::setNames(lapply(vs, function(v) eval(as.name(v), data, nl$env)), vs)
    bad <- vs[!lengths(nl$data) %in% c(1L, nrow(data))]
    if (length(bad)) stop_input(sprintf("The covariate%s %s of the nonlinear formula must have one value or %d values (one per row of data).",
                                        if (length(bad) > 1) "s" else "", paste(bad, collapse = ", "), nrow(data)))
  }
  for (k in setdiff(fam$extra$fixed_unless_formula, attr(fl, "given")))
    if (is.null(fixed[[k]])) fixed[[k]] <- fam$pars[[k]]$default
  # Parameters that take no formula, no data and are not estimated keep
  # their family default.
  for (k in setdiff(names(fam$pars), c(names(fl), nl$par, names(fixed), data_pars(fam)))) {
    pk <- fam$pars[[k]]
    if (is.null(pk$estimate) && !is.null(pk$default)) fixed[[k]] <- pk$default
  }
  for (k in names(fixed)) {
    if (!k %in% names(fam$pars)) stop_input(sprintf("Unknown fixed parameter '%s'.", k))
    fl[[k]] <- NULL
  }
  rpars <- names(fl)
  n <- nrow(data)
  y <- if (response) get_response(formula, data) else NULL
  ylev <- if (!is.null(prev)) prev$ylev else if (is.factor(y)) levels(y)
  if (is.factor(y)) y <- factor_response(y, fam, ylev, !is.null(prev))
  cmat <- NULL
  if (inherits(y, "regcomp_cens")) { cmat <- unclass(y); y <- cmat[, "lower"] }
  if (is.null(prev) && is.null(cmat) && fam$variate_type == "scalar") check_response(y, fam)
  # Links.
  lk <- list()
  for (k in nl$pars) lk[[k]] <- link("identity")
  for (k in c(setdiff(rpars, nl$pars), nl$par)) {
    nm <- if (!is.null(links[[k]])) links[[k]] else fam$pars[[k]]$link
    if (!nm %in% fam$pars[[k]]$links)
      stop_input(sprintf("Link '%s' is not allowed for parameter '%s' (%s).", nm, k,
                         paste(fam$pars[[k]]$links, collapse = ", ")))
    lk[[k]] <- link(nm)
  }
  preds <- list()
  for (k in rpars) {
    fk <- if (!is.null(prev$preds[[k]]$dot_formula)) prev$preds[[k]]$dot_formula else
      expand_dot(rhs_of(fl[[k]]), data, c(all.vars(formula[[2]]), unlist(lapply(dargs, all.vars)), all.vars(weights)))
    preds[[k]] <- build_predictor(fk, data, k,
      xlev = if (!is.null(prev)) prev$preds[[k]]$xlev,
      prev_pred = if (!is.null(prev)) prev$preds[[k]],
      smooth_specs = if (!is.null(prev)) prev$preds[[k]]$smooths,
      re_levels = if (!is.null(prev)) { pb <- prev$preds[[k]]$blocks; lab <- vapply(pb, `[[`, "", "label")
        c(stats::setNames(lapply(pb, `[[`, "levels"), lab), stats::setNames(lapply(pb, `[[`, "coords"), paste0(".coords.", lab))) })
    preds[[k]]$dot_formula <- fk
    if (length(preds[[k]]$aliased) && !allow_aliased)
      stop_input(sprintf("The design for '%s' is rank deficient; aliased columns: %s. Remove them from the formula, or use a shrinkage penalty (ridge, lasso and similar) with the ml engine.",
                         k, paste(preds[[k]]$aliased, collapse = ", ")))
    if (any(preds[[k]]$offset != 0) && isFALSE(fam$pars[[k]]$offset_allowed) &&
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
  fixed_data <- fixed[intersect(names(fixed), data_pars(fam))]
  fixed <- fixed[setdiff(names(fixed), data_pars(fam))]
  w <- if (is.null(weights)) rep(1, n) else {
    wv <- eval(weights, data, env)
    if (is.null(wv)) rep(1, n) else as.numeric(wv)
  }
  for (k in intersect(fam$extra$no_intercept, rpars)) {
    X <- preds[[k]]$X; keep <- colnames(X) != "(Intercept)"
    preds[[k]]$X <- X[, keep, drop = FALSE]
  }
  # Coefficient layout.
  beta_idx <- list(); off <- 0
  for (k in rpars) {
    p <- ncol(preds[[k]]$X)
    beta_idx[[k]] <- off + seq_len(p); off <- off + p
  }
  beta_names <- unlist(lapply(rpars, function(k) { cn <- colnames(preds[[k]]$X)
    if (length(cn)) paste0(k, ":", cn) else character() }))
  # Estimated vector parameters that are not regression parameters (ordinal
  # thresholds and similar); their size comes from the response.
  aux <- if (!is.null(prev)) prev$aux else {
    a <- list()
    for (k in names(fam$pars)) { sp <- fam$pars[[k]]
      if (is.null(sp$estimate) || k %in% names(fixed)) next
      K <- if (!is.null(ylev)) length(ylev) else max(y, na.rm = TRUE)
      a[[k]] <- list(type = sp$estimate, size = K - 1L,
                     labels = if (!is.null(ylev)) paste0(ylev[-K], "|", ylev[-1]) else paste0(seq_len(K - 1), "|", seq_len(K - 1) + 1))
    }
    a
  }
  for (k in names(aux)) {
    beta_idx[[k]] <- off + seq_len(aux[[k]]$size); off <- off + aux[[k]]$size
    beta_names <- c(beta_names, paste0(k, ":", aux[[k]]$labels))
  }
  # Random effect layout.
  blocks <- unlist(lapply(rpars, function(k) preds[[k]]$blocks), recursive = FALSE)
  Z <- list(); theta_idx <- list(); u_idx <- list(); toff <- 0; uoff <- 0
  Q <- sum(vapply(blocks, function(b) b$q * b$L, 1))
  for (j in seq_along(blocks)) {
    b <- blocks[[j]]
    nt <- block_ntheta(b)
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
  if (is.null(prev) && length(fixed)) {
    # Fixed values must satisfy the parameter domains and the family constraint.
    fx <- fixed[!vapply(fixed, is.null, TRUE)]
    for (k in names(fx)) if (!is.null(fam$pars[[k]]) && !all(check_par_domain(fam$pars[[k]], fx[[k]]), na.rm = TRUE))
      stop_domain(sprintf("Fixed parameter '%s' of '%s' is outside its %s domain.", k, fam$id, fam$pars[[k]]$domain))
  }
  spec <- structure(list(
    fam = fam, formulas = fl, rpars = rpars, links = lk, preds = preds,
    y = y, cmat = cmat, n = n, weights = w, data_pars = dp, fixed = fixed,
    beta_idx = beta_idx, beta_names = beta_names, blocks = blocks,
    theta_idx = theta_idx, u_idx = u_idx, Z = Z, Q = Q, ntheta = toff,
    data_args = dargs, formula = formula, fixed_data = fixed_data, response_name = paste(deparse(formula[[2]]), collapse = ""),
    aux = aux, ylev = ylev, nl = nl
  ), class = "regcomp_spec")
  if (is.null(prev) && length(fixed) && !is.null(fam$check) && length(rpars)) {
    # The family constraint with the fixed values and the starting values of
    # the free parameters (the check needs every parameter).
    pr <- tryCatch(eta_to_pars(spec, compute_eta(spec, start_values(spec, regcomp_control()))),
                   error = function(e) NULL)
    msg <- if (!is.null(pr)) tryCatch(fam$check(pr), error = function(e) NULL)
    if (is.character(msg) && length(msg) && nzchar(msg[1]))
      stop_domain(paste("With the fixed parameter values:", msg[1]))
  }
  # Dense copies of small random effect designs (faster Laplace iterations).
  if (isTRUE(spec$Q > 0 && spec$Q <= 300 && spec$Q * spec$n <= 5e6))
    spec$Zd <- lapply(spec$Z, function(z) if (!is.null(z)) as.matrix(z))
  spec
}

# Number of covariance parameters of a block, starting values and names.
block_ntheta <- function(b) switch(b$struct, us = b$q + b$q * (b$q - 1) / 2, diag = b$q,
  homdiag = 1, cs = b$q + 1, homcs = 2, 2)
block_theta_start <- function(b) {
  rng <- if (!is.null(b$coords) && nrow(b$coords) > 1) log(stats::median(stats::dist(b$coords))) else 0
  switch(b$struct, us = c(rep(log(0.5), b$q), rep(0, b$q * (b$q - 1) / 2)), diag = rep(log(0.5), b$q),
    homdiag = log(0.5), cs = c(rep(log(0.5), b$q), 0), homcs = c(log(0.5), 0), ar1 = c(log(0.5), 0),
    c(log(0.5), rng))
}
block_theta_names <- function(b) {
  pre <- paste0(b$par, ":", b$label, ":")
  sdn <- paste0(pre, "log_sd(", b$terms, ")")
  switch(b$struct,
    us = { idx <- which(lower.tri(diag(b$q)), arr.ind = TRUE)
           c(sdn, paste0(pre, "chol(", b$terms[idx[, 1]], ",", b$terms[idx[, 2]], ")")) },
    diag = sdn, homdiag = paste0(pre, "log_sd"), cs = c(sdn, paste0(pre, "rho_raw")),
    homcs = paste0(pre, c("log_sd", "rho_raw")), ar1 = paste0(pre, c("log_sd", "atanh_rho")),
    paste0(pre, c("log_sd", "log_range")))
}

# Matrices that receive taped values must be taped themselves (TMB engine).
ad_like <- function(M, th) if (inherits(th, "advector")) RTMB::AD(M) else M

# Covariance matrix of one block from its unconstrained parameters.
block_cov <- function(b, th) {
  q <- b$q
  if (is.null(b$struct)) b$struct <- if (b$corr) "us" else "diag"
  switch(b$struct,
    diag = return(diag(exp(2 * th[seq_len(q)]), q)),
    homdiag = return(diag(exp(2 * th[1]), q)),
    us = {
      sd <- exp(th[seq_len(q)])
      L <- ad_like(diag(q), th)
      L[lower.tri(L)] <- th[-seq_len(q)]
      L <- L / sqrt(drop((L * L) %*% rep(1, q)))
      D <- ad_like(diag(q), th); diag(D) <- sd
      return(D %*% L %*% t(L) %*% D)
    },
    cs = , homcs = {
      sd <- if (b$struct == "cs") exp(th[seq_len(q)]) else rep(exp(th[1]), q)
      a <- 1 / (q - 1); rho <- (1 + a) / (1 + exp(-th[length(th)])) - a   # plogis, written so RTMB can tape it
      R <- ad_like(matrix(0, q, q), th) + rho; diag(R) <- 1
      D <- ad_like(diag(q), th); diag(D) <- sd
      return(D %*% R %*% D)
    },
    ar1 = {
      R <- tanh(th[2])^abs(outer(seq_len(q), seq_len(q), `-`))
      return(exp(2 * th[1]) * R)
    })
  d <- as.matrix(stats::dist(b$coords)) / exp(th[2])
  R <- switch(b$struct, ou = , exp = exp(-d), gau = exp(-d^2),
              matern32 = (1 + sqrt(3) * d) * exp(-sqrt(3) * d))
  exp(2 * th[1]) * (R + diag(1e-10, q))
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
  for (k in names(spec$aux)) for (j in spec$beta_idx[[k]]) {
    h <- 1e-6 * max(1, abs(beta[j])); bp <- beta; bm <- beta
    bp[j] <- bp[j] + h; bm[j] <- bm[j] - h
    out[j] <- (ml_objective(spec, bp, u) - ml_objective(spec, bm, u)) / (2 * h)
  }
  out
}

# Inner Laplace problem: mode of the random effects for given beta, theta.
laplace_inner <- function(spec, beta, theta, u0, control, retry = TRUE) {
  Lam <- re_precision(spec, theta)
  u <- if (is.null(u0)) rep(0, spec$Q) else u0
  rk <- spec$rpars[!vapply(spec$rpars, function(k) is.null(spec$Z[[k]]), TRUE)]
  # Few random effects: dense base matrices avoid the sparse S4 overhead.
  dense <- spec$Q <= 300
  if (dense) { Lam <- as.matrix(Lam); Zs <- if (!is.null(spec$Zd)) spec$Zd[rk] else lapply(spec$Z[rk], as.matrix) } else Zs <- spec$Z[rk]
  fval <- function(u) {
    eta <- compute_eta(spec, beta, u)
    sum(spec$weights * loglik_eta(spec, eta)) - 0.5 * sum(u * as.vector(Lam %*% u))
  }
  # Gradient and negative Hessian in u; derivatives are needed only for the
  # parameters that carry random effects.
  derivs <- function(u) {
    d <- eta_derivs(spec, compute_eta(spec, beta, u), second = TRUE, which = rk)
    grad <- -as.vector(Lam %*% u); H <- Lam
    for (a in seq_along(rk)) {
      Za <- Zs[[a]]
      grad <- grad + as.vector(crossprod(Za, spec$weights * d$g[, a]))
      for (b in seq_along(rk))
        H <- H - crossprod(Za, Zs[[b]] * (spec$weights * d$W[, a, b]))
    }
    list(grad = grad, negH = if (dense) (H + t(H)) / 2 else Matrix::forceSymmetric(H))
  }
  chol1 <- function(M) {
    fail <- function(e) NULL
    if (dense) tryCatch(chol(M), error = fail, warning = fail)
    else tryCatch(Matrix::Cholesky(M, perm = TRUE, LDL = FALSE), error = fail, warning = fail)
  }
  factor <- function(negH) {
    ch <- chol1(negH)
    lam <- 0
    while (is.null(ch) && lam <= 1e8) {
      lam <- if (lam == 0) 1e-6 else lam * 10
      ch <- chol1(negH + lam * (if (dense) diag(spec$Q) else Matrix::Diagonal(spec$Q)))
    }
    list(ch = ch, lam = lam)
  }
  csolve <- function(ch, g) if (dense) backsolve(ch, backsolve(ch, g, transpose = TRUE)) else as.vector(Matrix::solve(ch, g))
  f <- fval(u)
  for (it in seq_len(control$inner_maxit)) {
    dd <- derivs(u)
    if (max(abs(dd$grad)) < control$inner_tol) break
    fc <- factor(dd$negH)
    step <- csolve(fc$ch, dd$grad)
    if (max(abs(step)) < 1e-9) break
    t <- 1
    repeat {
      un <- u + t * step
      fn <- fval(un)
      if (is.finite(fn) && fn >= f - 1e-10) break
      t <- t / 2
      if (t < 1e-8) { un <- u; fn <- f; break }
    }
    moved <- max(abs(un - u))
    u <- un; f <- fn
    if (moved == 0) break
  }
  # Curvature at the accepted point. The Laplace approximation needs a
  # maximum: a positive definite negative Hessian without damping and a
  # small score. Otherwise retry once from zero, then give up.
  dd <- derivs(u)
  fc <- factor(dd$negH)
  ok <- fc$lam == 0 && !is.null(fc$ch) && is.finite(f) &&
    max(abs(dd$grad)) < sqrt(control$inner_tol) * (1 + max(abs(u)))
  if (!ok && retry && !is.null(u0) && any(u0 != 0))
    return(laplace_inner(spec, beta, theta, NULL, control, retry = FALSE))
  ld <- if (!ok) NA_real_ else if (dense) 2 * sum(log(diag(fc$ch))) else
    2 * as.numeric(Matrix::determinant(fc$ch, logarithm = TRUE, sqrt = TRUE)$modulus)
  list(u = u, f = f, logdet = ld, negH = dd$negH, converged = ok)
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
  for (k in spec$nl$pars) {
    v <- control$start[[k]]
    if (is.null(v)) stop_input(sprintf("Nonlinear parameters need starting values, as in regcomp_control(start = list(%s)).",
                                       paste0(spec$nl$pars, " = 1", collapse = ", ")))
    X <- spec$preds[[k]]$X; idx <- spec$beta_idx[[k]]
    if (length(v) == length(idx)) beta[idx] <- v else beta[idx[colnames(X) == "(Intercept)"]] <- v[1]
  }
  for (k in setdiff(spec$rpars, spec$nl$pars)) {
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
  for (k in names(spec$aux)) {
    K <- spec$aux[[k]]$size + 1
    cp <- cumsum(tabulate(y[ok], K))[-K] / sum(ok)
    v <- stats::qlogis(pmin(pmax(cp, 0.02), 0.98))
    if (spec$aux[[k]]$type == "ordered") v <- sort(v + seq_along(v) * 1e-3)
    beta[spec$beta_idx[[k]]] <- aux_untransform(spec$aux[[k]]$type, v)
    if (!is.null(control$start[[k]])) beta[spec$beta_idx[[k]]] <- control$start[[k]]
  }
  beta
}

# Estimated vector parameters are optimized on their natural scale so that
# coefficients and covariances need no back transformation; ordering of
# thresholds is enforced by the likelihood (out of order values give -Inf).
aux_transform <- function(type, b) b
aux_untransform <- function(type, v) v

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
  bs <- coef_scale(spec, beta0)
  if (!has_re) {
    obj <- function(b) {
      v <- ml_objective(spec, b) - pen_val(b)
      if (!is.finite(v)) 1e300 else -v
    }
    gr <- function(b) -(ml_gradient(spec, b) - pen_grad(b))
    opt <- robust_nlminb(beta0, obj, gr, control, bs)
    beta <- opt$par
    if (!is.null(pen) && pen$type == "flic") beta <- flic_adjust(spec, beta, control)
    theta <- numeric(0); u <- numeric(0)
    cov <- NULL
    if (control$hessian && !is.null(pen) && pen$type %in% c("firth", "flic")) {
      cov <- safe_inverse(info_matrix(spec, beta))
    } else if (control$hessian) {
      H <- num_jacobian(function(b) -gr(b), beta, bs)
      H <- (H + t(H)) / 2
      cov <- safe_inverse(-H, bs)
    }
    ll <- ml_objective(spec, beta)
    approx <- "none"
  } else {
    np <- length(beta0)
    # Stage 1: fixed effects only for good starting values.
    obj0 <- function(b) { v <- ml_objective(spec, b) - pen_val(b); if (!is.finite(v)) 1e300 else -v }
    gr0 <- function(b) -(ml_gradient(spec, b) - pen_grad(b))
    b1 <- tryCatch(robust_nlminb(beta0, obj0, gr0, control, bs)$par, error = function(e) beta0)
    theta0 <- unlist(lapply(spec$blocks, block_theta_start))
    state <- new.env(); state$u <- rep(0, spec$Q)
    obj_re <- function(par) {
      # A singular random effect covariance (a correlation at its boundary)
      # is outside the parameter space.
      v <- tryCatch(laplace_objective(spec, par[seq_len(np)], par[-seq_len(np)], state, control),
                    error = function(e) -Inf) - pen_val(par[seq_len(np)])
      if (!is.finite(v)) 1e300 else -v
    }
    ps <- c(bs, rep(1, length(theta0)))
    opt <- robust_nlminb(c(b1, theta0), obj_re, NULL, control, ps)
    beta <- opt$par[seq_len(np)]; theta <- opt$par[-seq_len(np)]
    inner <- laplace_inner(spec, beta, theta, state$u, control)
    if (!isTRUE(inner$converged))
      warning("The random effects mode at the estimate is not a well defined maximum; ",
              "the Laplace approximation is unreliable.", call. = FALSE)
    u <- inner$u
    cov <- NULL
    if (control$hessian) {
      H <- num_hessian(function(p) -obj_re(p), opt$par, s = ps)
      cov <- safe_inverse(-H, ps)
    }
    ll <- -opt$objective + pen_val(beta)
    approx <- "laplace"
  }
  new_fit(spec, beta, theta, u, cov, ll, opt, approx, pen)
}

# Typical size of each coefficient: one, or the size of a nonlinear
# parameter's coefficients (a rate of 0.001 needs finite difference steps
# and optimizer scaling of that size).
coef_scale <- function(spec, beta) {
  bs <- rep(1, length(beta))
  for (k in spec$nl$pars) { i <- spec$beta_idx[[k]]; bs[i] <- max(abs(beta[i])) }
  bs[bs == 0] <- 1
  bs
}

robust_nlminb <- function(start, obj, gr, control, scale = 1) {
  best <- NULL
  starts <- list(start)
  if (control$restarts > 0) {
    seed <- if (exists(".Random.seed", globalenv())) get(".Random.seed", globalenv())
    on.exit(if (is.null(seed)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", seed, globalenv()))
    set.seed(control$seed)
    for (r in seq_len(control$restarts)) starts[[r + 1]] <- start + stats::rnorm(length(start), 0, 0.1 * r) * scale
  }
  # Keep the best point evaluated, so that an error late in a run (a
  # singular matrix near a boundary, say) returns that point flagged as not
  # converged instead of discarding the run.
  seen <- new.env(); seen$val <- Inf
  obj2 <- function(p) { v <- obj(p); if (is.finite(v) && v < seen$val) { seen$val <- v; seen$par <- p }; v }
  for (s in starts) {
    o <- tryCatch(stats::nlminb(s, obj2, gr, scale = 1 / scale, control = list(
      iter.max = control$maxit, eval.max = 2 * control$maxit, rel.tol = control$rel_tol,
      trace = as.integer(control$trace))), error = function(e) if (is.finite(seen$val))
        list(par = seen$par, objective = seen$val, convergence = 1L, iterations = NA_integer_,
             message = paste("stopped after an error:", conditionMessage(e))))
    if (is.null(o)) next
    if (is.null(best) || o$objective < best$objective) best <- o
    if (o$convergence == 0 && o$objective < 1e299) break
  }
  if (is.null(best) || best$objective >= 1e299)
    stop_numerical("The optimizer could not find a finite log likelihood from the starting values.")
  # "False convergence" with numerical gradients often stops at the optimum:
  # restart once; no further decrease confirms the solution.
  if (best$convergence != 0 && grepl("false convergence", best$message)) {
    o <- tryCatch(stats::nlminb(best$par, obj, gr, scale = 1 / scale, control = list(iter.max = control$maxit,
      eval.max = 2 * control$maxit, rel.tol = control$rel_tol)), error = function(e) NULL)
    if (!is.null(o)) {
      if (abs(best$objective - o$objective) <= 1e-6 * (1 + abs(best$objective))) {
        if (o$objective < best$objective) best <- o
        best$convergence <- 0L; best$message <- paste(best$message, "(confirmed by a restart)")
      } else if (o$objective < best$objective) best <- o
    }
  }
  best
}

# Forward difference Jacobian of a vector function.
num_jacobian <- function(f, x, s = 1) {
  f0 <- f(x)
  J <- matrix(0, length(f0), length(x))
  for (j in seq_along(x)) {
    h <- 1e-5 * max(rep_len(s, length(x))[j], abs(x[j]))
    xp <- x; xp[j] <- xp[j] + h
    xm <- x; xm[j] <- xm[j] - h
    J[, j] <- (f(xp) - f(xm)) / (2 * h)
  }
  J
}

# Central difference Hessian of a scalar function with one Richardson
# extrapolation step, which tolerates the small noise left by the inner
# Laplace solver (steps much below 1e-3 amplify that noise).
num_hessian <- function(f, x, rel = 1e-2, s = 1) {
  d <- length(x)
  f0 <- f(x)
  hess_h <- function(h) {
    H <- matrix(0, d, d)
    for (i in seq_len(d)) for (j in i:d) {
      if (i == j) {
        xp <- x; xp[i] <- xp[i] + h[i]; xm <- x; xm[i] <- xm[i] - h[i]
        H[i, i] <- (f(xp) - 2 * f0 + f(xm)) / h[i]^2
      } else {
        e <- function(si, sj) { z <- x; z[i] <- z[i] + si * h[i]; z[j] <- z[j] + sj * h[j]; f(z) }
        H[i, j] <- H[j, i] <- (e(1, 1) - e(1, -1) - e(-1, 1) + e(-1, -1)) / (4 * h[i] * h[j])
      }
    }
    H
  }
  h <- rel * pmax(s, abs(x))
  (4 * hess_h(h / 2) - hess_h(h)) / 3
}

# Inverse of an information matrix, which must be positive definite (to a
# relative tolerance); otherwise an NA matrix flagged as singular.
safe_inverse <- function(A, s = 1) {
  # Inverted on the scale of the coefficients' typical sizes s, so that
  # parameters of very different sizes are not mistaken for singularity.
  if (length(s) > 1) { S <- outer(s, s); out <- safe_inverse(A * S); return(if (isTRUE(attr(out, "singular"))) out else out * S) }
  A <- (A + t(A)) / 2
  ev <- tryCatch(eigen(A, symmetric = TRUE), error = function(e) NULL)
  if (is.null(ev) || any(!is.finite(ev$values)) || min(ev$values) <= 1e-10 * max(abs(ev$values), 1e-300))
    return(structure(matrix(NA_real_, nrow(A), ncol(A)), singular = TRUE))
  ev$vectors %*% diag(1 / ev$values, length(ev$values)) %*% t(ev$vectors)
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
    vcov = cov, logLik = ll, df = edf, nobs = obs_count(spec),
    optimizer = list(convergence = opt$convergence, message = opt$message,
                     iterations = opt$iterations, objective = opt$objective, kkt = opt$kkt,
                     controls = opt$controls),
    approximation = approx, penalty = pen, engine = engine,
    inference = "frequentist"
  ), class = "regcomp_fit")
  fit$diagnostics <- fit_diagnostics(fit)
  fit
}

theta_names <- function(spec) unlist(lapply(spec$blocks, block_theta_names))

fit_diagnostics <- function(fit) {
  spec <- fit$spec
  out <- list(converged = fit$optimizer$convergence == 0)
  if (!is.null(fit$optimizer$kkt)) {
    # Nonsmooth penalties: the proximal stationarity residual, since the
    # gradient of a lasso type objective does not vanish at its minimum.
    out$kkt_residual <- fit$optimizer$kkt
  } else if (fit$approximation == "none" && length(fit$coefficients)) {
    g <- tryCatch(ml_gradient(spec, unname(fit$coefficients)), error = function(e) NA)
    # Penalized fits: the gradient of the penalized objective.
    if (is.function(fit$penalty$grad)) g <- g - fit$penalty$grad(unname(fit$coefficients))
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
  if (!isTRUE(out$converged)) warning("The optimizer did not report convergence.", call. = FALSE)
  if (identical(out$hessian_positive_definite, FALSE))
    warning("The Hessian is not positive definite; standard errors are unavailable or unreliable.", call. = FALSE)
  out
}

# Factor responses: level codes for categorical and ordinal families, 0/1
# for a two level factor with a binary family.
factor_response <- function(y, fam, ylev, refit) {
  if (refit) y <- factor(as.character(y), levels = ylev)
  if (grepl("ordinal|categorical", fam$category)) return(as.integer(y))
  s <- tryCatch(fam$support(lapply(fam$pars, function(q) q$default %||% 1)), error = function(e) NULL)
  if (nlevels(y) == 2 && fam$discrete && isTRUE(s$lower == 0 && s$upper == 1)) return(as.integer(y) - 1L)
  stop_input(sprintf("A factor response is only allowed for categorical and ordinal families (and two levels with a binary family); '%s' needs a numeric response.", fam$id))
}

# Responses outside the family support: integrality for discrete laws and
# the support bounds that do not depend on the parameters.
check_response <- function(y, fam) {
  g <- lattice_step(fam)
  if (fam$discrete && g > 0 && any(bad <- abs(y / g - round(y / g)) > 1e-8))
    stop_input(sprintf("'%s' needs %s; rows %s are not.", fam$id,
                       if (g == 1) "integer responses" else sprintf("responses on its lattice of multiples of %s", format(g)),
                       paste(utils::head(which(bad), 5), collapse = ", ")))
  sup <- function(scale) tryCatch(fam$support(lapply(fam$pars, function(q) {
    v <- q$default %||% generic_start(q, y)
    if (q$domain == "unit") v / scale else if (q$domain == "real") v + scale - 1 else v * scale
  })), error = function(e) NULL)
  s1 <- sup(1); s2 <- sup(1.7)
  if (is.null(s1) || is.null(s2)) return(invisible())
  lo <- if (length(s1$lower) == 1 && isTRUE(s1$lower == s2$lower)) s1$lower else -Inf
  hi <- if (length(s1$upper) == 1 && isTRUE(s1$upper == s2$upper)) s1$upper else Inf
  bad <- which(y < lo | y > hi)
  if (length(bad))
    stop_input(sprintf("Responses outside the support [%s, %s] of '%s' in rows %s.", format(lo), format(hi), fam$id,
                       paste(utils::head(bad, 5), collapse = ", ")))
  invisible()
}

# Observations: rows (not cells) for vector valued responses.
obs_count <- function(spec) {
  y <- spec$y
  ok <- if (is.matrix(y)) stats::complete.cases(y) else !is.na(y)
  sum(rep_len(spec$weights, length(ok)) > 0 & ok)
}
