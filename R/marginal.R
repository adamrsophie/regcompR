# Marginal (population averaged) inference: cluster robust covariances for
# likelihood fits, and generalized estimating equations with quasi
# likelihood variance functions.

# Score contributions of each observation (rows) for every coefficient.
obs_scores <- function(fit) {
  spec <- fit$spec
  beta <- unname(fit$coefficients)
  eta <- compute_eta(spec, beta)
  d <- eta_derivs(spec, eta, second = FALSE)
  S <- matrix(0, spec$n, length(beta))
  for (a in seq_along(spec$rpars)) {
    k <- spec$rpars[a]
    S[, spec$beta_idx[[k]]] <- spec$preds[[k]]$X * (spec$weights * d$g[, a])
  }
  for (k in names(spec$aux)) for (j in spec$beta_idx[[k]]) {
    h <- 1e-6 * max(1, abs(beta[j])); bp <- beta; bm <- beta
    bp[j] <- bp[j] + h; bm[j] <- bm[j] - h
    S[, j] <- spec$weights * (loglik_eta(spec, compute_eta(spec, bp)) - loglik_eta(spec, compute_eta(spec, bm))) / (2 * h)
  }
  S[!is.finite(S)] <- 0
  S
}

#' Covariance of the estimates
#'
#' `type = "model"` gives the inverse observed information (or the
#' posterior covariance of a Bayesian fit). `type = "sandwich"` gives the
#' robust covariance `V M V` of Huber and White, where `V` is the model
#' covariance and `M` the summed outer products of score contributions;
#' with `cluster` the scores are first summed within clusters (Liang and
#' Zeger 1986) and `M` is multiplied by `G / (G - 1)` for `G` clusters.
#' This is the independence working model of generalized estimating
#' equations for any family; see [gee()] for other working correlations.
#'
#' @param object A fitted model.
#' @param type `"model"` or `"sandwich"`.
#' @param cluster Optional cluster identifier: a vector, or a one sided
#'   formula evaluated in the model data.
#' @param ... Unused.
#' @return A covariance matrix of the coefficients.
#' @export
#' @examples
#' set.seed(1)
#' d <- data.frame(g = rep(1:30, each = 5), x = rnorm(150))
#' d$y <- rpois(150, exp(0.3 + 0.5 * d$x + rnorm(30, 0, 0.5)[d$g]))
#' fit <- regcomp(y ~ x, family = "poisson", data = d)
#' sqrt(diag(vcov(fit, type = "sandwich", cluster = ~ g)))
vcov.regcomp_fit <- function(object, type = c("model", "sandwich"), cluster = NULL, ...) {
  type <- match.arg(type)
  bn <- object$spec$beta_names
  if (object$inference == "bayesian") {
    if (type == "sandwich") stop_unsupported("A sandwich covariance for a Bayesian fit")
    return(stats::cov(object$draws[, bn, drop = FALSE]))
  }
  v <- object$vcov
  if (is.null(v)) stop_unsupported("A covariance matrix for this fit (penalized or hessian = FALSE)")
  V <- v[bn, bn, drop = FALSE]
  if (type == "model") return(V)
  if (object$spec$Q > 0) stop_unsupported("A sandwich covariance for models with random effects or smooths")
  S <- obs_scores(object)
  if (!is.null(cluster)) {
    if (inherits(cluster, "formula"))
      cluster <- eval(cluster[[2]], eval(object$call$data, environment(object$spec$formulas[[1]])), environment(cluster))
    if (length(cluster) != nrow(S) || anyNA(cluster)) stop_input("cluster must have one nonmissing value per observation.")
    S <- rowsum(S, as.character(cluster), reorder = FALSE)
    G <- nrow(S)
    if (G < 2) stop_input("A clustered sandwich covariance needs at least two clusters.")
    M <- crossprod(S) * G / (G - 1)
  } else M <- crossprod(S)
  out <- V %*% M %*% V
  dimnames(out) <- list(bn, bn)
  out
}

#' Generalized estimating equations
#'
#' Fits a marginal mean model `g(E(y)) = X beta` with variance `phi V(mu)`
#' and a working correlation within clusters by the method of Liang and
#' Zeger (1986), with moment estimates of the correlation and scale and a
#' robust (sandwich) covariance. Quasi likelihood models (McCullagh and
#' Nelder 1989, ch. 9) are the independence case: `gee(y ~ x, data, family
#' = quasipoisson())` gives the quasi Poisson estimates with a Pearson
#' scale.
#'
#' @param formula Mean model formula; may contain an offset.
#' @param data Data frame.
#' @param id Cluster identifier (unquoted column name or vector). Rows of a
#'   cluster are taken in data order unless `waves` is given.
#' @param family A `stats` family object such as `gaussian()`, `poisson()`,
#'   `binomial()`, `Gamma()`, `inverse.gaussian()`, `quasi()`,
#'   `quasipoisson()` or `quasibinomial()`.
#' @param corstr Working correlation: `"independence"`, `"exchangeable"`,
#'   `"ar1"` or `"unstructured"`.
#' @param waves Optional within cluster time index (integers) for `"ar1"`
#'   and `"unstructured"`.
#' @param weights Optional prior weights.
#' @param maxit,tol Iteration limit and convergence tolerance.
#' @return An object of class `regcomp_gee` with `coefficients`, robust
#'   covariance `vcov`, naive (model based) covariance `vcov_naive`,
#'   `alpha` (working correlation parameters), `scale` and `fitted.values`.
#' @export
#' @examples
#' set.seed(1)
#' d <- data.frame(g = rep(1:40, each = 4), x = rnorm(160))
#' d$y <- rpois(160, exp(0.2 + 0.4 * d$x + rnorm(40, 0, 0.4)[d$g]))
#' gee(y ~ x, data = d, id = g, family = poisson(), corstr = "exchangeable")
gee <- function(formula, data, id, family = stats::gaussian(),
                corstr = c("independence", "exchangeable", "ar1", "unstructured"),
                waves = NULL, weights = NULL, maxit = 50, tol = 1e-8) {
  corstr <- match.arg(corstr)
  if (is.character(family)) family <- get(family, mode = "function")()
  if (is.function(family)) family <- family()
  mf <- stats::model.frame(formula, data, na.action = stats::na.fail)
  y <- stats::model.response(mf)
  X <- stats::model.matrix(stats::terms(mf), mf)
  off <- stats::model.offset(mf); if (is.null(off)) off <- rep(0, nrow(X))
  if (is.matrix(y)) { w0 <- rowSums(y); y <- y[, 1] / w0 } else w0 <- rep(1, length(y))
  # Binomial factor responses: first level is failure, as in glm().
  if (is.factor(y)) y <- as.numeric(y != levels(y)[1])
  if (is.logical(y)) y <- as.numeric(y)
  idv <- eval(substitute(id), data, parent.frame())
  if (length(idv) != length(y) || anyNA(idv)) stop_input("id must give a cluster for every row (no missing values).")
  wv <- eval(substitute(waves), data, parent.frame())
  pw <- eval(substitute(weights), data, parent.frame())
  pw <- w0 * if (is.null(pw)) 1 else pw
  n <- length(y); p <- ncol(X)
  cl <- split(seq_len(n), factor(idv, levels = unique(idv)))
  if (corstr != "independence" && all(lengths(cl) < 2)) {
    warning("No cluster has more than one observation; the working correlation is not estimable, so independence is used.", call. = FALSE)
    corstr <- "independence"
  }
  if (is.null(wv)) wv <- unlist(lapply(cl, seq_along))[order(unlist(cl))]
  Tmax <- max(wv)
  # Start from the independence (GLM) fit.
  beta <- stats::glm.fit(X, y, weights = pw, offset = off, family = family)$coefficients
  beta[is.na(beta)] <- 0
  Rfun <- function(alpha, t) {
    m <- length(t)
    switch(corstr, independence = diag(m),
      exchangeable = { R <- matrix(alpha, m, m); diag(R) <- 1; R },
      ar1 = alpha^abs(outer(t, t, `-`)),
      unstructured = { R <- alpha[t, t, drop = FALSE]; diag(R) <- 1; R })
  }
  alpha <- switch(corstr, independence = 0, exchangeable = 0, ar1 = 0, unstructured = diag(Tmax))
  for (it in seq_len(maxit)) {
    eta <- drop(X %*% beta) + off
    mu <- family$linkinv(eta); dmu <- family$mu.eta(eta)
    v <- family$variance(mu) / pw
    r <- (y - mu) / sqrt(v)
    phi <- sum(r^2) / (n - p)
    # Working correlation by least squares on the pairwise products of
    # standardized residuals (the second estimating equation of geepack).
    if (corstr == "exchangeable") {
      # Sum of pair products per cluster: ((sum r)^2 - sum r^2) / 2.
      sp <- 0; np <- 0
      for (i in cl) if (length(i) > 1) { ri <- r[i]; sp <- sp + (sum(ri)^2 - sum(ri^2)) / 2; np <- np + length(i) * (length(i) - 1) / 2 }
      mmax <- max(lengths(cl))
      alpha <- if (np > 0) sp / np / (sum(r^2) / n) else 0
      # Keep the working correlation positive definite for every cluster.
      if (mmax > 1) alpha <- min(max(alpha, -1 / (mmax - 1) + 1e-6), 1 - 1e-6)
    } else if (corstr == "ar1") {
      z <- NULL; lag <- NULL
      for (i in cl) if (length(i) > 1) {
        pr <- utils::combn(length(i), 2)
        z <- c(z, r[i][pr[1, ]] * r[i][pr[2, ]]); lag <- c(lag, abs(wv[i][pr[1, ]] - wv[i][pr[2, ]]))
      }
      z <- z / (sum(r^2) / n)
      alpha <- stats::optimize(function(a) sum((z - sign(a)^lag * abs(a)^lag)^2), c(-0.999, 0.999), tol = 1e-10)$minimum
    } else if (corstr == "unstructured") {
      S <- matrix(0, Tmax, Tmax); N <- matrix(0, Tmax, Tmax)
      for (i in cl) { t <- wv[i]; S[t, t] <- S[t, t] + tcrossprod(r[i]); N[t, t] <- N[t, t] + 1 }
      alpha <- S / (phi * pmax(N - p / Tmax, 1)); diag(alpha) <- 1
    }
    A <- matrix(0, p, p); U <- numeric(p)
    for (i in cl) {
      D <- X[i, , drop = FALSE] * dmu[i]
      Vi <- sqrt(v[i]) * t(sqrt(v[i]) * Rfun(alpha, wv[i]))
      Wi <- solve(Vi, D)
      A <- A + crossprod(D, Wi); U <- U + drop(crossprod(Wi, y[i] - mu[i]))
    }
    step <- solve(A, U)
    beta <- beta + step
    if (max(abs(step)) < tol * (1 + max(abs(beta)))) break
  }
  eta <- drop(X %*% beta) + off
  mu <- family$linkinv(eta); dmu <- family$mu.eta(eta); v <- family$variance(mu) / pw
  A <- matrix(0, p, p); B <- matrix(0, p, p)
  for (i in cl) {
    D <- X[i, , drop = FALSE] * dmu[i]
    Vi <- phi * sqrt(v[i]) * t(sqrt(v[i]) * Rfun(alpha, wv[i]))
    Wi <- solve(Vi, D)
    A <- A + crossprod(D, Wi)
    s <- drop(crossprod(Wi, y[i] - mu[i])); B <- B + tcrossprod(s)
  }
  Ai <- solve(A)
  names(beta) <- colnames(X)
  vr <- Ai %*% B %*% Ai; dimnames(vr) <- dimnames(Ai) <- list(colnames(X), colnames(X))
  structure(list(coefficients = beta, vcov = vr, vcov_naive = Ai, alpha = alpha, scale = phi,
                 corstr = corstr, family = family, fitted.values = mu, iterations = it,
                 converged = it < maxit, nclusters = length(cl), call = match.call()),
            class = "regcomp_gee")
}

#' @export
print.regcomp_gee <- function(x, ...) {
  cat(sprintf("<regcomp_gee> %s(%s), %s working correlation, %d clusters\n",
              x$family$family, x$family$link, x$corstr, x$nclusters))
  se <- sqrt(diag(x$vcov))
  print(cbind(estimate = x$coefficients, robust_se = se, z = x$coefficients / se))
  if (x$corstr %in% c("exchangeable", "ar1")) cat(sprintf("  correlation %.4f", x$alpha))
  cat(sprintf("  scale %.4f\n", x$scale))
  invisible(x)
}

#' @export
coef.regcomp_gee <- function(object, ...) object$coefficients

#' @export
vcov.regcomp_gee <- function(object, type = c("robust", "naive"), ...)
  if (match.arg(type) == "robust") object$vcov else object$vcov_naive
