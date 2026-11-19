# Stan program generation from a model specification.

stan_fill <- function(template, y, parmap) {
  out <- gsub("{y}", y, template, fixed = TRUE)
  for (k in names(parmap)) out <- gsub(sprintf("{%s}", k), parmap[[k]], out, fixed = TRUE)
  out
}

stan_ident <- function(x) gsub("[^A-Za-z0-9_]", "_", x)

# Family Stan functions plus helpers they need, without duplicates.
stan_family_functions <- function(fam) {
  s <- fam$stan
  helpers <- unique(c(s$needs, unlist(lapply(fam$extra$components, function(f) f$stan$needs))))
  paste(c(vapply(helpers, function(h) stan_helpers[[h]], ""), s$functions), collapse = "\n")
}

cens_type <- function(cmat) {
  if (is.null(cmat)) return(NULL)
  lo <- cmat[, "lower"]; up <- cmat[, "upper"]
  ifelse(lo == up, 0L, ifelse(is.infinite(up), 1L, ifelse(is.infinite(lo), 2L, 3L)))
}

#' Generate the Stan program of a model
#'
#' Returns the 'Stan' code that the `"bayes"` engine runs, together with the
#' data list. The program is generated from the same specification as the
#' frequentist engines: each distribution parameter has a linear predictor
#' with its link, random effects use a non-centered parameterization with
#' Cholesky factored correlations, smooth terms enter through their mixed
#' model representation, and the likelihood uses native 'Stan' functions
#' when the parameterization is identical or package generated functions
#' copied from the registry otherwise. Censored and truncated observations
#' use log CDF and log survival functions.
#'
#' @param formula,...,family,data,weights,links,fixed,parameterization As
#'   in [regcomp()].
#' @param priors Priors from [prior()].
#' @return A list with `code` and `data`.
#' @export
#' @examples
#' d <- data.frame(y = rpois(20, 3), x = rnorm(20))
#' cat(stan_code(y ~ x, family = "poisson", data = d)$code)
stan_code <- function(formula, ..., family, data, weights = NULL, links = NULL,
                      fixed = NULL, priors = NULL, parameterization = NULL) {
  fam <- model_family(family, parameterization, formula, as.data.frame(data))
  dots <- eval(substitute(alist(...)))
  spec <- build_spec(formula, dots, fam, as.data.frame(data), weights = substitute(weights),
                     links = links, fixed = fixed, env = parent.frame())
  list(code = stan_program(spec, priors), data = stan_data(spec))
}

stan_program <- function(spec, priors = NULL) {
  fam <- spec$fam
  if (is.null(fam$stan)) stop_unsupported("A Stan implementation", fam$id)
  if (!fam$variate_type %in% c("scalar", "categorical") && is.null(fam$stan$vector))
    stop_unsupported("Stan code generation for this vector valued family", fam$id)
  ct <- cens_type(spec$cmat)
  censored <- !is.null(ct) && any(ct > 0)
  truncated <- !is.null(spec$cmat) && any(is.finite(spec$cmat[, c("tlower", "tupper")]))
  if ((censored || truncated) && (is.null(fam$stan$lcdf) || is.null(fam$stan$lccdf)))
    stop_unsupported("Censored or truncated responses in Stan (no Stan log CDF)", fam$id)
  ytype <- if (fam$discrete) "array[N] int" else "vector[N]"
  d <- c("  int<lower=1> N;", sprintf("  %s y;", ytype), "  vector[N] w;")
  if (censored || truncated) {
    bt <- if (fam$discrete) "array[N] int" else "vector[N]"
    d <- c(d, "  array[N] int cens;", sprintf("  %s y_upper;", bt),
           "  array[N] int has_tl;", "  array[N] int has_tu;",
           sprintf("  %s tl;", bt), sprintf("  %s tu;", bt))
  }
  for (k in names(spec$data_pars)) {
    intd <- fam$pars[[k]]$domain %in% c("positive_integer", "nonneg_integer", "integer")
    d <- c(d, sprintf("  %s %s;", if (intd) "array[N] int" else "vector[N]", k))
  }
  for (k in names(spec$fixed)) d <- c(d, sprintf("  real fix_%s;", k))
  par <- character(); tpar <- character(); mdl <- character()
  pri <- character()
  for (k in spec$rpars) {
    P <- ncol(spec$preds[[k]]$X)
    d <- c(d, sprintf("  int<lower=1> P_%s;", k), sprintf("  matrix[N, P_%s] X_%s;", k, k),
           sprintf("  vector[N] off_%s;", k))
    par <- c(par, sprintf("  vector[P_%s] b_%s;", k, k))
    cn <- colnames(spec$preds[[k]]$X)
    for (j in seq_along(cn)) {
      cls <- if (cn[j] == "(Intercept)") "Intercept" else "b"
      pri <- c(pri, "  " %+% prior_stan(find_prior(priors, cls, k, cn[j]), sprintf("b_%s[%d]", k, j)))
    }
  }
  for (k in names(spec$aux)) {
    d <- c(d, sprintf("  int<lower=1> Kc_%s;", k))
    par <- c(par, sprintf("  %s[Kc_%s] %s;", if (spec$aux[[k]]$type == "ordered") "ordered" else "vector", k, k))
    pri <- c(pri, "  " %+% prior_stan(find_prior(priors, "cuts", k), k))
  }
  for (j in seq_along(spec$blocks)) {
    b <- spec$blocks[[j]]
    if (b$type != "smooth" && !b$struct %in% c("us", "diag")) {
      st <- stan_struct_block(b, j, priors)
      d <- c(d, st$data); par <- c(par, st$par); tpar <- c(tpar, st$tpar); pri <- c(pri, st$pri)
      next
    }
    if (b$type == "smooth") {
      d <- c(d, sprintf("  int<lower=1> L_%d;", j), sprintf("  matrix[N, L_%d] Zs_%d;", j, j))
      par <- c(par, sprintf("  real<lower=0> sd_%d;", j), sprintf("  vector[L_%d] z_%d;", j, j))
      pri <- c(pri, sprintf("  z_%d ~ std_normal();", j),
               "  " %+% prior_stan(find_prior(priors, "sd", b$par, b$label), sprintf("sd_%d", j)))
    } else {
      d <- c(d, sprintf("  int<lower=1> L_%d;", j), sprintf("  int<lower=1> q_%d;", j),
             if (is.null(b$W)) sprintf("  array[N] int<lower=1, upper=L_%d> g_%d;", j, j) else sprintf("  matrix[N, L_%d] W_%d;", j, j),
             sprintf("  matrix[N, q_%d] Xr_%d;", j, j))
      par <- c(par, sprintf("  vector<lower=0>[q_%d] sd_%d;", j, j),
               sprintf("  matrix[q_%d, L_%d] z_%d;", j, j, j))
      if (b$corr) {
        par <- c(par, sprintf("  cholesky_factor_corr[q_%d] Lc_%d;", j, j))
        tpar <- c(tpar, sprintf("  matrix[q_%d, L_%d] r_%d = diag_pre_multiply(sd_%d, Lc_%d) * z_%d;", j, j, j, j, j, j))
        pri <- c(pri, "  " %+% prior_stan(find_prior(priors, "cor", b$par, b$label), sprintf("Lc_%d", j)))
      } else {
        tpar <- c(tpar, sprintf("  matrix[q_%d, L_%d] r_%d = diag_pre_multiply(sd_%d, z_%d);", j, j, j, j, j))
      }
      pri <- c(pri, sprintf("  to_vector(z_%d) ~ std_normal();", j),
               "  " %+% prior_stan(find_prior(priors, "sd", b$par, b$label), sprintf("sd_%d", j)))
    }
  }
  # Linear predictors and parameters.
  for (k in spec$rpars) {
    mdl <- c(mdl, sprintf("  vector[N] eta_%s = X_%s * b_%s + off_%s;", k, k, k, k))
    for (j in which(vapply(spec$blocks, `[[`, "", "par") == k)) {
      b <- spec$blocks[[j]]
      if (b$type == "smooth") mdl <- c(mdl, sprintf("  eta_%s += Zs_%d * (sd_%d * z_%d);", k, j, j, j))
      else if (!is.null(b$W)) mdl <- c(mdl, sprintf("  for (c in 1:q_%d) eta_%s += Xr_%d[:, c] .* (W_%d * (r_%d[c])');", j, k, j, j, j))
      else mdl <- c(mdl, sprintf("  for (i in 1:N) eta_%s[i] += Xr_%d[i] * col(r_%d, g_%d[i]);", k, j, j, j))
    }
  }
  for (k in spec$rpars)
    mdl <- c(mdl, sprintf("  vector[N] %s = %s;", k, spec$links[[k]]$stan_inv(sprintf("eta_%s", k))))
  for (k in names(spec$fixed)) mdl <- c(mdl, sprintf("  vector[N] %s = rep_vector(fix_%s, N);", k, k))
  pm <- stats::setNames(lapply(names(fam$pars), function(k)
    if (k %in% names(spec$aux) || !is.null(fam$pars[[k]]$estimate)) k else sprintf("%s[i]", k)), names(fam$pars))
  for (k in intersect(names(spec$fixed), names(Filter(function(p) !is.null(p$estimate), fam$pars))))
    stop_unsupported("Fixed ordinal thresholds in Stan")
  lp <- stan_fill(fam$stan$lpdf, "y[i]", pm)
  lik <- if (!(censored || truncated)) {
    c("  for (i in 1:N) target += w[i] * " %+% lp %+% ";")
  } else {
    m1 <- if (fam$discrete) " - 1" else ""
    lcdf <- function(v) stan_fill(fam$stan$lcdf, v, pm)
    lccdf <- function(v) stan_fill(fam$stan$lccdf, v, pm)
    c("  for (i in 1:N) {",
      "    real lpi;",
      sprintf("    if (cens[i] == 0) lpi = %s;", lp),
      sprintf("    else if (cens[i] == 1) lpi = %s;", lccdf("y[i]" %+% m1)),
      sprintf("    else if (cens[i] == 2) lpi = %s;", lcdf("y_upper[i]")),
      sprintf("    else lpi = log_diff_exp(%s, %s);", lcdf("y_upper[i]"), lcdf("y[i]" %+% m1)),
      "    if (has_tl[i] == 1 && has_tu[i] == 1)",
      sprintf("      lpi -= log_diff_exp(%s, %s);", lcdf("tu[i]"), lcdf("tl[i]" %+% m1)),
      "    else if (has_tl[i] == 1)",
      sprintf("      lpi -= %s;", lccdf("tl[i]" %+% m1)),
      "    else if (has_tu[i] == 1)",
      sprintf("      lpi -= %s;", lcdf("tu[i]")),
      "    target += w[i] * lpi;",
      "  }")
  }
  code <- c("functions {", stan_family_functions(fam), "}",
            "data {", d, "}",
            "parameters {", par, "}",
            if (length(tpar)) c("transformed parameters {", tpar, "}"),
            "model {", mdl, pri, lik, "}")
  paste(code, collapse = "\n")
}

`%+%` <- function(a, b) paste0(a, b)

stan_data <- function(spec) {
  fam <- spec$fam
  out <- list(N = spec$n, w = spec$weights)
  y <- spec$y
  ct <- cens_type(spec$cmat)
  big <- 2147483647
  toint <- function(v) { v[is.infinite(v)] <- ifelse(v[is.infinite(v)] > 0, big - 1, -big + 1); as.integer(round(v)) }
  out$y <- if (fam$discrete) toint(y) else as.numeric(y)
  if (!is.null(ct) && (any(ct > 0) || any(is.finite(spec$cmat[, c("tlower", "tupper")])))) {
    cm <- spec$cmat
    out$cens <- as.integer(ct)
    yu <- cm[, "upper"]; yu[!is.finite(yu)] <- 0
    out$y_upper <- if (fam$discrete) toint(yu) else yu
    out$has_tl <- as.integer(is.finite(cm[, "tlower"]))
    out$has_tu <- as.integer(is.finite(cm[, "tupper"]))
    tl <- cm[, "tlower"]; tl[!is.finite(tl)] <- 0
    tu <- cm[, "tupper"]; tu[!is.finite(tu)] <- 0
    out$tl <- if (fam$discrete) toint(tl) else tl
    out$tu <- if (fam$discrete) toint(tu) else tu
    y2 <- y; y2[!is.finite(y2)] <- 0
    out$y <- if (fam$discrete) toint(y2) else y2
  }
  for (k in names(spec$data_pars)) {
    intd <- fam$pars[[k]]$domain %in% c("positive_integer", "nonneg_integer", "integer")
    out[[k]] <- if (intd) as.integer(spec$data_pars[[k]]) else as.numeric(spec$data_pars[[k]])
  }
  for (k in names(spec$fixed)) out[[paste0("fix_", k)]] <- spec$fixed[[k]]
  for (k in spec$rpars) {
    X <- spec$preds[[k]]$X
    out[[paste0("P_", k)]] <- ncol(X)
    out[[paste0("X_", k)]] <- unname(X)
    out[[paste0("off_", k)]] <- as.numeric(spec$preds[[k]]$offset)
  }
  for (k in names(spec$aux)) out[[paste0("Kc_", k)]] <- spec$aux[[k]]$size
  for (j in seq_along(spec$blocks)) {
    b <- spec$blocks[[j]]
    out[[paste0("L_", j)]] <- b$L
    if (b$type == "smooth") out[[paste0("Zs_", j)]] <- unname(b$Xr)
    else {
      out[[paste0("q_", j)]] <- b$q
      if (!is.null(b$W)) out[[paste0("W_", j)]] <- as.matrix(b$W)
      else {
        if (anyNA(b$gi)) stop_input("Grouping factors must not contain missing values.")
        out[[paste0("g_", j)]] <- as.integer(b$gi)
      }
      out[[paste0("Xr_", j)]] <- unname(b$Xr)
      if (!is.null(b$coords)) {
        D <- as.matrix(stats::dist(b$coords))
        out[[paste0("D_", j)]] <- unname(D)
        out[[paste0("Dmed_", j)]] <- stats::median(D[upper.tri(D)])
      }
    }
  }
  out
}

# Stan code of a structured random effect block: r_j = diag(sd) chol(R) z_j.
stan_struct_block <- function(b, j, priors) {
  q <- sprintf("q_%d", j); L <- sprintf("L_%d", j)
  data <- c(sprintf("  int<lower=1> L_%d;", j), sprintf("  int<lower=1> q_%d;", j),
            if (is.null(b$W)) sprintf("  array[N] int<lower=1, upper=L_%d> g_%d;", j, j) else sprintf("  matrix[N, L_%d] W_%d;", j, j),
            sprintf("  matrix[N, q_%d] Xr_%d;", j, j))
  hom <- b$struct %in% c("homdiag", "homcs", "ar1", "ou", "exp", "gau", "matern32")
  par <- c(if (hom) sprintf("  real<lower=0> sdh_%d;", j) else sprintf("  vector<lower=0>[%s] sd_%d;", q, j),
           sprintf("  matrix[%s, %s] z_%d;", q, L, j))
  tpar <- if (hom) sprintf("  vector[%s] sd_%d = rep_vector(sdh_%d, %s);", q, j, j, q)
  pri <- c(sprintf("  to_vector(z_%d) ~ std_normal();", j),
           "  " %+% prior_stan(find_prior(priors, "sd", b$par, b$label), if (hom) sprintf("sdh_%d", j) else sprintf("sd_%d", j)))
  fill <- function(expr) c(sprintf("  matrix[%s, %s] Rc_%d;", q, q, j),
    sprintf("  for (a in 1:%s) for (c in 1:%s) Rc_%d[a, c] = a == c ? 1 + 1e-9 : %s;", q, q, j, expr))
  if (b$struct %in% c("cs", "homcs")) {
    par <- c(par, sprintf("  real<lower=-1.0 / (%s - 1), upper=1> rho_%d;", q, j))
    tpar <- c(tpar, fill(sprintf("rho_%d", j)))
  } else if (b$struct == "ar1") {
    par <- c(par, sprintf("  real<lower=-1, upper=1> rho_%d;", j))
    tpar <- c(tpar, fill(sprintf("rho_%d^abs(a - c)", j)))
  } else if (b$struct == "homdiag") {
    tpar <- c(tpar, fill("0"))
  } else {
    data <- c(data, sprintf("  matrix[%s, %s] D_%d;", q, q, j), sprintf("  real<lower=0> Dmed_%d;", j))
    par <- c(par, sprintf("  real<lower=0> range_%d;", j))
    dd <- sprintf("D_%d[a, c] / range_%d", j, j)
    tpar <- c(tpar, fill(switch(b$struct, ou = , exp = sprintf("exp(-%s)", dd), gau = sprintf("exp(-square(%s))", dd),
      matern32 = sprintf("(1 + sqrt(3) * %s) * exp(-sqrt(3) * %s)", dd, dd))))
    pri <- c(pri, sprintf("  range_%d ~ lognormal(log(Dmed_%d), 1.5);", j, j))
  }
  tpar <- c(tpar, sprintf("  matrix[%s, %s] r_%d = diag_pre_multiply(sd_%d, cholesky_decompose(Rc_%d)) * z_%d;", q, L, j, j, j, j))
  list(data = data, par = par, tpar = tpar, pri = pri)
}
