# Model specification: per parameter formulas, fixed effects, random effects
# (lme4 style bars), smooth terms (mgcv bases) and offsets.

#' Censored, truncated and interval observations
#'
#' Wrap the response in `cens()` to describe the observation process.
#' A row with `lower == upper` is observed exactly. `upper = Inf` gives right
#' censoring at `lower`, `lower = -Inf` left censoring at `upper`, and finite
#' unequal bounds interval censoring. For discrete responses an interval
#' `[lower, upper]` includes both end points. `trunc_lower` and
#' `trunc_upper` give per observation truncation (delayed entry), with the
#' retained region `(trunc_lower, trunc_upper]` for continuous responses and
#' `[trunc_lower, trunc_upper]` for discrete ones. `surv()` is a shortcut for
#' right censored survival times with an event indicator.
#'
#' Likelihood contributions: exact `log f(y)`, left censored `log F(u)`,
#' right censored `log S(l)`, interval `log(F(u) - F(l))`, each minus
#' `log(F(tu) - F(tl))` under truncation.
#'
#' @param lower,upper Lower and upper bounds of the observed value.
#' @param trunc_lower,trunc_upper Truncation bounds.
#' @param time Event or censoring time.
#' @param event 1 for an observed event, 0 for right censoring.
#' @param entry Optional delayed entry (left truncation) time.
#' @return An object of class `regcomp_cens`.
#' @export
#' @examples
#' cens(c(1, 2, 3), c(1, Inf, 4))
#' surv(c(5, 8), c(1, 0))
cens <- function(lower, upper = lower, trunc_lower = -Inf, trunc_upper = Inf) {
  n <- max(length(lower), length(upper))
  out <- cbind(lower = rep_len(as.numeric(lower), n),
               upper = rep_len(as.numeric(upper), n),
               tlower = rep_len(as.numeric(trunc_lower), n),
               tupper = rep_len(as.numeric(trunc_upper), n))
  if (any(out[, "upper"] < out[, "lower"], na.rm = TRUE))
    stop_input("cens(): upper must not be smaller than lower.")
  structure(out, class = c("regcomp_cens", "matrix"))
}

#' @rdname cens
#' @export
surv <- function(time, event, entry = -Inf) {
  event <- as.numeric(event)
  if (any(!event %in% c(0, 1)))
    stop_input("surv(): event must be 0 (censored) or 1 (event).")
  cens(time, ifelse(event == 1, time, Inf), trunc_lower = entry)
}

#' @export
print.regcomp_cens <- function(x, ...) {
  lab <- ifelse(x[, "lower"] == x[, "upper"], format(x[, "lower"]),
         ifelse(is.infinite(x[, "upper"]), paste0(format(x[, "lower"]), "+"),
         ifelse(is.infinite(x[, "lower"]), paste0(format(x[, "upper"]), "-"),
                paste0("[", format(x[, "lower"]), ", ", format(x[, "upper"]), "]"))))
  print(noquote(lab))
  invisible(x)
}

# ---------------------------------------------------------------------------
# Random effect bars: (expr | group) and (expr || group).

find_bars <- function(term) {
  if (is.name(term) || !is.language(term)) return(NULL)
  if (term[[1]] == as.name("(")) return(find_bars(term[[2]]))
  if (term[[1]] == as.name("|") || term[[1]] == as.name("||")) return(list(term))
  if (is_struct_call(term)) return(list(term))
  if (length(term) == 2) return(find_bars(term[[2]]))
  c(find_bars(term[[2]]), find_bars(term[[3]]))
}

drop_bars <- function(term) {
  if (is.name(term) || !is.language(term)) return(term)
  if (term[[1]] == as.name("(") && is.call(term[[2]]) &&
      (term[[2]][[1]] == as.name("|") || term[[2]][[1]] == as.name("||"))) return(NULL)
  if (term[[1]] == as.name("|") || term[[1]] == as.name("||")) return(NULL)
  if (length(term) == 2) {
    inner <- drop_bars(term[[2]])
    if (is.null(inner)) return(NULL)
    term[[2]] <- inner
    return(term)
  }
  if (length(term) == 3 && (term[[1]] == as.name("+") || term[[1]] == as.name("-") ||
                            term[[1]] == as.name("~"))) {
    l <- drop_bars(term[[2]]); r <- drop_bars(term[[3]])
    if (term[[1]] == as.name("~")) {
      term[[3]] <- if (is.null(r)) 1 else r
      return(term)
    }
    if (is.null(l)) return(r)
    if (is.null(r)) return(l)
    term[[2]] <- l; term[[3]] <- r
    return(term)
  }
  term
}

# Covariance structures of random effect blocks (glmmTMB style wrappers):
# us (unstructured), diag, homdiag, cs (heterogeneous compound symmetry),
# homcs, ar1 (homogeneous first order autoregressive over the columns), and
# the distance based ou / exp (exponential), gau (Gaussian) and Matern
# with smoothness fixed at 3/2, matern32) whose left side gives the coordinates.
re_structs <- c("us", "diag", "homdiag", "cs", "homcs", "ar1", "ou", "exp", "gau", "matern32")
coord_structs <- c("ou", "exp", "gau", "matern32")
is_struct_call <- function(x) is.call(x) && is.name(x[[1]]) && as.character(x[[1]]) %in% re_structs &&
  length(x) == 2 && length(find_bars_plain(x[[2]])) == 1
find_bars_plain <- function(x) {
  while (is.call(x) && identical(x[[1]], as.name("("))) x <- x[[2]]
  if (is.call(x) && (identical(x[[1]], as.name("|")) || identical(x[[1]], as.name("||")))) list(x) else NULL
}

# Nested grouping: (x | a/b) is (x | a) + (x | a:b).
expand_nested <- function(bars) unlist(lapply(bars, function(b) {
  g <- b[[3]]
  if (is.call(g) && identical(g[[1]], as.name("/"))) {
    inner <- b; inner[[3]] <- call(":", g[[2]], g[[3]]); outer <- b; outer[[3]] <- g[[2]]
    return(c(expand_nested(list(outer)), list(inner)))
  }
  list(b)
}), recursive = FALSE)

# Grouping factor; a:b is the interaction of a and b, not a sequence.
eval_group <- function(gexpr, data, env) {
  if (is.call(gexpr) && identical(gexpr[[1]], as.name(":")))
    return(interaction(eval_group(gexpr[[2]], data, env), eval_group(gexpr[[3]], data, env), drop = TRUE, sep = ":"))
  g <- eval(gexpr, data, env)
  if (is.list(g)) g <- interaction(g, drop = TRUE)
  g
}

is_smooth_call <- function(x) is.call(x) &&
  as.character(x[[1]]) %in% c("s", "te", "ti", "t2")

find_smooths <- function(term) {
  if (is.name(term) || !is.language(term)) return(NULL)
  if (is_smooth_call(term)) return(list(term))
  if (length(term) == 2) return(find_smooths(term[[2]]))
  if (length(term) == 3) return(c(find_smooths(term[[2]]), find_smooths(term[[3]])))
  NULL
}

drop_smooths <- function(term) {
  if (is_smooth_call(term)) return(NULL)
  if (is.name(term) || !is.language(term)) return(term)
  if (length(term) == 3 && (term[[1]] == as.name("+") || term[[1]] == as.name("~"))) {
    l <- drop_smooths(term[[2]]); r <- drop_smooths(term[[3]])
    if (term[[1]] == as.name("~")) { term[[3]] <- if (is.null(r)) 1 else r; return(term) }
    if (is.null(l)) return(r)
    if (is.null(r)) return(l)
    term[[2]] <- l; term[[3]] <- r
    return(term)
  }
  term
}

# Build the design of one parameter: fixed matrix, offset, random blocks.
# `rhs` is a one sided formula.
build_predictor <- function(rhs, data, par, xlev = NULL, smooth_specs = NULL,
                            re_levels = NULL) {
  env <- environment(rhs)
  f_fixed <- drop_smooths(drop_bars(rhs))
  if (is.null(f_fixed)) f_fixed <- ~1
  f_fixed <- stats::as.formula(f_fixed, env = env)
  mf <- stats::model.frame(f_fixed, data, na.action = stats::na.pass,
                           xlev = xlev, drop.unused.levels = is.null(xlev))
  tt <- stats::terms(mf)
  X <- stats::model.matrix(tt, mf)
  off <- stats::model.offset(mf)
  if (is.null(off)) off <- rep(0, nrow(X))
  # Random effects.
  bars <- find_bars(rhs[[length(rhs)]])
  bars <- unlist(lapply(bars, function(b) {
    st <- if (is_struct_call(b)) as.character(b[[1]]) else NULL
    inner <- if (is.null(st)) b else find_bars_plain(b[[2]])[[1]]
    lapply(expand_nested(list(inner)), function(x) list(bar = x, struct = st,
      key = if (is.null(st)) x else call(st, x)))
  }), recursive = FALSE)
  blocks <- list()
  for (bb in bars) {
    b <- bb$bar
    struct <- if (!is.null(bb$struct)) bb$struct else if (b[[1]] == as.name("|")) "us" else "diag"
    coords <- NULL
    if (struct %in% coord_structs) {
      vars <- all.vars(b[[2]])
      C <- do.call(cbind, lapply(vars, function(v) as.numeric(eval(as.name(v), data, env))))
      ck <- do.call(paste, c(lapply(seq_len(ncol(C)), function(k) C[, k]), sep = "\r"))
      coords <- if (!is.null(re_levels[[paste0(".coords.", paste(deparse(bb$key, width.cutoff = 500L), collapse = ""))]]))
        re_levels[[paste0(".coords.", paste(deparse(bb$key, width.cutoff = 500L), collapse = ""))]] else
        unique(C[order(ck), , drop = FALSE])
      uk <- do.call(paste, c(lapply(seq_len(ncol(coords)), function(k) coords[, k]), sep = "\r"))
      Xr <- outer(match(ck, uk), seq_along(uk), `==`) * 1
      colnames(Xr) <- paste0("(", apply(coords, 1, paste, collapse = ","), ")")
    } else {
      lhs <- stats::as.formula(call("~", b[[2]]), env = env)
      Xr <- stats::model.matrix(lhs, stats::model.frame(lhs, data, na.action = stats::na.pass))
    }
    if (ncol(Xr) == 1 && struct %in% c("us", "cs", "homcs")) struct <- "diag"
    corr <- struct == "us"
    gexpr <- b[[3]]
    key <- paste(deparse(bb$key, width.cutoff = 500L), collapse = "")
    W <- NULL
    if (is.call(gexpr) && identical(gexpr[[1]], as.name("mm"))) {
      # Multiple membership: mm(g1, g2, ..., weights = cbind(w1, w2, ...)).
      args <- as.list(gexpr)[-1]
      wexpr <- args$weights; args$weights <- NULL
      G <- lapply(args, function(a) as.character(eval(a, data, env)))
      Wm <- if (is.null(wexpr)) matrix(1 / length(G), length(G[[1]]), length(G)) else as.matrix(eval(wexpr, data, env))
      lev <- if (!is.null(re_levels[[key]])) re_levels[[key]] else sort(unique(stats::na.omit(unlist(G))))
      ii <- rep(seq_along(G[[1]]), length(G)); ll <- match(unlist(G), lev); ww <- as.vector(Wm)
      keep <- !is.na(ll) & ww != 0
      W <- Matrix::sparseMatrix(i = ii[keep], j = ll[keep], x = ww[keep], dims = c(length(G[[1]]), length(lev)))
      gi <- rep(NA_integer_, length(G[[1]]))
    } else {
      g <- eval_group(gexpr, data, env)
      lev <- if (!is.null(re_levels[[key]])) re_levels[[key]] else
        levels(droplevels(as.factor(g)))
      gi <- match(as.character(g), lev)
    }
    blocks[[length(blocks) + 1L]] <- list(
      label = key, group = paste(deparse(gexpr), collapse = ""),
      terms = colnames(Xr), Xr = Xr, gi = gi, levels = lev,
      q = ncol(Xr), L = length(lev), corr = corr && ncol(Xr) > 1,
      struct = struct, coords = coords, W = W, type = "re", par = par)
  }
  # Smooth terms through mgcv's mixed model representation.
  sm <- find_smooths(rhs[[length(rhs)]])
  smooth_info <- list()
  for (j in seq_along(sm)) {
    if (!requireNamespace("mgcv", quietly = TRUE))
      stop_unsupported("Smooth terms need the 'mgcv' package; it")
    sobj <- eval(sm[[j]], envir = asNamespace("mgcv"))
    key <- paste(deparse(sm[[j]], width.cutoff = 500L), collapse = "")
    if (is.null(smooth_specs[[key]])) {
      sc <- mgcv::smoothCon(sobj, data = as.data.frame(data), absorb.cons = TRUE,
                            scale.penalty = TRUE)[[1]]
      rp <- mgcv::smooth2random(sc, names(data), type = 2)
      spec <- list(sc = sc, rp = rp)
    } else spec <- smooth_specs[[key]]
    Xs <- mgcv::PredictMat(spec$sc, as.data.frame(data))
    # Transform to the random effect parameterization.
    Xs <- Xs %*% spec$rp$trans.D
    if (!is.null(spec$rp$trans.U)) Xs <- Xs %*% spec$rp$trans.U
    nr <- spec$rp$rind
    fixed_cols <- setdiff(seq_len(ncol(Xs)), unlist(nr))
    if (length(fixed_cols)) {
      Xf <- Xs[, fixed_cols, drop = FALSE]
      colnames(Xf) <- paste0(key, ".fx", seq_along(fixed_cols))
      X <- cbind(X, Xf)
    }
    for (r in seq_along(nr)) {
      Zr <- Xs[, nr[[r]], drop = FALSE]
      blocks[[length(blocks) + 1L]] <- list(
        label = paste0(key, if (length(nr) > 1) paste0(".", r)), group = key,
        terms = "smooth", Xr = Zr, gi = NULL, levels = paste0("b", seq_along(nr[[r]])),
        q = 1, L = length(nr[[r]]), corr = FALSE, struct = "diag", type = "smooth", par = par)
    }
    smooth_info[[key]] <- spec
  }
  list(X = X, offset = off, terms = tt, xlev = stats::.getXlevels(tt, mf),
       blocks = blocks, smooths = smooth_info, formula = rhs)
}

# Sparse Z matrix of one block: n x (q * L), columns ordered level major
# (all q effects of level 1, then level 2, ...).
block_Z <- function(b, n) {
  if (b$type == "smooth") return(Matrix::Matrix(b$Xr, sparse = TRUE))
  if (!is.null(b$W)) {
    tw <- Matrix::summary(b$W)
    return(Matrix::sparseMatrix(i = rep(tw$i, b$q), j = as.vector(outer((tw$j - 1L) * b$q, seq_len(b$q), `+`)),
      x = as.vector(tw$x * b$Xr[tw$i, , drop = FALSE]), dims = c(n, b$q * b$L)))
  }
  ok <- !is.na(b$gi)
  i <- rep(which(ok), b$q)
  j <- as.vector(outer((b$gi[ok] - 1L) * b$q, seq_len(b$q), `+`))
  x <- as.vector(b$Xr[ok, , drop = FALSE])
  Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(n, b$q * b$L))
}

# Parse the user's formulas into one formula per regression parameter.
parse_formulas <- function(formula, dots, fam) {
  rp <- regression_pars(fam)
  if (!length(rp)) stop_unsupported("Regression", fam$id)
  fl <- list()
  if (!inherits(formula, "formula") || length(formula) != 3)
    stop_input("The first formula must be two sided: response ~ predictors.")
  fl[[rp[1]]] <- formula
  for (f in dots) {
    if (length(f) != 3 || !is.name(f[[2]]))
      stop_input("Formulas for further parameters must look like 'sigma ~ x'.")
    nm <- as.character(f[[2]])
    if (!nm %in% names(fam$pars))
      stop_input(sprintf("'%s' is not a parameter of '%s'. Parameters: %s.", nm, fam$id,
                         paste(names(fam$pars), collapse = ", ")))
    if (!nm %in% rp)
      stop_input(sprintf("Parameter '%s' of '%s' cannot have a predictor.", nm, fam$id))
    fl[[nm]] <- f
  }
  for (k in rp) if (is.null(fl[[k]])) fl[[k]] <- if (isTRUE(fam$extra$share_formula))
    stats::as.formula(call("~", formula[[3]]), env = environment(formula)) else
    stats::as.formula("~ 1", env = environment(formula))
  attr(fl, "given") <- c(rp[1], names(fl)[vapply(names(fl), function(k) k != rp[1] && k %in% vapply(dots, function(f) as.character(f[[2]]), ""), TRUE)])
  structure(fl[rp], given = attr(fl, "given"))
}

# Expand `.` in a one sided formula to every data column except the
# response variables.
expand_dot <- function(f, data, response_vars) {
  if (!"." %in% all.vars(f) || is.null(data)) return(f)
  vars <- setdiff(names(data), response_vars)
  rhs <- paste(deparse(f[[length(f)]], width.cutoff = 500L), collapse = " ")
  rhs <- gsub("(^|[^A-Za-z0-9_.])\\.($|[^A-Za-z0-9_.])", paste0("\\1(", paste(sprintf("`%s`", vars), collapse = " + "), ")\\2"), rhs)
  stats::as.formula(paste("~", rhs), env = environment(f))
}

rhs_of <- function(f) {
  if (length(f) == 3) {
    out <- f[-2]
  } else out <- f
  out
}

# Response extraction (numeric vector, cens matrix or matrix for vector
# valued families).
get_response <- function(formula, data) {
  lhs <- formula[[2]]
  eval(lhs, data, environment(formula))
}
