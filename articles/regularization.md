# Penalized and bias reduced estimation

``` r

library(regcompR)
```

[`penalty()`](https://adamrsophie.github.io/regcompR/reference/penalty.md)
builds a penalty on the regression coefficients for the `"ml"` engine;
the estimate maximizes `loglik(beta) - P(beta)`. Intercepts are never
penalized except by `"firth"`. `pars =` restricts the penalty to some
distribution parameters, and `standardize = TRUE` (the default for the
shrinkage penalties) scales each column by its standard deviation.

## Bias reduction and separation

Firth’s penalty (`"firth"`, also `"jeffreys"`) removes the first order
bias of the maximum likelihood estimator for canonical links and gives
finite estimates under separation. `"flic"` re-estimates the intercept
so that the mean predicted probability equals the observed proportion,
and `"flac"` fits maximum likelihood on data augmented with Firth’s
pseudo observations (Puhr et al. 2017). `"log_f"` (log F(m, m) priors)
and `"cauchy"` (independent Cauchy priors with scale `scale`) are weakly
informative alternatives.

``` r

sep <- data.frame(x = c(-2, -1, -0.5, 0.5, 1, 2), y = c(0, 0, 0, 1, 1, 1))
coef(suppressWarnings(regcomp(y ~ x, family = "bernoulli", data = sep)))   # diverges
#> prob:(Intercept)           prob:x 
#>    -1.084897e-04     1.168769e+02
fits <- lapply(c("firth", "flic", "flac", "log_f", "cauchy"), function(t)
  coef(regcomp(y ~ x, family = "bernoulli", data = sep, penalty = penalty(t))))
do.call(rbind, setNames(fits, c("firth", "flic", "flac", "log_f", "cauchy")))
#>        prob:(Intercept)   prob:x
#> firth     -7.071069e-12 1.583008
#> flic      -7.071069e-12 1.583008
#> flac      -4.957817e-12 1.583008
#> log_f     -1.449355e-16 2.436251
#> cauchy     3.662057e-16 2.504253
```

Grouped binomial data give the same Firth, FLIC and FLAC estimates as
the equivalent Bernoulli rows:

``` r

g <- data.frame(x = c(-2, -1, 0, 1, 2), y = c(0, 1, 2, 6, 3), n = c(3, 7, 5, 10, 4))
coef(regcomp(y ~ x, family = "binomial", size = n, data = g, penalty = penalty("flac")))
#> prob:(Intercept)           prob:x 
#>       -0.6370961        0.9367984
```

For `"firth"`, `"log_f"` and `"cauchy"` the inverse penalized Hessian is
the usual large sample covariance.

## Shrinkage

- `"ridge"`: `lambda / 2 * sum(beta^2)`;
- `"lasso"`: `lambda * sum(|beta|)`;
- `"elastic_net"`:
  `lambda * (alpha * sum(|beta|) + (1 - alpha) / 2 * sum(beta^2))`;
- `"adaptive_lasso"`: weights `1 / |b|^gamma` from an initial fit, or
  `weights =`;
- `"group_lasso"`: `lambda * sum_g sqrt(p_g) ||beta_g||` over `groups =`
  (by default the columns of one model term form a group);
- `"scad"` and `"mcp"`: folded concave penalties fitted by local linear
  approximation.

`lambda` is on the scale of the summed log likelihood:
`lambda = n * lambda_glmnet` reproduces glmnet and ncvreg.

``` r

n <- 100
X <- matrix(rnorm(n * 8), n, 8, dimnames = list(NULL, paste0("x", 1:8)))
d <- data.frame(X); d$y <- 1 + 2 * d$x1 - d$x2 + rnorm(n)
types <- c("ridge", "lasso", "elastic_net", "adaptive_lasso", "scad", "mcp")
est <- sapply(types, function(t) coef(regcomp(y ~ ., family = "normal", data = d,
  fixed = list(sigma = 1), penalty = penalty(t, lambda = 20, alpha = if (t == "elastic_net") 0.5))))
round(est, 3)
#>                 ridge  lasso elastic_net adaptive_lasso   scad    mcp
#> mu:(Intercept)  1.053  1.031       1.035          1.022  0.999  0.998
#> mu:x1           1.741  1.860       1.802          1.964  2.084  2.086
#> mu:x2          -0.824 -0.787      -0.805         -0.778 -0.998 -0.997
#> mu:x3          -0.008  0.000       0.000          0.000  0.000  0.000
#> mu:x4           0.214  0.041       0.132          0.000  0.040  0.060
#> mu:x5           0.077  0.000       0.000          0.000  0.000  0.000
#> mu:x6           0.144  0.000       0.059          0.000  0.000  0.000
#> mu:x7           0.047  0.000       0.000          0.000  0.000  0.000
#> mu:x8          -0.029  0.000       0.000          0.000  0.000  0.000
fg <- regcomp(y ~ ., family = "normal", data = d, fixed = list(sigma = 1),
              penalty = penalty("group_lasso", lambda = 30, groups = rep(c("a", "b"), each = 4)))
round(coef(fg), 3)
#> mu:(Intercept)          mu:x1          mu:x2          mu:x3          mu:x4 
#>          1.068          1.493         -0.707          0.003          0.166 
#>          mu:x5          mu:x6          mu:x7          mu:x8 
#>          0.000          0.000          0.000          0.000
```

Penalties apply to any distribution parameter; here only the dispersion
model of a negative binomial regression is shrunk:

``` r

dn <- data.frame(X); dn$y <- rnbinom(n, mu = exp(1 + 0.5 * dn$x1), size = 2)
fpar <- regcomp(y ~ x1, phi ~ ., family = "nb2", data = dn,
                penalty = penalty("lasso", lambda = 10, pars = "phi"))
round(coef(fpar), 3)
#>  mu:(Intercept)           mu:x1 phi:(Intercept)          phi:x1          phi:x2 
#>           0.995           0.406           0.456           0.000           0.000 
#>          phi:x3          phi:x4          phi:x5          phi:x6          phi:x7 
#>           0.000           0.000           0.000           0.000           0.000 
#>          phi:x8 
#>           0.000
```

### Diagnostics and controls

Lasso type penalties are fitted by accelerated proximal gradient
descent. Their optimum is not a zero of the gradient, so the fit reports
the proximal stationarity residual, which is zero at the minimizer. The
iteration limit and tolerance of this solver are `prox_maxit` and
`prox_tol` in
[`regcomp_control()`](https://adamrsophie.github.io/regcompR/reference/regcomp_control.md).

``` r

fl <- regcomp(y ~ ., family = "normal", data = d, fixed = list(sigma = 1),
              penalty = penalty("lasso", lambda = 20),
              control = regcomp_control(prox_maxit = 10000, prox_tol = 1e-10))
fl$diagnostics$kkt_residual
#> [1] 2.569308e-08
fl$optimizer$objective + as.numeric(logLik(fl))    # the stored penalty value
#> [1] 49.3062
fl$penalty$value(unname(coef(fl)))
#> [1] 49.3062
```

### More coefficients than observations

A shrinkage penalty identifies designs with more columns than rows or
with aliased columns; unpenalized fits refuse them.

``` r

wide <- data.frame(y = rnorm(10), matrix(rnorm(10 * 20), 10, 20))
fw <- regcomp(y ~ ., family = "normal", data = wide, fixed = list(sigma = 1),
              penalty = penalty("ridge", lambda = 2))
length(coef(fw))
#> [1] 21
try(regcomp(y ~ ., family = "normal", data = wide, fixed = list(sigma = 1)))
#> Error in build_spec(formula, mr$dots, fam, mr$data, weights = mr$weights,  : 
#>   The design for 'mu' is rank deficient; aliased columns: X10, X11, X12, X13, X14, X15, X16, X17, X18, X19, X20. Remove them from the formula, or use a shrinkage penalty (ridge, lasso and similar) with the ml engine.
```

### Choosing lambda

Choose `lambda` by cross validation:

``` r

lams <- c(5, 20, 60)
cvs <- sapply(lams, function(l) cv_regcomp(regcomp(y ~ ., family = "normal", data = d,
  fixed = list(sigma = 1), penalty = penalty("lasso", lambda = l)), K = 5)$estimate)
setNames(cvs, lams)
#>        5       20       60 
#> 1.519753 1.597449 2.118530
```

Penalized estimates are not maximum likelihood estimates; for the
shrinkage penalties p values are suppressed. Smooth penalties (ridge,
Firth, log F, Cauchy) also work with random effects and smooth terms;
the lasso type penalties need models without them.
