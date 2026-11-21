# Penalized quasi-likelihood for generalized linear mixed models

Fits a GLMM by penalized quasi-likelihood (Breslow and Clayton 1993):
repeated linear mixed model fits, by maximum likelihood, to the working
response `z = eta + (y - mu) deta/dmu` with variance weights
`(dmu/deta)^2 / V(mu)`, as
[`MASS::glmmPQL()`](https://rdrr.io/pkg/MASS/man/glmmPQL.html). PQL is
not a likelihood method (no log likelihood, AIC or likelihood ratio
tests); it is biased for binary data with few observations per cluster,
where the Laplace fit of
[`regcomp()`](https://adamrsophie.github.io/regcompR/reference/regcomp.md)
is preferable. It handles quasi families.

## Usage

``` r
pql_regcomp(formula, family, data, niter = 10, tol = 1e-06)
```

## Arguments

- formula:

  A formula with lme4 style random effects, `y ~ x + (1 | g)`.

- family:

  A [`stats::family()`](https://rdrr.io/r/stats/family.html) object
  (such as [`binomial()`](https://rdrr.io/r/stats/family.html),
  [`poisson()`](https://rdrr.io/r/stats/family.html),
  [`quasipoisson()`](https://rdrr.io/r/stats/family.html)).

- data:

  A data frame.

- niter:

  Maximum number of iterations.

- tol:

  Convergence tolerance on the linear predictor.

## Value

The final working normal mixed model fit (a `regcomp_fit`), whose `mu`
coefficients are the fixed effects on the link scale and whose `theta`
are the random effect parameters, with class `regcomp_pql` added and the
number of iterations in `iterations`.

## Examples

``` r
set.seed(1)
d <- data.frame(g = factor(rep(1:20, each = 8)), x = rnorm(160))
d$y <- rpois(160, exp(0.3 + 0.4 * d$x + rnorm(20, 0, 0.5)[d$g]))
pql_regcomp(y ~ x + (1 | g), family = poisson(), data = d)
#> <regcomp_pql> poisson (log link) GLMM by penalized quasi-likelihood, 5 iterations, n = 160
#>             estimate std_error      z
#> (Intercept)   0.5573    0.1309 4.2560
#> x             0.3785    0.0616 6.1499
#> Random effect parameters:
#> mu:1 | g:log_sd((Intercept)) 
#>                      -0.6645 
```
