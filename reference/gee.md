# Generalized estimating equations

Fits a marginal mean model `g(E(y)) = X beta` with variance `phi V(mu)`
and a working correlation within clusters by the method of Liang and
Zeger (1986), with moment estimates of the correlation and scale and a
robust (sandwich) covariance. Quasi likelihood models (McCullagh and
Nelder 1989, ch. 9) are the independence case:
`gee(y ~ x, data, family = quasipoisson())` gives the quasi Poisson
estimates with a Pearson scale.

## Usage

``` r
gee(
  formula,
  data,
  id,
  family = stats::gaussian(),
  corstr = c("independence", "exchangeable", "ar1", "unstructured"),
  waves = NULL,
  weights = NULL,
  maxit = 50,
  tol = 1e-08
)
```

## Arguments

- formula:

  Mean model formula; may contain an offset.

- data:

  Data frame.

- id:

  Cluster identifier (unquoted column name or vector). Rows of a cluster
  are taken in data order unless `waves` is given.

- family:

  A `stats` family object such as
  [`gaussian()`](https://rdrr.io/r/stats/family.html),
  [`poisson()`](https://rdrr.io/r/stats/family.html),
  [`binomial()`](https://rdrr.io/r/stats/family.html),
  [`Gamma()`](https://rdrr.io/r/stats/family.html),
  [`inverse.gaussian()`](https://rdrr.io/r/stats/family.html),
  [`quasi()`](https://rdrr.io/r/stats/family.html),
  [`quasipoisson()`](https://rdrr.io/r/stats/family.html) or
  [`quasibinomial()`](https://rdrr.io/r/stats/family.html).

- corstr:

  Working correlation: `"independence"`, `"exchangeable"`, `"ar1"` or
  `"unstructured"`.

- waves:

  Optional within cluster time index (integers) for `"ar1"` and
  `"unstructured"`.

- weights:

  Optional prior weights.

- maxit, tol:

  Iteration limit and convergence tolerance.

## Value

An object of class `regcomp_gee` with `coefficients`, robust covariance
`vcov`, naive (model based) covariance `vcov_naive`, `alpha` (working
correlation parameters), `scale` and `fitted.values`.

## Examples

``` r
set.seed(1)
d <- data.frame(g = rep(1:40, each = 4), x = rnorm(160))
d$y <- rpois(160, exp(0.2 + 0.4 * d$x + rnorm(40, 0, 0.4)[d$g]))
gee(y ~ x, data = d, id = g, family = poisson(), corstr = "exchangeable")
#> <regcomp_gee> poisson(log), exchangeable working correlation, 40 clusters
#>              estimate  robust_se        z
#> (Intercept) 0.3373722 0.07980235 4.227597
#> x           0.5769809 0.09642471 5.983745
#>   correlation 0.0899  scale 1.3128
```
