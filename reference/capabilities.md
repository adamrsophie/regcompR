# Capability queries

`capabilities()` reports, for a family or distribution, how each
operation is provided: `"analytic"` (a method specific to the family: a
closed form, or a dedicated algorithm such as a series, recursion or
quadrature written for that law), `"numeric"` (the generic numerical
fallback: inversion, summation or integration of the density, with
tolerance checks) or `"unavailable"`. Ordinal and nominal families draw
random values from their class probabilities. It also reports whether
the family can be used as a regression response and which fitting
engines support it. `supports()` answers a single question with `TRUE`
or `FALSE`.

## Usage

``` r
capabilities(x)

supports(x, what = NULL, backend = NULL)
```

## Arguments

- x:

  A family, distribution, identifier or fitted model.

- what:

  Operation name: `"log_prob"`, `"cdf"`, `"quantile"`, `"rng"`,
  `"mean"`, `"variance"`, `"skewness"`, `"kurtosis"`, `"score"`,
  `"hessian"` or `"regression"`.

- backend:

  Engine name: `"ml"`, `"tmb"`, `"stan"` (also `"bayes"`) or `"bugs"`
  (also `"jags"`, `"openbugs"`).

## Value

`capabilities()`: a named list. `supports()`: a logical scalar.

## Examples

``` r
capabilities("conway_maxwell_poisson")
#> $log_prob
#> [1] "analytic"
#> 
#> $cdf
#> [1] "analytic"
#> 
#> $quantile
#> [1] "numeric"
#> 
#> $rng
#> [1] "numeric"
#> 
#> $mean
#> [1] "analytic"
#> 
#> $variance
#> [1] "analytic"
#> 
#> $skewness
#> [1] "numeric"
#> 
#> $kurtosis
#> [1] "numeric"
#> 
#> $score
#> [1] "numeric"
#> 
#> $hessian
#> [1] "numeric"
#> 
#> $regression
#> [1] TRUE
#> 
#> $ml
#> [1] TRUE
#> 
#> $stan
#> [1] "generated"
#> 
#> $tmb
#> [1] "unavailable"
#> 
#> $bugs
#> [1] "unavailable"
#> 
#> $status
#> [1] "validated"
#> 
supports("cauchy", "mean")
#> [1] TRUE
supports("poisson", backend = "stan")
#> [1] TRUE
```
