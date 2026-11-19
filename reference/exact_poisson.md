# Exact conditional Poisson regression

Conditional inference for one coefficient of a Poisson log linear model
(Cox 1970; Hirji 2006), as in Stata's `expoisson` discussed by Hilbe
(2014, sec. 9.4). The counts are conditioned on their totals within the
strata, which removes the intercept and the stratum effects; the
conditional distribution of `T = sum(x * y)` is computed exactly by
polynomial convolution. Returns the conditional maximum likelihood
estimate (which equals the unconditional Poisson estimate with stratum
intercepts when it is finite), the median unbiased estimate when `T` is
at an end of its support, the exact two sided p value of `beta = 0`
(probability method), and the exact confidence interval.

## Usage

``` r
exact_poisson(y, x, strata = NULL, exposure = NULL, level = 0.95)
```

## Arguments

- y:

  Counts.

- x:

  Integer covariate whose coefficient is tested.

- strata:

  Optional factor of strata (nuisance effects).

- exposure:

  Optional exposure (offset on the rate scale).

- level:

  Confidence level.

## Value

A list with `estimate`, `type` ("CMLE" or "MUE"), `p_value`, `conf_int`,
`statistic` (observed `T`) and the conditional distribution `support`
and `counts` (log weights).

## Examples

``` r
exact_poisson(y = c(0, 1, 3, 4), x = c(0, 0, 1, 1), exposure = c(10, 12, 9, 11))
#> $estimate
#> [1] 2.04122
#> 
#> $type
#> [1] "CMLE"
#> 
#> $p_value
#> [1] 0.0315777
#> 
#> $conf_int
#> [1] -0.01082803  5.84941620
#> 
#> $level
#> [1] 0.95
#> 
#> $statistic
#> [1] 7
#> 
#> $support
#> [1] 0 1 2 3 4 5 6 7 8
#> 
#> $log_weights
#> [1] -5.173017 -3.188886 -2.031433 -1.433596 -1.305763 -1.624217 -2.412674
#> [8] -3.760747 -5.935499
#> 
```
