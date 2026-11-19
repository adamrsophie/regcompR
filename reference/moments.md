# Moments of a distribution object

[`mean()`](https://rdrr.io/r/base/mean.html), `variance()`, `skewness()`
and `kurtosis()` (excess kurtosis) return `NA` with a
`regcomp_undefined_moment` warning where the moment does not exist for
the parameter values. Without closed forms, moments are computed by
quadrature or summation when the family allows it.

## Usage

``` r
# S3 method for class 'regcomp_dist'
mean(x, ...)

variance(x, ...)

# S3 method for class 'regcomp_dist'
variance(x, ...)

skewness(x, ...)

# S3 method for class 'regcomp_dist'
skewness(x, ...)

kurtosis(x, ...)

# S3 method for class 'regcomp_dist'
kurtosis(x, ...)
```

## Arguments

- x:

  A distribution object.

- ...:

  Unused.

## Value

Numeric vector.

## Examples

``` r
d <- distribution("student_t", mu = 0, sigma = 1, nu = c(1.5, 3, 10))
suppressWarnings(variance(d))
#> [1]   NA 3.00 1.25
```
