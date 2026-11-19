# Stable numerical building blocks

Log scale helpers used throughout the package. `log1mexp(x)` computes
`log(1 - exp(-x))` for `x > 0` (Maechler 2012); `log1pexp(x)` computes
`log(1 + exp(x))`; `log_sum_exp()` and `log_diff_exp()` combine log
scale values without overflow.

## Usage

``` r
log1mexp(x)

log1pexp(x)

log_sum_exp(x)

log_diff_exp(a, b)
```

## Arguments

- x, a, b:

  Numeric vectors.

## Value

Numeric vector.

## Examples

``` r
log1mexp(1e-20)
#> [1] -46.0517
log1pexp(800)
#> [1] 800
log_sum_exp(c(1000, 1000))
#> [1] 1000.693
log_diff_exp(log(3), log(2))
#> [1] 2.220446e-16
```
