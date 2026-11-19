# Probability interface

Evaluate a distribution object. `log_prob()` returns the log density (or
log mass for discrete laws), `prob()` the density or mass, `cdf()` the
distribution function, `survival()` the survival function, and `rng()`
random draws. Values outside the support give log density `-Inf`.
Methods without a closed form use a declared numerical fallback
(quadrature, summation or inversion with tolerance checks); when none is
justified a `regcomp_unsupported` error is raised.

## Usage

``` r
log_prob(x, y, ...)

prob(x, y, ...)

cdf(x, q, ...)

# S3 method for class 'regcomp_dist'
cdf(x, q, lower_tail = TRUE, log_p = FALSE, ...)

survival(x, q, ...)

# S3 method for class 'regcomp_dist'
survival(x, q, log_p = FALSE, ...)

rng(x, n, ...)

# S3 method for class 'regcomp_dist'
rng(x, n = NULL, ...)
```

## Arguments

- x:

  A distribution object.

- y, q:

  Response values or quantiles.

- ...:

  Unused.

- lower_tail:

  Logical; `FALSE` gives upper tail probabilities.

- log_p:

  Logical; return log probabilities.

- n:

  Number of draws. With vectorized parameters of length `m > 1`, `n`
  must equal `m` or be omitted.

## Value

Numeric vector.

## Examples

``` r
d <- distribution("gamma", mu = 2, sigma = 0.5)
log_prob(d, 1:3)
#> [1] -1.0191707 -0.9397292 -1.7233339
prob(d, 1)
#> [1] 0.3608941
cdf(d, 2)
#> [1] 0.5665299
survival(d, 2, log_p = TRUE)
#> [1] -0.8359324
set.seed(1); rng(d, 5)
#> [1] 1.213062 3.214959 3.142633 2.159352 3.104412
```
