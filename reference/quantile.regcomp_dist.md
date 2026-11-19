# Quantiles of a distribution object

For discrete laws returns the smallest support point `y` with
`F(y) >= p`. Without a closed form the CDF is inverted numerically with
bracketing; failures raise `regcomp_numerical_failure`.

## Usage

``` r
# S3 method for class 'regcomp_dist'
quantile(
  x,
  probs = c(0.025, 0.5, 0.975),
  lower_tail = TRUE,
  log_p = FALSE,
  ...
)
```

## Arguments

- x:

  A distribution object.

- probs:

  Probabilities.

- lower_tail, log_p:

  As in
  [`cdf()`](https://adamrsophie.github.io/regcompR/reference/probability-interface.md).

- ...:

  Unused.

## Value

Numeric vector.

## Examples

``` r
quantile(distribution("weibull", shape = 2, scale = 1), c(0.25, 0.5))
#> [1] 0.5363600 0.8325546
```
