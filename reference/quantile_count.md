# Quantile regression for counts

The jittering method of Machado and Santos Silva (2005), discussed by
Hilbe (2014, sec. 9.5): `z = y + U` with uniform noise, the transform
`log(z - tau)` (with a small floor below `tau`) is fitted by linear
quantile regression, and estimates are averaged over `m` jittered
samples. The conditional quantile of the count is
`ceiling(tau + exp(x' beta)) - 1`.

## Usage

``` r
quantile_count(formula, data, tau = 0.5, m = 100, seed = 1)
```

## Arguments

- formula:

  Formula `count ~ predictors`.

- data:

  Data frame.

- tau:

  Quantile level in (0, 1).

- m:

  Number of jittered samples.

- seed:

  Seed for the jittering.

## Value

A list with `coefficients`, their `vcov` (between plus mean within
jitter variance), `tau` and a `predict(newdata)` function giving the
count quantiles.

## Examples

``` r
if (requireNamespace("quantreg", quietly = TRUE)) {
  set.seed(1)
  d <- data.frame(x = rnorm(300))
  d$y <- rpois(300, exp(1 + 0.4 * d$x))
  q <- quantile_count(y ~ x, d, tau = 0.5, m = 20)
  q$coefficients
}
#> (Intercept)           x 
#>   0.9511749   0.4049839 
```
