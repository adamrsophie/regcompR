# Variance inflation factors

Generalized variance inflation factors (Fox and Monette 1992) of the
terms of each distribution parameter's regression, from the correlations
of the coefficient estimates other than the intercept (as in
`car::vif()`). Values above about 5 or 10 indicate strong collinearity.
`gvif_adj` is `GVIF^(1 / (2 df))`, comparable across terms with several
columns.

## Usage

``` r
vif(object)
```

## Arguments

- object:

  A fitted model.

## Value

A data frame with one row per term and parameter.

## Examples

``` r
d <- data.frame(x = rnorm(100)); d$z <- d$x + rnorm(100, 0, 0.3); d$y <- rnorm(100)
vif(regcomp(y ~ x + z, family = "normal", data = d))
#>   parameter term df     gvif gvif_adj
#> 1        mu    x  1 11.05067 3.324255
#> 2        mu    z  1 11.05067 3.324255
```
