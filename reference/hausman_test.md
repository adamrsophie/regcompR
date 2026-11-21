# Hausman test

Compares a consistent estimator (fixed effects, instrumental variables)
with an estimator that is efficient under the null hypothesis (random
effects, least squares) over their common coefficients:
`(b1 - b0)' (V1 - V0)^- (b1 - b0)`, chi squared with the rank of
`V1 - V0` degrees of freedom. Equal to `plm::phtest()` given the same
fits. Works with regcompR fits (the first parameter's coefficients, by
term),
[`panel_within()`](https://adamrsophie.github.io/regcompR/reference/panel_within.md)
and
[`iv_regcomp()`](https://adamrsophie.github.io/regcompR/reference/iv_regcomp.md)
results and other fits with
[`coef()`](https://rdrr.io/r/stats/coef.html) and
[`vcov()`](https://rdrr.io/r/stats/vcov.html) methods.

## Usage

``` r
hausman_test(consistent, efficient)
```

## Arguments

- consistent, efficient:

  Fitted models.

## Value

A one row data frame with statistic, df and p value.

## Examples

``` r
set.seed(1)
d <- data.frame(id = factor(rep(1:50, each = 5)), x = rnorm(250))
d$y <- 1 + 0.5 * d$x + rnorm(50)[d$id] + rnorm(250)
fe <- panel_within(y ~ x, data = d, id = "id")
re <- regcomp(y ~ x + (1 | id), family = "normal", data = d)
hausman_test(fe, re)
#>   statistic df   p_value
#> 1 0.3925989  1 0.5309365
```
