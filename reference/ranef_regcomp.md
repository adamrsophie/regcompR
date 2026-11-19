# Random effect estimates

Returns the conditional modes (frequentist fits) or posterior means
(Bayesian fits) of the random effects, one data frame per term, and the
estimated standard deviations and correlations.

## Usage

``` r
ranef_regcomp(object, ...)
```

## Arguments

- object:

  A fitted model.

- ...:

  Unused.

## Value

A list with `effects` (list of data frames) and `varcorr`.

## Examples

``` r
set.seed(3)
g <- gl(20, 10)
d <- data.frame(g, x = rnorm(200))
d$y <- rpois(200, exp(0.3 + 0.4 * d$x + rnorm(20, 0, 0.5)[g]))
fit <- regcomp(y ~ x + (1 | g), family = "poisson", data = d)
ranef_regcomp(fit)$varcorr
#> $`1 | g`
#> $`1 | g`$sd
#> [1] 0.5321327
#> 
#> $`1 | g`$cor
#>      [,1]
#> [1,]    1
#> 
#> 
```
