# Simple slopes under moderation

The slope of a numeric covariate `x` on the linear predictor of one
distribution parameter at chosen values of a moderator, averaged over
the analysis rows with the moderator set to each value (the average
marginal effect on the link scale), with delta method standard errors
and Wald tests. With a linear interaction `x + z + x:z` this is
`b_x + b_xz z`. Random effects do not enter the slope.

## Usage

``` r
simple_slopes(object, x, moderator, at = NULL, par = NULL, vcov = "model", ...)
```

## Arguments

- object:

  A frequentist fit.

- x:

  Name of the focal numeric covariate.

- moderator:

  Name of the moderator.

- at:

  Values of the moderator; by default its mean and the mean plus and
  minus one standard deviation (its levels for a factor).

- par:

  Distribution parameter; the first regression parameter by default.

- vcov, ...:

  Covariance choice passed to
  [`vcov.regcomp_fit()`](https://adamrsophie.github.io/regcompR/reference/vcov.regcomp_fit.md).

## Value

A data frame with one row per moderator value.

## Examples

``` r
d <- data.frame(x = rnorm(300), z = rnorm(300))
d$y <- rpois(300, exp(0.5 + 0.3 * d$x + 0.2 * d$x * d$z))
simple_slopes(regcomp(y ~ x * z, family = "poisson", data = d), "x", "z")
#>             z     slope  std_error statistic      p_value
#> 1 -1.13100134 0.2142504 0.06307503  3.396754 6.819015e-04
#> 2 -0.08128138 0.3122955 0.04111748  7.595201 3.073155e-14
#> 3  0.96843858 0.4103406 0.05216571  7.866099 3.658718e-15
```
