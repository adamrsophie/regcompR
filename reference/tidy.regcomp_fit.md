# Tidy summaries of a fitted model

Methods for the generics `tidy()`, `glance()` and `augment()` of the
generics package (also used by broom, modelsummary and gtsummary).
`tidy()` returns one row per coefficient with the distribution parameter
it belongs to (`component`), and one row per random effect variance
parameter (`effect = "ran_pars"`). Random effect standard deviations are
reported on their natural scale (`sd(...)`). For Bayesian fits the
estimate is the posterior mean, `std.error` the posterior standard
deviation, the interval a credible interval, and `rhat` and `ess` are
added.

## Usage

``` r
# S3 method for class 'regcomp_fit'
tidy(
  x,
  conf.int = TRUE,
  conf.level = 0.95,
  exponentiate = FALSE,
  vcov = "model",
  standardize = FALSE,
  ...
)

# S3 method for class 'regcomp_fit'
glance(x, ...)

# S3 method for class 'regcomp_fit'
augment(x, data = NULL, newdata = NULL, ...)
```

## Arguments

- x:

  A fitted model from
  [`regcomp()`](https://adamrsophie.github.io/regcompR/reference/regcomp.md).

- conf.int:

  Add confidence (credible) intervals.

- conf.level:

  Interval level.

- exponentiate:

  Exponentiate coefficients and their intervals (odds ratios, rate
  ratios); random effect parameters are not changed.

- vcov:

  `"model"`, `"sandwich"` or `"quasi"` standard errors (see
  [`summary.regcomp_fit()`](https://adamrsophie.github.io/regcompR/reference/summary.regcomp_fit.md)).

- standardize:

  Multiply each coefficient (and its interval) by the standard deviation
  of its design column: the change of the linear predictor per standard
  deviation of the covariate (partially standardized coefficients;
  meaningful for any family, unlike scaling the response). Applied
  before `exponentiate`.

- ...:

  Unused.

- data, newdata:

  For `augment()`: data to add columns to; defaults to the analysis rows
  of the fit.

## Value

A data frame.

## Examples

``` r
d <- data.frame(x = rnorm(100)); d$y <- rpois(100, exp(0.5 + 0.3 * d$x))
f <- regcomp(y ~ x, family = "poisson", data = d)
if (requireNamespace("generics", quietly = TRUE)) {
  print(generics::tidy(f, exponentiate = TRUE))
  print(generics::glance(f))
}
#>          term component effect estimate std.error statistic      p.value
#> 1 (Intercept)    lambda  fixed 1.614615 0.1324938  5.838433 5.269420e-09
#> 2           x    lambda  fixed 1.521487 0.1276337  5.002996 5.644616e-07
#>   conf.low conf.high
#> 1 1.374739  1.896347
#> 2 1.290813  1.793385
#>    family engine nobs df    logLik      AIC      BIC converged
#> 1 Poisson     ml  100  2 -158.1277 320.2554 325.4657      TRUE
```
