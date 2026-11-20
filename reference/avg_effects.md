# Average marginal effects and contrasts

For a numeric covariate, the average over the rows of the derivative of
the predicted quantity with respect to the covariate (the average
marginal effect, or average slope); for a factor, logical or character
covariate, the average over the rows of the difference between the
predictions with every row set to each level and with every row set to
the reference level (g-computation; the average treatment effect under
the usual assumptions). These are the `avg_slopes()` and
`avg_comparisons()` of marginaleffects and the average marginal effects
of Stata's `margins, dydx()`.

## Usage

``` r
avg_effects(
  object,
  variables = NULL,
  type = c("mean", "parameter", "link"),
  par = NULL,
  newdata = NULL,
  re = c("conditional", "zero"),
  level = 0.95,
  vcov = "model",
  ...
)
```

## Arguments

- object:

  A fitted model.

- variables:

  Covariates; by default every covariate of the model.

- type:

  `"mean"` (the predicted mean), `"parameter"` (a distribution
  parameter) or `"link"` (its linear predictor).

- par:

  Distribution parameter for `type = "parameter"` or `"link"`.

- newdata:

  Rows to average over; the analysis rows by default.

- re:

  `"conditional"` (with the estimated random effects) or `"zero"`
  (random effects set to zero).

- level:

  Interval level.

- vcov, ...:

  Covariance choice for frequentist fits (see
  [`vcov.regcomp_fit()`](https://adamrsophie.github.io/regcompR/reference/vcov.regcomp_fit.md)).

## Value

A data frame with one row per covariate (or contrast).

## Details

Frequentist fits use the delta method with the chosen covariance (random
effects are held at their estimates); Bayesian fits summarize the effect
over the posterior draws.

## Examples

``` r
d <- data.frame(x = rnorm(300), g = factor(sample(c("a", "b", "c"), 300, TRUE)))
d$y <- rbinom(300, 1, plogis(-0.5 + 0.8 * d$x + (d$g == "b")))
f <- regcomp(y ~ x + g, family = "bernoulli", data = d)
avg_effects(f)
#>   term contrast   estimate  std_error statistic      p_value      lower
#> 1    x    dY/dX  0.1572576 0.02594746  6.060616 1.356008e-09  0.1064015
#> 2    g    b - a  0.2472136 0.06423042  3.848856 1.186709e-04  0.1213243
#> 3    g    c - a -0.1056021 0.06346644 -1.663904 9.613154e-02 -0.2299940
#>        upper
#> 1 0.20811367
#> 2 0.37310294
#> 3 0.01878985
```
