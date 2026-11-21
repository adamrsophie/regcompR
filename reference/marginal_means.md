# Marginal means (average predictions)

The prediction averaged over rows, overall or with every row set to each
value of a covariate (counterfactual or g-computation marginal means, as
`avg_predictions(variables = )` of marginaleffects). Unlike the
estimated marginal means of emmeans, which average over a balanced
reference grid, these average over the observed rows.

## Usage

``` r
marginal_means(
  object,
  variable = NULL,
  at = NULL,
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

- variable:

  Covariate to set, or `NULL` for the overall average.

- at:

  Values of the covariate; its levels (or mean and the mean plus and
  minus one standard deviation for a numeric covariate) by default.

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

A data frame with one row per value.

## Examples

``` r
d <- data.frame(x = rnorm(300), g = factor(sample(c("a", "b"), 300, TRUE)))
d$y <- rpois(300, exp(0.2 + 0.5 * d$x + 0.4 * (d$g == "b")))
marginal_means(regcomp(y ~ x + g, family = "poisson", data = d), "g")
#>   g estimate std_error statistic      p_value    lower    upper
#> 1 a 1.619548 0.1055833  15.33905 4.192763e-53 1.412608 1.826487
#> 2 b 2.212646 0.1205646  18.35237 3.160093e-75 1.976344 2.448948
```
