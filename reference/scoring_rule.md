# Proper scoring rules for predictive distributions

Scores each observation under the fitted predictive distribution
(Gneiting and Raftery 2007). Smaller is better for every rule.

- `"log"`: negative log predictive density or mass.

- `"crps"`: continuous ranked probability score
  `int (F(z) - 1(y <= z))^2 dz`, by quadrature for continuous laws and
  by summation for counts.

- `"dss"`: Dawid-Sebastiani score `(y - m)^2 / v + log v`.

- `"brier"` and `"rps"`: Brier and ranked probability scores for ordinal
  and nominal families.

## Usage

``` r
scoring_rule(
  object,
  newdata = NULL,
  rule = c("log", "crps", "dss", "brier", "rps"),
  re = c("conditional", "zero")
)
```

## Arguments

- object:

  A fitted model.

- newdata:

  Optional data frame with the response; defaults to the fitting data.

- rule:

  Scoring rule.

- re:

  Random effects used for prediction, as in
  [`predict()`](https://rdrr.io/r/stats/predict.html): `"conditional"`
  uses the estimated effects of known groups, `"zero"` sets them to zero
  (new groups).

## Value

Numeric vector of scores, one per observation.

## Examples

``` r
set.seed(1)
d <- data.frame(x = rnorm(100)); d$y <- rnorm(100, 1 + d$x)
fit <- regcomp(y ~ x, family = "normal", data = d)
mean(scoring_rule(fit, rule = "crps"))
#> [1] 0.539673
```
