# Predictions from a fitted model

Predictions from a fitted model

## Usage

``` r
# S3 method for class 'regcomp_fit'
predict(
  object,
  newdata = NULL,
  type = c("mean", "parameter", "link", "variance", "sd", "median", "quantile",
    "probability", "density", "cdf", "survival"),
  probs = c(0.1, 0.5, 0.9),
  at = 0,
  interval = c("none", "confidence"),
  level = 0.95,
  re = c("conditional", "zero", "marginal"),
  ndraws = 1000,
  nmarginal = 200,
  seed = NULL,
  ...
)
```

## Arguments

- object:

  A fitted model.

- newdata:

  Optional data frame.

- type:

  What to predict: `"parameter"` (all distribution parameters), `"link"`
  (linear predictors), `"mean"`, `"variance"`, `"sd"`, `"median"`,
  `"quantile"` (at `probs`), `"probability"` (mass or density at `at`),
  `"density"`, `"cdf"` or `"survival"` (at `at`). For ordinal and
  nominal families `"mean"` and `"probability"` give the matrix of
  category probabilities (all categories by default).

- probs:

  Probabilities for `type = "quantile"`.

- at:

  Response values for `"probability"`, `"density"`, `"cdf"` and
  `"survival"`.

- interval:

  `"none"`, or `"confidence"` for uncertainty of the predicted quantity.
  Frequentist fits use parametric simulation from the asymptotic normal
  distribution of the estimates (`ndraws` draws); Bayesian fits use
  posterior draws and return credible intervals.

- level:

  Interval level.

- re:

  `"conditional"`, `"zero"` or `"marginal"`. `"marginal"` averages the
  quantity over the random effect distribution by Monte Carlo
  integration (`nmarginal` draws), giving population averaged
  predictions.

- ndraws:

  Number of simulation draws for confidence intervals.

- nmarginal:

  Number of random effect draws for `re = "marginal"`.

- seed:

  Optional seed for the simulation.

- ...:

  Unused.

## Value

A numeric vector, a matrix (several `probs` or `at` values, or
`type = "parameter"`), or a data frame with `estimate`, `lower` and
`upper` when an interval is requested.
