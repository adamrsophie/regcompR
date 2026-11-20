# Cross validation of a fitted model

Refits the model on training folds and scores the held out observations
with a proper scoring rule. With `group` the folds leave out whole
groups (leave group out cross validation), and held out groups are
scored as new groups: by default with the marginal predictive law, the
mixture over draws of their random effects.

## Usage

``` r
cv_regcomp(
  object,
  K = 10,
  rule = "log",
  group = NULL,
  folds = NULL,
  seed = 1,
  re = NULL,
  nmarginal = 200
)
```

## Arguments

- object:

  A fitted model from
  [`regcomp()`](https://adamrsophie.github.io/regcompR/reference/regcomp.md).

- K:

  Number of folds.

- rule:

  Scoring rule passed to
  [`scoring_rule()`](https://adamrsophie.github.io/regcompR/reference/scoring_rule.md).

- group:

  Optional grouping variable name; folds then hold out whole groups.

- folds:

  Optional integer vector assigning each analysis row (the
  [`nobs()`](https://rdrr.io/r/stats/nobs.html) rows used in the fit) to
  a fold; a vector with one entry per row of the original data is also
  accepted and is restricted to the analysis rows.

- seed:

  Seed for the fold assignment.

- re:

  Random effects for scoring held out rows, as in
  [`scoring_rule()`](https://adamrsophie.github.io/regcompR/reference/scoring_rule.md).
  By default `"marginal"` with `group` (held out groups are new groups)
  and `"conditional"` otherwise.

- nmarginal:

  Number of random effect draws for `re = "marginal"`.

## Value

A list with `scores` (per observation), `estimate` (mean score), `se`
and `folds`. With `group` the standard error is clustered by group
(scores of one group are dependent); otherwise rows are treated as
independent.

## Examples

``` r
set.seed(1)
d <- data.frame(x = rnorm(120)); d$y <- rpois(120, exp(0.5 + 0.3 * d$x))
fit <- regcomp(y ~ x, family = "poisson", data = d)
cv_regcomp(fit, K = 5)$estimate
#> [1] 1.690031
```
