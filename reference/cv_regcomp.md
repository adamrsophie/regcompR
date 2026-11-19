# Cross validation of a fitted model

Refits the model on training folds and scores the held out observations
with a proper scoring rule. With `group` the folds leave out whole
groups (leave group out cross validation), and held out groups are
predicted with their random effects at zero.

## Usage

``` r
cv_regcomp(object, K = 10, rule = "log", group = NULL, folds = NULL, seed = 1)
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

  Optional integer vector assigning each row to a fold.

- seed:

  Seed for the fold assignment.

## Value

A list with `scores` (per observation), `estimate` (mean score), `se`
and `folds`.

## Examples

``` r
set.seed(1)
d <- data.frame(x = rnorm(120)); d$y <- rpois(120, exp(0.5 + 0.3 * d$x))
fit <- regcomp(y ~ x, family = "poisson", data = d)
cv_regcomp(fit, K = 5)$estimate
#> [1] 1.690031
```
