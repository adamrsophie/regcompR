# Simulate responses from a fitted model

Draws new responses from the fitted conditional distributions. Bayesian
fits draw from the posterior predictive distribution (one posterior draw
per simulated data set).

## Usage

``` r
# S3 method for class 'regcomp_fit'
simulate(
  object,
  nsim = 1,
  seed = NULL,
  newdata = NULL,
  re = c("conditional", "new"),
  ...
)
```

## Arguments

- object:

  A fitted model.

- nsim:

  Number of simulated data sets.

- seed:

  Optional seed.

- newdata:

  Optional data frame.

- re:

  `"conditional"` keeps the estimated random effects; `"new"` draws new
  random effects from their estimated distribution.

- ...:

  Unused.

## Value

A matrix with one column per simulation.
