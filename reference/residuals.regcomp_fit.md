# Residuals of a fitted model

`type = "quantile"` gives randomized quantile residuals (Dunn and Smyth
1996, <doi:10.1080/10618600.1996.10474708>): the PIT value of each
observation under its fitted distribution, randomized uniformly over the
probability mass of discrete or censored observations, mapped to the
standard normal scale. Under a correct model they are independent
standard normal. `type = "pit"` returns the PIT values themselves,
`"response"` the difference between the observation and the fitted mean,
and `"pearson"` that difference divided by the fitted standard
deviation.

## Usage

``` r
# S3 method for class 'regcomp_fit'
residuals(
  object,
  type = c("quantile", "pit", "response", "pearson"),
  seed = NULL,
  ...
)
```

## Arguments

- object:

  A fitted model.

- type:

  Residual type.

- seed:

  Optional seed for the randomization.

- ...:

  Unused.

## Value

Numeric vector.
