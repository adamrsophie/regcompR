# Pointwise log likelihood

Returns the log likelihood contribution of each observation (including
censoring and truncation). For Bayesian fits returns a draws by
observations matrix suitable for 'loo'.

## Usage

``` r
log_lik(object, newdata = NULL, re = c("conditional", "zero"), ...)
```

## Arguments

- object:

  A fitted model.

- newdata:

  Optional data frame containing the response.

- re:

  `"conditional"` uses the fitted random effects, `"zero"` sets them to
  zero.

- ...:

  Unused.

## Value

Numeric vector or matrix.
