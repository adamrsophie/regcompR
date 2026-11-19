# Pointwise log likelihood

Returns the log likelihood contribution of each observation (including
censoring and truncation). For Bayesian fits returns a draws by
observations matrix suitable for 'loo'.

## Usage

``` r
log_lik(object, newdata = NULL, ...)
```

## Arguments

- object:

  A fitted model.

- newdata:

  Optional data frame containing the response.

- ...:

  Unused.

## Value

Numeric vector or matrix.
