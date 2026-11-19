# Validate parameters and responses

`validate_parameters()` signals a `regcomp_domain_error` when a value
lies outside its domain or violates a family constraint.
`validate_response()` returns a logical vector that is `TRUE` where a
response value lies in the support.

## Usage

``` r
validate_parameters(x, pars)

validate_response(x, y)
```

## Arguments

- x:

  A family, a distribution or an identifier.

- pars:

  Named list of parameter values.

- y:

  Response values.

## Value

`validate_parameters()` returns `pars` invisibly; `validate_response()`
returns a logical vector.

## Examples

``` r
validate_parameters("poisson", list(lambda = 2))
validate_response(distribution("poisson", lambda = 2), c(-1, 0, 1.5, 3))
#> [1] FALSE  TRUE FALSE  TRUE
```
