# Generate the BUGS program of a model

Returns the BUGS language program run by the `"bayes"` engine with the
`"jags"` or `"openbugs"` runners, and its data list. Families with a
native BUGS distribution of the same parameterization use it directly;
other families with a closed form log density use the zeros trick (a
Poisson observation of zero with mean `C - log f(y)`), which gives the
same posterior. Families whose density needs series, special functions
outside the BUGS language or numerical integration are not available in
BUGS and raise `regcomp_unsupported`; use a 'Stan' runner for them.

## Usage

``` r
bugs_code(
  formula,
  ...,
  family,
  data,
  weights = NULL,
  links = NULL,
  fixed = NULL,
  priors = NULL,
  parameterization = NULL,
  dialect = c("jags", "openbugs")
)
```

## Arguments

- formula, ..., family, data, weights, links, fixed, parameterization:

  As in
  [`regcomp()`](https://adamrsophie.github.io/regcompR/reference/regcomp.md).

- priors:

  Priors from
  [`prior()`](https://adamrsophie.github.io/regcompR/reference/prior.md).

- dialect:

  `"jags"` or `"openbugs"` (they differ in truncation syntax).

## Value

A list with `code` and `data`.

## Examples

``` r
d <- data.frame(y = rpois(20, 3), x = rnorm(20))
cat(bugs_code(y ~ x, family = "poisson", data = d)$code)
#> model {
#>   for (i in 1:N) {
#>     eta_lambda[i] <- inprod(X_lambda[i, ], b_lambda[]) + off_lambda[i]
#>     lambda[i] <- exp(eta_lambda[i])
#>     y[i] ~ dpois(lambda[i])
#>   }
#>   b_lambda[1] ~ dnorm(0, 0.04)
#>   b_lambda[2] ~ dnorm(0, 0.04)
#> }
```
