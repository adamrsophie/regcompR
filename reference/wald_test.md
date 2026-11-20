# Wald tests of linear hypotheses

Tests `L beta = rhs` with the statistic
`(L b - rhs)' (L V L')^-1 (L b - rhs)`, chi square with `nrow(L)`
degrees of freedom, where `V` is the chosen covariance.

## Usage

``` r
wald_test(object, L, rhs = 0, vcov = "model", ...)
```

## Arguments

- object:

  A frequentist fit.

- L:

  The hypothesis: a character vector of coefficient names (all zero), a
  named numeric vector (one linear combination, for example
  `c("mu:x" = 1, "mu:z" = -1)` for equal coefficients), or a matrix with
  one row per combination and coefficient names as column names.

- rhs:

  Right hand side, recycled to the number of rows of `L`.

- vcov:

  `"model"`, `"sandwich"` or `"quasi"` (see
  [`vcov.regcomp_fit()`](https://adamrsophie.github.io/regcompR/reference/vcov.regcomp_fit.md)).

- ...:

  Passed to
  [`vcov.regcomp_fit()`](https://adamrsophie.github.io/regcompR/reference/vcov.regcomp_fit.md)
  (such as `cluster`).

## Value

A data frame with the estimates of the combinations, the statistic, its
degrees of freedom and the p value.

## Examples

``` r
d <- data.frame(x = rnorm(200), z = rnorm(200)); d$y <- rpois(200, exp(0.5 + 0.3 * d$x + 0.3 * d$z))
f <- regcomp(y ~ x + z, family = "poisson", data = d)
wald_test(f, c("lambda:x" = 1, "lambda:z" = -1))
#>   estimate      chisq df   p_value
#> 1 -0.01019 0.02043228  1 0.8863364
wald_test(f, c("lambda:x", "lambda:z"))
#>        estimate    chisq df      p_value
#> 1 0.319, 0.3292 71.99382  2 2.326696e-16
```
