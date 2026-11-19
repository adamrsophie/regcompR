# Priors for Bayesian fits

Priors live on the coefficients of each linear predictor (on the link
scale), on random effect standard deviations and on random effect
correlations. A prior is written as a distribution call such as
`"normal(0, 2.5)"`; the same specification is translated to 'Stan' and
to the BUGS language used by 'JAGS' and 'OpenBUGS'.

## Usage

``` r
prior(dist, par = NULL, class = c("b", "Intercept", "sd", "cor"), coef = NULL)
```

## Arguments

- dist:

  Prior distribution as a string.

- par:

  Distribution parameter the prior applies to (`NULL` for all).

- class:

  `"b"` (coefficients other than the intercept), `"Intercept"`, `"sd"`
  or `"cor"`.

- coef:

  Optional coefficient name (for `class = "b"`) or random effect term
  label (for `"sd"` and `"cor"`).

## Value

An object of class `regcomp_prior`; combine several with
[`c()`](https://rdrr.io/r/base/c.html).

## Details

Supported prior distributions: `normal(mean, sd)`,
`student_t(df, location, scale)`, `cauchy(location, scale)`,
`logistic(location, scale)`, `exponential(rate)`, `gamma(shape, rate)`,
`uniform(lower, upper)` and, for correlations, `lkj(eta)`. Standard
deviation priors are truncated to positive values.

Defaults: `normal(0, 5)` for intercepts and coefficients,
`student_t(3, 0, 2.5)` for standard deviations and `lkj(2)` for
correlations ('Stan'; BUGS engines use an inverse Wishart with identity
scale for correlated random effects, which is noted in the fit).

## Examples

``` r
p <- c(prior("normal(0, 1)", par = "mu", class = "b"),
       prior("student_t(3, 0, 1)", class = "sd"))
p
#> normal(0, 1)                 class b         par mu     coef (all)
#> student_t(3, 0, 1)           class sd        par (all)  coef (all)
```
