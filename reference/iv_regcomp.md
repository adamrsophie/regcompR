# Instrumental variable regression

Two stage least squares for a linear model with endogenous regressors,
with the usual diagnostics: a weak instrument F test for each endogenous
regressor (the first stage F test of the excluded instruments), the
Wu-Hausman test of endogeneity and the Sargan test of overidentifying
restrictions. The results equal `ivreg::ivreg()` (and `AER::ivreg()`).

## Usage

``` r
iv_regcomp(formula, data)
```

## Arguments

- formula:

  A three part formula `y ~ regressors | instruments`, as in ivreg:
  exogenous regressors appear on both sides, endogenous regressors only
  on the left of `|`.

- data:

  A data frame.

## Value

An object of class `regcomp_lm` with coefficients, covariance
(`sigma^2 (X'P_Z X)^-1` with `sigma^2 = RSS/(n - p)`), residuals and a
`diagnostics` table.

## Examples

``` r
set.seed(1)
n <- 500; z1 <- rnorm(n); z2 <- rnorm(n); u <- rnorm(n)
x <- z1 + 0.5 * z2 + u + rnorm(n)
d <- data.frame(y = 1 + 2 * x + u, x, z1, z2)
iv_regcomp(y ~ x | z1 + z2, data = d)
#> <regcomp_lm> instrumental variables (2SLS), n = 500, residual df = 498
#>             estimate std_error       t p_value
#> (Intercept)   0.9953    0.0476 20.9157       0
#> x             1.9470    0.0420 46.3696       0
#> 
#> Diagnostic tests:
#>                      df1 df2 statistic p_value
#> Weak instruments (x)   2 497  148.7356  0.0000
#> Wu-Hausman             1 497  202.4646  0.0000
#> Sargan                 1  NA    1.5290  0.2163
```
