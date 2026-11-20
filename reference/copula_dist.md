# Bivariate copula models

Joint law of a pair `(y1, y2)` with margins `m1`, `m2` and a copula:
`"gaussian"` (correlation `rho`), `"clayton"` (`theta > 0`), `"gumbel"`
(`theta > 1`, parameterized as `1 + exp(eta)`) or `"frank"` (`theta`
real, not 0). For continuous margins the log density is
`log c(F1(y1), F2(y2)) + log f1(y1) + log f2(y2)`; for discrete margins
the rectangle probability
`C(F1, F2) - C(F1-, F2) - C(F1, F2-) + C(F1-, F2-)` is used (Trivedi and
Zimmer 2007). Margins may differ in type, and continuous margins may
have point masses: at a jump of one margin and a continuous value of the
other the density is `(h(F1, F2) - h(F1-, F2)) f2(y2)` with `h = dC/dv`.
Margin parameters get suffixes `.1` and `.2`; every parameter can take a
formula, giving copula regression. The response is `cbind(y1, y2)`.

## Usage

``` r
copula_dist(m1, m2, copula = c("gaussian", "clayton", "gumbel", "frank"))
```

## Arguments

- m1, m2:

  Scalar margin families (continuous, discrete, or continuous with point
  masses).

- copula:

  Copula name.

## Value

A `regcomp_family`.

## Examples

``` r
cp <- copula_dist("normal", "gamma", "clayton")
names(cp$pars)
#> [1] "mu.1"    "sigma.1" "mu.2"    "sigma.2" "theta"  
```
