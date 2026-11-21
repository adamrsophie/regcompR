# Mean, mode, geometric and harmonic mean parameterizations

Replace a parameter of a continuous family by a summary of the law, so
that the regression is on that summary: the mean, the mode, the
geometric mean `exp(E log Y)` or the harmonic mean `1 / E(1 / Y)`. For
given values of the other parameters the original parameter is found by
root finding, as in
[`quantile_param()`](https://adamrsophie.github.io/regcompR/reference/quantile_param.md)
(which covers the median and other quantiles). Expectations use 40 point
Gauss-Hermite quadrature over `Y = Q(Phi(z))`, exact for log normal
laws; the mode is the root of the derivative of the log density (or the
lower end of the support).

## Usage

``` r
summary_param(
  base,
  par = names(get_family(base)$pars)[1],
  type = c("mean", "mode", "geometric_mean", "harmonic_mean")
)
```

## Arguments

- base:

  Base family.

- par:

  Parameter to replace (default: the first one).

- type:

  `"mean"`, `"mode"`, `"geometric_mean"` or `"harmonic_mean"`.

## Value

A `regcomp_family` whose first parameter `q` is the summary.

## Examples

``` r
gm <- summary_param("gamma", type = "mode")
# \donttest{
d <- data.frame(x = runif(200)); d$y <- rgamma(200, shape = 3, scale = exp(1 + d$x) / 2)
coef(regcomp(y ~ x, family = gm, data = d))
#>     q:(Intercept)               q:x sigma:(Intercept) 
#>         1.0546571         1.0990234        -0.6016625 
# }
```
