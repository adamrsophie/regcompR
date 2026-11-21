# Predictive shrinkage of regression coefficients

Shrinks the slopes of one parameter by a common factor `c` to correct
the overfitting of predictions in new data: the heuristic factor
`c = (LR - df) / LR` of van Houwelingen and le Cessie (1990), with `LR`
the likelihood ratio statistic against the model without the parameter's
covariates and `df` its degrees of freedom, or a factor given by the
user. The intercept and the other parameters (precision of a beta
regression, say) are then re-estimated with the shrunken slopes held
fixed. Works for any family (beta regression, GLMs, distributional
models) without random effects.

## Usage

``` r
shrink_regcomp(object, par = NULL, factor = NULL)
```

## Arguments

- object:

  A fitted model (frequentist, no random effects).

- par:

  Parameter whose slopes are shrunk (default the first).

- factor:

  Optional shrinkage factor in (0, 1\]; default heuristic.

## Value

The fit with shrunken coefficients (without standard errors), with the
factor as attribute `"shrinkage"`.

## References

van Houwelingen and le Cessie (1990), Predictive value of statistical
models. Stat. Med., 9, 1303-1325.

## Examples

``` r
d <- data.frame(x1 = rnorm(60), x2 = rnorm(60), x3 = rnorm(60))
d$y <- rbeta(60, plogis(0.3 * d$x1) * 10, (1 - plogis(0.3 * d$x1)) * 10)
f <- regcomp(y ~ x1 + x2 + x3, family = "beta", data = d)
attr(shrink_regcomp(f), "shrinkage")
#> [1] 0.8551958
```
