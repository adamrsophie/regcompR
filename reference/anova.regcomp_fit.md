# Likelihood ratio tests of nested models

Compares fits of the same response on the same rows, ordered by their
number of parameters, each against the previous one. The models must be
nested; that is not checked. When the smaller model sets a variance
component to zero the null value is on the boundary and the chi square p
value is conservative (about twice the correct one for one variance;
Self and Liang 1987).

## Usage

``` r
# S3 method for class 'regcomp_fit'
anova(object, ..., boundary = FALSE)
```

## Arguments

- object, ...:

  Fitted models.

- boundary:

  If `TRUE`, p values for an added variance component use the 50:50
  mixture of chi squares with `df - 1` and `df` degrees of freedom (Self
  and Liang 1987), as the variance lies on its boundary.

## Value

A data frame of class `anova` with degrees of freedom, log likelihood,
AIC, the likelihood ratio statistic, its degrees of freedom and p value.

## Examples

``` r
d <- data.frame(x = rnorm(200), z = rnorm(200)); d$y <- rpois(200, exp(0.5 + 0.3 * d$x))
f1 <- regcomp(y ~ x, family = "poisson", data = d)
f2 <- regcomp(y ~ x + z, family = "poisson", data = d)
anova(f1, f2)
#> Likelihood ratio tests
#> 
#>    df  logLik    AIC  Chisq Chi Df Pr(>Chisq)
#> f1  2 -312.23 628.45                         
#> f2  3 -312.02 630.04 0.4059      1     0.5241
```
