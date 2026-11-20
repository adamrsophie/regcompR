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
anova(object, ...)
```

## Arguments

- object, ...:

  Fitted models.

## Value

A data frame of class `anova` with degrees of freedom, log likelihood,
AIC, the likelihood ratio statistic, its degrees of freedom and p value.

## Examples

``` r
d <- data.frame(x = rnorm(200), z = rnorm(200)); d$y <- rpois(200, exp(0.5 + 0.3 * d$x))
anova(regcomp(y ~ x, family = "poisson", data = d), regcomp(y ~ x + z, family = "poisson", data = d))
#> Likelihood ratio tests
#> 
#>                                                  df  logLik    AIC  Chisq
#> regcomp(y ~ x, family = "poisson", data = d)      2 -312.23 628.45       
#> regcomp(y ~ x + z, family = "poisson", data = d)  3 -312.02 630.04 0.4059
#>                                                  Chi Df Pr(>Chisq)
#> regcomp(y ~ x, family = "poisson", data = d)                      
#> regcomp(y ~ x + z, family = "poisson", data = d)      1     0.5241
```
