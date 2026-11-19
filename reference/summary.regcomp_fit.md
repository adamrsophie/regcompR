# Summarize a fitted model

Frequentist fits report estimates, Wald standard errors, z statistics
and p values (labeled as Wald; penalized fits show no p values).
Bayesian fits report posterior means, standard deviations, quantiles,
R-hat and effective sample sizes. The two are different objects and are
never presented as interchangeable.

## Usage

``` r
# S3 method for class 'regcomp_fit'
summary(object, level = 0.95, ...)
```

## Arguments

- object:

  A fitted model.

- level:

  Interval level.

- ...:

  Unused.

## Value

An object of class `summary.regcomp_fit`.
