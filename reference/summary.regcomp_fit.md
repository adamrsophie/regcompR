# Summarize a fitted model

Frequentist fits report estimates, Wald standard errors, z statistics
and p values (labeled as Wald; penalized fits show no p values).
Bayesian fits report posterior means, standard deviations, quantiles,
R-hat and effective sample sizes. The two are different objects and are
never presented as interchangeable.

## Usage

``` r
# S3 method for class 'regcomp_fit'
summary(
  object,
  level = 0.95,
  vcov = c("model", "sandwich", "quasi"),
  cluster = NULL,
  ...
)
```

## Arguments

- object:

  A fitted model.

- level:

  Interval level.

- vcov:

  Covariance for the Wald statistics: `"model"`, `"sandwich"` or
  `"quasi"` (see
  [`vcov.regcomp_fit()`](https://adamrsophie.github.io/regcompR/reference/vcov.regcomp_fit.md));
  the last two cover the regression coefficients only.

- cluster:

  Cluster identifier for `vcov = "sandwich"`.

- ...:

  Unused.

## Value

An object of class `summary.regcomp_fit`.
