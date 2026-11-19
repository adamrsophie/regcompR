# Compare fitted models

Tabulates log likelihood, number of parameters, AIC and BIC, and
optionally a cross validated score, ordered by AIC.

## Usage

``` r
compare_models(..., cv = NULL, K = 10)
```

## Arguments

- ...:

  Fitted models, optionally named.

- cv:

  Optional scoring rule name; when given,
  [`cv_regcomp()`](https://adamrsophie.github.io/regcompR/reference/cv_regcomp.md)
  with `K` folds is added.

- K:

  Number of folds for cross validation.

## Value

A data frame.

## Examples

``` r
set.seed(1)
d <- data.frame(x = rnorm(200)); d$y <- rnbinom(200, mu = exp(1 + 0.4 * d$x), size = 2)
compare_models(poisson = regcomp(y ~ x, family = "poisson", data = d),
               nb2 = regcomp(y ~ x, family = "nb2", data = d))
#>     model                family    logLik df      AIC      BIC
#> 1     nb2 negative_binomial_nb2 -431.7319  3 869.4639 879.3589
#> 2 poisson               poisson -482.1963  2 968.3925 974.9892
```
