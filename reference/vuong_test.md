# Vuong test of non nested models

Compares two models fitted to the same observations through their
pointwise log likelihood ratios `m_i`: `z = sqrt(n) mean(m) / sd(m)`,
with the AIC and BIC corrected versions that subtract `(k1 - k2)` and
`(k1 - k2) log(n) / 2` from the summed ratio (Vuong 1989). Positive
values favor the first model. The raw statistic equals `pscl::vuong()`;
the corrections count every estimated parameter (pscl counts only the
regression coefficients).

## Usage

``` r
vuong_test(object1, object2)
```

## Arguments

- object1, object2:

  Fitted models on the same rows.

## Value

A data frame with the three statistics and one sided p values.

## Examples

``` r
d <- data.frame(x = rnorm(200)); d$y <- rnbinom(200, mu = exp(1 + 0.3 * d$x), size = 2)
vuong_test(regcomp(y ~ x, family = "negative_binomial_nb2", data = d),
           regcomp(y ~ x, family = "poisson", data = d))
#>               statistic      p_value preferred
#> Raw            4.240469 1.115264e-05   model 1
#> AIC-corrected  4.135843 1.768269e-05   model 1
#> BIC-corrected  3.963298 3.696073e-05   model 1
```
