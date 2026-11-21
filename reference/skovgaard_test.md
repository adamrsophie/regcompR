# Skovgaard's modified directed likelihood ratio

Higher order inference on one coefficient: the signed root of the
likelihood ratio statistic `r` and Skovgaard's (1996) modification
`r* = r + log(u / r) / r + log(C) / r`, whose normal approximation has
relative error of third order. The sample space derivatives are
approximated by covariances of the score and of the log likelihood ratio
over `R` responses simulated from the fit (Skovgaard 2001), as in
`likelihoodAsy::rstar()`. For models without random effects.

## Usage

``` r
skovgaard_test(object, term, value = 0, R = 1000, seed = 1)
```

## Arguments

- object:

  A fitted model (frequentist, no random effects).

- term:

  Coefficient name as in `coef(object)`.

- value:

  Value of the coefficient under the null hypothesis.

- R:

  Number of simulated samples.

- seed:

  Random seed.

## Value

A data frame with the Wald statistic, `r` and `r*` and their two sided
normal p values.

## References

Skovgaard (1996), An explicit large-deviation approximation to
one-parameter tests. Bernoulli, 2, 145-165. Skovgaard (2001), Likelihood
asymptotics. Scand. J. Statist., 28, 3-32.

## Examples

``` r
d <- data.frame(x = rnorm(15)); d$y <- rpois(15, exp(0.5 + 0.4 * d$x))
f <- regcomp(y ~ x, family = "poisson", data = d)
skovgaard_test(f, "lambda:x", R = 500)
#>      statistic     p_value
#> Wald  2.559293 0.010488532
#> r     2.735483 0.006228873
#> r*    2.711145 0.006705129
```
