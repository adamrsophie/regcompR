# Censored distributions

Censors a continuous law at `lower` and/or `upper`: values beyond a
bound are recorded at the bound, which then carries the tail probability
as a point mass. `censor("normal", lower = 0)` is the Tobit model and
the rectified Gaussian; `censor("beta", ...)` style limits give two
limit models. Data at the bounds contribute log tail probabilities.
Count laws can be censored too (`censor("poisson", upper = 3)` records
every count of 3 or more as 3), as in censored Poisson and negative
binomial regression.

## Usage

``` r
censor(base, lower = -Inf, upper = Inf)
```

## Arguments

- base:

  Base family (continuous or count).

- lower, upper:

  Censoring points.

## Value

A `regcomp_family` with mixed support.

## References

Tobin (1958), Estimation of relationships for limited dependent
variables. Econometrica, 26, 24-36.

## Examples

``` r
tobit <- censor("normal", lower = 0)
d <- distribution(tobit, mu = 0.5, sigma = 1)
prob(d, c(0, 1))
#> [1] 0.3085375 0.3520653
```
