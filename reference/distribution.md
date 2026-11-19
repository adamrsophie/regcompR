# Create a distribution object

Combines a registered family (or one built by an operator such as
[`inflate()`](https://adamrsophie.github.io/regcompR/reference/inflate.md))
with parameter values. Parameters are vectorized: a distribution object
with vectors of length `n` represents `n` distributions, one per
position, as in regression predictions.

## Usage

``` r
distribution(family, ..., parameterization = NULL)
```

## Arguments

- family:

  Identifier, alias or `regcomp_family` object.

- ...:

  Named parameter values.

- parameterization:

  Optional parameterization name, see
  [`get_family()`](https://adamrsophie.github.io/regcompR/reference/get_family.md).

## Value

An object of class `regcomp_dist`.

## See also

[`log_prob()`](https://adamrsophie.github.io/regcompR/reference/probability-interface.md),
[`cdf()`](https://adamrsophie.github.io/regcompR/reference/probability-interface.md),
[`quantile.regcomp_dist()`](https://adamrsophie.github.io/regcompR/reference/quantile.regcomp_dist.md),
[`rng()`](https://adamrsophie.github.io/regcompR/reference/probability-interface.md),
[`capabilities()`](https://adamrsophie.github.io/regcompR/reference/capabilities.md),
[`plot.regcomp_dist()`](https://adamrsophie.github.io/regcompR/reference/plot.regcomp_dist.md).

## Examples

``` r
d <- distribution("negative_binomial", parameterization = "nb2",
                  mu = 4, phi = 2)
d
#> <regcomp_dist> Negative binomial (NB2)
#>   mu = 4
#>   phi = 2
prob(d, 0:5)
#> [1] 0.1111111 0.1481481 0.1481481 0.1316872 0.1097394 0.0877915
cdf(d, 3)
#> [1] 0.5390947
quantile(d, c(0.1, 0.5, 0.9))
#> [1] 0 3 9
mean(d); variance(d)
#> [1] 4
#> [1] 12
distribution("conway_maxwell_poisson_mean", mu = 3, nu = 1.5)
#> <regcomp_dist> Conway-Maxwell-Poisson (mean parameterization)
#>   mu = 3
#>   nu = 1.5
```
