# Inflated distributions

Adds point masses to a base law: with weights `pi_j` at points `c_j`,
`P(Y = c_j) = pi_j + (1 - sum(pi)) P0(Y = c_j)` and elsewhere the base
law is scaled by `1 - sum(pi)`. For a continuous base law the result is
a mixed distribution whose log likelihood uses the point mass at `c_j`
and the scaled density elsewhere. One inflation parameter per point is
added (named `zi`, `oi` or `pi1`, `pi2`, ...), each in (0, 1) with a
logit link; with several points the parameters are mapped to weights
through a multinomial logit (stick breaking is not used), so
`pi_j = exp(eta_j) / (1 + sum(exp(eta)))`.

## Usage

``` r
inflate(base, at = 0, names = NULL)
```

## Arguments

- base:

  Base family (identifier or object).

- at:

  Inflation points.

- names:

  Optional parameter names.

## Value

A `regcomp_family`.

## Details

Zero inflation can only add zeros; compare
[`hurdle()`](https://adamrsophie.github.io/regcompR/reference/hurdle.md).

## Examples

``` r
zip <- inflate("poisson", at = 0)
d <- distribution(zip, lambda = 2, zi = 0.3)
prob(d, 0:3)
#> [1] 0.3947347 0.1894694 0.1894694 0.1263129
mean(d)
#> [1] 1.4
zoib <- inflate("beta", at = c(0, 1))
zoib
#> <regcomp_family inflate(beta, at = 0,1)> Beta (mean and precision, Ferrari and Cribari-Neto) inflated at 0, 1
#>   category: bounded, variate: scalar, continuous, status: implemented
#>   parameters:
#>     mu       unit              link logit     [mean]
#>     phi      positive          link log       
#>     pi1      unit              link logit     
#>     pi2      unit              link logit     
#>   formula: P(Y = c) = pi + (1 - pi) P0(Y = c); elsewhere (1 - pi) f0(y) 
#>   built by operator: inflate 
```
