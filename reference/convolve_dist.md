# Convolution of two independent laws

The law of `X1 + X2` for independent `X1` and `X2`: the density is
`int f1(y - x) f2(x) dx` by adaptive quadrature for continuous laws and
a finite sum for counts. Parameters get suffixes `.1` and `.2`. Examples
include the exponentially modified laws, the normal plus Laplace law of
Reed and Jorgensen and Poisson plus negative binomial counts.

## Usage

``` r
convolve_dist(a, b)
```

## Arguments

- a, b:

  Families, both continuous or both discrete.

## Value

A `regcomp_family`.

## Examples

``` r
cv <- convolve_dist("normal", "exponential")
d <- distribution(cv, mu.1 = 0, sigma.1 = 1, lambda.2 = 2)
prob(d, 1)
#> [1] 0.3173105
```
