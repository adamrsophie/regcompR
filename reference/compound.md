# Compound (mixed) distributions

Integrates one parameter of a base law over a mixing law:
`f(y) = int f(y | par = t) g(t) dt`. The mixing family's parameters are
appended (suffixed with `_mix`). Covers mixed Poisson laws
(`compound("poisson", "lambda", "gamma_rate")` is the negative
binomial), Rayleigh mixtures, normal variance mixtures and other
compound laws.

## Usage

``` r
compound(base, par, mixing)
```

## Arguments

- base:

  Base family.

- par:

  Parameter of `base` that is integrated out.

- mixing:

  Family of the mixing law, on the domain of `par`.

## Value

A `regcomp_family`.

## References

Karlis and Xekalaki (2005), Mixed Poisson distributions. International
Statistical Review, 73, 35-58.

## Examples

``` r
nb <- compound("poisson", "lambda", "gamma_rate")
```
