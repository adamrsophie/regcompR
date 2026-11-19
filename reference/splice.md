# Spliced (composite) distributions

Joins a body law below a threshold `u` and a tail law above it:
`f(y) = w f_b(y) / F_b(u)` for `y <= u` and `(1 - w) f_t(y) / S_t(u)`
for `y > u` (Cooray and Ananda 2005; Scollnik 2007). The threshold is a
parameter (`threshold`, identity link unless the body is positive) and
`w` the body weight (logit link). Parameters of the body get suffix `.b`
and of the tail `.t`.

## Usage

``` r
splice(body, tail)
```

## Arguments

- body, tail:

  Families.

## Value

A `regcomp_family`.

## Examples

``` r
sp <- splice("lognormal", "lomax")
names(sp$pars)
#> [1] "mu.b"      "sigma.b"   "scale.t"   "alpha.t"   "threshold" "w"        
```
