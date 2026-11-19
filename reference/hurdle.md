# Hurdle distributions

Replaces the base probability at `at` by a separate hurdle probability:
`P(Y = at) = pi` and, for `y != at`,
`P(Y = y) = (1 - pi) f0(y) / (1 - P0(Y = at))`. Unlike inflation, a
hurdle can produce fewer as well as more values at the hurdle than the
base law. For a continuous base law the conditioning term is one, giving
a two part (zero adjusted) model such as the zero adjusted gamma.

## Usage

``` r
hurdle(base, at = 0, name = "hu")
```

## Arguments

- base:

  Base family.

- at:

  Hurdle point.

- name:

  Name of the hurdle probability parameter.

## Value

A `regcomp_family`.

## Examples

``` r
hnb <- hurdle("negative_binomial_nb2", at = 0)
d <- distribution(hnb, mu = 3, phi = 1, hu = 0.1)
prob(d, 0:3)
#> [1] 0.1000000 0.2250000 0.1687500 0.1265625
zaga <- hurdle("gamma", at = 0)
```
