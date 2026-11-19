# Discretized and rounded distributions

`discretize()` turns a continuous law into a discrete one on the
integers by interval probabilities
`P(Y = k) = F(k + 1 - offset) - F(k - offset)` (the default `offset = 0`
maps `[k, k + 1)` to `k`; `offset = 0.5` rounds to the nearest integer).
`round_dist()` gives the law of a value rounded to `digits` decimals,
`P(Y = y) = F(y + h/2) - F(y - h/2)` with `h = 10^-digits`, on the grid
of multiples of `h`.

## Usage

``` r
discretize(base, offset = 0)

round_dist(base, digits = 0)
```

## Arguments

- base:

  Continuous base family.

- offset:

  Shift of the interval boundaries.

- digits:

  Number of decimals.

## Value

A `regcomp_family`.

## Examples

``` r
dg <- discretize("gamma")
d <- distribution(dg, mu = 3, sigma = 0.5)
sum(prob(d, 0:200))
#> [1] 1
r1 <- round_dist("normal", digits = 1)
```
