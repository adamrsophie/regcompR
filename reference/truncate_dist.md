# Truncated distributions

Restricts a base law to `[lower, upper]` (discrete; both end points
kept) or `(lower, upper]` (continuous) and renormalizes:
`f_T(y) = f(y) / (F(upper) - F(lower-))`, computed with log scale
differences. For observation specific truncation use
[`cens()`](https://adamrsophie.github.io/regcompR/reference/cens.md).

## Usage

``` r
truncate_dist(base, lower = -Inf, upper = Inf)
```

## Arguments

- base:

  Base family.

- lower, upper:

  Truncation bounds.

## Value

A `regcomp_family`.

## Examples

``` r
ztp <- truncate_dist("poisson", lower = 1)
d <- distribution(ztp, lambda = 1.5)
prob(d, 0:3); mean(d)
#> [1] 0.0000000 0.4308254 0.3231190 0.1615595
#> [1] 1.930825
tn <- truncate_dist("normal", lower = 0)
```
