# Wrapped distributions

Wraps a law on the real line (or the positive half line) onto the circle
of period `2 pi`: `f_W(y) = sum_k f(y + 2 pi k)`. The sum is truncated
once terms fall below `1e-14` of the running total (at most 2000 terms
on each side). Gives the wrapped normal, Cauchy, Laplace, asymmetric
Laplace, exponential, Levy, Lindley, xgamma and any other wrapped law.

## Usage

``` r
wrap(base)
```

## Arguments

- base:

  Base family on the real or positive line.

## Value

A circular `regcomp_family`.

## References

Mardia and Jupp (2000), Directional Statistics, sec. 3.5.

## Examples

``` r
we <- wrap("exponential")
d <- distribution(we, lambda = 0.8)
prob(d, c(-1, 0.5, 2))
#> [1] 0.01175932 0.53979788 0.16258400
```
