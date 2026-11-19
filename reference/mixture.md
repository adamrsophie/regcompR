# Finite mixtures

Mixture of `K` component families with weights `w_k`:
`f(y) = sum_k w_k f_k(y)`. Component parameters are renamed with
suffixes `.1`, `.2`, ... (for example `mu.1`). Weights are parameters
`w2, ..., wK` on the multinomial logit scale relative to the first
component (`w1` is implied), each with an identity link. Likelihoods are
marginalized over the component labels; label switching is possible and
is reported by
[`identifiability_check()`](https://adamrsophie.github.io/regcompR/reference/identifiability_check.md).

## Usage

``` r
mixture(...)
```

## Arguments

- ...:

  Component families.

## Value

A `regcomp_family`.

## Examples

``` r
mx <- mixture("normal", "normal")
d <- distribution(mx, mu.1 = -2, sigma.1 = 1, mu.2 = 2, sigma.2 = 0.5, w2 = 0)
plot(d)
```
