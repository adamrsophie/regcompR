# Skew-symmetric distributions

Azzalini type perturbation of a symmetric location scale law:
`f(y) = (2 / sigma) f0(z) G(alpha z)`, `z = (y - mu) / sigma`, with `G`
a symmetric CDF (`"normal"`, `"logistic"`, `"cauchy"`, `"laplace"` or
`"t"` with `nu` degrees of freedom). This covers the skew elliptical and
skew symmetric classes (skew normal, skew t, skew logistic, ...).

## Usage

``` r
skew_symmetric(
  base,
  skew_cdf = c("normal", "logistic", "cauchy", "laplace", "t"),
  nu = 4
)
```

## Arguments

- base:

  Symmetric family with parameters `mu` and `sigma`.

- skew_cdf:

  Symmetric skewing CDF.

- nu:

  Degrees of freedom for `skew_cdf = "t"`.

## Value

A `regcomp_family` with added parameter `alpha`.

## References

Azzalini (1985), A class of distributions which includes the normal
ones. Scand. J. Stat., 12, 171-178. Wang, Boyer and Genton (2004), A
skew-symmetric representation of multivariate distributions. Statistica
Sinica, 14, 1259-1270.

## Examples

``` r
sl <- skew_symmetric("logistic", "logistic")
```
