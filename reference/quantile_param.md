# Quantile parameterized families

Replace a parameter of a family by the `tau` quantile of the law, so a
regression on that parameter is a parametric quantile regression. For
given values of the remaining parameters, the original parameter is
found by root finding on its link scale so that `Q(tau) = q`.

## Usage

``` r
quantile_param(base, par = names(get_family(base)$pars)[1], tau = 0.5)
```

## Arguments

- base:

  Base family.

- par:

  Parameter to replace (default: the first one).

- tau:

  Quantile level in (0, 1).

## Value

A `regcomp_family` whose first parameter `q` is the quantile.

## References

Noufaily and Jones (2013), Parametric quantile regression based on the
generalized gamma distribution. JRSS C, 62, 723-740.

## Examples

``` r
qw <- quantile_param("weibull", par = "scale", tau = 0.5)
```
