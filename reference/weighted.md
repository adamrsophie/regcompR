# Weighted and contaminated distributions

`weighted()` gives the weighted law `f_w(y) = w(y) f(y) / E[w(Y)]` with
`w(y) = y^k` (length biased for `k = 1`, size biased in general). The
normalizing moment is computed by quadrature (or summation for a
discrete base).

## Usage

``` r
weighted(base, k = 1)

contaminated(base, par)
```

## Arguments

- base:

  Base family.

- k:

  Weight exponent.

- par:

  Name of the parameter inflated in the contaminating component.

## Value

A `regcomp_family`.

## Details

`contaminated()` gives `delta f(y; theta) + (1 - delta) f(y; theta*)`,
where `theta*` multiplies the parameter `par` by an inflation factor
`eta > 1` (Punzo and coauthors' contaminated models). Parameters `delta`
(logit link) and `eta` (log of `eta - 1`) are added.

## References

Patil and Rao (1978), Weighted distributions and size-biased sampling.
Biometrics, 34, 179-189. Punzo, Mazza and Maruotti (2018), Fitting
insurance and economic data with outliers: a flexible approach based on
finite mixtures of contaminated gamma distributions. J. Appl. Stat., 45,
2563-2584.

## Examples

``` r
lb <- weighted("birnbaum_saunders")
cg <- contaminated("gamma", par = "sigma")
```
