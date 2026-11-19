# Sine-skewed circular distributions

`f(y) = f0(y) (1 + lambda sin(y - mu))` for a symmetric circular base
law `f0` with location `mu` (Abe and Pewsey 2011), `lambda` in (-1, 1).

## Usage

``` r
sine_skew(base)
```

## Arguments

- base:

  Symmetric circular family with a parameter named `mu`.

## Value

A `regcomp_family`.

## References

Abe and Pewsey (2011), Sine-skewed circular distributions. Statistical
Papers, 52, 683-707.

## Examples

``` r
ssvm <- sine_skew("von_mises")
```
