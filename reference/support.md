# Support of a distribution object

Support of a distribution object

## Usage

``` r
support(x, ...)
```

## Arguments

- x:

  A distribution object.

- ...:

  Unused.

## Value

A list with `lower` and `upper` bounds (vectors when the support depends
on parameters) and `discrete`.

## Examples

``` r
support(distribution("generalized_pareto", mu = 0, sigma = 1, xi = -0.5))
#> $lower
#> [1] 0
#> 
#> $upper
#> [1] 2
#> 
#> $discrete
#> [1] FALSE
#> 
```
