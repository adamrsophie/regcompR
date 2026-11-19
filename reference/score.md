# Score and Hessian of the log likelihood

Derivatives of `log_prob(x, y)` with respect to the natural parameters,
per observation. Analytic expressions are used when a family provides
them; otherwise central finite differences with a step scaled to each
parameter.

## Usage

``` r
score(x, y, ...)

hessian(x, y, ...)
```

## Arguments

- x:

  A distribution object.

- y:

  Response values.

- ...:

  Unused.

## Value

`score()`: a matrix with one column per parameter. `hessian()`: an array
`n x k x k`.

## Examples

``` r
d <- distribution("normal", mu = 0, sigma = 1)
score(d, c(-1, 0, 2))
#>      mu         sigma
#> [1,] -1  1.554312e-10
#> [2,]  0 -1.000000e+00
#> [3,]  2  3.000000e+00
```
