# Monotone transformations

The law of `Y = g(X)` for a monotone, differentiable transformation `g`
of a continuous base law:
`log f_Y(y) = log f_X(g^{-1}(y)) + log |d g^{-1}(y) / dy|`. Named
transformations: `"log"` (`Y = log X`), `"exp"` (`Y = exp X`), `"logit"`
(`Y = logit X`), `"inv_logit"`, `"reciprocal"` (`Y = 1 / X` on positive
support) and `"sqrt"`; or give `forward`, `inverse` and `log_jacobian`
functions yourself.

## Usage

``` r
transform_dist(
  base,
  transform = c("log", "exp", "logit", "inv_logit", "reciprocal", "sqrt"),
  forward = NULL,
  inverse = NULL,
  log_jacobian = NULL,
  increasing = TRUE,
  support = NULL
)
```

## Arguments

- base:

  Base family.

- transform:

  Name of a built in transformation, or `NULL` when the functions are
  given.

- forward, inverse:

  Functions `g` and `g^{-1}`.

- log_jacobian:

  Function of `y` returning `log |d g^{-1}(y) / dy|`.

- increasing:

  Whether `g` is increasing.

- support:

  Optional function of the parameters returning the new support bounds.

## Value

A `regcomp_family`.

## Examples

``` r
ln <- transform_dist("normal", "exp")   # a lognormal law
d <- distribution(ln, mu = 0, sigma = 1)
prob(d, 1); dlnorm(1)
#> [1] 0.3989423
#> [1] 0.3989423
```
