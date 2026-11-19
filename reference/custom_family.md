# User defined families

Builds a family from a log density (or log mass) so that every engine
path that needs only the likelihood (maximum likelihood, random effects,
penalties, prediction by numerical inversion, scoring) works with it.
The normalization is checked numerically at `check` parameter values.

## Usage

``` r
custom_family(
  id,
  pars,
  lpdf,
  discrete = FALSE,
  support = NULL,
  cdf = NULL,
  quantile = NULL,
  rng = NULL,
  mean = NULL,
  variance = NULL,
  stan = NULL,
  check = NULL,
  ...
)
```

## Arguments

- id:

  Identifier.

- pars:

  List of parameters from
  [`par_spec()`](https://adamrsophie.github.io/regcompR/reference/par_spec.md).

- lpdf:

  Function `(y, p)` returning log densities, vectorized over `y` and the
  parameter vectors in `p`.

- discrete:

  Whether the law is discrete.

- support:

  Optional function of the parameters returning `list(lower, upper)`.

- cdf, quantile, rng, mean, variance:

  Optional functions.

- stan:

  Optional Stan template from the registry helpers.

- check:

  Optional named list of parameter values at which total mass is
  verified.

- ...:

  Further arguments to the family constructor.

## Value

A `regcomp_family`.

## Examples

``` r
f <- custom_family("my_exp", list(par_spec("rate", "positive", role = "rate")),
  lpdf = function(y, p) dexp(y, p$rate, log = TRUE),
  support = function(p) list(lower = 0, upper = Inf), check = list(rate = 2))
set.seed(1)
d <- data.frame(y = rexp(100, 2))
coef(regcomp(y ~ 1, family = f, data = d))
#> rate:(Intercept) 
#>        0.6629319 
```
