# Plot a distribution

Draws the density (continuous) or probability mass (discrete), the
distribution function, the survival function, the hazard or the quantile
function of a distribution object. Vectorized distribution objects are
overlaid (up to `max_curves`).

## Usage

``` r
# S3 method for class 'regcomp_dist'
plot(
  x,
  type = c("density", "cdf", "survival", "hazard", "quantile"),
  xlim = NULL,
  n = 401,
  max_curves = 8,
  col = NULL,
  add = FALSE,
  legend = TRUE,
  main = NULL,
  xlab = NULL,
  ylab = NULL,
  ...
)
```

## Arguments

- x:

  A distribution object.

- type:

  `"density"`, `"cdf"`, `"survival"`, `"hazard"` or `"quantile"`.

- xlim:

  Optional range of the response.

- n:

  Number of evaluation points for continuous laws.

- max_curves:

  Maximum number of overlaid distributions.

- col:

  Colors.

- add:

  Add to an existing plot.

- legend:

  Show a legend when several curves are drawn.

- main, xlab, ylab:

  Titles.

- ...:

  Passed to the plotting functions.

## Value

The input, invisibly.

## Examples

``` r
plot(distribution("gamma", mu = c(1, 2, 4), sigma = 0.5))

plot(distribution("negative_binomial_nb2", mu = 5, phi = 2))

plot(distribution("weibull", shape = 1.5, scale = 2), type = "hazard")
```
