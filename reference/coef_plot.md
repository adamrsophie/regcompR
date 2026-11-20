# Coefficient plot

Estimates with confidence (credible) intervals of one or several fits,
one row per coefficient, grouped by distribution parameter; several fits
are drawn side by side in different colors.

## Usage

``` r
coef_plot(
  ...,
  conf.level = 0.95,
  exponentiate = FALSE,
  intercept = FALSE,
  standardize = FALSE,
  vcov = "model"
)
```

## Arguments

- ...:

  Fitted models, optionally named.

- conf.level:

  Interval level.

- exponentiate:

  Exponentiate the estimates (then the reference line is at 1).

- intercept:

  Show intercepts.

- standardize:

  Standardize coefficients (see
  [`tidy.regcomp_fit()`](https://adamrsophie.github.io/regcompR/reference/tidy.regcomp_fit.md)).

- vcov:

  Covariance choice for frequentist fits.

## Value

The plotted data, invisibly.

## Examples

``` r
d <- data.frame(x = rnorm(200), z = rnorm(200)); d$y <- rnbinom(200, mu = exp(1 + 0.4 * d$x), size = 2)
coef_plot(Poisson = regcomp(y ~ x + z, family = "poisson", data = d),
          NB = regcomp(y ~ x + z, family = "nb2", data = d))
```
