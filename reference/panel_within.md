# Fixed effects (within) estimator for panel data

Least squares after subtracting the unit means (and, for two way
effects, after also projecting out time dummies), as
`plm(model = "within")`. Equal to least squares with unit (and time)
dummies, without building them. A three part formula `y ~ x | z` gives
the within (fixed effects) instrumental variable estimator: two stage
least squares on the demeaned data, as
`plm(y ~ x | z, model = "within")`.

## Usage

``` r
panel_within(
  formula,
  data,
  id,
  time = NULL,
  effect = c("individual", "twoways")
)
```

## Arguments

- formula:

  A formula `y ~ x1 + x2` without an intercept term (absorbed by the
  unit effects), or `y ~ regressors | instruments` for instrumental
  variables.

- data:

  A data frame.

- id, time:

  Names of the unit and time variables in `data`.

- effect:

  `"individual"` or `"twoways"`.

## Value

An object of class `regcomp_lm`.

## Examples

``` r
set.seed(1)
d <- data.frame(id = rep(1:50, each = 5), t = rep(1:5, 50), x = rnorm(250))
a <- rnorm(50); d$x <- d$x + a[d$id]; d$y <- 1 + 0.5 * d$x + a[d$id] + rnorm(250)
panel_within(y ~ x, data = d, id = "id")
#> <regcomp_lm> within estimator (individual effects), n = 250, residual df = 199
#>   estimate std_error      t p_value
#> x   0.5026    0.0773 6.5003       0
```
