# Power by simulation

Simulates responses from the fitted model, optionally with chosen
coefficient values (the effect size to detect) and a new design
(`newdata`, for another sample size), refits, and reports the share of
simulations whose Wald test of `term` rejects at `alpha`, with a
binomial confidence interval. This is the approach of the simr package.

## Usage

``` r
power_sim(
  object,
  term,
  nsim = 100,
  alpha = 0.05,
  coefficients = NULL,
  newdata = NULL,
  seed = 1
)
```

## Arguments

- object:

  A fitted model (frequentist).

- term:

  Coefficient name as in `coef(object)` (for example `"mu:x"`).

- nsim:

  Number of simulations.

- alpha:

  Test level.

- coefficients:

  Optional named vector of coefficient values to simulate from (other
  coefficients keep their estimates).

- newdata:

  Optional design to simulate on (must contain every covariate and the
  response column, whose values are replaced).

- seed:

  Random seed.

## Value

A one row data frame with power, its 95% interval, the number of
successful refits and of failures.

## Examples

``` r
d <- data.frame(x = rnorm(50)); d$y <- rpois(50, exp(0.5 + 0.3 * d$x))
f <- regcomp(y ~ x, family = "poisson", data = d)
power_sim(f, "lambda:x", nsim = 20, coefficients = c("lambda:x" = 0.2))
#>       term power     lower     upper nsim failed
#> 1 lambda:x   0.5 0.2719578 0.7280422   20      0
```
