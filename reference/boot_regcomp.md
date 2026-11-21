# Bootstrap a fitted model

The parametric bootstrap simulates responses from the fitted model (with
new random effects, as
[`lme4::bootMer()`](https://rdrr.io/pkg/lme4/man/bootMer.html) by
default) and refits; the case bootstrap resamples rows, or whole
clusters with `cluster`. Returns the replicates of `statistic` with
bootstrap standard errors and percentile intervals.

## Usage

``` r
boot_regcomp(
  object,
  R = 199,
  type = c("parametric", "case"),
  statistic = stats::coef,
  cluster = NULL,
  level = 0.95,
  seed = 1
)
```

## Arguments

- object:

  A fitted model (frequentist).

- R:

  Number of bootstrap replicates.

- type:

  `"parametric"` or `"case"`.

- statistic:

  Function of a fit returning a numeric vector; default the
  coefficients.

- cluster:

  For the case bootstrap, a one sided formula or vector naming the
  clusters to resample.

- level:

  Interval level.

- seed:

  Random seed.

## Value

An object of class `regcomp_boot` with `t0`, the replicate matrix `t`
and a `table` of estimates, bootstrap standard errors and percentile
intervals. Failed refits are dropped and counted.

## Examples

``` r
d <- data.frame(x = rnorm(60)); d$y <- rpois(60, exp(0.5 + 0.3 * d$x))
f <- regcomp(y ~ x, family = "poisson", data = d)
boot_regcomp(f, R = 20)
#> <regcomp_boot> parametric bootstrap, 20 replicates
#>                    estimate std_error  lower  upper
#> lambda:(Intercept)   0.4645    0.1249 0.1848 0.6006
#> lambda:x             0.2735    0.1508 0.0496 0.4871
```
