# Plot a fitted model

Diagnostic and interpretation plots:

- `"qq"`: normal QQ plot of randomized quantile residuals with a
  simulation envelope.

- `"pit"`: histogram of PIT values (uniform under a correct model).

- `"residuals"`: quantile residuals against the fitted mean.

- `"rootogram"`: hanging rootogram of observed and expected frequencies
  (discrete responses).

- `"ppc"`: observed response distribution against simulated replicates
  (posterior predictive check for Bayesian fits).

- `"effects"`: a distribution parameter or the mean against one
  covariate, other covariates at their mean or most frequent level, with
  a confidence or credible band.

- `"fitted"`: fitted distributions of a few observations.

- `"trace"`: trace plots of posterior draws (Bayesian fits).

- `"ranef"`: caterpillar plot of the random effects of each grouping
  term with approximate 95% intervals (conditional standard deviations
  from the Laplace Hessian, or posterior quantiles).

## Usage

``` r
# S3 method for class 'regcomp_fit'
plot(
  x,
  which = c("qq", "pit", "residuals", "rootogram"),
  variable = NULL,
  par = "mean",
  nsim = 50,
  rows = 1:4,
  ask = FALSE,
  ...
)
```

## Arguments

- x:

  A fitted model.

- which:

  Plot type(s).

- variable:

  Covariate for `"effects"`. Pass `by =` (a moderator) and optionally
  `at =` (its values; mean and mean plus and minus one standard
  deviation by default) to draw one curve per moderator value (an
  interaction plot).

- par:

  Distribution parameter for `"effects"`, or `"mean"`.

- nsim:

  Number of simulated replicates for `"qq"` and `"ppc"`.

- rows:

  Observations for `"fitted"`.

- ask:

  Ask before each new page when several plots are drawn.

- ...:

  Passed to plotting functions.

## Value

The input, invisibly.

## Details

For ordinal and nominal families `"effects"` draws the category
probabilities against the covariate.

## Examples

``` r
set.seed(1)
d <- data.frame(x = runif(200))
d$y <- rpois(200, exp(1 + d$x))
fit <- regcomp(y ~ x, family = "poisson", data = d)
plot(fit, which = "rootogram")

plot(fit, which = "effects", variable = "x")
```
