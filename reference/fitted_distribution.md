# Fitted distributions

Returns the conditional distribution of the response for each row of
`newdata` (or the data used for fitting) as a vectorized distribution
object, so every method of
[`distribution()`](https://adamrsophie.github.io/regcompR/reference/distribution.md)
(quantiles, probabilities, moments, plots) applies to a fitted model.
For Bayesian fits the parameters are posterior means of each
distribution parameter; use
[`predict.regcomp_fit()`](https://adamrsophie.github.io/regcompR/reference/predict.regcomp_fit.md)
for full posterior summaries.

## Usage

``` r
fitted_distribution(fit, newdata = NULL, re = c("conditional", "zero"))
```

## Arguments

- fit:

  A fitted model.

- newdata:

  Optional data frame.

- re:

  Random effects: `"conditional"` uses the estimated random effects (new
  levels get zero), `"zero"` sets all random effects to zero.

## Value

A `regcomp_dist` object.

## Examples

``` r
set.seed(1)
d <- data.frame(x = runif(100))
d$y <- rgamma(100, shape = 4, rate = 4 / exp(1 + d$x))
fit <- regcomp(y ~ x, family = "gamma", data = d)
fd <- fitted_distribution(fit, data.frame(x = c(0, 1)))
quantile(fd, 0.9)
#> [1]  4.32943 12.24137
```
