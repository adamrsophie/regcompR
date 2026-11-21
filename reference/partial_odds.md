# Partial proportional odds family

A cumulative link model whose thresholds depend on covariates:
`P(Y <= k) = F(theta_k + z' gamma_k - mu)`, where `mu` is the
proportional (parallel) part of the main formula and `nonparallel` gives
the covariates `z` with category specific effects. This is the nominal
effects model of `ordinal::clm(nominal = )` and the partial proportional
odds model of Peterson and Harrell (1990). In
[`regcomp()`](https://adamrsophie.github.io/regcompR/reference/regcomp.md)
the family adds parameters `d1`, ..., `d(K-1)` (the threshold shifts
`z' gamma_k`), each with the `nonparallel` formula unless given its own.

## Usage

``` r
partial_odds(
  nonparallel,
  link = c("logit", "probit", "cloglog", "loglog", "cauchit")
)
```

## Arguments

- nonparallel:

  One sided formula of the covariates with category specific effects. A
  covariate should not also appear in the main formula.

- link:

  `"logit"`, `"probit"`, `"cloglog"`, `"loglog"` or `"cauchit"`.

## Value

A family object for
[`regcomp()`](https://adamrsophie.github.io/regcompR/reference/regcomp.md).

## Examples

``` r
set.seed(1)
d <- data.frame(x = rnorm(400), z = rnorm(400))
eta <- 0.8 * d$x
u <- runif(400)
d$y <- factor(1 + (u > plogis(-1 - eta)) + (u > plogis(1 + 0.5 * d$z - eta)))
f <- regcomp(y ~ x, family = partial_odds(~ z), data = d)
coef(f)
#>        mu:x        d1:z        d2:z    cuts:1|2    cuts:2|3 
#>  0.86086728  0.07016505  0.66578264 -0.96785350  0.96581296 
```
