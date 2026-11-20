# Stepwise selection of terms

Forward, backward or bidirectional selection of the terms of one
distribution parameter's formula by AIC or BIC, as in
[`stats::step()`](https://rdrr.io/r/stats/step.html) and StepReg.
Selection makes the reported standard errors, intervals and p values of
the final model too optimistic (Harrell 2015, sec. 4.3); prefer
penalized estimation (see
[`penalty()`](https://adamrsophie.github.io/regcompR/reference/penalty.md))
or cross validation for prediction, and treat stepwise results as
exploratory.

## Usage

``` r
step_regcomp(
  object,
  scope = NULL,
  direction = c("both", "backward", "forward"),
  criterion = c("AIC", "BIC"),
  par = NULL,
  trace = FALSE
)
```

## Arguments

- object:

  A frequentist fit: the starting model.

- scope:

  A one sided formula with the largest set of terms considered for
  `par`; by default the terms of the starting model (backward only).

- direction:

  `"both"`, `"backward"` or `"forward"`.

- criterion:

  `"AIC"` or `"BIC"`.

- par:

  Distribution parameter whose formula is searched; the main one by
  default. Random effect terms are kept.

- trace:

  Print each step.

## Value

The selected fit, with the path of steps in `attr(, "steps")`.

## Examples

``` r
d <- data.frame(x1 = rnorm(200), x2 = rnorm(200), x3 = rnorm(200))
d$y <- rpois(200, exp(0.5 + 0.4 * d$x1))
f0 <- regcomp(y ~ 1, family = "poisson", data = d)
fs <- step_regcomp(f0, scope = ~ x1 + x2 + x3, direction = "forward", criterion = "BIC")
formula(fs)
#> $lambda
#> y ~ x1
#> <environment: 0x5893b22b25f8>
#> 
#> attr(,"given")
#> [1] "lambda"
attr(fs, "steps")
#>       step      BIC
#> 1    start 690.8945
#> + x1  + x1 646.6848
```
