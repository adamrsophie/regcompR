# Best subset selection

Fits every subset of the fixed effect terms of one parameter (keeping
random effects) and ranks them by AIC, AICc, BIC or the Hannan-Quinn
criterion. Practical for up to about 15 candidate terms.

## Usage

``` r
best_subset(
  object,
  scope = NULL,
  criterion = c("AIC", "AICc", "BIC", "HQ"),
  par = NULL,
  top = 5
)
```

## Arguments

- object:

  A fitted model.

- scope:

  Optional one sided formula of candidate terms (default the terms of
  the fit).

- criterion:

  `"AIC"`, `"AICc"`, `"BIC"` or `"HQ"`.

- par:

  Parameter whose formula is searched (default the first).

- top:

  Number of best models to return.

## Value

A data frame of the best models (terms, criterion, df), with the best
fit as attribute `"fit"`.

## Examples

``` r
d <- data.frame(x1 = rnorm(80), x2 = rnorm(80), x3 = rnorm(80))
d$y <- rpois(80, exp(0.5 + 0.4 * d$x1))
best_subset(regcomp(y ~ x1 + x2 + x3, family = "poisson", data = d), criterion = "BIC")
#>          terms      BIC df
#> 1           x1 257.7691  2
#> 2      x1 + x2 260.2847  3
#> 3      x1 + x3 261.6825  3
#> 4 x1 + x2 + x3 264.2270  4
#> 5       (none) 278.3815  1
```
