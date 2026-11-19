# Censored, truncated and interval observations

Wrap the response in `cens()` to describe the observation process. A row
with `lower == upper` is observed exactly. `upper = Inf` gives right
censoring at `lower`, `lower = -Inf` left censoring at `upper`, and
finite unequal bounds interval censoring. For discrete responses an
interval `[lower, upper]` includes both end points. `trunc_lower` and
`trunc_upper` give per observation truncation (delayed entry), with the
retained region `(trunc_lower, trunc_upper]` for continuous responses
and `[trunc_lower, trunc_upper]` for discrete ones. `surv()` is a
shortcut for right censored survival times with an event indicator.

## Usage

``` r
cens(lower, upper = lower, trunc_lower = -Inf, trunc_upper = Inf)

surv(time, event, entry = -Inf)
```

## Arguments

- lower, upper:

  Lower and upper bounds of the observed value.

- trunc_lower, trunc_upper:

  Truncation bounds.

- time:

  Event or censoring time.

- event:

  1 for an observed event, 0 for right censoring.

- entry:

  Optional delayed entry (left truncation) time.

## Value

An object of class `regcomp_cens`.

## Details

Likelihood contributions: exact `log f(y)`, left censored `log F(u)`,
right censored `log S(l)`, interval `log(F(u) - F(l))`, each minus
`log(F(tu) - F(tl))` under truncation.

## Examples

``` r
cens(c(1, 2, 3), c(1, Inf, 4))
#> [1] 1        2+       [3,   4]
surv(c(5, 8), c(1, 0))
#> [1] 5  8+
```
