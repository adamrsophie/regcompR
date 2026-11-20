# Person period data for discrete time hazard models

Expands one row per subject into one row per subject and period at risk
(Allison 1982; Singer and Willett 2003; Goldstein 2011, ch. 9). The new
column `event` is 1 in the period of the event and 0 otherwise, and
`period` counts the periods. A binary regression of `event` on period
effects and covariates is the discrete time hazard model: the `logit`
link gives the proportional odds model and the `cloglog` link the
grouped proportional hazards model; random effects such as `(1 | id)`
give frailty and multilevel survival models.

## Usage

``` r
person_period(data, time, status, start = 1, cuts = NULL)
```

## Arguments

- data:

  Data frame with one row per subject.

- time:

  Name of the column with the last observed period (a positive integer),
  or the observed time when `cuts` is given.

- status:

  Name of the column with 1 for an event in period `time` and 0 for
  censoring after it.

- start:

  First period (default 1).

- cuts:

  `NULL` for integer periods, `"events"` for the distinct event times,
  or a sorted numeric vector of period end points. Subjects whose time
  is before the first cut are dropped, as in the Cox risk sets. The
  column `period_time` holds the cut of each row.

## Value

The expanded data frame with columns `period` and `event`.

## Details

With `cuts = "events"` the time is continuous and the periods are the
distinct event times: each subject has one row per event time at or
before its own time, the risk set of that time. A Poisson regression of
`event` on `0 + factor(period)` and covariates then gives exactly the
Cox proportional hazards estimates with Breslow ties and their standard
errors (Whitehead 1980; Laird and Olivier 1981); random effects give Cox
frailty models, and `engine = "bayes"` a Bayesian Cox model.

## Examples

``` r
d <- data.frame(id = 1:3, t = c(2, 3, 1), d = c(1, 0, 1), x = c(0.5, -1, 2))
pp <- person_period(d, "t", "d")
pp
#>   id t d    x period event
#> 1  1 2 1  0.5      1     0
#> 2  1 2 1  0.5      2     1
#> 3  2 3 0 -1.0      1     0
#> 4  2 3 0 -1.0      2     0
#> 5  2 3 0 -1.0      3     0
#> 6  3 1 1  2.0      1     1
fit <- regcomp(event ~ factor(period) + x, family = "bernoulli",
               links = list(prob = "cloglog"), data = person_period(
                 data.frame(t = rpois(200, 3) + 1, d = rbinom(200, 1, 0.7),
                            x = rnorm(200)), "t", "d"))
#> Warning: The Hessian is not positive definite; standard errors are unavailable or unreliable.
```
