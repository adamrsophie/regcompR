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
person_period(data, time, status, start = 1)
```

## Arguments

- data:

  Data frame with one row per subject.

- time:

  Name of the column with the last observed period (a positive integer).

- status:

  Name of the column with 1 for an event in period `time` and 0 for
  censoring after it.

- start:

  First period (default 1).

## Value

The expanded data frame with columns `period` and `event`.

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
```
