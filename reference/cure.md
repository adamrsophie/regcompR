# Mixture cure models

A fraction `cure` of subjects never has the event:
`S(t) = cure + (1 - cure) S0(t)` and `f(t) = (1 - cure) f0(t)` (Boag
1949; Berkson and Gage 1952; Farewell 1982). The cure fraction is a
parameter in (0, 1) with a logit link, so `cure ~ z` gives the incidence
model and the base parameters give the latency model. Use with right
censored responses
([`surv()`](https://adamrsophie.github.io/regcompR/reference/cens.md));
the law is improper (its distribution function tends to `1 - cure`).

## Usage

``` r
cure(base)
```

## Arguments

- base:

  A family for the event time of uncured subjects.

## Value

A `regcomp_family`.

## Examples

``` r
cw <- cure("weibull")
set.seed(1)
n <- 400; cured <- runif(n) < 0.3
t <- ifelse(cured, Inf, rweibull(n, 1.5, 2)); cens <- runif(n, 0, 6)
d <- data.frame(time = pmin(t, cens), event = as.numeric(t <= cens))
fit <- regcomp(surv(time, event) ~ 1, family = cw, data = d)
plogis(coef(fit)["cure:(Intercept)"])
#> cure:(Intercept) 
#>        0.2884821 
```
