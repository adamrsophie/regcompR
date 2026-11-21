# Censored data and survival models

``` r

library(regcompR)
```

## Right censoring and parametric survival

`surv(time, event)` marks right censored responses (with optional
delayed entry); `cens(lower, upper)` gives interval censoring and
truncation for any family, discrete or continuous.

``` r

n <- 400
d <- data.frame(x = rbinom(n, 1, 0.5), age = rnorm(n))
tt <- rweibull(n, shape = 1.4, scale = exp(1 + 0.6 * d$x - 0.2 * d$age))
cc <- runif(n, 0, 8)
d$time <- pmin(tt, cc); d$event <- as.numeric(tt <= cc)
fw <- regcomp(surv(time, event) ~ x + age, family = "weibull", data = d)
summary(fw)
#> <regcomp_fit> Weibull regression, engine 'ml' (frequentist)
#>   shape (log link): ~1
#>   scale (log link): ~x + age
#>   log likelihood -520.057 on 4 parameters, n = 400
#> 
#> Wald estimates with 95% confidence intervals:
#>                   estimate std_error   lower   upper       z p_value
#> shape:(Intercept)   0.3137    0.0510  0.2137  0.4136  6.1516       0
#> scale:(Intercept)   1.0871    0.0627  0.9643  1.2099 17.3480       0
#> scale:x             0.4561    0.0981  0.2639  0.6483  4.6510       0
#> scale:age          -0.2479    0.0467 -0.3396 -0.1563 -5.3038       0
#> max |gradient| 1.66e-03; Hessian positive definite
```

The main formula models the Weibull scale, so the estimates match
[`survival::survreg()`](https://rdrr.io/pkg/survival/man/survreg.html)
on the log time scale. Other choices include the lognormal, log
logistic, gamma, generalized gamma and Gompertz families; any parameter
can take a formula.

``` r

fl <- regcomp(surv(time, event) ~ x + age, family = "lognormal", data = d)
AIC(fw, fl)
#>    df      AIC
#> fw  4 1048.114
#> fl  4 1053.875
nd <- data.frame(x = c(0, 1), age = 0)
predict(fw, nd, type = "survival", at = c(1, 3, 5))
#>   survival(1) survival(3) survival(5)
#> 1   0.7977997   0.3620883   0.1295408
#> 2   0.8860239   0.5803113   0.3345988
predict(fw, nd, type = "median")
#> [1] 2.268868 3.580171
```

## Cure models

When a fraction of subjects never has the event,
[`cure()`](https://adamrsophie.github.io/regcompR/reference/cure.md)
gives the mixture cure model `S(t) = cure + (1 - cure) S0(t)`. A formula
for `cure` is the incidence model; the results match flexsurvcure.

``` r

cured <- runif(n) < plogis(-0.5 + d$x)
tc <- ifelse(cured, Inf, rweibull(n, 1.4, 2))
d$time2 <- pmin(tc, cc); d$event2 <- as.numeric(tc <= cc)
fc <- regcomp(surv(time2, event2) ~ 1, cure ~ x, family = cure("weibull"), data = d)
coef(fc)
#> shape:(Intercept) scale:(Intercept)  cure:(Intercept)            cure:x 
#>         0.4036345         0.5431813        -0.3634276         1.1080392
```

## Interval censoring

Event times known only to the year, with right censored subjects:

``` r

d$lo <- ifelse(d$event == 1, floor(d$time), d$time)
d$hi <- ifelse(d$event == 1, floor(d$time) + 1, Inf)
fi <- regcomp(cens(lo, hi) ~ x + age, family = "weibull", data = d)
coef(fi)
#> shape:(Intercept) scale:(Intercept)           scale:x         scale:age 
#>         0.3042646         1.1021716         0.4452103        -0.2499371
```

## Truncation and delayed entry

`surv(time, event, entry)` handles delayed entry (left truncation):
subjects enter observation at `entry` and contribute the conditional
likelihood given survival to that time.
`cens(lower, upper, trunc_lower, trunc_upper)` gives truncation limits
for any family, together with interval censoring.

``` r

d$entry <- runif(n, 0, 1)
keep <- d$time > d$entry
fe <- regcomp(surv(time, event, entry) ~ x + age, family = "weibull", data = d[keep, ])
coef(fe)
#> shape:(Intercept) scale:(Intercept)           scale:x         scale:age 
#>         0.1950867         1.0196206         0.5200054        -0.2714023
# Event times observed only when below 6 (right truncation).
tt2 <- rweibull(n, shape = 1.4, scale = exp(1 + 0.6 * d$x))
dt <- data.frame(x = d$x, tr = tt2)[tt2 < 6, ]
ft <- regcomp(cens(tr, tr, trunc_upper = 6) ~ x, family = "weibull", data = dt)
coef(ft)
#> shape:(Intercept) scale:(Intercept)           scale:x 
#>         0.2439162         0.9628687         0.9734912
```

The
[`censor()`](https://adamrsophie.github.io/regcompR/reference/censor.md)
operator gives the law of a censored variable itself, with point masses
at the censoring limits (a Tobit model when the base is normal):

``` r

tob <- data.frame(x = rnorm(200)); tob$y <- pmax(0, 0.5 + tob$x + rnorm(200))
coef(regcomp(y ~ x, family = censor("normal", lower = 0), data = tob))
#>    mu:(Intercept)              mu:x sigma:(Intercept) 
#>         0.4277192         0.9548162         0.1095207
```

## Cox proportional hazards and frailty

`person_period(cuts = "events")` builds the risk set of every distinct
event time. A Poisson regression of `event` on one intercept per event
time and the covariates gives exactly the Cox estimates with Breslow
ties and their standard errors (Whitehead 1980). Random effects then
give Cox frailty models (the `"tmb"` engine is fastest here), and
`engine = "bayes"` a Bayesian Cox model.

``` r

lung <- na.omit(survival::lung[, c("time", "status", "age", "sex", "inst")])
lung$status <- lung$status - 1
pp <- person_period(lung, "time", "status", cuts = "events")
fc <- regcomp(event ~ 0 + factor(period) + age + sex, family = "poisson", data = pp)
rbind(regcompR = coef(fc)[c("lambda:age", "lambda:sex")],
      coxph = coef(survival::coxph(survival::Surv(time, status) ~ age + sex, lung, ties = "breslow")))
#>          lambda:age lambda:sex
#> regcompR 0.01700002 -0.5109971
#> coxph    0.01700005 -0.5109966
fr <- regcomp(event ~ 0 + factor(period) + age + sex + (1 | inst), family = "poisson",
              data = transform(pp, inst = factor(inst)), engine = "tmb")
st <- summary(fr)$table
st[!grepl("period", rownames(st)), c("estimate", "std_error")]
#>                                        estimate    std_error
#> lambda:age                           0.01700007 9.222569e-03
#> lambda:sex                          -0.51099724 1.676826e-01
#> lambda:1 | inst:log_sd((Intercept)) -9.78357380 1.901896e+03
```

## Discrete time and frailty

See
[`vignette("hierarchical-longitudinal")`](https://adamrsophie.github.io/regcompR/articles/hierarchical-longitudinal.md)
for discrete time hazard models built with
[`person_period()`](https://adamrsophie.github.io/regcompR/reference/person_period.md)
and for shared frailty through random effects.
