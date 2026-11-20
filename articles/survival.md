# Censored data and survival models

``` r

library(regcompR)
#> 
#> Attaching package: 'regcompR'
#> The following object is masked from 'package:base':
#> 
#>     capabilities
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

## Discrete time and frailty

See
[`vignette("hierarchical-longitudinal")`](https://adamrsophie.github.io/regcompR/articles/hierarchical-longitudinal.md)
for discrete time hazard models built with
[`person_period()`](https://adamrsophie.github.io/regcompR/reference/person_period.md)
and for shared frailty through random effects.
