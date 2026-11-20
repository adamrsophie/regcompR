# Ordinal and categorical responses

``` r

library(regcompR)
#> 
#> Attaching package: 'regcompR'
#> The following object is masked from 'package:base':
#> 
#>     capabilities
families("ordinal")
#>  [1] "adjacent_category_logit"    "continuation_ratio_cloglog"
#>  [3] "continuation_ratio_logit"   "continuation_ratio_probit" 
#>  [5] "cumulative_cauchit"         "cumulative_cloglog"        
#>  [7] "cumulative_logit"           "cumulative_loglog"         
#>  [9] "cumulative_probit"          "stereotype"
```

## Cumulative link models

An ordered factor response with a cumulative family gives the
proportional odds model (logit), the ordered probit, and the cloglog,
loglog and cauchit variants. The thresholds are estimated on their
natural scale; the main formula has no intercept.

``` r

n <- 400
d <- data.frame(x = rnorm(n), z = rnorm(n))
d$g <- factor(rep(1:40, each = 10))
d$y <- cut(0.8 * d$x - 0.5 * d$z + rnorm(40, 0, 0.7)[d$g] + rlogis(n), c(-Inf, -1, 0.5, 2, Inf),
           labels = c("low", "mid", "high", "top"), ordered_result = TRUE)
fit <- regcomp(y ~ x + z, family = "cumulative_logit", data = d)
summary(fit)
#> <regcomp_fit> Cumulative logit (ordinal) regression, engine 'ml' (frequentist)
#>   mu (identity link): ~x + z
#>   log likelihood -516.412 on 5 parameters, n = 400
#> 
#> Wald estimates with 95% confidence intervals:
#>               estimate std_error   lower   upper       z p_value
#> mu:x            0.6056    0.0998  0.4099  0.8013  6.0653       0
#> mu:z           -0.3720    0.0863 -0.5412 -0.2028 -4.3091       0
#> cuts:low|mid   -0.7929    0.1116 -1.0116 -0.5743 -7.1064       0
#> cuts:mid|high   0.4777    0.1074  0.2671  0.6883  4.4460       0
#> cuts:high|top   1.8334    0.1444  1.5505  2.1164 12.7003       0
#> max |gradient| 1.63e-04; Hessian positive definite
```

The estimates and standard errors equal those of
[`MASS::polr()`](https://rdrr.io/pkg/MASS/man/polr.html). A formula for
`scale` gives the location scale model of
[`ordinal::clm()`](https://rdrr.io/pkg/ordinal/man/clm.html):

``` r

fit_ls <- regcomp(y ~ x + z, scale ~ x, family = "cumulative_logit", data = d)
AIC(fit, fit_ls)
#>        df      AIC
#> fit     5 1042.824
#> fit_ls  6 1043.989
```

Category probabilities:

``` r

predict(fit, data.frame(x = c(-1, 0, 1), z = 0), type = "probability")
#>      probability(low) probability(mid) probability(high) probability(top)
#> [1,]        0.4532963        0.2938149         0.1726431        0.0802457
#> [2,]        0.3115365        0.3056633         0.2449699        0.1378303
#> [3,]        0.1980494        0.2700178         0.3053738        0.2265590
plot(fit, which = "effects", variable = "x")
```

![](ordinal-categorical_files/figure-html/unnamed-chunk-5-1.png)

Continuation ratio (`continuation_ratio_logit`, `_probit`, `_cloglog`),
adjacent category (`adjacent_category_logit`) and stereotype models are
also available, and random effects work as for any other family:

``` r

fit_m <- regcomp(y ~ x + z + (1 | g), family = "cumulative_logit", data = d)
exp(fit_m$theta)   # random intercept standard deviation (true value 0.7)
#> mu:1 | g:log_sd((Intercept)) 
#>                    0.5912834
```

## Nominal responses

`multinomial_logit` fits the baseline category logit model of
[`nnet::multinom()`](https://rdrr.io/pkg/nnet/man/multinom.html); the
main formula is shared by every non-baseline category.

``` r

d$k <- factor(sample(c("a", "b", "c"), n, TRUE, prob = c(0.5, 0.3, 0.2)))
fm <- regcomp(k ~ x, family = "multinomial_logit", data = d)
coef(fm)
#> eta2:(Intercept)           eta2:x eta3:(Intercept)           eta3:x 
#>       -0.5299411       -0.1745779       -1.0110616        0.3332306
```

## Bivariate binary responses and selection

`bivariate_binary_odds_ratio` models two binary outcomes with an odds
ratio regression (as in VGAM’s `binom2.or`), and `heckman_selection`
fits the sample selection model with response
`cbind(selected, outcome)`.
