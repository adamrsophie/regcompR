# Ordinal and categorical responses

``` r

library(regcompR)
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

Scoring rules for categorical responses use the class probabilities:

``` r

c(brier = mean(scoring_rule(fit, rule = "brier")), rps = mean(scoring_rule(fit, rule = "rps")))
#>     brier       rps 
#> 0.7037549 0.5429779
```

## Bivariate binary responses

`bivariate_binary_odds_ratio` models two binary outcomes with a
regression for each margin (`p1`, `p2`) and for their log odds ratio
(`psi`), as in VGAM’s `binom2.or`. The response is `cbind(y1, y2)`.

``` r

d$y1 <- rbinom(n, 1, plogis(d$x)); d$y2 <- rbinom(n, 1, plogis(0.5 * d$x + 0.4 * d$y1))
fb <- regcomp(cbind(y1, y2) ~ x, p2 ~ x, psi ~ 1, family = "bivariate_binary_odds_ratio", data = d)
coef(fb)
#>  p1:(Intercept)            p1:x  p2:(Intercept)            p2:x psi:(Intercept) 
#>       0.0679298       0.9583180       0.3143582       0.6379752       0.3864641
```

## Sample selection

`heckman_selection` fits the sample selection model by maximum
likelihood. The response is `cbind(selected, outcome)` (the outcome
value of unselected rows is ignored); the main formula is the outcome
equation and `gamma ~` gives the probit selection equation; `rho` is the
error correlation on the Fisher z scale.

``` r

z <- rnorm(n); e <- MASS::mvrnorm(n, c(0, 0), matrix(c(1, 0.5, 0.5, 1), 2))
sel <- as.numeric(0.3 + z + e[, 1] > 0)
hs <- data.frame(sel, out = ifelse(sel == 1, 1 + d$x + e[, 2], 0), x = d$x, z)
fh <- regcomp(cbind(sel, out) ~ x, gamma ~ z, family = "heckman_selection", data = hs)
coef(fh)
#>    mu:(Intercept)              mu:x gamma:(Intercept)           gamma:z 
#>        1.15521959        1.10947346        0.21387059        0.89044899 
#> sigma:(Intercept)   rho:(Intercept) 
#>       -0.03359496        0.33915035
tanh(coef(fh)["rho:(Intercept)"])     # the error correlation
#> rho:(Intercept) 
#>       0.3267187
```

## Partial proportional odds and multinomial probit

`partial_odds(~ z)` relaxes proportional odds for the covariates `z`:
their effects shift each threshold separately (parameters `d1`, `d2`,
…), as the nominal effects of
[`ordinal::clm()`](https://rdrr.io/pkg/ordinal/man/clm.html).

``` r

po <- data.frame(x = rnorm(400), z = rnorm(400)); u <- runif(400)
po$y <- factor(1 + (u > plogis(-1 - 0.8 * po$x)) + (u > plogis(1 + 0.5 * po$z - 0.8 * po$x)))
fpo <- regcomp(y ~ x, family = partial_odds(~ z), data = po)
coef(fpo)
#>       mu:x       d1:z       d2:z   cuts:1|2   cuts:2|3 
#>  0.8858083 -0.2667663  0.2123071 -1.0016558  0.9846142
anova(regcomp(y ~ x + z, family = "cumulative_logit", data = po), fpo)
#> Likelihood ratio tests
#> 
#>                                                            df  logLik    AIC
#> regcomp(y ~ x + z, family = "cumulative_logit", data = po)  4 -391.99 791.99
#> fpo                                                         5 -385.07 780.15
#>                                                             Chisq Chi Df
#> regcomp(y ~ x + z, family = "cumulative_logit", data = po)              
#> fpo                                                        13.842      1
#>                                                            Pr(>Chisq)    
#> regcomp(y ~ x + z, family = "cumulative_logit", data = po)               
#> fpo                                                         0.0001989 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

`multinomial_probit` is the nominal model with independent normal
utility errors; like the multinomial logit, it gives each category but
the first its own predictor.

``` r

xm <- rnorm(300)
ym <- max.col(cbind(0, 0.5 + xm, -0.3 - xm) + matrix(rnorm(900), 300))
coef(regcomp(ym ~ xm, family = "multinomial_probit", data = data.frame(xm, ym)))
#> eta2:(Intercept)          eta2:xm eta3:(Intercept)          eta3:xm 
#>        0.6907868        1.2721510       -0.3542442       -0.9623227
```
