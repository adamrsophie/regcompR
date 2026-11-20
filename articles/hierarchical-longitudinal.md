# Hierarchical and longitudinal models

``` r

library(regcompR)
#> 
#> Attaching package: 'regcompR'
#> The following object is masked from 'package:base':
#> 
#>     capabilities
```

Random effects can enter any distribution parameter. Models are fitted
by the Laplace approximation (`engine = "ml"` or the faster
`engine = "tmb"`) or by MCMC (`engine = "bayes"`). Marginal models are
fitted by generalized estimating equations or with cluster robust
covariances.

## Random intercepts and slopes

The syntax follows lme4: `(1 | g)` for a random intercept, `(1 + x | g)`
for correlated intercepts and slopes, `(1 + x || g)` for independent
ones, `(1 | a/b)` for nested and `(1 | a) + (1 | b)` for crossed
factors.

``` r

G <- 30; m <- 8
dat <- data.frame(g = factor(rep(1:G, each = m)), x = rnorm(G * m))
u <- rnorm(G, 0, 0.6); v <- rnorm(G, 0, 0.3)
dat$y <- rpois(G * m, exp(0.5 + (0.4 + v[dat$g]) * dat$x + u[dat$g]))
fit <- regcomp(y ~ x + (1 + x | g), family = "poisson", data = dat)
summary(fit)
#> <regcomp_fit> Poisson regression, engine 'ml' (frequentist)
#>   marginal likelihood approximation: laplace 
#>   lambda (log link): ~x + (1 + x | g)
#>   log likelihood -417.032 on 5 parameters, n = 240
#> 
#> Wald estimates with 95% confidence intervals:
#>                                      estimate std_error   lower   upper       z
#> lambda:(Intercept)                     0.5495    0.0885  0.3760  0.7230  6.2086
#> lambda:x                               0.5157    0.0799  0.3591  0.6724  6.4524
#> lambda:1 + x | g:log_sd((Intercept))  -0.9795    0.2038 -1.3789 -0.5801 -4.8069
#> lambda:1 + x | g:log_sd(x)            -1.2338    0.2799 -1.7824 -0.6853 -4.4086
#> lambda:1 + x | g:chol(x,(Intercept))  -0.0504    0.3680 -0.7717  0.6709 -0.1370
#>                                      p_value
#> lambda:(Intercept)                    0.0000
#> lambda:x                              0.0000
#> lambda:1 + x | g:log_sd((Intercept))  0.0000
#> lambda:1 + x | g:log_sd(x)            0.0000
#> lambda:1 + x | g:chol(x,(Intercept))  0.8911
#> Hessian positive definite
head(ranef_regcomp(fit)[[1]])
#> $`1 + x | g`
#>     (Intercept)            x
#> 1   0.583929099 -0.180972803
#> 2   0.478751753  0.280042955
#> 3  -0.285983708  0.001543361
#> 4  -0.210814667  0.376985749
#> 5   0.335876984  0.149176474
#> 6  -0.473163656  0.010990201
#> 7  -0.113403071 -0.406488130
#> 8  -0.062558484 -0.051810670
#> 9  -0.096987898 -0.410290910
#> 10  0.302147000  0.137122424
#> 11  0.013581798 -0.181108213
#> 12 -0.041378676  0.035682969
#> 13  0.124956522 -0.049013906
#> 14 -0.017131363 -0.309312740
#> 15  0.004747919 -0.252650319
#> 16  0.414047746 -0.040981870
#> 17 -0.526162353  0.176744631
#> 18 -0.033395476  0.150885710
#> 19  0.081493809  0.070933571
#> 20 -0.273605249  0.069168185
#> 21  0.321316680 -0.208833315
#> 22  0.026052690 -0.137350104
#> 23  0.021885560  0.045963124
#> 24 -0.041586699  0.002324960
#> 25  0.834461443  0.370682108
#> 26  0.028317080  0.108106409
#> 27 -0.019611769  0.256640288
#> 28 -0.327077449  0.049022255
#> 29 -0.024232429 -0.174961272
#> 30 -0.328376876  0.162451922
plot(fit, which = "ranef")
```

![](hierarchical-longitudinal_files/figure-html/unnamed-chunk-3-1.png)![](hierarchical-longitudinal_files/figure-html/unnamed-chunk-3-2.png)

The estimates, standard errors and log likelihood agree with lme4 and
glmmTMB for the families they share.

## Random effects in other parameters

A random effect in the dispersion gives a model in which the clusters
differ in variability as well as in level.

``` r

dat$z <- rnorm(G * m, 1 + dat$x + u[dat$g], exp(-0.5 + v[dat$g]))
fit2 <- regcomp(z ~ x + (1 | g), sigma ~ 1 + (1 | g), family = "normal", data = dat)
fit2$theta
#>    mu:1 | g:log_sd((Intercept)) sigma:1 | g:log_sd((Intercept)) 
#>                      -0.6445014                      -1.2631436
```

## Covariance structures for longitudinal data

Structured random effects describe serial correlation within subjects:

- `ar1(0 + time | id)`: first order autoregressive over a factor of
  times;
- `cs(0 + time | id)`, `homcs()`: compound symmetry;
- [`diag()`](https://rdrr.io/r/base/diag.html), `homdiag()`, `us()`:
  diagonal and unstructured;
- `ou(t | id)`, `exp(t | id)`, `gau(t | id)`, `matern32(t | id)`:
  continuous time (Ornstein Uhlenbeck) and spatial kernels on numeric
  times or coordinates.

``` r

Tn <- 6; N <- 40
long <- data.frame(id = factor(rep(1:N, each = Tn)), time = factor(rep(1:Tn, N)),
                   t = rep(1:Tn, N), x = rnorm(N * Tn))
e <- unlist(lapply(1:N, function(i) as.numeric(arima.sim(list(ar = 0.6), Tn, sd = 0.8))))
long$y <- 1 + 0.5 * long$x + e + rnorm(N * Tn, 0, 0.3)
f_ar1 <- regcomp(y ~ x + ar1(0 + time | id), family = "normal", data = long, engine = "tmb")
f_ou <- regcomp(y ~ x + ou(t | id), family = "normal", data = long, engine = "tmb")
c(ar1 = AIC(f_ar1), ou = AIC(f_ou))
#>      ar1       ou 
#> 665.8609 665.8609
```

With equally spaced times the two models coincide
(`rho = exp(-1 / range)`), so their AIC values agree; the OU form also
handles irregular times.

The TMB engine uses automatic differentiation and is much faster for
large multilevel models; it gives the same estimates as the default
engine.

## Multiple membership

`mm(g1, g2, weights = cbind(w1, w2))` lets an observation belong to
several groups with given weights, as in school mobility models.

``` r

mmd <- data.frame(s1 = factor(sample(1:15, 200, TRUE)), s2 = factor(sample(1:15, 200, TRUE)),
                  w1 = 0.6, w2 = 0.4, x = rnorm(200))
eff <- rnorm(15, 0, 0.7)
mmd$y <- rnorm(200, 1 + mmd$x + 0.6 * eff[mmd$s1] + 0.4 * eff[mmd$s2], 1)
fmm <- regcomp(y ~ x + (1 | mm(s1, s2, weights = cbind(w1, w2))), family = "normal", data = mmd)
fmm$theta
#> mu:1 | mm(s1, s2, weights = cbind(w1, w2)):log_sd((Intercept)) 
#>                                                     -0.6498441
```

## Marginal models

Population averaged effects come from generalized estimating equations
with independence, exchangeable, AR(1) or unstructured working
correlations, with results equal to geepack.

``` r

gee(y ~ x, data = dat, id = g, family = poisson(), corstr = "exchangeable")
#> <regcomp_gee> poisson(log), exchangeable working correlation, 30 clusters
#>              estimate  robust_se        z
#> (Intercept) 0.6425348 0.08650442 7.427768
#> x           0.5016570 0.07146660 7.019460
#>   correlation 0.1899  scale 1.5815
```

Any likelihood fit gets cluster robust standard errors:

``` r

f0 <- regcomp(y ~ x, family = "poisson", data = dat)
sqrt(diag(vcov(f0, type = "sandwich", cluster = ~ g)))
#> lambda:(Intercept)           lambda:x 
#>         0.08985036         0.06947197
```

## Discrete time survival and frailty

[`person_period()`](https://adamrsophie.github.io/regcompR/reference/person_period.md)
expands survival data to one row per subject and period; a binary model
with a random intercept is a frailty model. The logit link gives the
proportional odds model; `links = list(prob = "cloglog")` gives grouped
proportional hazards.

``` r

n <- 400
sd0 <- data.frame(id = 1:n, cl = factor(rep(1:40, 10)), x = rnorm(n))
tev <- rgeom(n, plogis(-1.5 + 0.5 * sd0$x + rnorm(40, 0, 0.5)[sd0$cl])) + 1
sd0$t <- pmin(tev, 8); sd0$d <- as.numeric(tev <= 8)
pp <- person_period(sd0, "t", "d")
fh <- regcomp(event ~ factor(period) + x + (1 | cl), family = "bernoulli", data = pp)
coef(fh)["prob:x"]
#>    prob:x 
#> 0.5543546
```
