# Distributions and operators

``` r

library(regcompR)
#> 
#> Attaching package: 'regcompR'
#> The following object is masked from 'package:base':
#> 
#>     capabilities
```

## The registry

[`registry_table()`](https://adamrsophie.github.io/regcompR/reference/registry_table.md)
lists every family with its category, parameters, status and the
backends that support it.

``` r

reg <- registry_table()
table(reg$category)
#> 
#>                bounded            categorical               circular 
#>                     35                      3                     10 
#>                  count            directional        finite_discrete 
#>                     51                      4                     16 
#>                  mixed           multivariate                ordinal 
#>                      2                     22                     10 
#>               positive quantile_parameterized                   real 
#>                     64                      5                     56
head(reg[, c("id", "category", "parameters", "status")])
#>                        id category      parameters    status
#> 1 adjacent_category_logit  ordinal mu, cuts, scale validated
#> 2                 arcsine  bounded            a, b validated
#> 3      asymmetric_laplace     real  mu, sigma, tau validated
#> 4         axial_von_mises circular       mu, kappa validated
#> 5         balding_nichols  bounded            F, p validated
#> 6                   bates  bounded               n validated
capabilities("conway_maxwell_poisson")
#> $log_prob
#> [1] "analytic"
#> 
#> $cdf
#> [1] "analytic"
#> 
#> $quantile
#> [1] "numeric"
#> 
#> $rng
#> [1] "numeric"
#> 
#> $mean
#> [1] "analytic"
#> 
#> $variance
#> [1] "analytic"
#> 
#> $skewness
#> [1] "numeric"
#> 
#> $kurtosis
#> [1] "numeric"
#> 
#> $score
#> [1] "numeric"
#> 
#> $hessian
#> [1] "numeric"
#> 
#> $regression
#> [1] TRUE
#> 
#> $ml
#> [1] TRUE
#> 
#> $stan
#> [1] "generated"
#> 
#> $tmb
#> [1] "unavailable"
#> 
#> $bugs
#> [1] "unavailable"
#> 
#> $status
#> [1] "validated"
```

Families are validated against reference implementations: densities
integrate or sum to one, distribution functions agree with numerical
integration, and random draws agree with the stated moments.

## Operators

Operators build new families from registered ones. The result is a
family like any other: it can be used in
[`distribution()`](https://adamrsophie.github.io/regcompR/reference/distribution.md)
and
[`regcomp()`](https://adamrsophie.github.io/regcompR/reference/regcomp.md).

``` r

zinb <- inflate("nb2")          # zero inflation
names(zinb$pars)
#> [1] "mu"  "phi" "zi"
hnb <- hurdle("poisson")                        # hurdle at zero
tn <- truncate_dist("normal", lower = 0)        # truncation
mx <- mixture("normal", "normal")               # finite mixture
d <- distribution(mx, mu.1 = -2, sigma.1 = 1, mu.2 = 2, sigma.2 = 0.7, w2 = 0.4)
plot(d)
```

![](distributions_files/figure-html/unnamed-chunk-4-1.png)

Further operators:

- [`censor()`](https://adamrsophie.github.io/regcompR/reference/censor.md),
  [`discretize()`](https://adamrsophie.github.io/regcompR/reference/discretize.md),
  [`round_dist()`](https://adamrsophie.github.io/regcompR/reference/discretize.md),
  [`transform_dist()`](https://adamrsophie.github.io/regcompR/reference/transform_dist.md),
  [`wrap()`](https://adamrsophie.github.io/regcompR/reference/wrap.md),
  [`weighted()`](https://adamrsophie.github.io/regcompR/reference/weighted.md);
- [`compound()`](https://adamrsophie.github.io/regcompR/reference/compound.md)
  for mixtures over a parameter, such as the Poisson lognormal law;
- [`cure()`](https://adamrsophie.github.io/regcompR/reference/cure.md)
  for mixture cure survival models;
- [`splice()`](https://adamrsophie.github.io/regcompR/reference/splice.md)
  for body and tail composite laws;
- [`convolve_dist()`](https://adamrsophie.github.io/regcompR/reference/convolve_dist.md)
  for sums of independent variables;
- [`copula_dist()`](https://adamrsophie.github.io/regcompR/reference/copula_dist.md)
  for bivariate copula models with any margins.

``` r

emg <- convolve_dist("normal", "exponential")
prob(distribution(emg, mu.1 = 0, sigma.1 = 1, lambda.2 = 2), c(0, 1, 2))
#> [1] 0.3362040 0.3173105 0.1353353
sp <- splice("lognormal", "lomax")
plot(distribution(sp, mu.b = 0, sigma.b = 0.6, scale.t = 2, alpha.t = 2.5,
                  threshold = 2, w = 0.7), xlim = c(0, 8))
```

![](distributions_files/figure-html/unnamed-chunk-5-1.png)

## Generators

Generators such as
[`kumaraswamy_g()`](https://adamrsophie.github.io/regcompR/reference/exponentiated.md),
[`beta_g()`](https://adamrsophie.github.io/regcompR/reference/exponentiated.md),
[`marshall_olkin()`](https://adamrsophie.github.io/regcompR/reference/exponentiated.md),
[`transmuted()`](https://adamrsophie.github.io/regcompR/reference/exponentiated.md)
and
[`exponentiated()`](https://adamrsophie.github.io/regcompR/reference/exponentiated.md)
produce the generalized families of the distribution literature from any
baseline.

``` r

kw <- kumaraswamy_g("weibull")
names(kw$pars)
#> [1] "shape" "scale" "a"     "b"
```

## Custom families

[`custom_family()`](https://adamrsophie.github.io/regcompR/reference/custom_family.md)
registers a law from its log density; regcompR checks that it is
normalized.

``` r

myexp <- custom_family("my_exponential", list(par_spec("rate", "positive", role = "rate")),
  lpdf = function(y, p) dexp(y, p$rate, log = TRUE),
  support = function(p) list(lower = 0, upper = Inf), check = list(rate = 2))
coef(regcomp(y ~ 1, family = myexp, data = data.frame(y = rexp(200, 2))))
#> rate:(Intercept) 
#>         0.693163
```
