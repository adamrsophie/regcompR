# Bayesian distributional regression

``` r

library(regcompR)
```

`engine = "bayes"` fits the same model by MCMC. regcompR writes the Stan
program (run by bayescompR, cmdstanr or rstan) or a BUGS program (run by
JAGS through R2jags or by OpenBUGS through R2OpenBUGS).

## The generated program

``` r

d <- data.frame(x = rnorm(100), g = factor(rep(1:10, each = 10)))
d$y <- rpois(100, exp(0.5 + 0.3 * d$x))
cat(stan_code(y ~ x + (1 | g), family = "poisson", data = d)$code)
#> functions {
#> 
#> }
#> data {
#>   int<lower=1> N;
#>   array[N] int y;
#>   vector[N] w;
#>   int<lower=0> P_lambda;
#>   matrix[N, P_lambda] X_lambda;
#>   vector[N] off_lambda;
#>   int<lower=1> L_1;
#>   int<lower=1> q_1;
#>   array[N] int<lower=1, upper=L_1> g_1;
#>   matrix[N, q_1] Xr_1;
#> }
#> parameters {
#>   vector[P_lambda] b_lambda;
#>   vector<lower=0>[q_1] sd_1;
#>   matrix[q_1, L_1] z_1;
#> }
#> transformed parameters {
#>   matrix[q_1, L_1] r_1 = diag_pre_multiply(sd_1, z_1);
#> }
#> model {
#>   vector[N] eta_lambda = X_lambda * b_lambda + off_lambda;
#>   for (i in 1:N) eta_lambda[i] += Xr_1[i] * col(r_1, g_1[i]);
#>   vector[N] lambda = exp(eta_lambda);
#>   b_lambda[1] ~ normal(0, 5);
#>   b_lambda[2] ~ normal(0, 5);
#>   to_vector(z_1) ~ std_normal();
#>   sd_1 ~ student_t(3, 0, 2.5);
#>   for (i in 1:N) target += w[i] * poisson_lpmf(y[i] | lambda[i]);
#> }
```

## Priors

Default priors are weakly informative.
[`prior()`](https://adamrsophie.github.io/regcompR/reference/prior.md)
sets priors by class (`b`, `Intercept`, `sd`, `cor`), by distribution
parameter and by coefficient.

``` r

pr <- c(prior("normal(0, 1)", class = "b"), prior("exponential(2)", class = "sd"))
pr
#> normal(0, 1)                 class b         par (all)  coef (all)
#> exponential(2)               class sd        par (all)  coef (all)
```

## Fitting

The examples below need a working Stan or JAGS installation.

``` r

fit <- regcomp(y ~ x + (1 | g), family = "poisson", data = d, engine = "bayes",
               priors = pr, control = regcomp_control(runner = "cmdstanr", chains = 4))
summary(fit)
plot(fit, which = c("trace", "ppc"))
predict(fit, type = "mean", interval = "confidence")

fit_jags <- regcomp(y ~ x + (1 | g), family = "poisson", data = d, engine = "bayes",
                    control = regcomp_control(runner = "jags"))
```

Every family with a Stan density, every operator built from such
families, ordinal and nominal models, all random effect covariance
structures and multiple membership terms are supported by the Stan
backends. The BUGS backends support the families with BUGS densities and
unstructured random effects;
[`capabilities()`](https://rdrr.io/r/base/capabilities.html) reports
what is available for a family.

``` r

capabilities("generalized_poisson")
#> named logical(0)
```
