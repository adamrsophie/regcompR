# Bayesian distributional regression

``` r

library(regcompR)
```

`engine = "bayes"` fits the same model by MCMC. regcompR writes a Stan
program, run by bayescompR (`runner = "bayescompR"`), cmdstanr
(`"cmdstanr"`) or rstan (`"rstan"`), or a BUGS program, run by JAGS
(`"jags"`, through R2jags) or OpenBUGS (`"openbugs"`, through
R2OpenBUGS).

## The generated programs

[`stan_code()`](https://adamrsophie.github.io/regcompR/reference/stan_code.md)
and
[`bugs_code()`](https://adamrsophie.github.io/regcompR/reference/bugs_code.md)
return the program and its data without running anything.

``` r

d <- data.frame(x = rnorm(100), g = factor(rep(1:10, each = 10)))
d$y <- rpois(100, exp(0.5 + 0.3 * d$x))
sc <- stan_code(y ~ x + (1 | g), family = "poisson", data = d)
cat(sc$code)
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
names(sc$data)
#>  [1] "N"          "w"          "y"          "P_lambda"   "X_lambda"  
#>  [6] "off_lambda" "L_1"        "q_1"        "g_1"        "Xr_1"
cat(bugs_code(y ~ x + (1 | g), family = "poisson", data = d)$code)
#> model {
#>   for (i in 1:N) {
#>     eta_lambda[i] <- inprod(X_lambda[i, ], b_lambda[]) + off_lambda[i] + inprod(Xr_1[i, ], r_1[g_1[i], ])
#>     lambda[i] <- exp(eta_lambda[i])
#>     y[i] ~ dpois(lambda[i])
#>   }
#>   b_lambda[1] ~ dnorm(0, 0.04)
#>   b_lambda[2] ~ dnorm(0, 0.04)
#>   for (e in 1:q_1) {
#>     for (l in 1:L_1) { r_1[l, e] ~ dnorm(0, 1 / pow(sd_1[e], 2)) }
#>     sd_1[e] ~ dt(0, 0.16, 3) T(0,)
#>   }
#> }
```

Each distribution parameter has a linear predictor with its link; random
effects use a non centered parameterization with Cholesky factored
correlations; smooth terms enter through their mixed model form;
censored and truncated observations use log CDF and log survival
functions. Families whose parameters have a joint constraint that no
link expresses (the generalized log series needs `theta * beta < 1`, for
example) reject draws outside the valid region:

``` r

gl <- stan_code(y ~ 1, family = "generalized_log_series", data = data.frame(y = c(1, 2, 3, 2)))$code
grep("negative_infinity", strsplit(gl, "\n")[[1]], value = TRUE)
#> [1] "  for (i in 1:N) if (!(theta[i] * beta[i] < 1)) target += negative_infinity();"
```

Composed families get helper functions named from their full settings,
so different truncations or hurdles in one model never clash:

``` r

mt <- stan_code(y ~ 1, family = mixture(truncate_dist("normal", 0, 1), truncate_dist("normal", 2, 3)),
                data = data.frame(y = c(0.3, 0.7, 2.3, 2.7)))$code
regmatches(mt, gregexpr("real rc_[A-Za-z0-9_]+", mt))[[1]]
#> [1] "real rc_trunc_normal_3ee7444f_lpdf"                           
#> [2] "real rc_trunc_normal_3ee74cd3_lpdf"                           
#> [3] "real rc_mix_truncate_normal_0_1_truncate_normal_34b9fc20_lpdf"
```

## Priors

Default priors are weakly informative: `normal(0, 5)` for intercepts and
coefficients, `student_t(3, 0, 2.5)` for standard deviations and
`lkj(2)` for correlation matrices.
[`prior()`](https://adamrsophie.github.io/regcompR/reference/prior.md)
sets priors by class (`b`, `Intercept`, `sd`, `cor`), by distribution
parameter (`par =`, using the parameter names of the family, such as
`mu.1` for a mixture) and by coefficient or random effect term
(`coef =`). The most specific matching prior wins.

``` r

pr <- c(prior("normal(0, 1)", class = "b"), prior("exponential(2)", class = "sd"),
        prior("normal(1, 0.5)", class = "b", coef = "x"))
pr
#> normal(0, 1)                 class b         par (all)  coef (all)
#> exponential(2)               class sd        par (all)  coef (all)
#> normal(1, 0.5)               class b         par (all)  coef x
code <- stan_code(y ~ x + (1 | g), family = "poisson", data = d, priors = pr)$code
grep("~ (normal|exponential)", strsplit(code, "\n")[[1]], value = TRUE)
#> [1] "  b_lambda[1] ~ normal(0, 5);"   "  b_lambda[2] ~ normal(1, 0.5);"
#> [3] "  sd_1 ~ exponential(2);"
m <- stan_code(y ~ 1, family = mixture("normal", "normal"), data = data.frame(y = c(-1, 0.2, 3, 4)),
               priors = prior("normal(-1, 1)", par = "mu.1", class = "Intercept"))$code
grep("b_mu_1", strsplit(m, "\n")[[1]], value = TRUE)
#> [1] "  vector[P_mu_1] b_mu_1;"                          
#> [2] "  vector[N] eta_mu_1 = X_mu_1 * b_mu_1 + off_mu_1;"
#> [3] "  b_mu_1[1] ~ normal(-1, 1);"
```

A prior that matches nothing in the model (a misspelled parameter or
coefficient, or a class the model does not have) is an error, so a prior
is never silently replaced by the default:

``` r

try(stan_code(y ~ x, family = "poisson", data = d, priors = prior("normal(0, 1)", par = "lambd")))
#> Error in check_priors_used(priors) : 
#>   These priors match nothing in the model: normal(0, 1) (class b, par lambd). Check par, class and coef.
```

For `ar1()` and `cs()` random effects, a `cor` prior applies to the
scalar correlation `rho`; `lkj` priors are for correlation matrices of
`us()` blocks.

``` r

lg <- data.frame(y = rnorm(60), t = factor(rep(1:3, 20)), id = factor(rep(1:20, each = 3)))
ar <- stan_code(y ~ 1 + ar1(0 + t | id), family = "normal", data = lg,
                priors = prior("uniform(0, 0.9)", class = "cor"))$code
grep("rho", strsplit(ar, "\n")[[1]], value = TRUE)
#> [1] "  real<lower=-1, upper=1> rho_1;"                                                      
#> [2] "  for (a in 1:q_1) for (c in 1:q_1) Rc_1[a, c] = a == c ? 1 + 1e-9 : rho_1^abs(a - c);"
#> [3] "  rho_1 ~ uniform(0, 0.9);"
```

A nonlinear formula (`y ~ nl(expr)`) computes its parameter row by row
from the nonlinear parameters. These have no natural scale, so each
needs an intercept prior; their starting values also start the chains.

``` r

dn <- data.frame(x = runif(60, 0, 5))
dn$y <- rpois(60, 10 * exp(-0.6 * dn$x))
pn <- c(prior("normal(2, 1)", par = "loga", class = "Intercept"),
        prior("normal(0.5, 0.5)", par = "b", class = "Intercept"))
nc <- stan_code(y ~ nl(loga - b * x), loga ~ 1, b ~ 1, family = "poisson", data = dn, priors = pn)$code
grep("lambda", strsplit(nc, "\n")[[1]], value = TRUE)
#> [1] "  vector[N] lambda;"                                              
#> [2] "  for (i in 1:N) lambda[i] = exp(loga[i] - b[i] * nlx_x[i]);"     
#> [3] "  for (i in 1:N) target += w[i] * poisson_lpmf(y[i] | lambda[i]);"
```

## Fitting

The examples below need a working Stan or JAGS installation.

``` r

fit <- regcomp(y ~ x + (1 | g), family = "poisson", data = d, engine = "bayes",
               priors = pr, control = regcomp_control(runner = "cmdstanr", chains = 4,
                                                      iter_warmup = 1000, iter_sampling = 1000))
summary(fit)                      # posterior means, credible intervals, R-hat, ESS
plot(fit, which = c("trace", "ppc"))
predict(fit, type = "mean", interval = "confidence")    # credible intervals
scoring_rule(fit, rule = "crps")  # posterior predictive CRPS
ll <- log_lik(fit)                # draws by observations, for loo::loo()

fit_rstan <- update(fit, control = regcomp_control(runner = "rstan"))
fit_jags <- regcomp(y ~ x + (1 | g), family = "poisson", data = d, engine = "bayes",
                    control = regcomp_control(runner = "jags"))
fit_ob <- regcomp(y ~ x + (1 | g), family = "poisson", data = d, engine = "bayes",
                  control = regcomp_control(runner = "openbugs"))
```

Mixture posteriors are multimodal: a chain started near the label
swapped mode stays there, which the R-hat warning reports. Informative
priors on the component parameters together with `init`, a list with one
named list of starting values per chain using the names of the generated
program, keep every chain in the intended mode:

``` r

ym <- c(rnorm(150, -2, 1), rnorm(150, 2, 0.7))
ini <- rep(list(list(b_mu_1 = array(-2), b_mu_2 = array(2))), 4)
fm <- regcomp(y ~ 1, family = mixture("normal", "normal"), data = data.frame(y = ym), engine = "bayes",
              priors = c(prior("normal(-2, 0.5)", par = "mu.1", class = "Intercept"),
                         prior("normal(2, 0.5)", par = "mu.2", class = "Intercept")),
              control = regcomp_control(runner = "cmdstanr", init = ini))
```

[`predict()`](https://rdrr.io/r/stats/predict.html) for a Bayesian fit
summarizes the posterior of the conditional quantity (for `"variance"`,
the posterior mean of the conditional variance);
[`scoring_rule()`](https://adamrsophie.github.io/regcompR/reference/scoring_rule.md)
uses the posterior predictive law.

## What each backend supports

Families with a Stan density, operators built from them, ordinal and
nominal models, all random effect covariance structures and multiple
membership terms are supported by the Stan backends, as are a single
inflation or hurdle point on a continuous law (zero inflated or zero
adjusted gamma, say). Other operators applied to a law that already has
point masses (truncation, censoring, copulas, several inflation points)
are fitted by the `"ml"` and `"tmb"` engines only. The BUGS backends
support families with BUGS densities and unstructured random effects.
[`family_capabilities()`](https://adamrsophie.github.io/regcompR/reference/family_capabilities.md)
reports what is available for a family.

``` r

family_capabilities("generalized_poisson")[c("stan", "tmb", "bugs")]
#> $stan
#> [1] "generated"
#> 
#> $tmb
#> [1] "unavailable"
#> 
#> $bugs
#> [1] "unavailable"
supports(hurdle("gamma"), backend = "stan")
#> [1] TRUE
```

## Marginal likelihoods and Bayes factors

[`bridge_regcomp()`](https://adamrsophie.github.io/regcompR/reference/bridge_regcomp.md)
estimates the log marginal likelihood of a Stan fit by bridge sampling
and
[`bayes_factor_regcomp()`](https://adamrsophie.github.io/regcompR/reference/bridge_regcomp.md)
compares two fits. The model is rewritten with all normalizing
constants, so the estimates match direct numerical integration (for
example, -51.941 against -51.940 for a Poisson regression with normal
priors). Use long chains.

``` r

bayes_factor_regcomp(fit_with_x, fit_without_x)
```
