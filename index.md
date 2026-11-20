# regcompR: Distributional Regression in R

***regcompR*** fits distributional regression models: every parameter of
a response distribution (location, scale, shape, zero inflation, cure
fraction, ordinal thresholds, copula dependence) can depend on
covariates through its own formula, and random effects can enter any of
them. One model specification runs on several engines:

- maximum or penalized likelihood, with a Laplace approximation for
  random effects and smooth terms;
- automatic differentiation and the Laplace approximation of
  [TMB](https://github.com/kaskr/adcomp), through
  [RTMB](https://CRAN.R-project.org/package=RTMB);
- Bayesian inference with generated [Stan](https://mc-stan.org) programs
  (run by bayescompR, cmdstanr or rstan) or BUGS programs (run by JAGS
  or OpenBUGS).

Distributions live in a machine readable registry of 278 families. Each
entry records its support, parameter domains and links, moments with the
conditions under which they exist, the methods available on each engine,
references, and its validation status. Operators build new families from
registered ones, and every family, built or registered, works in every
part of the package.

## Installation

``` r

# Development version
# install.packages("remotes")
remotes::install_github("adamrsophie/regcompR")
```

regcompR needs R 4.3.0 or later. It contains no compiled code of its
own. Optional engines need their own packages: `RTMB` for
`engine = "tmb"`, a Stan runner (`cmdstanr`, `rstan` or `bayescompR`),
or JAGS with `R2jags` for `engine = "bayes"`.

## Quick start

``` r

library(regcompR)
set.seed(2026)

# Distribution objects
length(families())
#> [1] 278
d <- distribution("nb2", mu = 4, phi = 2)
prob(d, 0:3)
#> [1] 0.1111111 0.1481481 0.1481481 0.1316872
quantile(d, c(0.1, 0.5, 0.9))
#> [1] 0 3 9
c(mean = mean(d), variance = variance(d))
#>     mean variance
#>        4       12

# A negative binomial model with a random intercept in the mean and a
# regression on the dispersion, fitted through TMB
n <- 500
dat <- data.frame(x = runif(n), g = factor(rep(1:50, each = 10)))
dat$y <- rnbinom(n, mu = exp(1 + dat$x + rnorm(50, 0, 0.4)[dat$g]),
                 size = exp(0.3 + 1.2 * dat$x))
fit <- regcomp(y ~ x + (1 | g), phi ~ x, family = "nb2", data = dat,
               engine = "tmb")
round(coef(fit), 3)
#>  mu:(Intercept)            mu:x phi:(Intercept)           phi:x
#>           0.940           1.104           0.393           0.890
round(exp(fit$theta), 3)       # random intercept standard deviation
#> mu:1 | g:log_sd((Intercept))
#>                        0.261

# Predictive quantiles and event probabilities for a new group
nd <- data.frame(x = c(0.2, 0.8))
predict(fit, nd, type = "quantile", probs = c(0.1, 0.9), re = "zero")
#>      q0.1 q0.9
#> [1,]    0    7
#> [2,]    1   12
predict(fit, nd, type = "probability", at = 0, re = "zero")
#>      probability(0)
#> [1,]     0.16123603
#> [2,]     0.03444051

# Operators give new families, used like any other
zinb <- regcomp(y ~ x, zi ~ 1, family = inflate("nb2"), data = dat)
AIC(fit, zinb)
#>      df      AIC
#> fit   5 2557.586
#> zinb  4 2571.487

# Diagnostics
plot(fit, which = c("rootogram", "qq", "ranef"))

# The same model in Stan
fit_b <- regcomp(y ~ x + (1 | g), phi ~ x, family = "nb2", data = dat,
                 engine = "bayes")
```

## What is covered

| Area | Highlights |
|----|----|
| Registry | 278 families: count, binary, finite discrete, bounded, ordinal, nominal, real line, positive, unit interval, simplex, circular and directional, extreme value, multivariate and matrix valued; density, distribution function, quantile, random generation, moments, score and Hessian for each |
| Count data | Poisson, NB1, NB2, NB-P, generalized Poisson, Conway Maxwell Poisson, Poisson inverse Gaussian, Poisson lognormal, Sichel, Delaporte, generalized Waring, generalized negative binomial, Faddy’s extended Poisson process models, multivariate counts; zero truncated, zero inflated and hurdle variants |
| Operators | [`inflate()`](https://adamrsophie.github.io/regcompR/reference/inflate.md), [`hurdle()`](https://adamrsophie.github.io/regcompR/reference/hurdle.md), [`truncate_dist()`](https://adamrsophie.github.io/regcompR/reference/truncate_dist.md), [`censor()`](https://adamrsophie.github.io/regcompR/reference/censor.md), [`transform_dist()`](https://adamrsophie.github.io/regcompR/reference/transform_dist.md), [`discretize()`](https://adamrsophie.github.io/regcompR/reference/discretize.md), [`round_dist()`](https://adamrsophie.github.io/regcompR/reference/discretize.md), [`mixture()`](https://adamrsophie.github.io/regcompR/reference/mixture.md), [`compound()`](https://adamrsophie.github.io/regcompR/reference/compound.md), [`convolve_dist()`](https://adamrsophie.github.io/regcompR/reference/convolve_dist.md), [`splice()`](https://adamrsophie.github.io/regcompR/reference/splice.md), [`cure()`](https://adamrsophie.github.io/regcompR/reference/cure.md), [`copula_dist()`](https://adamrsophie.github.io/regcompR/reference/copula_dist.md), [`wrap()`](https://adamrsophie.github.io/regcompR/reference/wrap.md), [`weighted()`](https://adamrsophie.github.io/regcompR/reference/weighted.md), skew symmetric constructions, and generators such as Kumaraswamy G, beta G, Marshall Olkin, transmuted and exponentiated families |
| Regression | A formula for every parameter, offsets, weights, custom links, fixed parameters, smooth terms through mgcv bases, penalties (ridge, lasso, elastic net, adaptive lasso, group lasso, SCAD, MCP) and bias reduction (Firth, FLIC, FLAC, log F) |
| Hierarchical models | Random intercepts and slopes in any parameter, nested and crossed factors, correlated or independent effects, multiple membership terms |
| Longitudinal models | AR(1), compound symmetry, unstructured, diagonal, Ornstein Uhlenbeck and spatial (exponential, Gaussian, Matern 3/2) covariance structures; GEE with independence, exchangeable, AR(1) and unstructured working correlations; cluster robust covariances |
| Ordinal and categorical | Cumulative (logit, probit, cloglog, loglog, cauchit), continuation ratio, adjacent category and stereotype models, location scale ordinal models, multinomial logit, bivariate binary odds ratio and sample selection models |
| Survival | Right, left and interval censoring, truncation and delayed entry for any family; parametric survival models; mixture cure models; discrete time hazards with [`person_period()`](https://adamrsophie.github.io/regcompR/reference/person_period.md) and frailty through random effects |
| Special methods | Exact conditional Poisson regression, quantile regression for counts, custom families from a log density |
| Inference and prediction | Wald and posterior summaries, sandwich covariances, predictions of the mean, any parameter, quantiles, probabilities and survival with confidence or credible intervals, conditional or marginal over random effects |
| Checking and comparison | Randomized quantile residuals, PIT histograms, rootograms, effect and caterpillar plots, posterior predictive checks, log, CRPS, Dawid Sebastiani, Brier and ranked probability scores, K fold cross validation, identifiability checks |
| Bayesian backends | Generated Stan code for every family with a Stan density, all random effect structures and multiple membership; BUGS code for JAGS and OpenBUGS; priors by class, parameter and coefficient |

[`family_capabilities()`](https://adamrsophie.github.io/regcompR/reference/family_capabilities.md)
reports what each engine supports for a family, and
[`registry_table()`](https://adamrsophie.github.io/regcompR/reference/registry_table.md)
returns the whole registry as a data frame.

## Validation

- Families are checked for normalization, for agreement between its
  density, distribution function and quantile function, and for moments
  against simulation; families with a reference implementation are
  compared with it (`stats`, `extraDistr`, `VGAM`, `actuar`, `evd`,
  `circular` and others).
- Model fits are tested against established packages on the same data:
  linear and generalized linear mixed models against `lme4` and
  `glmmTMB`, covariance structures against `glmmTMB`, ordinal models
  against `MASS`, `ordinal` and `VGAM`, multinomial models against
  `nnet`, GEE against `geepack`, sandwich covariances against
  `sandwich`, parametric survival against `survival`, cure models
  against `flexsurvcure`, copulas against `copula`, and scoring rules
  against `scoringRules`.
- The TMB engine reproduces the default engine’s estimates and standard
  errors; Bayesian fits are checked against the likelihood fits.
- `R CMD check --as-cran` passes with no errors or warnings.

## Documentation

The [package website](https://adamrsophie.github.io/regcompR/) has the
full reference and these articles:

- [Getting
  started](https://adamrsophie.github.io/regcompR/articles/regcompR.html)
- [Distributions and
  operators](https://adamrsophie.github.io/regcompR/articles/distributions.html)
- [Count data
  models](https://adamrsophie.github.io/regcompR/articles/count-models.html)
- [Hierarchical and longitudinal
  models](https://adamrsophie.github.io/regcompR/articles/hierarchical-longitudinal.html)
- [Ordinal and categorical
  responses](https://adamrsophie.github.io/regcompR/articles/ordinal-categorical.html)
- [Censored data and survival
  models](https://adamrsophie.github.io/regcompR/articles/survival.html)
- [Model checking and
  comparison](https://adamrsophie.github.io/regcompR/articles/model-checking.html)
- [Bayesian distributional
  regression](https://adamrsophie.github.io/regcompR/articles/bayesian.html)

## Citing regcompR

``` r

citation("regcompR")
```

The repository ships machine readable citation metadata:

- [`CITATION.cff`](https://adamrsophie.github.io/regcompR/CITATION.cff):
  Citation File Format (GitHub “Cite this repository”, Zenodo).
- [`.zenodo.json`](https://adamrsophie.github.io/regcompR/.zenodo.json):
  Zenodo deposit metadata.
- `inst/CITATION`: the R level citation (`citation("regcompR")`).

Once the first GitHub release is archived on Zenodo, its concept DOI
will resolve to the latest release.

## Contributing

Bug reports, feature requests, and pull requests are welcome; see
[`CONTRIBUTING.md`](https://adamrsophie.github.io/regcompR/CONTRIBUTING.md).
A wrong answer is the most useful bug report: include a reproducible
example and, if you can, what a reference package returns for the same
model. Please follow the [Code of
Conduct](https://adamrsophie.github.io/regcompR/CODE_OF_CONDUCT.md), and
report security issues privately as described in
[`SECURITY.md`](https://adamrsophie.github.io/regcompR/SECURITY.md).

## Author

- Adam R Sophie (author, maintainer, copyright holder),
  <adamrsophie@gmail.com>

## Acknowledgments

regcompR builds on the distributional regression framework of Rigby and
Stasinopoulos (2005), on TMB (Kristensen et al. 2016) through RTMB, and
on Stan, JAGS and OpenBUGS for Bayesian inference. The packages used as
references in its tests are gratefully acknowledged.

## License

GPL-3. See
[`LICENSE.md`](https://adamrsophie.github.io/regcompR/LICENSE.md) or
<https://www.gnu.org/licenses/gpl-3.0> for the full license text.
