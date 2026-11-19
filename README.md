# regcompR: Distributional Regression with a Verified Registry of Distributions <img src="man/figures/logo.png" align="right" height="139" alt="regcompR logo" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/adamrsophie/regcompR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/adamrsophie/regcompR/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/adamrsophie/regcompR/actions/workflows/pkgdown.yaml/badge.svg)](https://adamrsophie.github.io/regcompR/)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

***regcompR*** fits distributional regression models: every parameter of a
response distribution (location, scale, shape, zero inflation, cure
fraction, thresholds, copula dependence) can depend on covariates through
its own formula, with random effects in any parameter. One model
specification runs on several engines:

* maximum likelihood with a Laplace approximation for random effects and smooths;
* TMB automatic differentiation through RTMB;
* Bayesian inference with generated Stan programs (bayescompR, cmdstanr, rstan) or BUGS programs (JAGS, OpenBUGS).

## Highlights

* **278 families** in a machine readable registry: count, binary, bounded,
  ordinal, nominal, real line, positive, unit interval, simplex, circular,
  extreme value, multivariate and matrix valued laws. Each records support,
  parameter domains, links, moments with existence conditions, references
  and validation status.
* **Count models** from Hilbe's *Modeling Count Data* and the wider
  literature: NB1, NB2, NB-P, generalized Poisson, COM-Poisson, PIG,
  Poisson lognormal, Sichel, Delaporte, generalized Waring, Faddy's EPPM,
  zero truncated, zero inflated and hurdle variants, exact Poisson
  regression and quantile regression for counts.
* **Operators** build new laws: inflation, hurdles, truncation, censoring,
  transformation, discretization, mixtures, compounding, convolution,
  splicing, mixture cure models and copulas, plus the generalized G
  families (Kumaraswamy G, beta G, Marshall Olkin, transmuted, ...).
* **Hierarchical and longitudinal models**: nested and crossed random
  effects, random slopes, AR(1), compound symmetry, unstructured, Ornstein
  Uhlenbeck and spatial covariance structures, multiple membership, GEE
  with working correlations and cluster robust covariances.
* **Ordinal and categorical models**: cumulative, continuation ratio,
  adjacent category and stereotype models, location scale ordinal models,
  multinomial logit, bivariate binary and sample selection models.
* **Survival**: right, left and interval censoring, truncation, delayed
  entry, cure models and discrete time hazards with frailty.
* **Checking and comparison**: randomized quantile residuals, PIT,
  rootograms, effect and caterpillar plots, scoring rules, cross
  validation and identifiability checks.

## Installation

```r
# install.packages("pak")
pak::pak("adamrsophie/regcompR")
```

## Example

```r
library(regcompR)

d <- data.frame(x = runif(400))
d$y <- rnbinom(400, mu = exp(1 + d$x), size = exp(0.5 + d$x))

fit <- regcomp(y ~ x, phi ~ x, family = "nb2", data = d)
summary(fit)
predict(fit, data.frame(x = c(0.2, 0.8)), type = "quantile", probs = c(0.1, 0.9))
plot(fit, which = "rootogram")

# random effects in any parameter, faster through TMB
d$g <- factor(rep(1:40, each = 10))
fit_re <- regcomp(y ~ x + (1 | g), family = "nb2", data = d, engine = "tmb")

# the same model in Stan
fit_b <- regcomp(y ~ x + (1 | g), family = "nb2", data = d, engine = "bayes")
```

See the [articles](https://adamrsophie.github.io/regcompR/articles/) for
count, hierarchical, ordinal, survival, Bayesian and model checking
workflows.

## Validation

Families are checked for normalization, agreement of distribution
functions and quantiles, and moments against simulation. Fits are tested
against lme4, glmmTMB, MASS, ordinal, nnet, VGAM, geepack, sandwich,
survival, flexsurvcure, copula and scoringRules.

## License

GPL-3
