# Fit a distributional regression model

Every regression parameter of the response distribution gets its own
linear predictor on its link scale. The first formula gives the response
and the predictor of the main parameter: the parameter equal to the mean
if there is one, else the location, scale or rate parameter, else the
first regression parameter (for example `scale` of the Weibull law, so
that survival models are accelerated failure time models); further
formulas such as `sigma ~ x` or `shape ~ x` give predictors for other
parameters, which otherwise are constant. Predictors may contain
offsets, factors and interactions, random effects in `lme4` syntax
(`(1 | g)`, `(1 + x | g)`, uncorrelated `(1 + x || g)`, nested
`(1 | a/b)`, crossed `(1 | a) + (1 | b)`), structured random effects for
longitudinal and spatial data written as glmmTMB style wrappers, and
smooth terms `s()` and `t2()` from 'mgcv', which are fitted in their
mixed model representation (as in
[`mgcv::gamm()`](https://rdrr.io/pkg/mgcv/man/gamm.html)). Tensor
products `te()` and `ti()` have no such representation; use `t2()`
instead.

## Usage

``` r
regcomp(
  formula,
  ...,
  family,
  data,
  engine = c("ml", "tmb", "bayes"),
  weights = NULL,
  links = NULL,
  fixed = NULL,
  penalty = NULL,
  priors = NULL,
  control = regcomp_control(),
  parameterization = NULL,
  subset = NULL
)
```

## Arguments

- formula:

  Two sided formula `response ~ predictors` for the first regression
  parameter. The response may be wrapped in
  [`cens()`](https://adamrsophie.github.io/regcompR/reference/cens.md)
  or
  [`surv()`](https://adamrsophie.github.io/regcompR/reference/cens.md).
  A nonlinear predictor is written `response ~ nl(expr)`, as in
  `y ~ nl(a * exp(-b * x))`: the names of `expr` that have their own
  formula in `...` (`a ~ 1 + (1 | g)`, `b ~ z`) are nonlinear
  parameters, each with a linear predictor that may hold covariates and
  random effects; the other names are covariates. `expr` gives the first
  parameter on its link scale, as in
  [`nls()`](https://rdrr.io/r/stats/nls.html), `nlme()` and `nlmer()`.
  Nonlinear parameters need starting values
  (`control = regcomp_control(start = list(a = 10, b = 0.5))`, also the
  initial values of the Stan chains) and are fitted by `engine = "ml"`
  or by `engine = "bayes"` with a Stan runner.

- ...:

  Further formulas for other parameters, such as `sigma ~ x`, and data
  arguments of the family (for example `trials = n` for the binomial),
  evaluated in `data`.

- family:

  A family identifier, family object or a family built by an operator
  such as
  [`inflate()`](https://adamrsophie.github.io/regcompR/reference/inflate.md).

- data:

  A data frame.

- engine:

  `"ml"`, `"tmb"` or `"bayes"`.

- weights:

  Optional prior weights multiplying each log likelihood contribution.

- links:

  Optional named list or character vector of link names per parameter.

- fixed:

  Optional named list of parameters held fixed at given values (not
  estimated).

- penalty:

  Optional penalty from
  [`penalty()`](https://adamrsophie.github.io/regcompR/reference/penalty.md)
  (engine `"ml"`).

- priors:

  Optional priors from
  [`prior()`](https://adamrsophie.github.io/regcompR/reference/prior.md)
  (engine `"bayes"`).

- control:

  List of control settings, see
  [`regcomp_control()`](https://adamrsophie.github.io/regcompR/reference/regcomp_control.md).

- parameterization:

  Optional parameterization name passed to
  [`get_family()`](https://adamrsophie.github.io/regcompR/reference/get_family.md).

- subset:

  Optional logical or index vector selecting rows of `data`.

## Value

An object of class `regcomp_fit`.

## Random effect structures

A bar term may be wrapped in a covariance structure:

- `us(x | g)`: unstructured (the default for `|`); `diag(x | g)` and
  `homdiag(x | g)`: independent effects with separate or common variance
  (the default for `||`).

- `cs(0 + time | g)` and `homcs(0 + time | g)`: compound symmetry with
  heterogeneous or common variances.

- `ar1(0 + time | g)`: first order autoregressive correlation
  `rho^|i - j|` across the columns (time must be a factor with equally
  spaced, ordered levels), with a common variance.

- `ou(t | g)` or `exp(t | g)`: continuous time autoregressive
  (Ornstein-Uhlenbeck) correlation `exp(-|t_i - t_j| / range)` for
  unequally spaced numeric times; `gau(x + y | g)` and
  `matern32(x + y | g)` give Gaussian and Matern (smoothness 3/2)
  correlation in the coordinates on the left of the bar, for spatial
  effects. One random effect is created per distinct coordinate.

- `relmat(id | 1, K)`: effects over the levels of `id` with a known
  relatedness (kinship, genomic or phylogenetic) matrix `K` whose rows
  are the levels (by row name, else in level order).

- `car(region | 1, W)`, `leroux(region | 1, W)` and
  `icar(region | 1, W)`: areal effects from a symmetric adjacency matrix
  `W`: the proper CAR model (precision `D - rho W`), the Leroux model
  (precision `(1 - rho) I + rho (D - W)`) and the intrinsic CAR (sum to
  zero, scaled to unit geometric mean variance). `icar()` plus
  `(1 | region)` is the BYM model.

- `spacetime(region + time | 1, W)`: a separable space time effect, the
  Kronecker product of a Leroux spatial and an AR(1) temporal
  correlation (time a factor of equally spaced periods).

- `spde(x + y | 1, mesh)`: a Matern field (smoothness 1) by the SPDE
  approach of Lindgren, Rue and Lindstrom (2011) on an fmesher
  triangulation, with parameters `log_sd` and `log_range`. The precision
  is sparse but the current engines use it densely, so meshes of up to a
  few hundred nodes are practical.

A group other than `1` on the right of the bar gives independent
replicates of the field (one per year, say).

Structured random effects are available in every parameter and for every
family, so a latent autoregressive process can enter a count, survival
or ordinal model. The TMB engine does not support `car()`, `leroux()`,
`spacetime()` and `spde()`; Stan does not support `spacetime()` and
`spde()`.

## Engines

- `"ml"`: maximum likelihood in R. Without random effects or smooths the
  log likelihood is maximized by
  [`nlminb()`](https://rdrr.io/r/stats/nlminb.html) with an analytic
  chain rule gradient. With random effects or smooths the marginal
  likelihood is approximated by the Laplace method (joint over all
  random effects, so crossed and nested factors and random effects in
  several parameters are allowed); the fitted object records
  `approximation = "laplace"`. Penalties (see
  [`penalty()`](https://adamrsophie.github.io/regcompR/reference/penalty.md))
  give penalized likelihood estimates.

- `"tmb"`: the same model, random effects and covariance structures
  through automatic differentiation and the Laplace approximation of
  'TMB' (via 'RTMB'), for the normal, Poisson, binomial, Bernoulli, NB1,
  NB2, Student t, gamma, lognormal, Weibull, beta and exponential
  families (see
  [`family_capabilities()`](https://adamrsophie.github.io/regcompR/reference/family_capabilities.md)).
  It gives the same estimates as `"ml"` and is much faster for large
  multilevel models; the taped density is checked against the registry
  density before fitting.

- `"bayes"`: Bayesian inference. A 'Stan' program is generated from the
  same specification (see
  [`stan_code()`](https://adamrsophie.github.io/regcompR/reference/stan_code.md))
  and run by 'bayescompR' (default, no compiler needed), 'cmdstanr' or
  'rstan'. Priors are set with
  [`prior()`](https://adamrsophie.github.io/regcompR/reference/prior.md).

## Examples

``` r
set.seed(1)
d <- data.frame(x = rnorm(300), z = rnorm(300))
d$y <- rnbinom(300, mu = exp(1 + 0.5 * d$x), size = exp(0.5 + 0.4 * d$z))
fit <- regcomp(y ~ x, phi ~ z, family = "nb2", data = d)
summary(fit)
#> <regcomp_fit> Negative binomial (NB2) regression, engine 'ml' (frequentist)
#>   mu (log link): ~x
#>   phi (log link): ~z
#>   log likelihood -637.197 on 4 parameters, n = 300
#> 
#> Wald estimates with 95% confidence intervals:
#>                 estimate std_error  lower  upper       z p_value
#> mu:(Intercept)    0.9818    0.0547 0.8745 1.0890 17.9431  0.0000
#> mu:x              0.6232    0.0593 0.5070 0.7395 10.5090  0.0000
#> phi:(Intercept)   0.7865    0.1651 0.4628 1.1101  4.7631  0.0000
#> phi:z             0.4709    0.1645 0.1486 0.7933  2.8635  0.0042
#> max |gradient| 2.84e-05; Hessian positive definite 
head(predict(fit, type = "mean"))
#>        1        2        3        4        5        6 
#> 1.806405 2.992792 1.585622 7.213584 3.277602 1.600675 
predict(fit, newdata = data.frame(x = 0, z = 0), type = "quantile",
        probs = c(0.1, 0.9))
#> [1] 0 6
```
