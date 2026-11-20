# Changelog

## regcompR 0.1.0

- First release.
- [`family_capabilities()`](https://adamrsophie.github.io/regcompR/reference/family_capabilities.md)
  (formerly
  [`capabilities()`](https://rdrr.io/r/base/capabilities.html), renamed
  so that attaching the package masks no base R function) labels family
  specific numerical methods `"algorithm"`, separately from closed forms
  (`"analytic"`) and the generic fallback (`"numeric"`).
- `scoring_rule(re = "marginal")` scores the predictive law of a new
  group, mixing over draws of its random effects. `cv_regcomp(group = )`
  now scores held out groups this way by default (previously with the
  effects set to zero; `re = "zero"` restores that) and clusters the
  standard error by group.
- The multivariate hypergeometric law takes the number of draws `size`.
- Finite difference scores and Hessians respect joint parameter
  constraints.
- Registry of 278 distribution families with densities, distribution
  functions, quantiles, random generation, moments and references.
- Distributional regression by maximum likelihood with a Laplace
  approximation for random effects and smooths, through TMB (RTMB), and
  by Bayesian inference in Stan, JAGS and OpenBUGS.
- Random effects in any parameter with nested, crossed, structured
  (AR(1), compound symmetry, unstructured, Ornstein Uhlenbeck, spatial
  kernels) and multiple membership terms.
- Ordinal, nominal, bivariate binary and sample selection models.
- Operators for inflation, hurdles, truncation, censoring, mixtures,
  compounding, convolution, splicing, cure models and copulas.
- Generalized estimating equations, cluster robust covariances, exact
  Poisson regression, quantile regression for counts and person period
  expansion for discrete time hazards.
- Residuals, scoring rules, cross validation, model comparison and
  plots.
