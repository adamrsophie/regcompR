# regcompR 0.1.0

* First release.
* `pql_regcomp()` (penalized quasi-likelihood, equal to MASS::glmmPQL),
  `re_nongauss()` (gamma, t, logistic and Laplace random intercepts by
  adaptive quadrature), `spde()` Matern fields on fmesher meshes, and
  `bridge_regcomp()` / `bayes_factor_regcomp()` (marginal likelihoods of
  Stan fits by bridge sampling, with all normalizing constants). Offsets
  are allowed in log scale parameters (variance weights).
* `panel_within(y ~ x | z)` gives the within instrumental variable
  estimator (equal to plm); `vcov(type = "HAC")` (Newey-West, equal to
  sandwich); `anova(boundary = TRUE)` for variance components on their
  boundary; `shrink_regcomp()` (predictive shrinkage, any family including
  beta regression); `mediate_regcomp()` (causal mediation, as mediation).
* `summary_param()` reparameterizes a family by its mean, mode, geometric
  or harmonic mean; `uplift()` and `qini()` (equal to tools4uplift) for
  uplift models; `skovgaard_test()` gives Skovgaard's modified directed
  likelihood ratio r* (equal to likelihoodAsy::rstar()).
* Inference: `vcov(type = "HC1")` to `"HC3"` (equal to sandwich::vcovHC),
  `boot_regcomp()` (parametric and case bootstrap, clustered),
  `power_sim()` (power by simulation), `vuong_test()` and `best_subset()`;
  `step_regcomp()` and `best_subset()` also use AICc and Hannan-Quinn.
* Random effect structures `relmat()` (known relatedness matrix),
  `car()`, `leroux()` and `icar()` (areal models; `icar()` plus an iid term
  is BYM) and `spacetime()` (separable Leroux by AR(1)); exact normal
  marginal likelihoods checked against direct computation.
* Econometrics: `iv_regcomp()` (two stage least squares with weak
  instrument, Wu-Hausman and Sargan tests; equal to ivreg),
  `panel_within()` (fixed effects within estimator; equal to plm) and
  `hausman_test()` (equal to plm::phtest).
* Families: `huber` (its ML estimate with fixed scale is the Huber M
  estimate of MASS::rlm), `stochastic_frontier` (normal minus half normal),
  `multinomial_probit` (independent errors) and `partial_odds()` (partial
  proportional odds; equal to ordinal::clm with nominal effects).
* Nonlinear formulas: `y ~ nl(expr)` with a formula for each nonlinear
  parameter (covariates and random effects allowed) for any family, fitted
  by `engine = "ml"` (equal to `nls()` and to `lme4::nlmer()`) or by Stan.
* Reporting: `regcomp_table()` builds regression tables of one or several
  fits (rows grouped by distribution parameter, random effect standard
  deviations and correlations, fit statistics, optional standardized or
  exponentiated coefficients and sandwich or quasi standard errors) that
  print in the console and render in R Markdown and Quarto (HTML, PDF,
  Word). `tidy()`, `glance()` and `augment()` methods work with broom and
  modelsummary. `export_table()` writes `.csv`, `.html`, `.tex`, `.md` and
  `.docx`; `export_plot()` saves any diagnostic plot; `coef_plot()` draws
  estimates of several fits.
* Inference: `anova()` likelihood ratio tests of nested fits,
  `wald_test()` for linear hypotheses, `simple_slopes()` and
  `johnson_neyman()` for moderation, interaction plots through
  `plot(which = "effects", by = )`, `vif()` generalized variance inflation
  factors, and `step_regcomp()` stepwise selection by AIC or BIC.
* Effects: `avg_effects()` (average marginal effects and contrasts, with
  delta method or posterior uncertainty), `marginal_means()` (average
  predictions at covariate values), `r2()` (McFadden, Cox and Snell,
  Nagelkerke, and Nakagawa marginal and conditional R squared) and `icc()`.
* `vcov(type = "quasi")` and `summary(vcov = )` give quasi-likelihood
  (Pearson dispersion) standard errors, matching `quasibinomial` and
  `quasipoisson`.
* Dirichlet regression (matrix response; one log linked concentration per
  component) and the shifted lognormal family.
* `person_period(cuts = "events")` builds Cox risk sets: a Poisson fit
  reproduces `coxph()` with Breslow ties exactly, and random effects give
  Cox frailty models.
* `family_capabilities()` (formerly `capabilities()`, renamed so that
  attaching the package masks no base R function) labels family specific
  numerical methods `"algorithm"`, separately from closed forms
  (`"analytic"`) and the generic fallback (`"numeric"`).
* `scoring_rule(re = "marginal")` scores the predictive law of a new group,
  mixing over draws of its random effects. `cv_regcomp(group = )` now
  scores held out groups this way by default (previously with the effects
  set to zero; `re = "zero"` restores that) and clusters the standard
  error by group.
* The multivariate hypergeometric law takes the number of draws `size`.
* Finite difference scores and Hessians respect joint parameter
  constraints.
* Truncation, censoring, rounding, monotone transformations, splicing,
  CDF generators, convolutions (including two laws that both have point
  masses) and copulas now accept continuous laws with point masses;
  copula margins may also differ in type (a count and a continuous
  amount, say).
* `regcomp_control(init = )` sets initial values for the Stan runners.
* Registry of 278 distribution families with densities, distribution
  functions, quantiles, random generation, moments and references.
* Distributional regression by maximum likelihood with a Laplace
  approximation for random effects and smooths, through TMB (RTMB), and
  by Bayesian inference in Stan, JAGS and OpenBUGS.
* Random effects in any parameter with nested, crossed, structured
  (AR(1), compound symmetry, unstructured, Ornstein Uhlenbeck, spatial
  kernels) and multiple membership terms.
* Ordinal, nominal, bivariate binary and sample selection models.
* Operators for inflation, hurdles, truncation, censoring, mixtures,
  compounding, convolution, splicing, cure models and copulas.
* Generalized estimating equations, cluster robust covariances, exact
  Poisson regression, quantile regression for counts and person period
  expansion for discrete time hazards.
* Residuals, scoring rules, cross validation, model comparison and plots.
