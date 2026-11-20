# Package index

## Fitting models

- [`regcomp()`](https://adamrsophie.github.io/regcompR/reference/regcomp.md)
  : Fit a distributional regression model
- [`regcomp_control()`](https://adamrsophie.github.io/regcompR/reference/regcomp_control.md)
  : Control settings for model fitting
- [`penalty()`](https://adamrsophie.github.io/regcompR/reference/penalty.md)
  : Penalties for regularized and bias reduced estimation
- [`prior()`](https://adamrsophie.github.io/regcompR/reference/prior.md)
  : Priors for Bayesian fits
- [`gee()`](https://adamrsophie.github.io/regcompR/reference/gee.md) :
  Generalized estimating equations
- [`exact_poisson()`](https://adamrsophie.github.io/regcompR/reference/exact_poisson.md)
  : Exact conditional Poisson regression
- [`quantile_count()`](https://adamrsophie.github.io/regcompR/reference/quantile_count.md)
  : Quantile regression for counts
- [`person_period()`](https://adamrsophie.github.io/regcompR/reference/person_period.md)
  : Person period data for discrete time hazard models

## Distributions

- [`distribution()`](https://adamrsophie.github.io/regcompR/reference/distribution.md)
  : Create a distribution object
- [`log_prob()`](https://adamrsophie.github.io/regcompR/reference/probability-interface.md)
  [`prob()`](https://adamrsophie.github.io/regcompR/reference/probability-interface.md)
  [`cdf()`](https://adamrsophie.github.io/regcompR/reference/probability-interface.md)
  [`survival()`](https://adamrsophie.github.io/regcompR/reference/probability-interface.md)
  [`rng()`](https://adamrsophie.github.io/regcompR/reference/probability-interface.md)
  : Probability interface
- [`mean(`*`<regcomp_dist>`*`)`](https://adamrsophie.github.io/regcompR/reference/moments.md)
  [`variance()`](https://adamrsophie.github.io/regcompR/reference/moments.md)
  [`skewness()`](https://adamrsophie.github.io/regcompR/reference/moments.md)
  [`kurtosis()`](https://adamrsophie.github.io/regcompR/reference/moments.md)
  : Moments of a distribution object
- [`quantile(`*`<regcomp_dist>`*`)`](https://adamrsophie.github.io/regcompR/reference/quantile.regcomp_dist.md)
  : Quantiles of a distribution object
- [`support()`](https://adamrsophie.github.io/regcompR/reference/support.md)
  : Support of a distribution object
- [`score()`](https://adamrsophie.github.io/regcompR/reference/score.md)
  [`hessian()`](https://adamrsophie.github.io/regcompR/reference/score.md)
  : Score and Hessian of the log likelihood
- [`plot(`*`<regcomp_dist>`*`)`](https://adamrsophie.github.io/regcompR/reference/plot.regcomp_dist.md)
  : Plot a distribution
- [`get_family()`](https://adamrsophie.github.io/regcompR/reference/get_family.md)
  [`families()`](https://adamrsophie.github.io/regcompR/reference/get_family.md)
  : Look up a registered family
- [`registry_table()`](https://adamrsophie.github.io/regcompR/reference/registry_table.md)
  : Registry table
- [`family_capabilities()`](https://adamrsophie.github.io/regcompR/reference/family_capabilities.md)
  [`supports()`](https://adamrsophie.github.io/regcompR/reference/family_capabilities.md)
  : Capability queries
- [`validate_parameters()`](https://adamrsophie.github.io/regcompR/reference/validate_parameters.md)
  [`validate_response()`](https://adamrsophie.github.io/regcompR/reference/validate_parameters.md)
  : Validate parameters and responses
- [`ordinal-families`](https://adamrsophie.github.io/regcompR/reference/ordinal-families.md)
  : Cumulative link ordinal families

## Building families

- [`par_spec()`](https://adamrsophie.github.io/regcompR/reference/par_spec.md)
  : Parameter specification
- [`custom_family()`](https://adamrsophie.github.io/regcompR/reference/custom_family.md)
  : User defined families
- [`link()`](https://adamrsophie.github.io/regcompR/reference/link.md)
  [`links()`](https://adamrsophie.github.io/regcompR/reference/link.md)
  : Link functions
- [`quantile_param()`](https://adamrsophie.github.io/regcompR/reference/quantile_param.md)
  : Quantile parameterized families

## Operators

- [`inflate()`](https://adamrsophie.github.io/regcompR/reference/inflate.md)
  : Inflated distributions
- [`hurdle()`](https://adamrsophie.github.io/regcompR/reference/hurdle.md)
  : Hurdle distributions
- [`truncate_dist()`](https://adamrsophie.github.io/regcompR/reference/truncate_dist.md)
  : Truncated distributions
- [`censor()`](https://adamrsophie.github.io/regcompR/reference/censor.md)
  : Censored distributions
- [`cens()`](https://adamrsophie.github.io/regcompR/reference/cens.md)
  [`surv()`](https://adamrsophie.github.io/regcompR/reference/cens.md) :
  Censored, truncated and interval observations
- [`transform_dist()`](https://adamrsophie.github.io/regcompR/reference/transform_dist.md)
  : Monotone transformations
- [`discretize()`](https://adamrsophie.github.io/regcompR/reference/discretize.md)
  [`round_dist()`](https://adamrsophie.github.io/regcompR/reference/discretize.md)
  : Discretized and rounded distributions
- [`mixture()`](https://adamrsophie.github.io/regcompR/reference/mixture.md)
  : Finite mixtures
- [`compound()`](https://adamrsophie.github.io/regcompR/reference/compound.md)
  : Compound (mixed) distributions
- [`convolve_dist()`](https://adamrsophie.github.io/regcompR/reference/convolve_dist.md)
  : Convolution of two independent laws
- [`splice()`](https://adamrsophie.github.io/regcompR/reference/splice.md)
  : Spliced (composite) distributions
- [`cure()`](https://adamrsophie.github.io/regcompR/reference/cure.md) :
  Mixture cure models
- [`copula_dist()`](https://adamrsophie.github.io/regcompR/reference/copula_dist.md)
  : Bivariate copula models
- [`wrap()`](https://adamrsophie.github.io/regcompR/reference/wrap.md) :
  Wrapped distributions
- [`weighted()`](https://adamrsophie.github.io/regcompR/reference/weighted.md)
  [`contaminated()`](https://adamrsophie.github.io/regcompR/reference/weighted.md)
  : Weighted and contaminated distributions
- [`skew_symmetric()`](https://adamrsophie.github.io/regcompR/reference/skew_symmetric.md)
  : Skew-symmetric distributions
- [`sine_skew()`](https://adamrsophie.github.io/regcompR/reference/sine_skew.md)
  : Sine-skewed circular distributions
- [`exponentiated()`](https://adamrsophie.github.io/regcompR/reference/exponentiated.md)
  [`kumaraswamy_g()`](https://adamrsophie.github.io/regcompR/reference/exponentiated.md)
  [`beta_g()`](https://adamrsophie.github.io/regcompR/reference/exponentiated.md)
  [`mcdonald_g()`](https://adamrsophie.github.io/regcompR/reference/exponentiated.md)
  [`exp_generalized_g()`](https://adamrsophie.github.io/regcompR/reference/exponentiated.md)
  [`marshall_olkin()`](https://adamrsophie.github.io/regcompR/reference/exponentiated.md)
  [`transmuted()`](https://adamrsophie.github.io/regcompR/reference/exponentiated.md)
  [`record_transmuted()`](https://adamrsophie.github.io/regcompR/reference/exponentiated.md)
  [`odd_log_logistic()`](https://adamrsophie.github.io/regcompR/reference/exponentiated.md)
  [`odd_burr()`](https://adamrsophie.github.io/regcompR/reference/exponentiated.md)
  [`weibull_g()`](https://adamrsophie.github.io/regcompR/reference/exponentiated.md)
  [`gamma_g()`](https://adamrsophie.github.io/regcompR/reference/exponentiated.md)
  [`topp_leone_g()`](https://adamrsophie.github.io/regcompR/reference/exponentiated.md)
  [`lindley_g()`](https://adamrsophie.github.io/regcompR/reference/exponentiated.md)
  [`harris_g()`](https://adamrsophie.github.io/regcompR/reference/exponentiated.md)
  [`sine_g()`](https://adamrsophie.github.io/regcompR/reference/exponentiated.md)
  [`dus()`](https://adamrsophie.github.io/regcompR/reference/exponentiated.md)
  [`km()`](https://adamrsophie.github.io/regcompR/reference/exponentiated.md)
  [`power_series_max()`](https://adamrsophie.github.io/regcompR/reference/exponentiated.md)
  [`power_series_min()`](https://adamrsophie.github.io/regcompR/reference/exponentiated.md)
  : Generated distribution families

## Tables and export

- [`regcomp_table()`](https://adamrsophie.github.io/regcompR/reference/regcomp_table.md)
  : Regression tables
- [`tidy(`*`<regcomp_fit>`*`)`](https://adamrsophie.github.io/regcompR/reference/tidy.regcomp_fit.md)
  [`glance(`*`<regcomp_fit>`*`)`](https://adamrsophie.github.io/regcompR/reference/tidy.regcomp_fit.md)
  [`augment(`*`<regcomp_fit>`*`)`](https://adamrsophie.github.io/regcompR/reference/tidy.regcomp_fit.md)
  : Tidy summaries of a fitted model
- [`export_table()`](https://adamrsophie.github.io/regcompR/reference/export_table.md)
  [`export_plot()`](https://adamrsophie.github.io/regcompR/reference/export_table.md)
  : Export tables and plots
- [`coef_plot()`](https://adamrsophie.github.io/regcompR/reference/coef_plot.md)
  : Coefficient plot

## Tests and effects

- [`avg_effects()`](https://adamrsophie.github.io/regcompR/reference/avg_effects.md)
  : Average marginal effects and contrasts
- [`marginal_means()`](https://adamrsophie.github.io/regcompR/reference/marginal_means.md)
  : Marginal means (average predictions)
- [`r2()`](https://adamrsophie.github.io/regcompR/reference/r2.md)
  [`icc()`](https://adamrsophie.github.io/regcompR/reference/r2.md) : R
  squared and intraclass correlation
- [`anova(`*`<regcomp_fit>`*`)`](https://adamrsophie.github.io/regcompR/reference/anova.regcomp_fit.md)
  : Likelihood ratio tests of nested models
- [`wald_test()`](https://adamrsophie.github.io/regcompR/reference/wald_test.md)
  : Wald tests of linear hypotheses
- [`simple_slopes()`](https://adamrsophie.github.io/regcompR/reference/simple_slopes.md)
  : Simple slopes under moderation
- [`johnson_neyman()`](https://adamrsophie.github.io/regcompR/reference/johnson_neyman.md)
  : Johnson-Neyman regions of significance
- [`vif()`](https://adamrsophie.github.io/regcompR/reference/vif.md) :
  Variance inflation factors
- [`step_regcomp()`](https://adamrsophie.github.io/regcompR/reference/step_regcomp.md)
  : Stepwise selection of terms

## Results and prediction

- [`summary(`*`<regcomp_fit>`*`)`](https://adamrsophie.github.io/regcompR/reference/summary.regcomp_fit.md)
  : Summarize a fitted model
- [`predict(`*`<regcomp_fit>`*`)`](https://adamrsophie.github.io/regcompR/reference/predict.regcomp_fit.md)
  : Predictions from a fitted model
- [`fitted_distribution()`](https://adamrsophie.github.io/regcompR/reference/fitted_distribution.md)
  : Fitted distributions
- [`simulate(`*`<regcomp_fit>`*`)`](https://adamrsophie.github.io/regcompR/reference/simulate.regcomp_fit.md)
  : Simulate responses from a fitted model
- [`residuals(`*`<regcomp_fit>`*`)`](https://adamrsophie.github.io/regcompR/reference/residuals.regcomp_fit.md)
  : Residuals of a fitted model
- [`vcov(`*`<regcomp_fit>`*`)`](https://adamrsophie.github.io/regcompR/reference/vcov.regcomp_fit.md)
  : Covariance of the estimates
- [`log_lik()`](https://adamrsophie.github.io/regcompR/reference/log_lik.md)
  : Pointwise log likelihood
- [`ranef_regcomp()`](https://adamrsophie.github.io/regcompR/reference/ranef_regcomp.md)
  : Random effect estimates
- [`plot(`*`<regcomp_fit>`*`)`](https://adamrsophie.github.io/regcompR/reference/plot.regcomp_fit.md)
  : Plot a fitted model

## Model checking and comparison

- [`scoring_rule()`](https://adamrsophie.github.io/regcompR/reference/scoring_rule.md)
  : Proper scoring rules for predictive distributions
- [`cv_regcomp()`](https://adamrsophie.github.io/regcompR/reference/cv_regcomp.md)
  : Cross validation of a fitted model
- [`compare_models()`](https://adamrsophie.github.io/regcompR/reference/compare_models.md)
  : Compare fitted models
- [`identifiability_check()`](https://adamrsophie.github.io/regcompR/reference/identifiability_check.md)
  : Identifiability diagnostics of a fit

## Bayesian code generation

- [`stan_code()`](https://adamrsophie.github.io/regcompR/reference/stan_code.md)
  : Generate the Stan program of a model
- [`bugs_code()`](https://adamrsophie.github.io/regcompR/reference/bugs_code.md)
  : Generate the BUGS program of a model

## Utilities

- [`log1mexp()`](https://adamrsophie.github.io/regcompR/reference/numerics.md)
  [`log1pexp()`](https://adamrsophie.github.io/regcompR/reference/numerics.md)
  [`log_sum_exp()`](https://adamrsophie.github.io/regcompR/reference/numerics.md)
  [`log_diff_exp()`](https://adamrsophie.github.io/regcompR/reference/numerics.md)
  : Stable numerical building blocks
- [`regcomp-conditions`](https://adamrsophie.github.io/regcompR/reference/regcomp-conditions.md)
  : Conditions raised by regcompR
- [`regcompR`](https://adamrsophie.github.io/regcompR/reference/regcompR-package.md)
  [`regcompR-package`](https://adamrsophie.github.io/regcompR/reference/regcompR-package.md)
  : regcompR: Distributional Regression with a Verified Registry of
  Probability Distributions
