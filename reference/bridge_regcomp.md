# Marginal likelihood and Bayes factors by bridge sampling

Estimates the log marginal likelihood of a Bayesian fit made with the
rstan or cmdstanr runner by bridge sampling (Meng and Wong 1996; Gronau
et al. 2017) with the bridgesampling package. The model is re-expressed
with all normalizing constants (sampling statements as `target +=`,
priors on standard deviations and correlations renormalized to their
support), so marginal likelihoods of different models are comparable.
`bayes_factor_regcomp()` gives the Bayes factor of two fits. Run long
chains (several thousand draws) for stable estimates.

## Usage

``` r
bridge_regcomp(object, ...)

bayes_factor_regcomp(object1, object2, ...)
```

## Arguments

- object, object1, object2:

  Bayesian fits (`engine = "bayes"`, Stan runner rstan or cmdstanr).

- ...:

  Passed to `bridgesampling::bridge_sampler()`.

## Value

`bridge_regcomp()`: a `bridge` object with the estimate `logml`;
`bayes_factor_regcomp()`: a data frame with the Bayes factor of the
first fit against the second and its logarithm.
