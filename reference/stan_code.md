# Generate the Stan program of a model

Returns the 'Stan' code that the `"bayes"` engine runs, together with
the data list. The program is generated from the same specification as
the frequentist engines: each distribution parameter has a linear
predictor with its link, random effects use a non-centered
parameterization with Cholesky factored correlations, smooth terms enter
through their mixed model representation, and the likelihood uses native
'Stan' functions when the parameterization is identical or package
generated functions copied from the registry otherwise. Censored and
truncated observations use log CDF and log survival functions.

## Usage

``` r
stan_code(
  formula,
  ...,
  family,
  data,
  weights = NULL,
  links = NULL,
  fixed = NULL,
  priors = NULL,
  parameterization = NULL
)
```

## Arguments

- formula, ..., family, data, weights, links, fixed, parameterization:

  As in
  [`regcomp()`](https://adamrsophie.github.io/regcompR/reference/regcomp.md).

- priors:

  Priors from
  [`prior()`](https://adamrsophie.github.io/regcompR/reference/prior.md).

## Value

A list with `code` and `data`.

## Examples

``` r
d <- data.frame(y = rpois(20, 3), x = rnorm(20))
cat(stan_code(y ~ x, family = "poisson", data = d)$code)
#> functions {
#> 
#> }
#> data {
#>   int<lower=1> N;
#>   array[N] int y;
#>   vector[N] w;
#>   int<lower=1> P_lambda;
#>   matrix[N, P_lambda] X_lambda;
#>   vector[N] off_lambda;
#> }
#> parameters {
#>   vector[P_lambda] b_lambda;
#> }
#> model {
#>   vector[N] eta_lambda = X_lambda * b_lambda + off_lambda;
#>   vector[N] lambda = exp(eta_lambda);
#>   b_lambda[1] ~ normal(0, 5);
#>   b_lambda[2] ~ normal(0, 5);
#>   for (i in 1:N) target += w[i] * poisson_lpmf(y[i] | lambda[i]);
#> }
```
