# Control settings for model fitting

Control settings for model fitting

## Usage

``` r
regcomp_control(
  maxit = 1000,
  rel_tol = 1e-10,
  inner_maxit = 50,
  inner_tol = 1e-08,
  restarts = 2,
  hessian = TRUE,
  trace = FALSE,
  chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  seed = 1,
  runner = c("bayescompR", "cmdstanr", "rstan", "jags", "openbugs"),
  cores = 1,
  adapt_delta = 0.8,
  start = NULL,
  prox_maxit = 5000,
  prox_tol = 1e-08
)
```

## Arguments

- maxit:

  Maximum outer iterations of
  [`nlminb()`](https://rdrr.io/r/stats/nlminb.html).

- rel_tol:

  Relative tolerance of the outer optimization.

- inner_maxit:

  Maximum Newton iterations of the Laplace inner problem.

- inner_tol:

  Convergence tolerance of the inner problem.

- restarts:

  Number of random restarts tried when the first optimization fails.

- hessian:

  Compute the covariance of the estimates.

- trace:

  Print progress.

- chains, iter_warmup, iter_sampling, seed, runner, cores:

  Bayesian settings: chains, iterations, seed, the runner
  (`"bayescompR"`, `"cmdstanr"` or `"rstan"` for Stan, `"jags"` through
  'R2jags' or `"openbugs"` through 'R2OpenBUGS') and parallel chains.

- adapt_delta:

  Target acceptance for NUTS.

- start:

  Optional named list of starting values for the coefficients of each
  parameter (vectors on the link scale).

- prox_maxit, prox_tol:

  Iteration limit and relative step tolerance of the proximal gradient
  (FISTA) solver used for lasso type penalties (`maxit` and `rel_tol`
  apply to [`nlminb()`](https://rdrr.io/r/stats/nlminb.html) only).

## Value

A list.

## Examples

``` r
regcomp_control(maxit = 500)
#> $maxit
#> [1] 500
#> 
#> $rel_tol
#> [1] 1e-10
#> 
#> $inner_maxit
#> [1] 50
#> 
#> $inner_tol
#> [1] 1e-08
#> 
#> $restarts
#> [1] 2
#> 
#> $hessian
#> [1] TRUE
#> 
#> $trace
#> [1] FALSE
#> 
#> $chains
#> [1] 4
#> 
#> $iter_warmup
#> [1] 1000
#> 
#> $iter_sampling
#> [1] 1000
#> 
#> $seed
#> [1] 1
#> 
#> $runner
#> [1] "bayescompR"
#> 
#> $cores
#> [1] 1
#> 
#> $adapt_delta
#> [1] 0.8
#> 
#> $start
#> NULL
#> 
#> $prox_maxit
#> [1] 5000
#> 
#> $prox_tol
#> [1] 1e-08
#> 
```
