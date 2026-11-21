# Random intercepts with non-normal distributions

Fits a model with one random intercept per level of `group` in the first
parameter's linear predictor, whose distribution is not normal: a
log-gamma effect (`exp(b)` gamma with mean 1 and variance `v`, the gamma
frailty of a Poisson or survival model), a Student t effect with scale
`s` and estimated degrees of freedom, a logistic or a Laplace effect
with scale `s`. The marginal likelihood is computed by adaptive
Gauss-Hermite quadrature over `b = F^-1(Phi(z))` with standard normal
`z`, so any family of the package can be used.

## Usage

``` r
re_nongauss(
  formula,
  group,
  family,
  data,
  re = c("gamma", "t", "logistic", "laplace"),
  nq = 20,
  ...
)
```

## Arguments

- formula:

  Model formula for the first parameter, without random effects; other
  parameters (`sigma ~ x`) go in `...`.

- group:

  A one sided formula naming the grouping factor, `~ g`.

- family:

  Family (any regcompR family).

- data:

  A data frame.

- re:

  `"gamma"`, `"t"`, `"logistic"` or `"laplace"`.

- nq:

  Number of quadrature nodes.

- ...:

  Formulas for the other parameters.

## Value

An object of class `regcomp_nongauss` with coefficients (including
`re:log_v` or `re:log_scale`, and `re:log_df_minus_2` for t effects),
their covariance, the log likelihood, AIC and the predicted random
effects (posterior means).

## Examples

``` r
set.seed(1)
d <- data.frame(g = factor(rep(1:30, each = 5)), x = rnorm(150))
d$y <- rpois(150, exp(0.3 + 0.4 * d$x) * rgamma(30, 2, 2)[d$g])
re_nongauss(y ~ x, ~ g, family = "poisson", data = d, re = "gamma")
#> <regcomp_nongauss> poisson model with a gamma random intercept, n = 150
#>                    estimate std_error       z
#> lambda:(Intercept)   0.3724    0.1444  2.5790
#> lambda:x             0.3374    0.0750  4.4978
#> re:log_v            -0.7327    0.3396 -2.1578
#> log likelihood -232.933, AIC 471.866
```
