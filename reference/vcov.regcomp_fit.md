# Covariance of the estimates

`type = "model"` gives the inverse observed information (or the
posterior covariance of a Bayesian fit). `type = "sandwich"` gives the
robust covariance `V M V` of Huber and White, where `V` is the model
covariance and `M` the summed outer products of score contributions;
with `cluster` the scores are first summed within clusters (Liang and
Zeger 1986) and `M` is multiplied by `G / (G - 1)` for `G` clusters.
This is the independence working model of generalized estimating
equations for any family; see
[`gee()`](https://adamrsophie.github.io/regcompR/reference/gee.md) for
other working correlations. `type = "quasi"` multiplies the model
covariance by the Pearson dispersion `sum(r^2) / (n - p)` of the Pearson
residuals `r`: the quasi-likelihood covariance of
`glm(family = quasibinomial)` and `quasipoisson` (McCullagh and Nelder
1989, ch. 9) when the family is binomial or Poisson, whose estimates are
the same. `"HC0"` is the sandwich; `"HC1"` multiplies it by
`n / (n - p)` with `p` the number of coefficients; `"HC2"` and `"HC3"`
divide each score by `sqrt(1 - h_i)` and `1 - h_i`, with leverages `h_i`
the diagonal of the weighted hat matrix (MacKinnon and White 1985). For
a one parameter GLM (Poisson, binomial) they equal
[`sandwich::vcovHC()`](https://zeileis.codeberg.page/sandwich/reference/vcovHC.html).

## Usage

``` r
# S3 method for class 'regcomp_fit'
vcov(
  object,
  type = c("model", "sandwich", "quasi", "HC0", "HC1", "HC2", "HC3"),
  cluster = NULL,
  ...
)
```

## Arguments

- object:

  A fitted model.

- type:

  `"model"`, `"sandwich"` (`"HC0"`), `"HC1"`, `"HC2"`, `"HC3"` or
  `"quasi"`.

- cluster:

  Optional cluster identifier: a vector, or a one sided formula
  evaluated in the model data.

- ...:

  Unused.

## Value

A covariance matrix of the coefficients.

## Examples

``` r
set.seed(1)
d <- data.frame(g = rep(1:30, each = 5), x = rnorm(150))
d$y <- rpois(150, exp(0.3 + 0.5 * d$x + rnorm(30, 0, 0.5)[d$g]))
fit <- regcomp(y ~ x, family = "poisson", data = d)
sqrt(diag(vcov(fit, type = "sandwich", cluster = ~ g)))
#> lambda:(Intercept)           lambda:x 
#>         0.10418495         0.08933988 
```
