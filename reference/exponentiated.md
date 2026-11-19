# Generated distribution families

Build a new family from a baseline family `G` through a transformation
of its CDF, `F(y) = T(G(y))`, where `T` is a CDF on `[0, 1]`. The
density is `f(y) = t(G(y)) g(y)` for a continuous baseline and
`P(Y = y) = T(G(y)) - T(G(y - 1))` for a discrete one. New shape
parameters are appended to the baseline parameters (renamed with a
trailing underscore on a clash). All of them can take regression
predictors. Quantiles use the baseline quantile function where `T` has a
closed form inverse.

## Usage

``` r
exponentiated(base)

kumaraswamy_g(base)

beta_g(base)

mcdonald_g(base)

exp_generalized_g(base)

marshall_olkin(base)

transmuted(base)

record_transmuted(base)

odd_log_logistic(base)

odd_burr(base)

weibull_g(base)

gamma_g(base, type = c("zografos", "ristic"))

topp_leone_g(base)

lindley_g(base)

harris_g(base)

sine_g(base)

dus(base)

km(base)

power_series_max(
  base,
  series = c("geometric", "poisson", "logarithmic", "binomial"),
  size = 2
)

power_series_min(
  base,
  series = c("geometric", "poisson", "logarithmic", "binomial"),
  size = 2
)
```

## Arguments

- base:

  Baseline family (identifier or family object).

- type:

  For `gamma_g()`: `"zografos"` or `"ristic"`.

- series:

  Power series type.

- size:

  Binomial series size.

## Value

A `regcomp_family`.

## Details

- `exponentiated()`: `T(u) = u^a` (Lehmann alternative,
  exponentiated-G).

- `kumaraswamy_g()`: `T(u) = 1 - (1 - u^a)^b` (Cordeiro and de Castro
  2011).

- `beta_g()`: `T(u) = I_u(a, b)` (Eugene, Lee and Famoye 2002).

- `mcdonald_g()`: `T(u) = I_{u^c}(a, b)` (generalized beta-G, GB1-G).

- `exp_generalized_g()`: `T(u) = (1 - (1 - u)^a)^b` (Cordeiro, Ortega
  and da Cunha 2013).

- `marshall_olkin()`: `T(u) = u / (alpha + (1 - alpha) u)` (Marshall and
  Olkin 1997).

- `transmuted()`: `T(u) = (1 + lambda) u - lambda u^2`, `|lambda| <= 1`
  (Shaw and Buckley 2007).

- `record_transmuted()`: `T(u) = u + p (1 - u) log(1 - u)` (Balakrishnan
  and He 2021).

- `odd_log_logistic()`: `T(u) = u^a / (u^a + (1 - u)^a)` (Gleaton and
  Lynch 2006).

- `odd_burr()`: `T(u) = 1 - (1 + (u / (1 - u))^a)^(-b)` (Alizadeh et al.
  2017).

- `weibull_g()`: `T(u) = 1 - exp(-a (u / (1 - u))^b)` (Bourguignon,
  Silva and Cordeiro 2014).

- `gamma_g()`: Zografos and Balakrishnan (2009)
  `T(u) = P(a, -log(1 - u))`, or Ristic and Balakrishnan (2012)
  `T(u) = 1 - P(a, -log u)`.

- `topp_leone_g()`: `T(u) = (u (2 - u))^a`.

- `lindley_g()`:
  `T(u) = u^theta (1 + theta - theta log u) / (1 + theta)`.

- `harris_g()`: survival `(theta S^k / (1 - (1 - theta) S^k))^(1/k)`
  (Aly and Benkherouf 2011).

- `sine_g()`: `T(u) = sin(pi u / 2)` (Kumar, Singh and Singh 2015).

- `dus()`: `T(u) = (e^u - 1) / (e - 1)` (Kumar, Singh and Singh 2015).

- `km()`: `T(u) = e (1 - e^-u) / (e - 1)` (Kavya and Manoharan 2021).

- `power_series_max()` and `power_series_min()`: the maximum or minimum
  of a zero truncated power series number of baseline draws, with series
  `"geometric"`, `"poisson"`, `"logarithmic"` or `"binomial"`;
  `T(u) = C(theta u) / C(theta)` and
  `T(u) = 1 - C(theta (1 - u)) / C(theta)`.

Generators nest, so for example `beta_g(kumaraswamy_g("weibull"))` gives
the beta Kumaraswamy Weibull law.

## References

Eugene, Lee and Famoye (2002), Beta-normal distribution and its
applications. Comm. Stat. Theory Methods, 31, 497-512. Cordeiro and de
Castro (2011), A new family of generalized distributions. J. Stat.
Comput. Simul., 81, 883-898. Marshall and Olkin (1997), A new method for
adding a parameter to a family of distributions. Biometrika, 84,
641-652. Shaw and Buckley (2007), The alchemy of probability
distributions. arXiv 0901.0434. Tahir and Cordeiro (2016), Compounding
of distributions: a survey and new generalized classes. J. Stat.
Distrib. Appl., 3, 13.

## Examples

``` r
bw <- beta_g("weibull")
d <- distribution(bw, shape = 1.5, scale = 2, a = 2, b = 0.5)
cdf(d, 1)
#> [1] 0.03725503
kwn <- kumaraswamy_g("normal")
mow <- marshall_olkin("weibull")
ewp <- power_series_max(exponentiated("weibull"), "poisson")
```
