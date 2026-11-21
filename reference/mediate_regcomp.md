# Causal mediation analysis

Average causal mediation effects (ACME), average direct effects (ADE),
the total effect and the proportion mediated from a mediator model and
an outcome model, by the quasi-Bayesian Monte Carlo method of Imai,
Keele and Tingley (2010), as `mediation::mediate()`: coefficients are
drawn from their asymptotic normal distribution, mediators are simulated
from the mediator model under treatment and control, and outcomes are
predicted under each combination. Any families (and distributional
models) without random effects.

## Usage

``` r
mediate_regcomp(
  med,
  out,
  treat,
  mediator,
  treated = 1,
  control = 0,
  sims = 1000,
  level = 0.95,
  seed = 1
)
```

## Arguments

- med:

  Fitted mediator model (with the treatment as covariate).

- out:

  Fitted outcome model (with treatment and mediator).

- treat:

  Name of the treatment variable.

- mediator:

  Name of the mediator variable (the response of `med`).

- treated, control:

  Treatment values compared.

- sims:

  Number of Monte Carlo draws.

- level:

  Interval level.

- seed:

  Random seed.

## Value

A data frame with estimates, intervals and two sided Monte Carlo p
values for ACME and ADE under control and treatment, their averages, the
total effect and the proportion mediated (the median of its draws).

## References

Imai, Keele and Tingley (2010), A general approach to causal mediation
analysis. Psychological Methods, 15, 309-334.

## Examples

``` r
set.seed(1); n <- 300
d <- data.frame(t = rbinom(n, 1, 0.5), x = rnorm(n))
d$m <- rnorm(n, 0.5 * d$t + 0.3 * d$x); d$y <- rpois(n, exp(0.2 + 0.3 * d$m + 0.2 * d$t))
mm <- regcomp(m ~ t + x, family = "normal", data = d)
my <- regcomp(y ~ t + m + x, family = "poisson", data = d)
mediate_regcomp(mm, my, "t", "m", sims = 200)
#>                 estimate      lower     upper p_value
#> ACME (control) 0.1984468 0.07616170 0.3416720    0.00
#> ACME (treated) 0.2524252 0.09527715 0.4505869    0.00
#> ADE (control)  0.3529925 0.05475094 0.6288538    0.01
#> ADE (treated)  0.4069709 0.06451547 0.7124198    0.01
#> Total effect   0.6054177 0.27320196 0.9423842    0.00
#> ACME (average) 0.2254360 0.08521905 0.3823918    0.00
#> ADE (average)  0.3799817 0.05963320 0.6707096    0.01
#> Prop. mediated 0.3716966 0.15453078 0.7811429    0.00
```
