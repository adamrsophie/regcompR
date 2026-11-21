# Individual uplift predictions

The predicted mean response with the treatment set to its treated value
minus the prediction with it set to the control value, for each row: the
individual uplift (conditional average treatment effect) of a model with
treatment by covariate interactions.

## Usage

``` r
uplift(object, treat, newdata = NULL, treated = NULL, control = NULL)
```

## Arguments

- object:

  A fitted model.

- treat:

  Name of the treatment variable (0/1, logical or a two level factor).

- newdata:

  Optional data (default the model data).

- treated, control:

  Values of `treat` for treatment and control (default 1 and 0, or the
  second and first factor levels).

## Value

A numeric vector of uplifts.

## Examples

``` r
set.seed(1); n <- 400
d <- data.frame(t = rbinom(n, 1, 0.5), x = rnorm(n))
d$y <- rbinom(n, 1, plogis(-0.5 + 0.3 * d$x + d$t * (0.2 + 0.8 * d$x)))
f <- regcomp(y ~ t * x, family = "bernoulli", data = d)
u <- uplift(f, "t")
qini(u, d$t, d$y)
#>  cum_per T_Y1 T_n C_Y1 C_n incremental_Y1 inc_uplift   uplift
#>      0.1   18  23   11  19         4.6842     2.5597  19.4373
#>      0.2   33  43   20  39        10.9487     5.9829  30.0000
#>      0.3   43  59   30  63        14.9048     8.1447  20.8333
#>      0.4   51  73   41  88        16.9886     9.2834  13.1429
#>      0.5   63  94   50 108        19.4815    10.6456  12.1429
#>      0.6   69 110   61 132        18.1667     9.9271  -8.3333
#>      0.7   75 130   67 152        17.6974     9.6707   0.0000
#>      0.8   80 151   75 171        13.7719     7.5256 -18.2957
#>      0.9   85 169   82 193        13.1969     7.2114  -4.0404
#>      1.0   87 183   87 217        13.6313     7.4488  -6.5476
#> Qini coefficient: 3.573
```
