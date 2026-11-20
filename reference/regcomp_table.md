# Regression tables

Builds a table of estimates for one or several fits, with rows grouped
by distribution parameter, random effect standard deviations and fit
statistics. It prints in the console and renders as a formatted table in
R Markdown and Quarto documents (HTML, PDF and Word output) through
[`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html);
[`export_table()`](https://adamrsophie.github.io/regcompR/reference/export_table.md)
writes it to a file.

## Usage

``` r
regcomp_table(
  ...,
  conf.level = 0.95,
  exponentiate = FALSE,
  show = c("ci", "se"),
  stars = FALSE,
  digits = 3,
  stats = c("nobs", "logLik", "AIC", "BIC"),
  vcov = "model",
  standardize = FALSE
)
```

## Arguments

- ...:

  Fitted models, optionally named (the names become column headers).

- conf.level:

  Interval level.

- exponentiate:

  Exponentiate coefficients (all of them, or a logical vector with one
  entry per model).

- show:

  `"ci"` for estimates with intervals, `"se"` for estimates with
  standard errors.

- stars:

  Add significance stars (`* p < 0.05`, `** p < 0.01`, `*** p < 0.001`)
  to Wald estimates.

- digits:

  Number of decimal places.

- stats:

  Fit statistics to add: any of `"nobs"`, `"logLik"`, `"AIC"`, `"BIC"`;
  Bayesian fits report draws and the largest R-hat instead.

- vcov:

  Standard errors for frequentist fits: `"model"`, `"sandwich"` or
  `"quasi"`.

- standardize:

  Standardized coefficients (see
  [`tidy.regcomp_fit()`](https://adamrsophie.github.io/regcompR/reference/tidy.regcomp_fit.md)).

## Value

A `regcomp_table`: a data frame with columns `parameter`, `term` and one
column per model.

## Examples

``` r
d <- data.frame(x = rnorm(200)); d$y <- rnbinom(200, mu = exp(1 + 0.4 * d$x), size = 2)
tab <- regcomp_table(Poisson = regcomp(y ~ x, family = "poisson", data = d),
                     NB = regcomp(y ~ x, family = "nb2", data = d), exponentiate = TRUE)
tab
#>  parameter term        Poisson              NB                  
#>  lambda    (Intercept) 2.953 [2.718, 3.208]                     
#>            x           1.404 [1.304, 1.513]                     
#>  mu        (Intercept)                      2.961 [2.620, 3.347]
#>            x                                1.388 [1.233, 1.562]
#>  phi       (Intercept)                      2.354 [1.619, 3.422]
#>            nobs        200                  200                 
#>            logLik      -480.856             -437.106            
#>            AIC         965.711              880.213             
#>            BIC         972.308              890.108             
#> 95% confidence intervals in brackets; coefficients exponentiated. 
export_table(tab, tempfile(fileext = ".html"))
```
