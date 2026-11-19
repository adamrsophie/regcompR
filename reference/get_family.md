# Look up a registered family

`get_family()` returns the family specification for a registered
distribution, by identifier or alias. `families()` lists identifiers,
optionally filtered by category or status.

## Usage

``` r
get_family(id, parameterization = NULL)

families(category = NULL, status = NULL)
```

## Arguments

- id:

  Identifier or alias, for example `"negative_binomial_nb2"` or `"nb2"`.

- parameterization:

  Optional parameterization name. When given, the identifier
  `paste(id, parameterization, sep = "_")` is tried as well, so
  `get_family("negative_binomial", "nb2")` works.

- category:

  Optional category filter.

- status:

  Optional status filter.

## Value

A `regcomp_family` object, or a character vector for `families()`.

## Examples

``` r
get_family("poisson")
#> <regcomp_family poisson> Poisson
#>   category: count, variate: scalar, discrete, status: validated
#>   parameters:
#>     lambda   positive          link log       [mean]
#>   formula: p(y) = exp(-lambda) lambda^y / y! 
get_family("negative_binomial", "nb2")
#> <regcomp_family negative_binomial_nb2> Negative binomial (NB2)
#>   category: count, variate: scalar, discrete, status: validated
#>   parameters:
#>     mu       positive          link log       [mean]
#>     phi      positive          link log       
#>   formula: p(y) = C(y + phi - 1, y) (mu / (mu + phi))^y (phi / (mu + phi))^phi 
head(families(category = "count"))
#> [1] "bell"                        "beta_negative_binomial"     
#> [3] "borel"                       "borel_tanner"               
#> [5] "conway_maxwell_poisson"      "conway_maxwell_poisson_mean"
```
