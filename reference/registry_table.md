# Registry table

One row per registered distribution with its category, variate type,
parameters, validation status and the methods available in R, the
maximum likelihood engines and Stan. Capabilities are computed from the
specification, never asserted separately.

## Usage

``` r
registry_table(category = NULL)
```

## Arguments

- category:

  Optional category filter.

## Value

A data frame.

## Examples

``` r
tab <- registry_table()
table(tab$category)
#> 
#>                bounded            categorical               circular 
#>                     35                      3                     10 
#>                  count            directional        finite_discrete 
#>                     51                      4                     16 
#>                  mixed           multivariate                ordinal 
#>                      2                     22                     10 
#>               positive quantile_parameterized                   real 
#>                     65                      5                     56 
head(tab[, c("id", "status", "cdf", "quantile", "stan")])
#>                        id    status         cdf    quantile        stan
#> 1 adjacent_category_logit validated unavailable unavailable   generated
#> 2                 arcsine validated    analytic    analytic unavailable
#> 3      asymmetric_laplace validated    analytic    analytic   generated
#> 4         axial_von_mises validated     numeric     numeric unavailable
#> 5         balding_nichols validated    analytic    analytic unavailable
#> 6                   bates validated    analytic     numeric unavailable
```
