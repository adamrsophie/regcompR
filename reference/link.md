# Link functions

Returns a link object by name. Links available: `identity`, `log`,
`softplus`, `sqrt`, `inverse`, `logit`, `probit`, `cloglog`, `cauchit`,
`loglog`, `tan_half` (for circular locations), `logshift1` and
`logshift2` (parameters constrained above 1 or 2).

## Usage

``` r
link(name)

links()
```

## Arguments

- name:

  Link name.

## Value

An object of class `regcomp_link` with elements `linkfun`, `linkinv`,
`mu_eta` and `valid`.

## Examples

``` r
l <- link("logit")
l$linkinv(0)
#> [1] 0.5
links()
#>  [1] "identity"  "log"       "softplus"  "sqrt"      "inverse"   "logit"    
#>  [7] "probit"    "cloglog"   "cauchit"   "loglog"    "tan_half"  "fisherz"  
#> [13] "logshift1" "logshift2"
```
