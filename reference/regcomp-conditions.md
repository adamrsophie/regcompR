# Conditions raised by regcompR

All errors raised on purpose by the package inherit from
`regcomp_condition` and one of the classes below, so they can be caught
selectively with [`tryCatch()`](https://rdrr.io/r/base/conditions.html).

## Details

- `regcomp_unsupported`: an operation (for example a quantile function
  or a backend) is not available for a distribution or model. The
  package never substitutes an unjustified approximation.

- `regcomp_domain_error`: a parameter or response value lies outside its
  domain or support.

- `regcomp_numerical_failure`: a series, root finder or integral did not
  reach its declared tolerance.

- `regcomp_input_error`: malformed arguments.

Warnings use the classes `regcomp_undefined_moment` (a moment that does
not exist for the given parameters; `NA` is returned) and
`regcomp_identifiability`.

## Examples

``` r
res <- tryCatch(mean(distribution("cauchy", mu = 0, sigma = 1)),
                regcomp_undefined_moment = function(w) "undefined")
res
#> [1] "undefined"
```
