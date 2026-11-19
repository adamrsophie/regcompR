# Parameter specification

Describes one parameter of a distribution: its domain, links, role and
how it enters regression. Used when defining distributions with
[`custom_family()`](https://adamrsophie.github.io/regcompR/reference/custom_family.md).

## Usage

``` r
par_spec(
  name,
  domain = "real",
  role = "shape",
  link = NULL,
  links = NULL,
  symbol = name,
  is_mean = FALSE,
  offset_allowed = FALSE,
  regression = TRUE,
  data = FALSE,
  interpretation = "",
  lower = NULL,
  upper = NULL,
  default = NULL,
  estimate = NULL
)
```

## Arguments

- name:

  Parameter name used in formulas and
  [`distribution()`](https://adamrsophie.github.io/regcompR/reference/distribution.md).

- domain:

  One of `"real"`, `"positive"`, `"nonneg"`, `"unit"`, `"circular"`,
  `"greater_than_one"`, `"greater_than_two"`, `"positive_integer"`,
  `"nonneg_integer"`, `"integer"`, `"unit_closed"` or `"any"`.

- role:

  Free text role such as `"location"`, `"scale"`, `"shape"`,
  `"dispersion"`, `"rate"`, `"skewness"`, `"tail"`, `"inflation"` or
  `"mixing"`.

- link:

  Default link name.

- links:

  Allowed link names.

- symbol:

  Mathematical symbol for documentation.

- is_mean:

  `TRUE` only when the parameter equals the expectation of the response
  for every parameter value. Never inferred from the name.

- offset_allowed:

  Whether an offset may enter this parameter's predictor.

- regression:

  Whether the parameter may have a predictor. Fixed data arguments such
  as binomial trials use `FALSE` together with `data = TRUE`.

- data:

  `TRUE` for data arguments (for example the number of trials).

- interpretation:

  Text shown in documentation.

- lower, upper:

  Optional numeric bounds overriding the domain.

- default:

  Optional value used by
  [`distribution()`](https://adamrsophie.github.io/regcompR/reference/distribution.md)
  when omitted.

- estimate:

  For a non regression parameter, `"ordered"` or `"free"` to estimate it
  as a vector whose length comes from the response (ordinal thresholds);
  `NULL` otherwise.

## Value

An object of class `regcomp_par`.

## Examples

``` r
par_spec("mu", "positive", role = "location", is_mean = TRUE)
#> <parameter mu> domain positive, role location, link log, equals the mean
```
