# Qini curve and coefficient

Ranks observations by predicted uplift and, for the top `k / groups`
share, computes the incremental number of responders
`R_T - R_C N_T / N_C` and the observed uplift in each group (Radcliffe
2007). The Qini coefficient is the area between the Qini curve
(incremental responders as a share of all treated, by share of treated)
and the line of random targeting. Equal to
`tools4uplift::PerformanceUplift()` and `QiniArea()`.

## Usage

``` r
qini(uplift, treat, outcome, groups = 10)
```

## Arguments

- uplift:

  Predicted uplifts (for example from
  [`uplift()`](https://adamrsophie.github.io/regcompR/reference/uplift.md)).

- treat:

  Treatment indicator (1 treated, 0 control).

- outcome:

  Binary outcome (1 response).

- groups:

  Number of groups.

## Value

An object of class `regcomp_qini`: a data frame with columns `cum_per`,
`T_Y1`, `T_n`, `C_Y1`, `C_n`, `incremental_Y1`, `inc_uplift` and
`uplift` (percent), with the Qini coefficient as attribute `"qini"`.
