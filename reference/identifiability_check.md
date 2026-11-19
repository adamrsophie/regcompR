# Identifiability diagnostics of a fit

Reports the conditioning of the observed information, the parameter
combinations that are nearly flat (eigenvectors of the smallest
eigenvalues of the information), random effect variances at the
boundary, and for mixtures the possibility of label switching.

## Usage

``` r
identifiability_check(object, tol = 1e-08)
```

## Arguments

- object:

  A frequentist fit.

- tol:

  Eigenvalue ratio below which a direction is reported.

## Value

A list with `condition_number`, `weak_directions` (a list of named
loading vectors), `boundary_variances` and `notes`.
