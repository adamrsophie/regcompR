## Submission

This is the first submission of regcompR.

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission. 'cmdstanr' is a suggested package outside the
  mainstream repositories, available from the `Additional_repositories`
  given in DESCRIPTION.

## Test environments

* local macOS (aarch64), R release
* GitHub Actions: ubuntu-latest (devel, release, oldrel-1), windows-latest (release)

## Notes for the reviewers

* The many suggested packages are reference implementations used only in
  tests (each guarded by `skip_if_not_installed()`), optional engines
  ('TMB', 'RTMB', 'rstan', 'cmdstanr', 'R2jags') and optional table and plot
  outputs. None are needed for the core functionality, which depends only
  on base R and 'Matrix'.
* 'cmdstanr' is available from the Stan r-universe, given in
  `Additional_repositories`.
* The optional Stan runner 'bayescompR' and the optional 'R2OpenBUGS'
  runner are not declared dependencies; they are called only when the user
  selects them and they are installed.
* Examples, tests and vignettes that would run long or need a Stan or BUGS
  installation are skipped on CRAN or wrapped in `\donttest{}`.
