# Contributing to regcompR

Thank you for considering a contribution to `regcompR`. The package is a
statistical tool, so a change is judged on two things: whether the R
code is sound, and whether its numbers agree with an independent
reference.

## Before You Start

- Use GitHub issues for bug reports, feature requests, documentation
  gaps, and requests for new distributions or models.
- For larger changes, such as a new engine or a family of distributions,
  open an issue first so the design can be discussed before
  implementation.
- Keep pull requests focused. A new family, a bug fix, and a refactor
  should usually be separate pull requests.

## Reporting a Wrong Answer

A wrong result is the most valuable bug report. Please include:

- a reproducible R snippet (simulated data is fine),
- what regcompR returns,
- what a reference implementation returns for the same model (for
  example `stats`, `lme4`, `glmmTMB`, `gamlss`, `VGAM`, `MASS`,
  `ordinal` or `survival`), if you can run it,
- [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html).

## Development Setup

Install the development dependencies from the repository root:

``` r

install.packages(c("remotes", "testthat", "devtools", "roxygen2"))
remotes::install_deps(dependencies = TRUE)
```

regcompR has no compiled code of its own. The TMB engine needs `RTMB`;
the Bayesian engines need a Stan runner (`cmdstanr`, `rstan` or
`bayescompR`), JAGS with `R2jags`, or OpenBUGS with `R2OpenBUGS`.

## Adding a Distribution

- Register the family with `family_def()` and `register_family()` in the
  matching `R/zfam-*.R` file: parameters with their domains and links,
  the log density, and the distribution function, quantile, random
  generation and moments where they exist, with the conditions under
  which moments exist.
- Cite a primary source for the density and its parameterization.
- Add a test that checks normalization and compares the density,
  distribution function and moments with a reference package or with
  simulation. Do not write an expected value by hand when a reference
  can produce it.
- Give the family a Stan density where one can be written, so the
  Bayesian engine supports it.

## Adding a Model or Engine Feature

- Validate estimates, standard errors and the log likelihood against an
  established package on simulated and real data, and keep that
  comparison as a test (skipped when the reference package is not
  installed).
- Keep results identical across engines where the model is the same.

## Branches and Pull Requests

- Use one branch per logical change, named after the change
  (`add-sichel-family`, `fix-ordinal-prediction`).
- Write commit subjects in the imperative (“Add the Sichel family”).
- Keep the implementation, its tests, the roxygen source, the
  regenerated documentation, and the `NEWS.md` bullet together in the
  same pull request.

## Code Style

- Follow the style already used in the surrounding R files.
- Prefer clear validation and explicit errors over silent coercion; use
  the condition helpers in `R/conditions.R`.
- Work on the log scale for densities and use the stable helpers
  ([`log_sum_exp()`](https://adamrsophie.github.io/regcompR/reference/numerics.md),
  [`log1mexp()`](https://adamrsophie.github.io/regcompR/reference/numerics.md),
  [`log1pexp()`](https://adamrsophie.github.io/regcompR/reference/numerics.md)).
- Do not add new package dependencies unless they are necessary and
  justified.
- Regenerate the `.Rd` files with `devtools::document()` after changing
  roxygen comments.

## Testing

Run the tests while developing:

``` r

devtools::test()
```

Before a pull request, run the full check:

``` sh
R CMD build .
R CMD check --as-cran regcompR_*.tar.gz
```

## Pull Request Checklist

Before opening a pull request, confirm:

- Tests relevant to the change pass, with expected values from a
  reference.
- Documentation is updated.
- `NEWS.md` has an entry for user-facing changes.
- No generated build artifacts, local caches, `.Rcheck` directories, or
  tarballs are committed.

## Code of Conduct

All contributors are expected to follow the project [Code of
Conduct](https://adamrsophie.github.io/regcompR/CODE_OF_CONDUCT.md).
