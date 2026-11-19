# Security Policy

## Supported Versions

Security fixes are considered for the active development branch and the
most recent release.

| Version   | Supported |
|-----------|-----------|
| `0.1.x`   | Yes       |
| `< 0.1.0` | No        |

## Reporting a Vulnerability

Please do not report security vulnerabilities in public GitHub issues.

Report suspected vulnerabilities by email to `adamrsophie@gmail.com`, or
use GitHub private vulnerability reporting if it is enabled for the
repository.

Include as much detail as you can:

- affected version or commit,
- operating system and R version,
- minimal reproduction steps,
- any proof-of-concept code.

## Response Process

The maintainer will acknowledge reports as soon as practical,
investigate the issue, and coordinate a fix if the report is confirmed.
Public disclosure should wait until a fix or mitigation is available,
unless there is an overriding user safety reason to disclose earlier.

## Scope

Security-sensitive areas include:

- generated code: Stan and BUGS programs are written from model formulas
  and family definitions and compiled or run by external tools; a way
  for a formula, factor level or variable name to inject code into a
  generated program would be a vulnerability,
- user supplied functions in
  [`custom_family()`](https://adamrsophie.github.io/regcompR/reference/custom_family.md)
  and
  [`transform_dist()`](https://adamrsophie.github.io/regcompR/reference/transform_dist.md)
  run as given; running code the user did not supply would be a
  vulnerability,
- arbitrary file writes or path traversal, for example in temporary
  files used by the Stan and BUGS runners,
- unbounded memory or time on small inputs, for example in numerical
  integration or series summation,
- exposure of private data in examples, tests, or logs.

A statistically wrong answer that does not create a software-security
vulnerability should be opened as a normal GitHub issue instead.
