# regcompR: Distributional Regression with a Verified Registry of Probability Distributions

A framework for distributional regression in which every parameter of a
response distribution can depend on covariates through its own formula
(Rigby and Stasinopoulos 2005
[doi:10.1111/j.1467-9876.2005.00510.x](https://doi.org/10.1111/j.1467-9876.2005.00510.x)
). Distributions are held in a machine readable registry that records
support, parameter domains, links, moments and their existence
conditions, available methods, references and validation status. Count,
binary, bounded, categorical, ordinal, real line, positive, unit
interval, simplex, circular, extreme value, multivariate and matrix
valued laws are provided, each with density or mass, distribution
function, quantile, random generation and moments where these exist.
Operators build new laws from old ones: inflation, hurdles, truncation,
censoring, monotone transformation, discretization, rounding, finite
mixtures, convolution, compounding and splicing. Models are fitted by
maximum or penalized likelihood with a Laplace approximation for random
effects and smooth terms, optionally through 'TMB' (Kristensen et al.
2016 [doi:10.18637/jss.v070.i05](https://doi.org/10.18637/jss.v070.i05)
), or by Bayesian inference through generated 'Stan' programs run by
'bayescompR', 'cmdstanr' or 'rstan'. Predictions, quantiles, event
probabilities, randomized quantile residuals, proper scoring rules,
cross validation, and plots for every distribution and fit are included.

## See also

Useful links:

- <https://github.com/adamrsophie/regcompR>

- <https://adamrsophie.github.io/regcompR/>

- Report bugs at <https://github.com/adamrsophie/regcompR/issues>

## Author

**Maintainer**: Adam R Sophie <adamrsophie@gmail.com> \[copyright
holder\]

Authors:

- Adam R Sophie <adamrsophie@gmail.com> \[copyright holder\]
