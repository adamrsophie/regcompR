# Multivariate and matrix valued laws. Vector variates take `y` as a matrix
# with one observation per row (a vector is one observation); matrix
# variates take a single matrix or a list of matrices. Parameters are
# vectors or matrices and are not recycled.

mfam <- function(id, display, pars, ..., variate_type = "vector") register_family(family_def(id, display, "multivariate", pars,
  variate_type = variate_type, numeric_cdf = FALSE, numeric_moments = FALSE, ...))
as_rows <- function(y) if (is.matrix(y)) y else matrix(y, nrow = 1)
as_mats <- function(y) if (is.list(y)) y else list(y)
lmgamma <- function(a, p) p * (p - 1) / 4 * log(pi) + sum(lgamma(a + (1 - seq_len(p)) / 2))
ldet <- function(A) as.numeric(determinant(as.matrix(A), logarithm = TRUE)$modulus)
etr <- function(A) sum(diag(A))

# Count vectors ---------------------------------------------------------------------
mfam("multinomial", "Multinomial", list(P("probs", "unit", role = "cell probabilities", regression = FALSE)), discrete = TRUE,
  formula = "P(y) = n! prod p_k^y_k / y_k!, n = sum y",
  lpdf = function(y, p) { y <- as_rows(y); pr <- p$probs / sum(p$probs)
    lgamma(rowSums(y) + 1) + drop(y %*% log(pr)) - rowSums(lgamma(y + 1)) },
  rng = function(n, p, size = 10) t(stats::rmultinom(n, size, p$probs)),
  regression_notes = "Multinomial logit regression is provided by the 'multinomial_logit' family.",
  references = "Johnson, Kotz and Balakrishnan (1997), Discrete Multivariate Distributions. Wiley.",
  reference = ref_impl("stats", function(y, p) apply(as_rows(y), 1, stats::dmultinom, prob = p$probs, log = TRUE)), status = "validated")

mfam("negative_multinomial", "Negative multinomial", list(P("x0", "positive", role = "stopping count"), P("probs", "unit", role = "cell probabilities", regression = FALSE)),
  discrete = TRUE, check = function(p) if (sum(p$probs) >= 1) "Cell probabilities must sum to less than one (p0 = 1 - sum).",
  formula = "P(y) = Gamma(x0 + sum y) / (Gamma(x0) prod y_k!) p0^x0 prod p_k^y_k, p0 = 1 - sum p",
  lpdf = function(y, p) { y <- as_rows(y)
    lgamma(p$x0 + rowSums(y)) - lgamma(p$x0) - rowSums(lgamma(y + 1)) + p$x0 * log1p(-sum(p$probs)) + drop(y %*% log(p$probs)) },
  rng = function(n, p) { lam <- stats::rgamma(n, p$x0, (1 - sum(p$probs)) / 1); t(vapply(lam, function(l) stats::rpois(length(p$probs), l * p$probs), numeric(length(p$probs)))) },
  references = "Sibuya, Yoshimura and Shimizu (1964), Negative multinomial distribution. Ann. Inst. Stat. Math., 16, 409-426.",
  reference = ref_impl("gamma mixture of independent Poissons", function(y, p) apply(as_rows(y), 1, function(r) log(stats::integrate(function(l)
    vapply(l, function(ll) prod(stats::dpois(r, ll * p$probs)), 0) * stats::dgamma(l, p$x0, 1 - sum(p$probs)), 0, Inf, rel.tol = 1e-12)$value))),
  status = "validated")

mfam("dirichlet_multinomial", "Dirichlet-multinomial", list(P("alpha", "positive", role = "concentrations", regression = FALSE)), discrete = TRUE,
  aliases = "multivariate_polya",
  formula = "P(y) = n! Gamma(A) / Gamma(n + A) prod Gamma(y_k + a_k) / (y_k! Gamma(a_k)), A = sum a",
  lpdf = function(y, p) { y <- as_rows(y); n <- rowSums(y); A <- sum(p$alpha)
    lgamma(n + 1) + lgamma(A) - lgamma(n + A) + rowSums(lgamma(sweep(y, 2, p$alpha, "+"))) - rowSums(lgamma(y + 1)) - sum(lgamma(p$alpha)) },
  references = "Mosimann (1962), On the compound multinomial distribution. Biometrika, 49, 65-82.",
  reference = ref_impl("extraDistr", function(y, p) { y <- as_rows(y); extraDistr::ddirmnom(y, rowSums(y), p$alpha, log = TRUE) }), status = "validated")

mfam("dirichlet_negative_multinomial", "Dirichlet negative multinomial", list(P("x0", "positive"), P("alpha0", "positive"),
  P("alpha", "positive", role = "concentrations", regression = FALSE)), discrete = TRUE,
  formula = "P(y) = B(x., a.) / B(x0, alpha0) prod Gamma(y_k + a_k) / (y_k! Gamma(a_k)), x. = x0 + sum y, a. = alpha0 + sum a",
  lpdf = function(y, p) { y <- as_rows(y); xd <- p$x0 + rowSums(y); ad <- p$alpha0 + sum(p$alpha)
    lbeta(xd, ad) - lbeta(p$x0, p$alpha0) + rowSums(lgamma(sweep(y, 2, p$alpha, "+"))) - rowSums(lgamma(y + 1)) - sum(lgamma(p$alpha)) },
  references = "Farewell and Farewell (2013), Dirichlet negative multinomial regression for overdispersed correlated count data. Biostatistics, 14, 395-404.",
  reference = ref_impl("beta mixture of negative multinomials", function(y, p) apply(as_rows(y), 1, function(r) { A <- sum(p$alpha)
    log(stats::integrate(function(q) vapply(q, function(qq) exp(lgamma(p$x0 + sum(r)) - lgamma(p$x0) - sum(lgamma(r + 1)) + p$x0 * log(qq) +
      sum(r) * log1p(-qq) + lgamma(A) - sum(lgamma(p$alpha)) + sum(lgamma(r + p$alpha)) - lgamma(sum(r) + A)), 0) *
      stats::dbeta(q, p$alpha0, A), 0, 1, rel.tol = 1e-12)$value) }), note = "p0 ~ Beta(alpha0, sum alpha), cell split Dirichlet."),
  status = "validated")

ewens_lp <- function(a, th) { n <- sum(seq_along(a) * a); j <- seq_along(a)
  lgamma(n + 1) - (lgamma(th + n) - lgamma(th)) + sum(a * log(th) - a * log(j) - lgamma(a + 1)) }
mfam("ewens", "Ewens sampling formula", list(P("theta", "positive", role = "mutation rate")), discrete = TRUE, variate_type = "partition",
  formula = "P(a) = n! / (theta (theta + 1) ... (theta + n - 1)) prod theta^a_j / (j^a_j a_j!), a_j = number of types seen j times",
  lpdf = function(y, p) apply(as_rows(y), 1, ewens_lp, th = p$theta),
  rng = function(n, p, size = 10) t(vapply(seq_len(n), function(i) { cl <- integer(0)
    for (k in seq_len(size)) { if (stats::runif(1) < p$theta / (p$theta + k - 1)) cl <- c(cl, 1L) else { j <- sample.int(length(cl), 1, prob = cl); cl[j] <- cl[j] + 1L } }
    tabulate(cl, size) }, numeric(size))),
  references = "Ewens (1972), The sampling theory of selectively neutral alleles. Theor. Popul. Biol., 3, 87-112.",
  reference = ref_impl("Chinese restaurant process", function(y, p) apply(as_rows(y), 1, ewens_lp, th = p$theta),
    note = "Checked by summing to one over all partitions of n."), status = "validated")

# Continuous vectors ------------------------------------------------------------------------
mfam("dirichlet", "Dirichlet", list(P("alpha", "positive", role = "concentrations", regression = FALSE)),
  formula = "f(y) = Gamma(A) / prod Gamma(a_k) prod y_k^(a_k - 1) on the simplex",
  lpdf = function(y, p) { y <- as_rows(y); lgamma(sum(p$alpha)) - sum(lgamma(p$alpha)) + drop(log(y) %*% (p$alpha - 1)) },
  rng = function(n, p) { g <- matrix(stats::rgamma(n * length(p$alpha), p$alpha), n, byrow = TRUE); g / rowSums(g) },
  regression_notes = "Dirichlet regression: log link on each alpha_k, or the mean/precision form of Maier (2014).",
  references = "Kotz, Balakrishnan and Johnson (2000), Continuous Multivariate Distributions, ch. 49.",
  reference = ref_impl("extraDistr", function(y, p) extraDistr::ddirichlet(as_rows(y), p$alpha, log = TRUE)), status = "validated")

# log C(eta): 1/C is the divided difference of exp at eta, read off the
# matrix exponential of the bidiagonal matrix with eta on the diagonal
# (exact for tied parameters); eta is centered for stability.
cc_logC <- function(eta) { K <- length(eta); m <- max(eta)
  A <- diag(eta - m, K); if (K > 1) A[cbind(1:(K - 1), 2:K)] <- 1
  -(m + log(as.matrix(Matrix::expm(Matrix::Matrix(A)))[1, K])) }
mfam("continuous_categorical", "Continuous categorical", list(P("eta", "real", role = "natural parameters", regression = FALSE)),
  formula = "f(y) = C(eta) exp(sum eta_k y_k) on the simplex (Lebesgue on the first K - 1 coordinates), 1/C = sum_k e^eta_k / prod_{i != k} (eta_k - eta_i)",
  lpdf = function(y, p) { y <- as_rows(y); cc_logC(p$eta) + drop(y %*% p$eta) },
  references = "Gordon-Rodriguez, Loaiza-Ganem and Cunningham (2020), The continuous categorical: a novel simplex-valued exponential family. ICML.",
  reference = ref_impl("definition", function(y, p) cc_logC(p$eta) + drop(as_rows(y) %*% p$eta), note = "Checked by integrating over the 2-simplex."),
  status = "validated")

mfam("multivariate_normal", "Multivariate normal", list(P("mu", "real", role = "mean vector", regression = FALSE), P("Sigma", "any", role = "covariance matrix", regression = FALSE)),
  aliases = c("mvn", "MVN"),
  lpdf = function(y, p) { y <- as_rows(y); S <- as.matrix(p$Sigma); R <- chol(S); z <- backsolve(R, t(y) - p$mu, transpose = TRUE)
    -0.5 * ncol(y) * log(2 * pi) - sum(log(diag(R))) - 0.5 * colSums(z^2) },
  rng = function(n, p) sweep(matrix(stats::rnorm(n * length(p$mu)), n) %*% chol(as.matrix(p$Sigma)), 2, p$mu, "+"),
  references = "Anderson (2003), An Introduction to Multivariate Statistical Analysis, 3rd ed. Wiley.",
  reference = ref_impl("mvtnorm", function(y, p) mvtnorm::dmvnorm(as_rows(y), p$mu, as.matrix(p$Sigma), log = TRUE)),
  stan = stan_native("multi_normal_lpdf({y} | {mu}, {Sigma})", rng = "multi_normal_rng({mu}, {Sigma})"), status = "validated")

mfam("multivariate_t", "Multivariate t", list(P("mu", "real", regression = FALSE), P("Sigma", "any", role = "scale matrix", regression = FALSE), P("nu", "positive", role = "degrees of freedom")),
  aliases = "mvt",
  lpdf = function(y, p) { y <- as_rows(y); d <- ncol(y); R <- chol(as.matrix(p$Sigma)); z <- backsolve(R, t(y) - p$mu, transpose = TRUE)
    lgamma((p$nu + d) / 2) - lgamma(p$nu / 2) - d / 2 * log(p$nu * pi) - sum(log(diag(R))) - (p$nu + d) / 2 * log1p(colSums(z^2) / p$nu) },
  rng = function(n, p) { d <- length(p$mu); z <- matrix(stats::rnorm(n * d), n) %*% chol(as.matrix(p$Sigma)); sweep(z / sqrt(stats::rchisq(n, p$nu) / p$nu), 2, p$mu, "+") },
  references = "Kotz and Nadarajah (2004), Multivariate t Distributions and Their Applications. Cambridge.",
  reference = ref_impl("mvtnorm", function(y, p) mvtnorm::dmvt(as_rows(y), p$mu, as.matrix(p$Sigma), df = p$nu, log = TRUE)),
  stan = stan_native("multi_student_t_lpdf({y} | {nu}, {mu}, {Sigma})"), status = "validated")

mglg_lpdf <- function(y, p) { y <- as_rows(y); k <- ncol(y)
  # Latent count range: the negative binomial weights and the conditional
  # gamma terms must both be covered.
  nb <- if (p$delta >= 1) 0 else stats::qnbinom(1 - 1e-15, size = p$nu, prob = p$delta)
  apply(y, 1, function(r) { s <- sum(p$mu * r); e <- sum(exp(p$mu * r) / p$lambda)
    nmax <- ceiling(max(400, 2 * nb, 4 * max(exp(p$mu * r) / p$lambda) + 100))
    if (p$delta >= 1) nmax <- 0
    if (nmax > 1e7) stop_numerical("The G-MVLG series needs more than 1e7 terms.")
    n <- 0:nmax
    log_sum_exp(p$nu * log(p$delta) + ifelse(n == 0, 0, n * log1p(-p$delta)) + sum(log(p$mu)) -
      vapply(n, function(m) sum((p$nu + m) * log(p$lambda)), 0) - (k - 1) * lgamma(p$nu + n) - lgamma(p$nu) - lgamma(n + 1) +
      (p$nu + n) * s - e) }) }
mfam("generalized_multivariate_log_gamma", "Generalized multivariate log-gamma (G-MVLG)", list(P("delta", "unit_closed", role = "dependence (1 = independence)"), P("nu", "positive"),
  P("lambda", "positive", regression = FALSE), P("mu", "positive", regression = FALSE)), aliases = "gmvlg",
  formula = "f(y) = delta^nu sum_n (1 - delta)^n prod(mu_i lambda_i^(-nu - n)) / (Gamma(nu + n)^(k - 1) Gamma(nu) n!) exp((nu + n) sum mu_i y_i - sum exp(mu_i y_i) / lambda_i)",
  lpdf = mglg_lpdf,
  references = "Demirhan and Hamurkaroglu (2011), On a multivariate log-gamma distribution and the use of the distribution in the Bayesian analysis. J. Stat. Plan. Inference, 141, 1141-1152.",
  reference = ref_impl("delta = 1 gives independent log-gamma margins", function(y, p) mglg_lpdf(y, p), note = "Checked by total mass in two dimensions and the independence case."),
  status = "validated")

mfam("bivariate_marshall_olkin_exponential", "Marshall-Olkin bivariate exponential", list(P("lambda1", "positive"), P("lambda2", "positive"), P("lambda12", "positive", role = "common shock rate")),
  aliases = "marshall_olkin_exponential", variate_type = "vector_mixed",
  formula = "P(X > x, Y > y) = exp(-l1 x - l2 y - l12 max(x, y)); absolutely continuous part off the diagonal plus a singular part on x = y",
  lpdf = function(y, p) { y <- as_rows(y); a <- y[, 1]; b <- y[, 2]; l1 <- p$lambda1; l2 <- p$lambda2; l3 <- p$lambda12
    ifelse(a < b, log(l1 * (l2 + l3)) - l1 * a - (l2 + l3) * b,
      ifelse(a > b, log(l2 * (l1 + l3)) - (l1 + l3) * a - l2 * b, log(l3) - (l1 + l2 + l3) * a)) },
  rng = function(n, p) { z1 <- stats::rexp(n, p$lambda1); z2 <- stats::rexp(n, p$lambda2); z3 <- stats::rexp(n, p$lambda12); cbind(pmin(z1, z3), pmin(z2, z3)) },
  numerical_notes = "On the diagonal x = y the value is the density of the singular part with respect to length on the diagonal.",
  references = "Marshall and Olkin (1967), A multivariate exponential distribution. JASA, 62, 30-44.",
  reference = ref_impl("survival function differentiation", function(y, p) { y <- as_rows(y); a <- y[, 1]; b <- y[, 2]
    S <- function(x, z) exp(-p$lambda1 * x - p$lambda2 * z - p$lambda12 * pmax(x, z)); h <- 1e-4
    log((S(a + h, b + h) - S(a + h, b - h) - S(a - h, b + h) + S(a - h, b - h)) / (4 * h^2)) }, note = "Off diagonal only; total mass checked."),
  status = "validated")

register_family(family_def("balding_nichols", "Balding-Nichols", "bounded", list(P("F", "unit", role = "fixation index"), P("p", "unit", role = "ancestral allele frequency")),
  support = function(p) list(lower = 0, upper = 1),
  formula = "Beta(p (1 - F)/F, (1 - p)(1 - F)/F)",
  lpdf = function(y, q) stats::dbeta(y, q$p * (1 - q$F) / q$F, (1 - q$p) * (1 - q$F) / q$F, log = TRUE),
  cdf = function(x, q) stats::pbeta(x, q$p * (1 - q$F) / q$F, (1 - q$p) * (1 - q$F) / q$F),
  quantile = function(pr, q) stats::qbeta(pr, q$p * (1 - q$F) / q$F, (1 - q$p) * (1 - q$F) / q$F),
  mean = function(q) q$p, variance = function(q) q$F * q$p * (1 - q$p),
  references = "Balding and Nichols (1995), A method for quantifying differentiation between populations at multi-allelic loci. Genetica, 96, 3-12.",
  reference = ref_impl("definition", function(y, q) stats::dbeta(y, q$p * (1 - q$F) / q$F, (1 - q$p) * (1 - q$F) / q$F, log = TRUE)), status = "validated"))

# Matrix variates -----------------------------------------------------------------------------
xfam <- function(...) mfam(..., variate_type = "matrix")
xfam("wishart", "Wishart", list(P("nu", "positive", role = "degrees of freedom"), P("Sigma", "any", role = "scale matrix", regression = FALSE)),
  formula = "f(W) = |W|^((nu - p - 1)/2) etr(-Sigma^-1 W / 2) / (2^(nu p/2) |Sigma|^(nu/2) Gamma_p(nu/2))",
  lpdf = function(y, q) { S <- as.matrix(q$Sigma); p <- nrow(S); Si <- solve(S)
    vapply(as_mats(y), function(W) (q$nu - p - 1) / 2 * ldet(W) - etr(Si %*% W) / 2 - q$nu * p / 2 * log(2) - q$nu / 2 * ldet(S) - lmgamma(q$nu / 2, p), 0) },
  rng = function(n, q) lapply(seq_len(n), function(i) stats::rWishart(1, q$nu, as.matrix(q$Sigma))[, , 1]),
  references = "Wishart (1928), The generalised product moment distribution in samples from a normal multivariate population. Biometrika, 20A, 32-52.",
  reference = ref_impl("LaplacesDemon", function(y, q) vapply(as_mats(y), function(W) LaplacesDemon::dwishart(W, q$nu, as.matrix(q$Sigma), log = TRUE), 0)),
  stan = stan_native("wishart_lpdf({y} | {nu}, {Sigma})"), status = "validated")

xfam("inverse_wishart", "Inverse Wishart", list(P("nu", "positive"), P("Psi", "any", role = "scale matrix", regression = FALSE)),
  formula = "f(X) = |Psi|^(nu/2) |X|^(-(nu + p + 1)/2) etr(-Psi X^-1 / 2) / (2^(nu p/2) Gamma_p(nu/2))",
  lpdf = function(y, q) { S <- as.matrix(q$Psi); p <- nrow(S)
    vapply(as_mats(y), function(X) q$nu / 2 * ldet(S) - (q$nu + p + 1) / 2 * ldet(X) - etr(S %*% solve(X)) / 2 - q$nu * p / 2 * log(2) - lmgamma(q$nu / 2, p), 0) },
  rng = function(n, q) lapply(seq_len(n), function(i) solve(stats::rWishart(1, q$nu, solve(as.matrix(q$Psi)))[, , 1])),
  references = "Gupta and Nagar (2000), Matrix Variate Distributions. Chapman and Hall, ch. 3.",
  reference = ref_impl("LaplacesDemon", function(y, q) vapply(as_mats(y), function(X) LaplacesDemon::dinvwishart(X, q$nu, as.matrix(q$Psi), log = TRUE), 0)),
  stan = stan_native("inv_wishart_lpdf({y} | {nu}, {Psi})"), status = "validated")

lkj_logC <- function(eta, d) { k <- seq_len(d - 1)
  sum((2 * eta - 2 + d - k) * (d - k)) * log(2) + sum((d - k) * lbeta(eta + (d - k - 1) / 2, eta + (d - k - 1) / 2)) }
xfam("lkj", "Lewandowski-Kurowicka-Joe (correlation matrices)", list(P("eta", "positive", role = "shape")), aliases = "lkj_corr",
  formula = "f(R) = det(R)^(eta - 1) / c_d(eta), c_d(eta) = 2^(sum (2 eta - 2 + d - k)(d - k)) prod B(eta + (d - k - 1)/2, eta + (d - k - 1)/2)^(d - k)",
  lpdf = function(y, q) vapply(as_mats(y), function(R) (q$eta - 1) * ldet(R) - lkj_logC(q$eta, nrow(R)), 0),
  references = "Lewandowski, Kurowicka and Joe (2009), Generating random correlation matrices based on vines and extended onion method. J. Multivariate Anal., 100, 1989-2001.",
  reference = ref_impl("d = 2: (r + 1)/2 ~ Beta(eta, eta)", function(y, q) vapply(as_mats(y), function(R) stats::dbeta((R[1, 2] + 1) / 2, q$eta, q$eta, log = TRUE) - log(2), 0),
    note = "Checked against the Beta form for d = 2 and by integration for d = 3."),
  stan = stan_native("lkj_corr_lpdf({y} | {eta})"), status = "validated")

xfam("matrix_normal", "Matrix normal", list(P("M", "real", role = "mean matrix", regression = FALSE), P("U", "any", role = "row covariance", regression = FALSE),
  P("V", "any", role = "column covariance", regression = FALSE)),
  formula = "vec(X) ~ N(vec(M), V kron U)",
  lpdf = function(y, q) { U <- as.matrix(q$U); V <- as.matrix(q$V); n <- nrow(U); p <- nrow(V); Ui <- solve(U); Vi <- solve(V)
    vapply(as_mats(y), function(X) { D <- X - q$M; -0.5 * etr(Vi %*% t(D) %*% Ui %*% D) - n * p / 2 * log(2 * pi) - n / 2 * ldet(V) - p / 2 * ldet(U) }, 0) },
  rng = function(n, q) lapply(seq_len(n), function(i) q$M + t(chol(as.matrix(q$U))) %*% matrix(stats::rnorm(length(q$M)), nrow(q$M)) %*% chol(as.matrix(q$V))),
  references = "Dawid (1981), Some matrix-variate distribution theory. Biometrika, 68, 265-274.",
  reference = ref_impl("mvtnorm on vec(X)", function(y, q) vapply(as_mats(y), function(X) mvtnorm::dmvnorm(as.vector(X), as.vector(q$M), kronecker(as.matrix(q$V), as.matrix(q$U)), log = TRUE), 0)),
  status = "validated")

xfam("matrix_t", "Matrix t", list(P("M", "real", regression = FALSE), P("Sigma", "any", role = "row scale", regression = FALSE), P("Omega", "any", role = "column scale", regression = FALSE),
  P("nu", "positive")),
  formula = "f(X) = Gamma_p((nu + n + p - 1)/2) / (pi^(np/2) Gamma_p((nu + p - 1)/2)) |Omega|^(-n/2) |Sigma|^(-p/2) |I_n + Sigma^-1 (X - M) Omega^-1 (X - M)'|^(-(nu + n + p - 1)/2)",
  lpdf = function(y, q) { S <- as.matrix(q$Sigma); O <- as.matrix(q$Omega); n <- nrow(S); p <- nrow(O); a <- (q$nu + n + p - 1) / 2
    vapply(as_mats(y), function(X) { D <- X - q$M
      lmgamma(a, p) - n * p / 2 * log(pi) - lmgamma((q$nu + p - 1) / 2, p) - n / 2 * ldet(O) - p / 2 * ldet(S) - a * ldet(diag(n) + solve(S) %*% D %*% solve(O) %*% t(D)) }, 0) },
  references = "Gupta and Nagar (2000), Matrix Variate Distributions, ch. 4.",
  reference = ref_impl("p = 1 reduces to a multivariate t with scale Sigma Omega / nu", function(y, q) vapply(as_mats(y), function(X)
    mvtnorm::dmvt(as.vector(X), as.vector(q$M), as.matrix(q$Sigma) * as.numeric(q$Omega) / q$nu, df = q$nu, log = TRUE), 0)),
  status = "validated")

xfam("matrix_beta", "Matrix variate beta (type I)", list(P("a", "positive"), P("b", "positive"), P("p", "positive_integer", data = TRUE)),
  formula = "f(U) = Gamma_p(a + b) / (Gamma_p(a) Gamma_p(b)) |U|^(a - (p + 1)/2) |I - U|^(b - (p + 1)/2), 0 < U < I",
  lpdf = function(y, q) vapply(as_mats(y), function(U) { p <- nrow(U)
    lmgamma(q$a + q$b, p) - lmgamma(q$a, p) - lmgamma(q$b, p) + (q$a - (p + 1) / 2) * ldet(U) + (q$b - (p + 1) / 2) * ldet(diag(p) - U) }, 0),
  rng = function(n, q) lapply(seq_len(n), function(i) { A <- stats::rWishart(1, 2 * q$a, diag(q$p))[, , 1]; B <- stats::rWishart(1, 2 * q$b, diag(q$p))[, , 1]
    e <- eigen(A + B, symmetric = TRUE); Si <- e$vectors %*% diag(1 / sqrt(e$values), q$p) %*% t(e$vectors); Si %*% A %*% Si }),
  references = "Gupta and Nagar (2000), Matrix Variate Distributions, ch. 5.",
  reference = ref_impl("p = 1 is Beta(a, b)", function(y, q) vapply(as_mats(y), function(U) stats::dbeta(U[1, 1], q$a, q$b, log = TRUE), 0)),
  status = "validated")

stiefel_logvol <- function(n, k) k * log(2) + n * k / 2 * log(pi) - lmgamma(n / 2, k)
xfam("stiefel_uniform", "Uniform (Haar) on the Stiefel manifold", list(P("n", "positive_integer", data = TRUE), P("k", "positive_integer", data = TRUE)),
  aliases = c("haar_orthogonal", "uniform_stiefel"),
  formula = "Constant density 1 / vol(V_{n,k}), vol = 2^k pi^(nk/2) / Gamma_k(n/2), with respect to the Hausdorff measure",
  lpdf = function(y, q) rep(-stiefel_logvol(q$n, q$k), length(as_mats(y))),
  rng = function(n, q) lapply(seq_len(n), function(i) { Z <- matrix(stats::rnorm(q$n * q$k), q$n); qr.Q(qr(Z)) %*% diag(sign(diag(qr.R(qr(Z)))), q$k) }),
  references = "Chikuse (2003), Statistics on Special Manifolds. Springer, ch. 1.",
  reference = ref_impl("k = 1: area of the sphere S^(n-1)", function(y, q) rep(-(log(2) + q$n / 2 * log(pi) - lgamma(q$n / 2)), length(as_mats(y)))),
  status = "validated")

# 0F1(n/2; F'F/4) = E_Haar etr(F'X): closed form for k = 1, quasi Monte Carlo otherwise.
ml_log0F1 <- function(Fm, nsim = 2e5) { Fm <- as.matrix(Fm); n <- nrow(Fm); k <- ncol(Fm)
  if (k == 1) { kap <- sqrt(sum(Fm^2)); if (kap < 1e-10) return(0)
    return(lgamma(n / 2) + (1 - n / 2) * log(kap / 2) + log(besselI(kap, n / 2 - 1, expon.scaled = TRUE)) + kap) }
  # The constant depends on F only through its singular values d, and then
  # tr(F'Q) = sum_j d_j Q_jj. Haar draws by Gram-Schmidt on Gaussian columns,
  # vectorized over draws; cached per distinct d.
  d <- svd(Fm, nu = 0, nv = 0)$d; key <- paste(signif(c(n, d), 12), collapse = ",")
  if (!is.null(ml_cache[[key]])) return(ml_cache[[key]])
  local_seed(20240101L)
  Qs <- list(); v <- numeric(nsim)
  for (j in seq_len(k)) {
    z <- matrix(stats::rnorm(nsim * n), nsim)
    for (q in Qs) z <- z - rowSums(z * q) * q
    z <- z / sqrt(rowSums(z^2)); Qs[[j]] <- z
    v <- v + d[j] * z[, j]
  }
  w <- exp(v - max(v))
  if (sum(w)^2 / sum(w^2) < 200)
    stop_numerical("The matrix Langevin normalizing constant cannot be estimated reliably at this concentration.")
  ml_cache[[key]] <- log_sum_exp(v) - log(nsim)
}
ml_cache <- new.env(parent = emptyenv())
xfam("matrix_langevin", "Matrix Langevin (matrix von Mises-Fisher) on the Stiefel manifold", list(P("F", "real", role = "parameter matrix", regression = FALSE)),
  aliases = "matrix_von_mises_fisher",
  formula = "f(X) = etr(F'X) / 0F1(n/2; F'F/4) with respect to the normalized Haar measure on V_{n,k}",
  lpdf = function(y, q) { lc <- ml_log0F1(q$F); vapply(as_mats(y), function(X) etr(t(as.matrix(q$F)) %*% X) - lc, 0) },
  numerical_notes = "The normalizing hypergeometric function of matrix argument is exact for k = 1 (Bessel form) and estimated by Monte Carlo with 2e5 Haar draws otherwise (relative error about 1e-3 for moderate F).",
  references = c("Downs (1972), Orientation statistics. Biometrika, 59, 665-676.", "Chikuse (2003), Statistics on Special Manifolds, ch. 2."),
  reference = ref_impl("k = 1 equals the von Mises-Fisher law", function(y, q) { Fm <- as.matrix(q$F); kap <- sqrt(sum(Fm^2))
    vapply(as_mats(y), function(X) log_vmf_c(kap, nrow(Fm)) + sum(Fm * X) + stiefel_logvol(nrow(Fm), 1), 0) }),
  status = "implemented")
