# Portable October 13 replication. Run from the repository root; see docs/SIMULATIONS.txt.
# Run: Rscript check_information_per_unit.R [candidate_state_space_em.R] [output_directory] [application_fit.rds]
# Only sources function definitions and reads application data. No original writes.
args <- commandArgs(trailingOnly = TRUE)
base <- getwd()
code <- if (length(args)) args[1] else file.path(base, 'R/state_space_em.R')
out <- if (length(args) >= 2) args[2] else file.path(base, 'results/simulation/validation/check_unit_information')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
source(code)
stopifnot(requireNamespace('numDeriv', quietly = TRUE))

# Independent dense Gaussian distribution, without calling a Kalman routine.
latent_moments <- function(th, x, Tt) {
  p <- length(x); a <- th[p + 1]; P <- th[p + 2]; q <- th[p + 3]
  ph <- tanh(sum(x * th[seq_len(p)]))
  v <- numeric(Tt); v[1] <- P
  for (t in 2:Tt) v[t] <- ph^2 * v[t - 1] + q
  V <- outer(seq_len(Tt), seq_len(Tt), function(j, k) ph^abs(j - k) * v[pmin(j, k)])
  list(mean = a * ph^(0:(Tt - 1)), V = V)
}
unit_ll <- function(th, yy, xx, mui, r) {
  z <- latent_moments(th, xx, length(yy)); S <- z$V + diag(r, length(yy))
  U <- chol(S); e <- yy - mui - z$mean
  -.5 * (length(yy) * log(2 * pi) + 2 * sum(log(diag(U))) + sum(backsolve(U, e, transpose = TRUE)^2))
}
frozen_posterior <- function(th, yy, xx, mui, r) {
  z <- latent_moments(th, xx, length(yy)); S <- z$V + diag(r, length(yy))
  K <- z$V %*% solve(S)
  list(m = as.vector(z$mean + K %*% (yy - mui - z$mean)), V = z$V - K %*% z$V)
}
unit_Q <- function(th, post, xx) {
  p <- length(xx); a <- th[p + 1]; P <- th[p + 2]; q <- th[p + 3]
  ph <- tanh(sum(xx * th[seq_len(p)])); Tt <- length(post$m)
  m <- post$m; V <- post$V
  S <- V[1, 1] + (m[1] - a)^2
  R <- sum((m[2:Tt] - ph * m[1:(Tt - 1)])^2 + diag(V)[2:Tt] +
             ph^2 * diag(V)[1:(Tt - 1)] - 2 * ph * V[cbind(2:Tt, 1:(Tt - 1))])
  -.5 * log(P) - S / (2 * P) - (Tt - 1) / 2 * log(q) - R / (2 * q)
}
make_fit <- function(th, yy, XX, mu, r) {
  p <- ncol(XX)
  list(beta = th[seq_len(p)], a1 = th[p + 1], P1 = th[p + 2], sigma_eta2 = th[p + 3], sigma_eps2 = r,
       kf = kalman_run(yy, tanh(as.vector(XX %*% th[seq_len(p)])), mu, th[p + 1], th[p + 2], th[p + 3], r))
}
scaled_error <- function(x, y) max(abs(x - y) / (1 + abs(y)))
rows <- list()
record <- function(case, check, error, tolerance) {
  rows[[length(rows) + 1]] <<- data.frame(case, check, error, tolerance, pass = is.finite(error) && error < tolerance)
}

set.seed(29092026)
for (k in 1:5) {
  N <- c(7, 5, 8, 1, 6)[k]; Tt <- c(8, 2, 14, 5, 6)[k]; p <- c(3, 2, 3, 1, 1)[k]
  XX <- cbind(1, matrix(rnorm(N * max(0, p - 1)), N, max(0, p - 1)))
  if (k != 5) colnames(XX) <- paste0('b', seq_len(p))
  beta <- c(.25, .8, -.65)[seq_len(p)]
  th <- c(beta, .6, 1.3, .8); r <- c(.45, 1.2, 2, .7, .5)[k]; mu <- seq(-.2, .2, length.out = N)
  yy <- matrix(rnorm(N * Tt, .2, 1.5), N, Tt) # Deliberately not an optimum.
  fit <- make_fit(th, yy, XX, mu, r)
  Ii <- compute_info_per_unit_fixed_eps2(yy, XX, mu, fit)
  expected <- array(0, dim(Ii)); score_num <- matrix(0, N, length(th))
  for (i in seq_len(N)) {
    ll <- function(z) unit_ll(z, yy[i, ], XX[i, ], mu[i], r)
    expected[i, , ] <- -numDeriv::hessian(ll, th)
    score_num[i, ] <- numDeriv::grad(ll, th)
  }
  case <- paste0('off_optimum_', k)
  record(case, 'unit_observed_hessians', scaled_error(Ii, expected), 2e-5)
  record(case, 'raw_unit_symmetry', max(abs(Ii - aperm(Ii, c(1, 3, 2)))), 1e-9)
  groups <- rep(seq_len(max(1, floor(N / 2))), length.out = N)
  for (g in unique(groups)) {
    ix <- which(groups == g)
    gg <- function(z) sum(vapply(ix, function(i) unit_ll(z, yy[i, ], XX[i, ], mu[i], r), numeric(1)))
    record(case, paste0('cluster_', g, '_hessian'),
           scaled_error(apply(Ii[ix, , , drop = FALSE], c(2, 3), sum), -numDeriv::hessian(gg, th)), 2e-5)
  }
  # SEs themselves are not meaningful at an arbitrary parameter; suppress only
  # square-root warnings from indefinite inverse Hessians while examining D/I.
  agg <- suppressWarnings(compute_SE_regression_mv_fixed_eps2(yy, XX, mu, fit))
  record(case, 'sum_equals_aggregate', max(abs(apply(Ii, c(2, 3), sum) - agg$I_obs)), 1e-9)
  posts <- lapply(seq_len(N), function(i) frozen_posterior(th, yy[i, ], XX[i, ], mu[i], r))
  Q <- function(z) sum(vapply(seq_len(N), function(i) unit_Q(z, posts[[i]], XX[i, ]), numeric(1)))
  record(case, 'frozen_Q_direct_hessian', scaled_error(agg$D, -numDeriv::hessian(Q, th)), 2e-5)
  cr <- suppressWarnings(compute_SE_regression_mv_cluster_robust(yy, XX, mu, fit, groups))
  record(case, 'unit_scores', scaled_error(cr$score, score_num), 1e-6)
}

# Fit a reproducible, independent synthetic dataset, then check the optimum.
N <- 80; Tt <- 12; XX <- cbind(Intercept = 1, x = rnorm(N))
th <- c(.8, .3, .7, 1.4, .65); r <- .4; mu <- 0
YY <- t(vapply(seq_len(N), function(i) {
  zz <- latent_moments(th, XX[i, ], Tt)
  zz$mean + as.vector(t(chol(zz$V + diag(r, Tt))) %*% rnorm(Tt))
}, numeric(Tt)))
fit <- em_regression_mv(YY, XX, mu, init = list(beta = c(.5, .1), a1 = .5, P1 = 1,
                         sigma_eta2 = 1, sigma_eps2 = r), fix_sigma_eps2 = r, thres = 1e-9, max_iter = 3000)
stopifnot(fit$converged)
theta <- c(fit$beta, fit$a1, fit$P1, fit$sigma_eta2)
Ii <- compute_info_per_unit_fixed_eps2(YY, XX, mu, fit)
agg <- compute_SE_regression_mv_fixed_eps2(YY, XX, mu, fit)
ll <- function(z) sum(vapply(seq_len(N), function(i) unit_ll(z, YY[i, ], XX[i, ], mu, r), numeric(1)))
record('fitted_synthetic', 'total_dense_hessian', scaled_error(agg$I_obs, -numDeriv::hessian(ll, theta)), 2e-5)
for (i in c(1, 17, 61)) record('fitted_synthetic', paste0('unit_', i, '_hessian'),
  scaled_error(Ii[i, , ], -numDeriv::hessian(function(z) unit_ll(z, YY[i, ], XX[i, ], mu, r), theta)), 2e-5)
record('fitted_synthetic', 'sum_equals_aggregate', max(abs(apply(Ii, c(2, 3), sum) - agg$I_obs)), 1e-9)
cat('Synthetic fit converged in', fit$iter, 'iterations.\n')

# Application is read only. Supply the third argument to audit another saved
# fit on the same application cohort.
app_path <- if (length(args) >= 3) args[3] else file.path(base, 'results/application/FINAL_fitted_models_final.rds')
f <- readRDS(app_path); YY <- f$Y; XX <- f$Xmat; fit <- f$fit; p <- ncol(XX)
theta <- c(fit$beta, fit$a1, fit$P1, fit$sigma_eta2); r <- fit$sigma_eps2
Ii <- compute_info_per_unit_fixed_eps2(YY, XX, 0, fit)
agg <- compute_SE_regression_mv_fixed_eps2(YY, XX, 0, fit)
record('application', 'sum_equals_aggregate', max(abs(apply(Ii, c(2, 3), sum) - agg$I_obs)), 1e-8)
record('application', 'raw_unit_symmetry', max(abs(Ii - aperm(Ii, c(1, 3, 2)))), 1e-8)
# Scale parameter coordinates for numerical differentiation only; beta unit
# scale and nuisance scales of current estimates avoid bad finite differences.
scale <- c(rep(1, p), max(abs(fit$a1), 1), fit$P1, fit$sigma_eta2)
Dscale <- diag(scale)
for (sid in unique(f$meta$subject)[c(1, 7, 33)]) {
  ix <- which(f$meta$subject == sid)
  ll_scaled <- function(z) sum(vapply(ix, function(i) unit_ll(z * scale, YY[i, ], XX[i, ], 0, r), numeric(1)))
  expected <- -numDeriv::hessian(ll_scaled, theta / scale)
  analytic <- Dscale %*% apply(Ii[ix, , , drop = FALSE], c(2, 3), sum) %*% Dscale
  record('application', paste0(sid, '_scaled_cluster_hessian'), scaled_error(analytic, expected), 5e-4)
}
cr0 <- compute_SE_regression_mv_cluster_robust(YY, XX, 0, fit, f$meta$subject)
md <- compute_SE_regression_mv_cluster_MD(YY, XX, 0, fit, f$meta$subject)
cr2 <- compute_SE_regression_mv_cluster_CR2(YY, XX, 0, fit, f$meta$subject)
cat('Application source:', app_path, '\nN=', nrow(YY), 'T=', ncol(YY), 'converged=', fit$converged, '\n')
cat('Application total score:', format(colSums(cr0$score), digits = 5), '\n')
results <- do.call(rbind, rows)
print(results, row.names = FALSE)
write.csv(results, file.path(out, 'validation_results.csv'), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out, 'sessionInfo.txt'))
stopifnot(all(results$pass))
cat('PASS:', nrow(results), 'checks.\n')
