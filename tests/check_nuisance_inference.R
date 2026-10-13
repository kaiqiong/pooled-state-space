# Run from repository root: Rscript tests/check_nuisance_inference.R [R_source] [output_directory] [application_fit.rds]
args <- commandArgs(trailingOnly = TRUE)
base <- normalizePath('.', mustWork=TRUE)
code <- if (length(args)) args[1] else file.path(base, 'R/state_space_em.R')
out <- if (length(args) >= 2) args[2] else file.path(base, 'results/validation', 'check_nuisance_inference')
app_path <- if (length(args) >= 3) args[3] else file.path(base, 'results/application/FINAL_fitted_models_final.rds')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
source(code)
stopifnot(requireNamespace('numDeriv', quietly = TRUE))
options(width = 150)
checks <- list()
err <- function(a, b) max(abs(a - b) / (1 + abs(b)))
check <- function(case, what, value, tolerance = 1e-9) {
  checks[[length(checks) + 1]] <<- data.frame(case, check = what, error = value, tolerance,
    pass = is.finite(value) && value < tolerance)
}
rejects <- function(expr) as.numeric(!inherits(try(force(expr), silent = TRUE), 'try-error'))

set.seed(29092027)
for (scenario in 1:4) {
  d <- c(7, 5, 3, 6)[scenario]; p <- c(3, 1, 3, 2)[scenario]; N <- 17
  pn <- paste0('theta', seq_len(d)); bn <- pn[seq_len(p)]
  M <- matrix(rnorm(d * d), d, d)
  H <- crossprod(M) + diag(.8, d); dimnames(H) <- list(pn, pn)
  score <- matrix(rnorm(N * d), N, d, dimnames = list(NULL, pn))
  score <- scale(score, center = TRUE, scale = FALSE)
  cl <- rep(c('C', 'A', 'B', 'Z'), c(1, 3, 5, 8))
  weights <- runif(N); weights <- weights / sum(weights)
  Ii <- array(0, c(N, d, d), dimnames = list(NULL, pn, pn))
  for (i in seq_len(N)) Ii[i, , ] <- H * weights[i]
  z <- reduce_regression_cluster_scores(H, score, cl, bn, Ii)
  inv <- solve(H); S <- rowsum(score, cl, reorder = FALSE); full <- inv %*% crossprod(S) %*% t(inv)
  label <- paste0('algebra_', scenario)
  check(label, 'full_sandwich', err(z$vcov, full[bn, bn, drop = FALSE]))
  check(label, 'full_model_block', err(z$vcov_model_based, inv[bn, bn, drop = FALSE]))
  check(label, 'full_influence', err(z$influence, (S %*% t(inv))[, bn, drop = FALSE]))
  check(label, 'projected_information_sum', err(apply(z$projected_info_per_unit, c(2, 3), sum), z$I_eff))
  check(label, 'cluster_projection_sum', err(apply(z$projected_info_per_cluster, c(2, 3), sum), z$I_eff))
  perm <- sample(pn)
  zz <- reduce_regression_cluster_scores(H[perm, perm], score[, perm, drop = FALSE], cl, bn)
  check(label, 'coordinate_permutation', err(zz$vcov, z$vcov))
  if (p < d) {
    nu <- pn[-seq_len(p)]; J <- diag(d)
    L <- matrix(rnorm(length(nu)^2), length(nu), length(nu)) + diag(3, length(nu))
    J[-seq_len(p), -seq_len(p)] <- L
    Hz <- t(J) %*% H %*% J; dimnames(Hz) <- list(pn, pn)
    sz <- score %*% J; colnames(sz) <- pn
    zr <- reduce_regression_cluster_scores(Hz, sz, cl, bn)
    check(label, 'nuisance_reparameterization_covariance', err(zr$vcov, z$vcov))
    check(label, 'nuisance_reparameterization_score', err(zr$score, z$score))
    naive <- solve(H[bn, bn, drop = FALSE]) %*% crossprod(S[, bn, drop = FALSE]) %*% solve(H[bn, bn, drop = FALSE])
    check(label, 'detect_nontrivial_nuisance_effect', as.numeric(max(abs(naive - z$vcov)) < 1e-5), .5)
  }
}
check('input', 'reject_missing_cluster', rejects(reduce_regression_cluster_scores(H, score, replace(cl, 1, NA), bn)), .5)
check('input', 'reject_single_cluster', rejects(reduce_regression_cluster_scores(H, score, rep('a', N), bn)), .5)
Hb <- H; Hb[1, 2] <- Hb[1, 2] + 1
check('input', 'reject_asymmetry', rejects(reduce_regression_cluster_scores(Hb, score, cl, bn)), .5)
Hb <- H; Hb[1, 1] <- -1
check('input', 'reject_indefinite_information', rejects(reduce_regression_cluster_scores(Hb, score, cl, bn)), .5)
check('input', 'reject_wrong_score_coordinates', rejects(reduce_regression_cluster_scores(H, score[, rev(pn)], cl, bn)), .5)

# Dense Gaussian likelihood independent of Kalman recursions.
latent <- function(th, x, Tt) {
  p <- length(x); ph <- tanh(sum(x * th[seq_len(p)])); a <- th[p + 1]; P <- th[p + 2]; q <- th[p + 3]
  v <- numeric(Tt); v[1] <- P
  for (t in 2:Tt) v[t] <- ph^2 * v[t - 1] + q
  list(mean = a * ph^(0:(Tt - 1)), V = outer(seq_len(Tt), seq_len(Tt), function(j, k) ph^abs(j-k) * v[pmin(j,k)]))
}
unit_ll <- function(th, yy, xx, r) {
  z <- latent(th, xx, length(yy)); U <- chol(z$V + diag(r, length(yy)))
  e <- yy - z$mean
  -.5 * (length(yy) * log(2*pi) + 2 * sum(log(diag(U))) + sum(backsolve(U, e, transpose = TRUE)^2))
}
N <- 40; Tt <- 8; XX <- cbind(Intercept = 1, x = rnorm(N)); cl <- rep(1:8, c(2,3,4,5,5,6,7,8))
th <- c(.7, -.2, 1.6, 1.3, .6); r <- .5; p <- ncol(XX)
y <- t(vapply(seq_len(N), function(i) {
  z <- latent(th, XX[i, ], Tt)
  z$mean + as.vector(t(chol(z$V + diag(r, Tt))) %*% rnorm(Tt))
}, numeric(Tt)))
fit <- em_regression_mv(y, XX, 0, init = list(beta = c(.5, -.1), a1 = 1, P1 = 1, sigma_eta2 = .8, sigma_eps2 = r),
                       fix_sigma_eps2 = r, thres = 1e-9, max_iter = 3000)
stopifnot(fit$converged, !fit$at_cap)
# Refine the optimum independently, in log-variance coordinates.
to_raw <- function(z) c(z[1:(p+1)], exp(z[(p+2):(p+3)]))
to_log <- function(z) c(z[1:(p+1)], log(z[(p+2):(p+3)]))
negll <- function(z, w = rep(1, N)) -sum(vapply(seq_len(N), function(i) w[i] * unit_ll(to_raw(z), y[i, ], XX[i, ], r), numeric(1)))
z0 <- to_log(c(fit$beta, fit$a1, fit$P1, fit$sigma_eta2))
opt <- optim(z0, negll, method = 'BFGS', control = list(reltol = 1e-13, maxit = 2000, ndeps = rep(1e-5, length(z0))))
stopifnot(opt$convergence == 0)
th <- to_raw(opt$par); fit$beta <- th[1:p]; fit$a1 <- th[p+1]; fit$P1 <- th[p+2]; fit$sigma_eta2 <- th[p+3]
fit$kf <- kalman_run(y, tanh(as.vector(XX %*% fit$beta)), 0, fit$a1, fit$P1, fit$sigma_eta2, r)
z <- compute_SE_regression_mv_cluster_beta(y, XX, 0, fit, cl, return_unit_info = TRUE)
ll <- function(a) sum(vapply(seq_len(N), function(i) unit_ll(a, y[i, ], XX[i, ], r), numeric(1)))
Hnum <- -numDeriv::hessian(ll, th); pn <- c(colnames(XX), 'a1', 'P1', 'sigma_eta2'); dimnames(Hnum) <- list(pn, pn)
Snum <- t(vapply(seq_len(N), function(i) numDeriv::grad(function(a) unit_ll(a, y[i, ], XX[i, ], r), th), numeric(length(th))))
colnames(Snum) <- pn
ref <- reduce_regression_cluster_scores(Hnum, Snum, cl, colnames(XX))
check('synthetic', 'independent_dense_covariance', err(z$vcov, ref$vcov), 2e-6)
check('synthetic', 'independent_dense_adjusted_scores', err(z$score, ref$score), 2e-6)
check('synthetic', 'full_CR0_equivalence', z$diagnostics$existing_full_CR0_error)
check('synthetic', 'projected_sum', z$diagnostics$projected_information_sum_error, 1e-7)
# Central differences from weighted JOINT refits, not fixed-nuisance refits.
eps <- .01
for (g in c(1, 4, 8)) {
  opts <- lapply(c(-1, 1), function(sign) {
    w <- 1 + sign * eps * (cl == g)
    optim(opt$par, negll, w = w, method = 'BFGS', control = list(reltol = 1e-13, maxit = 2000, ndeps = rep(1e-5, length(z0))))
  })
  stopifnot(all(vapply(opts, function(o) o$convergence == 0, logical(1))))
  fd <- (opts[[2]]$par[1:p] - opts[[1]]$par[1:p]) / (2 * eps)
  check('synthetic', paste0('weighted_joint_refit_', g), err(z$influence[as.character(g), ], fd), 2e-4)
  cat('Weighted-refit derivative, cluster', g, ': analytic', z$influence[as.character(g), ], 'numeric', fd, '\n')
}
# Profile curvature from nuisance-only reoptimization with beta fixed.
profile <- function(b) {
  obj <- function(nu) negll(c(b, nu))
  op <- optim(opt$par[-seq_len(p)], obj, method = 'BFGS',
              control = list(reltol = 1e-13, maxit = 2000, ndeps = rep(1e-5, 3)))
  stopifnot(op$convergence == 0)
  op$value
}
h <- .002; b <- opt$par[1:p]; f0 <- profile(b); Hp <- matrix(0, p, p)
for (j in 1:p) for (k in 1:p) {
  ej <- diag(p)[,j] * h; ek <- diag(p)[,k] * h
  if (j == k) Hp[j,k] <- (profile(b+ej) - 2*f0 + profile(b-ej))/h^2
  else Hp[j,k] <- (profile(b+ej+ek)-profile(b+ej-ek)-profile(b-ej+ek)+profile(b-ej-ek))/(4*h^2)
}
check('synthetic', 'profile_curvature', err(z$I_eff, Hp), 1e-3)
# Contrast semantics and wrapper failure cases.
ci <- wald_regression_cluster_contrast(z, c(x=1))
check('contrast', 'named_estimate', abs(ci$estimate - fit$beta[2]))
check('contrast', 'variance', abs(ci$se^2 - z$vcov['x','x']))
ct <- wald_regression_cluster_contrast(z, c(x=1), reference='t', df=7)
check('contrast', 'explicit_t_reference', abs(ct$upper - (fit$beta[2] + qt(.975,7)*ct$se)))
check('input', 'reject_unspecified_t_df', rejects(wald_regression_cluster_contrast(z, c(x=1), reference='t')), .5)
check('input', 'reject_unknown_contrast', rejects(wald_regression_cluster_contrast(z, c(bad=1))), .5)
ff <- fit; ff$sigma_eps2_fixed <- FALSE
check('input', 'reject_free_variance', rejects(compute_SE_regression_mv_cluster_beta(y,XX,0,ff,cl)), .5)
ff <- fit; ff$sigma_eps2_fixed <- NULL
check('input', 'reject_unspecified_variance_status', rejects(compute_SE_regression_mv_cluster_beta(y,XX,0,ff,cl)), .5)
ff <- fit; ff$converged <- NULL
check('input', 'reject_unspecified_convergence', rejects(compute_SE_regression_mv_cluster_beta(y,XX,0,ff,cl)), .5)
ff <- fit; ff$converged <- FALSE
check('input', 'reject_nonconvergence', rejects(compute_SE_regression_mv_cluster_beta(y,XX,0,ff,cl)), .5)
ff <- fit; ff$at_cap <- TRUE
check('input', 'reject_active_cap', rejects(compute_SE_regression_mv_cluster_beta(y,XX,0,ff,cl)), .5)
ybad <- y; ybad[1,1] <- NA
check('input', 'reject_missing_observation', rejects(compute_SE_regression_mv_cluster_beta(ybad,XX,0,fit,cl)), .5)
Xu <- unname(XX); zu <- compute_SE_regression_mv_cluster_beta(y, Xu, 0, fit, cl)
check('synthetic', 'unnamed_covariates', err(zu$vcov, z$vcov))

if (length(args) >= 3) {
# Corrected application, read only.
f <- readRDS(app_path); y <- f$Y; XX <- f$Xmat; fit <- f$fit; cl <- f$meta$subject
app <- compute_SE_regression_mv_cluster_beta(y, XX, 0, fit, cl, return_unit_info=TRUE)
full <- compute_SE_regression_mv_cluster_robust(y, XX, 0, fit, cl)
legacy <- compute_SE_regression_mv_cluster_CR2(y, XX, 0, fit, cl)
bn <- colnames(XX); H <- solve(full$bread); S <- rowsum(full$score, cl, reorder=FALSE)
Ib_inv <- solve(H[bn,bn]); naive_v <- Ib_inv %*% crossprod(S[,bn]) %*% t(Ib_inv)
comparison <- data.frame(parameter=bn, estimate=as.vector(fit$beta),
  model_nuisance_fixed=sqrt(diag(Ib_inv)), model_nuisance_adjusted=app$se_model_based,
  CR0_nuisance_fixed=sqrt(diag(naive_v)), CR0_nuisance_adjusted=app$se,
  full_joint_CR0=full$se[bn], legacy_beta_only_CR2=legacy$se)
print(comparison, digits=9)
write.csv(comparison, file.path(out,'application_comparison.csv'), row.names=FALSE)
check('application','full_sandwich',app$diagnostics$full_sandwich_error,1e-8)
check('application','existing_full_CR0',app$diagnostics$existing_full_CR0_error,1e-8)
check('application','model_block',app$diagnostics$full_model_block_error,1e-8)
check('application','influence_block',app$diagnostics$full_influence_error,1e-8)
check('application','projected_information_sum',app$diagnostics$projected_information_sum_error,1e-7)
check('application','legacy_flag',as.numeric(legacy$nuisance_adjusted || legacy$validated_CR2),.5)
unit_eigen <- vapply(seq_len(dim(app$projected_info_per_unit)[1]),function(i) min(eigen(app$projected_info_per_unit[i,,],symmetric=TRUE,only.values=TRUE)$values),numeric(1))
cluster_eigen <- vapply(seq_len(dim(app$projected_info_per_cluster)[1]),function(i) min(eigen(app$projected_info_per_cluster[i,,],symmetric=TRUE,only.values=TRUE)$values),numeric(1))
influence <- data.frame(subject=rownames(app$influence),app$influence,check.names=FALSE)
write.csv(influence,file.path(out,'application_cluster_influence.csv'),row.names=FALSE)
contrast <- rbind(wald_regression_cluster_contrast(app,c(Diabetic=1)),
                  wald_regression_cluster_contrast(app,c(Diabetic=1),reference='t',df=app$n_clusters-1))
print(contrast,digits=9); write.csv(contrast,file.path(out,'diabetes_reference_comparison.csv'),row.names=FALSE)
cat('Application:',app_path,'N',nrow(y),'G',app$n_clusters,'\n')
cat('Minimum projected information eigenvalues: units',min(unit_eigen),'clusters',min(cluster_eigen),'\n')
cat('Negative projected eigenvalues (< -1e-8):',sum(unit_eigen < -1e-8),'units;',sum(cluster_eigen < -1e-8),'clusters\n')
print(app$diagnostics)

}
results <- do.call(rbind,checks); print(results,row.names=FALSE)
write.csv(results,file.path(out,'validation_results.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
stopifnot(all(results$pass))
cat('PASS:',nrow(results),'checks.\n')
