# Portable October 13 replication. Run from the repository root; see docs/SIMULATIONS.txt.
## Simulation: phi_i = tanh(beta0 + beta_age*Age_i + beta_sex*Sex_i + beta_bio*Biomarker_i)
## A more realistic covariate structure than the single-binary-covariate case:
##   - Age and Biomarker are continuous and correlated (age-related biomarker drift)
##   - Sex is binary but drawn at 80% female, matching the hEDS*Omics target population
##     (hEDS has strong female predominance) rather than a 50/50 split
##   - N = 150, T = 30 chosen to resemble a plausible pilot scale / a month of
##     daily-resolution wearable summaries (see the All of Us Fitbit-summary caveat
##     from the pilot-study planning note)
## See SIMULATION_LOG.md for the full write-up of this run.

source("R/state_space_em.R")

generate_covariates <- function(N, sex_prob = 0.8, age_biomarker_cor = 0.3) {
  Sex <- rbinom(N, 1, sex_prob)
  Sigma <- matrix(c(1, age_biomarker_cor, age_biomarker_cor, 1), 2, 2)
  Z <- matrix(rnorm(N * 2), N, 2) %*% chol(Sigma)
  cbind(Intercept = 1, Age = Z[, 1], Sex = Sex, Biomarker = Z[, 2])
}

generate_data_mv <- function(Xmat, beta, Tt, mu, a1, P1, sigma_eta2, sigma_eps2) {
  N <- nrow(Xmat)
  phi <- tanh(as.vector(Xmat %*% beta))
  alpha <- matrix(0, N, Tt); y <- matrix(0, N, Tt)
  alpha[, 1] <- rnorm(N, a1, sqrt(P1))
  y[, 1] <- mu + alpha[, 1] + rnorm(N, 0, sqrt(sigma_eps2))
  for (t in 2:Tt) {
    alpha[, t] <- phi * alpha[, t - 1] + rnorm(N, 0, sqrt(sigma_eta2))
    y[, t] <- mu + alpha[, t] + rnorm(N, 0, sqrt(sigma_eps2))
  }
  list(y = y, alpha = alpha, phi = phi)
}

dir.create("results/simulation/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/simulation/figures", recursive = TRUE, showWarnings = FALSE)

N <- 150; Tt <- 30; mu <- 0
truth <- c(Intercept = 0.3, Age = 0.15, Sex = -0.2, Biomarker = 0.25)
a1_true <- 7; P1_true <- 15; sigma_eta2_true <- 1; sigma_eps2_true <- 1

cat("=======================================================================\n")
cat("PART E1: multi-covariate model -- single-dataset correctness check\n")
cat("=======================================================================\n")

set.seed(11)
Xmat <- generate_covariates(N)
dat <- generate_data_mv(Xmat, truth, Tt, mu, a1_true, P1_true, sigma_eta2_true, sigma_eps2_true)

fit <- em_regression_mv(dat$y, Xmat, mu = mu,
                         init = list(beta = c(0, 0, 0, 0), a1 = mean(dat$y[, 1]),
                                     P1 = var(dat$y[, 1]), sigma_eta2 = 0.5, sigma_eps2 = 0.5),
                         thres = 1e-7, max_iter = 3000)
cat(sprintf("converged: %s in %d iterations\n", fit$converged, fit$iter))
beta_hat <- setNames(fit$beta, colnames(Xmat))
cat("beta estimates:\n"); print(beta_hat)
cat(sprintf("a1=%.3f P1=%.3f sigma_eta2=%.3f sigma_eps2=%.3f\n",
            fit$a1, fit$P1, fit$sigma_eta2, fit$sigma_eps2))

se_fit <- compute_SE_regression_mv(dat$y, Xmat, mu, fit)

theta_hat <- c(fit$beta, fit$a1, fit$P1, fit$sigma_eta2, fit$sigma_eps2)
names(theta_hat) <- rownames(se_fit$I_obs)
num <- numeric_I_obs(function(th, th_star) Q_value_regression_mv(th, th_star, dat$y, Xmat, mu), theta_hat)
dimnames(num$D) <- dimnames(se_fit$D); dimnames(num$X) <- dimnames(se_fit$X); dimnames(num$I_obs) <- dimnames(se_fit$I_obs)

cat("\nmax|D_analytic - D_numeric| =", max(abs(se_fit$D - num$D)), "\n")
cat("max|X_analytic - X_numeric| =", max(abs(se_fit$X - num$X)), "\n")
cat("max|I_obs_analytic - I_obs_numeric| =", max(abs(se_fit$I_obs - num$I_obs)), "\n")

se_num <- sqrt(diag(solve(num$I_obs)))
se_compare <- data.frame(parameter = names(theta_hat), estimate = theta_hat,
                          se_analytic = se_fit$se, se_numeric = se_num,
                          abs_diff = abs(se_fit$se - se_num))
cat("\nSE comparison (analytic vs numeric Hessian):\n")
print(se_compare, row.names = FALSE, digits = 5)

write.csv(se_compare, "results/simulation/tables/part_e1_mv_se_analytic_vs_numeric.csv", row.names = FALSE)

cat("\n=======================================================================\n")
cat("PART E2: multi-covariate model -- coverage simulation (M replicates)\n")
cat("=======================================================================\n")

set.seed(12)
M <- 150
pn <- c(colnames(Xmat), "a1", "P1", "sigma_eta2", "sigma_eps2")
truth_full <- c(truth, a1 = a1_true, P1 = P1_true, sigma_eta2 = sigma_eta2_true, sigma_eps2 = sigma_eps2_true)

est_mat <- matrix(NA, M, length(pn), dimnames = list(NULL, pn))
se_mat <- matrix(NA, M, length(pn), dimnames = list(NULL, pn))
cover_mat <- matrix(NA, M, length(pn), dimnames = list(NULL, pn))
n_iter_vec <- integer(M)

t0 <- Sys.time()
for (m in seq_len(M)) {
  Xm <- generate_covariates(N)
  dm <- generate_data_mv(Xm, truth, Tt, mu, a1_true, P1_true, sigma_eta2_true, sigma_eps2_true)
  fm <- em_regression_mv(dm$y, Xm, mu = mu,
                          init = list(beta = c(0, 0, 0, 0), a1 = mean(dm$y[, 1]),
                                      P1 = var(dm$y[, 1]), sigma_eta2 = 0.5, sigma_eps2 = 0.5),
                          thres = 1e-6, max_iter = 2000)
  sem <- compute_SE_regression_mv(dm$y, Xm, mu, fm)
  est <- c(fm$beta, fm$a1, fm$P1, fm$sigma_eta2, fm$sigma_eps2)
  names(est) <- pn
  est_mat[m, ] <- est
  se_mat[m, ] <- sem$se[pn]
  n_iter_vec[m] <- fm$iter
  lo <- est - 1.96 * sem$se[pn]; hi <- est + 1.96 * sem$se[pn]
  cover_mat[m, ] <- (truth_full >= lo) & (truth_full <= hi)
}
cat(sprintf("M = %d replicates done in %.1f s\n\n", M, as.numeric(Sys.time() - t0, units = "secs")))

summary_mv <- data.frame(
  parameter = pn, truth = truth_full,
  mean_est = colMeans(est_mat), empirical_sd = apply(est_mat, 2, sd),
  mean_analytic_se = colMeans(se_mat), coverage_95 = colMeans(cover_mat)
)
print(summary_mv, row.names = FALSE, digits = 4)

write.csv(summary_mv, "results/simulation/tables/part_e2_mv_coverage_summary.csv", row.names = FALSE)
saveRDS(list(est = est_mat, se = se_mat, cover = cover_mat, truth = truth_full,
             n_iter = n_iter_vec, M = M, N = N, Tt = Tt),
        "results/simulation/tables/part_e2_mv_coverage_raw.rds")

mc_se <- sqrt(0.95 * 0.05 / M)
png("results/simulation/figures/part_e2_mv_coverage.png", width = 1000, height = 550, res = 130)
par(mar = c(7, 4, 3, 1))
bp <- barplot(summary_mv$coverage_95, ylim = c(0, 1), names.arg = rep("", nrow(summary_mv)),
              col = "seagreen", ylab = "empirical coverage of 95% Wald CI",
              main = sprintf("Multi-covariate model, M = %d", M))
text(bp, par("usr")[3] - 0.03, labels = summary_mv$parameter, srt = 30, adj = c(1, 1), xpd = TRUE, cex = 0.85)
abline(h = 0.95, lty = 2, col = "red")
rect(min(bp) - 1, 0.95 - 1.96 * mc_se, max(bp) + 1, 0.95 + 1.96 * mc_se, col = rgb(1, 0, 0, 0.1), border = NA)
dev.off()

png("results/simulation/figures/part_e2_mv_estimates.png", width = 1300, height = 650, res = 130)
par(mfrow = c(2, 4), mar = c(3, 3, 2, 1))
for (j in seq_along(pn)) {
  boxplot(est_mat[, j], main = pn[j], ylab = "")
  abline(h = truth_full[j], col = "red", lty = 2)
}
dev.off()

png("results/simulation/figures/part_e2_mv_se_calibration.png", width = 700, height = 550, res = 130)
plot(summary_mv$empirical_sd, summary_mv$mean_analytic_se, pch = 19, col = "seagreen",
     xlab = "empirical SD across replicates", ylab = "mean analytic SE",
     main = "Multi-covariate model: SE calibration")
abline(0, 1, lty = 2, col = "grey40")
text(summary_mv$empirical_sd, summary_mv$mean_analytic_se, summary_mv$parameter, pos = 3, cex = 0.7)
dev.off()

cat("\nSaved: results/simulation/tables/part_e1_mv_se_analytic_vs_numeric.csv\n")
cat("Saved: results/simulation/tables/part_e2_mv_coverage_summary.csv\n")
cat("Saved: results/simulation/tables/part_e2_mv_coverage_raw.rds\n")
cat("Saved: results/simulation/figures/part_e2_mv_coverage.png\n")
cat("Saved: results/simulation/figures/part_e2_mv_estimates.png\n")
cat("Saved: results/simulation/figures/part_e2_mv_se_calibration.png\n")
