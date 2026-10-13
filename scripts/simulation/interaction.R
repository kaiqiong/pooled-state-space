# Portable October 13 replication. Run from the repository root; see docs/SIMULATIONS.txt.
## Simulation: phi_i = tanh(beta0 + beta_age*Age_i + beta_sex*Sex_i
##                          + beta_bio*Biomarker_i + beta_sexbio*(Sex_i*Biomarker_i))
##
## Motivation: the Sex:Biomarker interaction directly represents "a molecular
## feature that drives recovery dynamics mainly/only in one sex" -- exactly the
## kind of sex-specific effect the grant proposal flags as scientifically
## important given hEDS's strong female predominance. No new derivation is
## needed: em_regression_mv/compute_SE_regression_mv/Q_value_regression_mv are
## already written for an arbitrary design matrix (SIMULATION_LOG.md, Part E),
## and an interaction term is just one more column of X. What's new here is
## validating that under real collinearity (Sex:Biomarker correlates with
## Biomarker whenever Sex is imbalanced) everything still holds up, and
## checking what happens if an analyst fits a main-effects-only model to data
## that actually has the interaction (a realistic misspecification scenario).
##
## Same N=150, T=30, Sex~Bernoulli(0.8), correlated Age/Biomarker as Part E.
## Truth chosen so the Biomarker effect is small for males (0.10) and much
## larger for females (0.10+0.30=0.40) -- a literal "sex-specific molecular
## driver" scenario.

source("R/state_space_em.R")

generate_covariates_interaction <- function(N, sex_prob = 0.8, age_biomarker_cor = 0.3) {
  Sex <- rbinom(N, 1, sex_prob)
  Sigma <- matrix(c(1, age_biomarker_cor, age_biomarker_cor, 1), 2, 2)
  Z <- matrix(rnorm(N * 2), N, 2) %*% chol(Sigma)
  Age <- Z[, 1]; Biomarker <- Z[, 2]
  cbind(Intercept = 1, Age = Age, Sex = Sex, Biomarker = Biomarker, `Sex:Biomarker` = Sex * Biomarker)
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
truth <- c(Intercept = 0.3, Age = 0.15, Sex = -0.2, Biomarker = 0.10, `Sex:Biomarker` = 0.30)
a1_true <- 7; P1_true <- 15; sigma_eta2_true <- 1; sigma_eps2_true <- 1

cat("=======================================================================\n")
cat("PART H1: interaction model -- single-dataset correctness check\n")
cat("=======================================================================\n")
cat("True Biomarker effect: male (Sex=0) =", truth["Biomarker"],
    " | female (Sex=1) =", truth["Biomarker"] + truth["Sex:Biomarker"], "\n\n")

set.seed(41)
Xmat <- generate_covariates_interaction(N)
dat <- generate_data_mv(Xmat, truth, Tt, mu, a1_true, P1_true, sigma_eta2_true, sigma_eps2_true)

cat("condition number of Xmat (collinearity check):", kappa(Xmat), "\n")
cat("cor(Biomarker, Sex:Biomarker) =", cor(Xmat[, "Biomarker"], Xmat[, "Sex:Biomarker"]), "\n\n")

fit <- em_regression_mv(dat$y, Xmat, mu = mu,
                         init = list(beta = c(0, 0, 0, 0, 0), a1 = mean(dat$y[, 1]),
                                     P1 = var(dat$y[, 1]), sigma_eta2 = 0.5, sigma_eps2 = 0.5),
                         thres = 1e-7, max_iter = 3000)
cat(sprintf("converged: %s in %d iterations\n", fit$converged, fit$iter))
beta_hat <- setNames(fit$beta, colnames(Xmat))
cat("beta estimates:\n"); print(beta_hat)

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
write.csv(se_compare, "results/simulation/tables/part_h1_interaction_se_analytic_vs_numeric.csv", row.names = FALSE)

cat("\n=======================================================================\n")
cat("PART H2: interaction model -- coverage simulation (correctly specified)\n")
cat("=======================================================================\n")

set.seed(42)
M <- 150
pn <- c(colnames(Xmat), "a1", "P1", "sigma_eta2", "sigma_eps2")
truth_full <- c(truth, a1 = a1_true, P1 = P1_true, sigma_eta2 = sigma_eta2_true, sigma_eps2 = sigma_eps2_true)

est_mat <- matrix(NA, M, length(pn), dimnames = list(NULL, pn))
se_mat <- matrix(NA, M, length(pn), dimnames = list(NULL, pn))
cover_mat <- matrix(NA, M, length(pn), dimnames = list(NULL, pn))

## PART H3 (bonus, same replicates): what happens if an analyst fits a
## main-effects-only model (omits Sex:Biomarker) to this interaction-generated
## data? Compare the resulting Biomarker coefficient to the true male/female-
## specific effects and to the prevalence-weighted average effect.
reduced_pn <- c("Intercept", "Age", "Sex", "Biomarker", "a1", "P1", "sigma_eta2", "sigma_eps2")
est_reduced <- matrix(NA, M, length(reduced_pn), dimnames = list(NULL, reduced_pn))

t0 <- Sys.time()
for (m in seq_len(M)) {
  Xm <- generate_covariates_interaction(N)
  dm <- generate_data_mv(Xm, truth, Tt, mu, a1_true, P1_true, sigma_eta2_true, sigma_eps2_true)

  fm <- em_regression_mv(dm$y, Xm, mu = mu,
                          init = list(beta = c(0, 0, 0, 0, 0), a1 = mean(dm$y[, 1]),
                                      P1 = var(dm$y[, 1]), sigma_eta2 = 0.5, sigma_eps2 = 0.5),
                          thres = 1e-6, max_iter = 2000)
  sem <- compute_SE_regression_mv(dm$y, Xm, mu, fm)
  est <- c(fm$beta, fm$a1, fm$P1, fm$sigma_eta2, fm$sigma_eps2)
  names(est) <- pn
  est_mat[m, ] <- est
  se_mat[m, ] <- sem$se[pn]
  lo <- est - 1.96 * sem$se[pn]; hi <- est + 1.96 * sem$se[pn]
  cover_mat[m, ] <- (truth_full >= lo) & (truth_full <= hi)

  Xm_reduced <- Xm[, c("Intercept", "Age", "Sex", "Biomarker")]
  fr <- em_regression_mv(dm$y, Xm_reduced, mu = mu,
                          init = list(beta = c(0, 0, 0, 0), a1 = mean(dm$y[, 1]),
                                      P1 = var(dm$y[, 1]), sigma_eta2 = 0.5, sigma_eps2 = 0.5),
                          thres = 1e-6, max_iter = 2000)
  est_reduced[m, ] <- c(fr$beta, fr$a1, fr$P1, fr$sigma_eta2, fr$sigma_eps2)
}
cat(sprintf("M = %d replicates (both full and reduced model) done in %.1f s\n\n",
            M, as.numeric(Sys.time() - t0, units = "secs")))

summary_h <- data.frame(
  parameter = pn, truth = truth_full,
  mean_est = colMeans(est_mat), empirical_sd = apply(est_mat, 2, sd),
  mean_analytic_se = colMeans(se_mat), coverage_95 = colMeans(cover_mat)
)
cat("Part H2 -- correctly specified model, coverage:\n")
print(summary_h, row.names = FALSE, digits = 4)

write.csv(summary_h, "results/simulation/tables/part_h2_interaction_coverage_summary.csv", row.names = FALSE)
saveRDS(list(est = est_mat, se = se_mat, cover = cover_mat, truth = truth_full, M = M, N = N, Tt = Tt),
        "results/simulation/tables/part_h2_interaction_coverage_raw.rds")

cat("\n=======================================================================\n")
cat("PART H3: consequence of omitting the interaction (main-effects-only fit)\n")
cat("=======================================================================\n")

male_effect <- truth["Biomarker"]
female_effect <- truth["Biomarker"] + truth["Sex:Biomarker"]
prevalence_weighted <- mean(Xmat[, "Sex"]) * female_effect + (1 - mean(Xmat[, "Sex"])) * male_effect

reduced_summary <- data.frame(
  parameter = reduced_pn,
  mean_est_reduced_model = colMeans(est_reduced),
  empirical_sd = apply(est_reduced, 2, sd)
)
cat("Reduced (main-effects-only) model, mean estimates across", M, "replicates:\n")
print(reduced_summary, row.names = FALSE, digits = 4)

cat(sprintf("\nTrue Biomarker effect: male = %.3f, female = %.3f\n", male_effect, female_effect))
cat(sprintf("Prevalence-weighted average (%.0f%% female) = %.3f\n", 100 * mean(Xmat[, "Sex"]), prevalence_weighted))
cat(sprintf("Reduced model's mean Biomarker estimate = %.3f\n",
            reduced_summary$mean_est_reduced_model[reduced_summary$parameter == "Biomarker"]))
cat(sprintf("Reduced model's mean Sex estimate = %.3f (truth, correctly specified model = %.3f)\n",
            reduced_summary$mean_est_reduced_model[reduced_summary$parameter == "Sex"], truth["Sex"]))

write.csv(reduced_summary, "results/simulation/tables/part_h3_reduced_model_summary.csv", row.names = FALSE)
saveRDS(est_reduced, "results/simulation/tables/part_h3_reduced_model_raw.rds")

## ---- plots ----
mc_se <- sqrt(0.95 * 0.05 / M)
png("results/simulation/figures/part_h2_interaction_coverage.png", width = 1100, height = 550, res = 130)
par(mar = c(7, 4, 3, 1))
bp <- barplot(summary_h$coverage_95, ylim = c(0, 1), names.arg = rep("", nrow(summary_h)),
              col = "darkorchid", ylab = "empirical coverage of 95% Wald CI",
              main = sprintf("Interaction model, M = %d", M))
text(bp, par("usr")[3] - 0.03, labels = summary_h$parameter, srt = 30, adj = c(1, 1), xpd = TRUE, cex = 0.8)
abline(h = 0.95, lty = 2, col = "red")
rect(min(bp) - 1, 0.95 - 1.96 * mc_se, max(bp) + 1, 0.95 + 1.96 * mc_se, col = rgb(1, 0, 0, 0.1), border = NA)
dev.off()

png("results/simulation/figures/part_h3_omitted_interaction.png", width = 700, height = 550, res = 130)
biomarker_ests <- est_reduced[, "Biomarker"]
hist(biomarker_ests, breaks = 20, col = "grey80", border = "white",
     main = "Reduced (main-effects-only) model's Biomarker estimate\nwhen truth has a Sex:Biomarker interaction",
     xlab = "estimated Biomarker coefficient")
abline(v = male_effect, col = "steelblue", lwd = 2, lty = 2)
abline(v = female_effect, col = "firebrick", lwd = 2, lty = 2)
abline(v = prevalence_weighted, col = "black", lwd = 2)
legend("topleft", legend = c("true male effect", "true female effect", "prevalence-weighted avg"),
       col = c("steelblue", "firebrick", "black"), lty = c(2, 2, 1), lwd = 2, bty = "n", cex = 0.8)
dev.off()

cat("\nSaved: results/simulation/tables/part_h1_interaction_se_analytic_vs_numeric.csv\n")
cat("Saved: results/simulation/tables/part_h2_interaction_coverage_summary.csv\n")
cat("Saved: results/simulation/tables/part_h2_interaction_coverage_raw.rds\n")
cat("Saved: results/simulation/tables/part_h3_reduced_model_summary.csv\n")
cat("Saved: results/simulation/tables/part_h3_reduced_model_raw.rds\n")
cat("Saved: results/simulation/figures/part_h2_interaction_coverage.png\n")
cat("Saved: results/simulation/figures/part_h3_omitted_interaction.png\n")
