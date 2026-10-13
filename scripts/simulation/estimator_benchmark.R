# Portable October 13 replication. Run from the repository root; see docs/SIMULATIONS.txt.
## Extended benchmark (Part G): adds the two competitors suggested at the end of
## Part F to the SAME non-interaction multi-covariate design used there:
##   1. per_subject_kalman_two_stage() -- each subject's own state-space model
##      fit separately by Kalman-filter MLE (correct measurement-error handling,
##      but no pooling of variance components or beta across subjects), then
##      Fisher-z + OLS on covariates, exactly like naive_two_stage()'s second
##      stage. Isolates "correcting for measurement error alone" from "also
##      pooling", relative to the fully joint EM.
##   2. stratified_gls_corAR1() -- nlme::gls with an AR(1) residual correlation,
##      fit separately by Sex. No latent-state/measurement-error split, and
##      structurally cannot use continuous covariates (Age, Biomarker) at all --
##      included to show what a "reasonable applied statistician's" marginal-
##      model baseline can and cannot do.
##
## Because gls only ever produces ONE number per stratum, the common comparison
## metric across all four methods is the Sex-stratum-average phi (not full beta
## vectors) for Part G2; Part G1 compares full beta vectors across the three
## methods that can estimate one (joint EM, naive two-stage, per-subject Kalman
## two-stage).
##
## M reduced to 20 (from Part F's 60) because adding 150 additional per-subject
## MLE fits per replicate materially raises the per-replicate cost; see
## SIMULATION_LOG.md for the runtime accounting.

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
  list(y = y, phi = phi)
}

dir.create("results/simulation/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/simulation/figures", recursive = TRUE, showWarnings = FALSE)

N <- 150; Tt <- 30; mu <- 0
truth <- c(Intercept = 0.3, Age = 0.15, Sex = -0.2, Biomarker = 0.25)
a1_true <- 7; P1_true <- 15; sigma_eta2_true <- 1
M <- 60
snr_settings <- c(low_noise = 0.25, moderate_noise = 1, high_noise = 4)

cat("=======================================================================\n")
cat("PART G: extended benchmark -- + per-subject Kalman MLE, + stratified gls\n")
cat("=======================================================================\n")

beta_rows <- list()
strat_rows <- list()
raw <- list()
set.seed(31)

for (setting in names(snr_settings)) {
  sigma_eps2_true <- snr_settings[[setting]]
  cat(sprintf("\n--- setting: %s (sigma_eps2 = %.2f, sigma_eta2 = %.2f) ---\n",
              setting, sigma_eps2_true, sigma_eta2_true))

  cov_names <- colnames(generate_covariates(2))
  est_joint <- matrix(NA, M, 4, dimnames = list(NULL, cov_names))
  est_naive <- matrix(NA, M, 4, dimnames = list(NULL, cov_names))
  est_kalman <- matrix(NA, M, 4, dimnames = list(NULL, cov_names))
  strat_true <- matrix(NA, M, 2, dimnames = list(NULL, c("0", "1")))
  strat_joint <- matrix(NA, M, 2, dimnames = list(NULL, c("0", "1")))
  strat_naive <- matrix(NA, M, 2, dimnames = list(NULL, c("0", "1")))
  strat_kalman <- matrix(NA, M, 2, dimnames = list(NULL, c("0", "1")))
  strat_gls <- matrix(NA, M, 2, dimnames = list(NULL, c("0", "1")))

  t0 <- Sys.time()
  for (m in seq_len(M)) {
    Xm <- generate_covariates(N)
    dm <- generate_data_mv(Xm, truth, Tt, mu, a1_true, P1_true, sigma_eta2_true, sigma_eps2_true)
    sex <- Xm[, "Sex"]

    fm <- em_regression_mv(dm$y, Xm, mu = mu,
                            init = list(beta = c(0, 0, 0, 0), a1 = mean(dm$y[, 1]),
                                        P1 = var(dm$y[, 1]), sigma_eta2 = 0.5, sigma_eps2 = 0.5),
                            thres = 1e-6, max_iter = 1500)
    est_joint[m, ] <- fm$beta
    phi_joint <- tanh(as.vector(Xm %*% fm$beta))

    nv <- naive_two_stage(dm$y, Xm)
    est_naive[m, ] <- nv$beta_hat

    ks <- per_subject_kalman_two_stage(dm$y, Xm, mu = mu)
    est_kalman[m, ] <- ks$beta_hat

    gl <- stratified_gls_corAR1(dm$y, sex)

    strat_true[m, ]   <- tapply(dm$phi, sex, mean)[c("0", "1")]
    strat_joint[m, ]  <- tapply(phi_joint, sex, mean)[c("0", "1")]
    strat_naive[m, ]  <- tapply(nv$phi_hat, sex, mean)[c("0", "1")]
    strat_kalman[m, ] <- tapply(ks$phi_hat, sex, mean)[c("0", "1")]
    strat_gls[m, ]    <- gl[c("0", "1")]
  }
  cat(sprintf("  M = %d replicates done in %.1f s\n", M, as.numeric(Sys.time() - t0, units = "secs")))

  raw[[setting]] <- list(est_joint = est_joint, est_naive = est_naive, est_kalman = est_kalman,
                          strat_true = strat_true, strat_joint = strat_joint, strat_naive = strat_naive,
                          strat_kalman = strat_kalman, strat_gls = strat_gls, sigma_eps2 = sigma_eps2_true)

  for (method_name in c("joint_EM", "naive_two_stage", "per_subject_kalman")) {
    est <- switch(method_name, joint_EM = est_joint, naive_two_stage = est_naive, per_subject_kalman = est_kalman)
    for (j in seq_along(truth)) {
      beta_rows[[length(beta_rows) + 1]] <- data.frame(
        setting = setting, sigma_eps2 = sigma_eps2_true, method = method_name,
        parameter = names(truth)[j], truth = truth[j],
        mean_est = mean(est[, j]), bias = mean(est[, j]) - truth[j],
        rel_bias_pct = 100 * (mean(est[, j]) - truth[j]) / truth[j],
        rmse = sqrt(mean((est[, j] - truth[j])^2))
      )
    }
  }

  for (method_name in c("joint_EM", "naive_two_stage", "per_subject_kalman", "stratified_gls")) {
    est <- switch(method_name, joint_EM = strat_joint, naive_two_stage = strat_naive,
                  per_subject_kalman = strat_kalman, stratified_gls = strat_gls)
    for (s in c("0", "1")) {
      truth_s <- mean(strat_true[, s])
      strat_rows[[length(strat_rows) + 1]] <- data.frame(
        setting = setting, sigma_eps2 = sigma_eps2_true, method = method_name, sex = s,
        truth_avg_phi = truth_s, mean_est = mean(est[, s]), bias = mean(est[, s]) - truth_s,
        rmse = sqrt(mean((est[, s] - truth_s)^2))
      )
    }
  }
}

beta_df <- do.call(rbind, beta_rows); rownames(beta_df) <- NULL
strat_df <- do.call(rbind, strat_rows); rownames(strat_df) <- NULL

cat("\nPart G1 -- regression coefficient benchmark (joint EM / naive / per-subject Kalman):\n")
print(beta_df, row.names = FALSE, digits = 4)

cat("\nPart G2 -- Sex-stratum-average phi benchmark (+ stratified gls):\n")
print(strat_df, row.names = FALSE, digits = 4)

write.csv(beta_df, "results/simulation/tables/part_g1_M60_beta_benchmark.csv", row.names = FALSE)
write.csv(strat_df, "results/simulation/tables/part_g2_M60_stratum_phi_benchmark.csv", row.names = FALSE)
saveRDS(list(raw = raw, beta_df = beta_df, strat_df = strat_df, truth = truth, N = N, Tt = Tt, M = M),
        "results/simulation/tables/part_g_M60_extended_benchmark_raw.rds")

## ---- plot G1: % relative bias vs sigma_eps2, one panel per covariate, 3 methods ----
method_cols <- c(joint_EM = "steelblue", naive_two_stage = "firebrick", per_subject_kalman = "darkorange")
method_pch  <- c(joint_EM = 19, naive_two_stage = 17, per_subject_kalman = 15)

png("results/simulation/figures/part_g1_M60_beta_bias.png", width = 1200, height = 900, res = 130)
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
for (p in names(truth)) {
  sub <- beta_df[beta_df$parameter == p, ]
  ylim <- range(sub$rel_bias_pct, 0)
  plot(NULL, xlim = range(sub$sigma_eps2), ylim = ylim,
       xlab = expression(sigma[epsilon]^2), ylab = "% relative bias", main = p)
  abline(h = 0, lty = 2, col = "grey50")
  for (meth in names(method_cols)) {
    s <- sub[sub$method == meth, ]; s <- s[order(s$sigma_eps2), ]
    lines(s$sigma_eps2, s$rel_bias_pct, type = "b", col = method_cols[meth], pch = method_pch[meth])
  }
  if (p == names(truth)[1]) legend("bottomleft", legend = names(method_cols), col = method_cols,
                                    pch = method_pch, bty = "n", cex = 0.75)
}
dev.off()

png("results/simulation/figures/part_g1_M60_beta_rmse.png", width = 1200, height = 900, res = 130)
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
for (p in names(truth)) {
  sub <- beta_df[beta_df$parameter == p, ]
  ylim <- range(sub$rmse, 0)
  plot(NULL, xlim = range(sub$sigma_eps2), ylim = ylim,
       xlab = expression(sigma[epsilon]^2), ylab = "RMSE", main = p)
  for (meth in names(method_cols)) {
    s <- sub[sub$method == meth, ]; s <- s[order(s$sigma_eps2), ]
    lines(s$sigma_eps2, s$rmse, type = "b", col = method_cols[meth], pch = method_pch[meth])
  }
  if (p == names(truth)[1]) legend("topleft", legend = names(method_cols), col = method_cols,
                                    pch = method_pch, bty = "n", cex = 0.75)
}
dev.off()

## ---- plot G2: stratum-average-phi bias, all 4 methods, both strata ----
method_cols4 <- c(joint_EM = "steelblue", naive_two_stage = "firebrick",
                   per_subject_kalman = "darkorange", stratified_gls = "purple4")
method_pch4  <- c(joint_EM = 19, naive_two_stage = 17, per_subject_kalman = 15, stratified_gls = 18)

png("results/simulation/figures/part_g2_M60_stratum_phi_bias.png", width = 1200, height = 550, res = 130)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
for (s in c("0", "1")) {
  sub <- strat_df[strat_df$sex == s, ]
  ylim <- range(sub$bias, 0)
  plot(NULL, xlim = range(sub$sigma_eps2), ylim = ylim, xlab = expression(sigma[epsilon]^2),
       ylab = "bias in stratum-average phi",
       main = paste("Sex =", s, ifelse(s == "1", "(female)", "(male)")))
  abline(h = 0, lty = 2, col = "grey50")
  for (meth in names(method_cols4)) {
    ss <- sub[sub$method == meth, ]; ss <- ss[order(ss$sigma_eps2), ]
    lines(ss$sigma_eps2, ss$bias, type = "b", col = method_cols4[meth], pch = method_pch4[meth])
  }
  if (s == "0") legend("bottomleft", legend = names(method_cols4), col = method_cols4,
                        pch = method_pch4, bty = "n", cex = 0.7)
}
dev.off()

cat("\nSaved: results/simulation/tables/part_g1_M60_beta_benchmark.csv\n")
cat("Saved: results/simulation/tables/part_g2_M60_stratum_phi_benchmark.csv\n")
cat("Saved: results/simulation/tables/part_g_M60_extended_benchmark_raw.rds\n")
cat("Saved: results/simulation/figures/part_g1_M60_beta_bias.png\n")
cat("Saved: results/simulation/figures/part_g1_M60_beta_rmse.png\n")
cat("Saved: results/simulation/figures/part_g2_M60_stratum_phi_bias.png\n")
