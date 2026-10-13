# Portable October 13 replication. Run from the repository root; see docs/SIMULATIONS.txt.
source("R/state_space_em.R")

generate_data_homogeneous <- function(N, Tt, phi, mu, a1, P1, sigma_eta2, sigma_eps2) {
  alpha <- matrix(0, N, Tt); y <- matrix(0, N, Tt)
  alpha[, 1] <- rnorm(N, a1, sqrt(P1))
  y[, 1] <- mu + alpha[, 1] + rnorm(N, 0, sqrt(sigma_eps2))
  for (t in 2:Tt) {
    alpha[, t] <- phi * alpha[, t - 1] + rnorm(N, 0, sqrt(sigma_eta2))
    y[, t] <- mu + alpha[, t] + rnorm(N, 0, sqrt(sigma_eps2))
  }
  list(y = y, alpha = alpha)
}

generate_data_regression <- function(N, Tt, beta0, beta1, mu, a1, P1, sigma_eta2, sigma_eps2) {
  x <- rbinom(N, 1, 0.5)
  phi <- tanh(beta0 + beta1 * x)
  alpha <- matrix(0, N, Tt); y <- matrix(0, N, Tt)
  alpha[, 1] <- rnorm(N, a1, sqrt(P1))
  y[, 1] <- mu + alpha[, 1] + rnorm(N, 0, sqrt(sigma_eps2))
  for (t in 2:Tt) {
    alpha[, t] <- phi * alpha[, t - 1] + rnorm(N, 0, sqrt(sigma_eta2))
    y[, t] <- mu + alpha[, t] + rnorm(N, 0, sqrt(sigma_eps2))
  }
  list(y = y, x = x, alpha = alpha, phi = phi)
}

## ---------------------------------------------------------------------------
## PART C: original (loop-over-i-and-t) E-step vs vectorized kalman_run
## Original em_e_step1 copied verbatim from Report_sep_3_26/ENURA26_code_final.R
## ---------------------------------------------------------------------------
em_e_step1_original <- function(init_params.v, x, y, mu, N, T) {
  beta0_init <- init_params.v[1]; beta1_init <- init_params.v[2]
  a1_init <- init_params.v[3]; log_P1_init <- init_params.v[4]
  log_sigma_eta2_init <- init_params.v[5]; log_sigma_eps2_init <- init_params.v[6]
  P1_init <- exp(log_P1_init); sigma_eta2_init <- exp(log_sigma_eta2_init)
  sigma_eps2_init <- exp(log_sigma_eps2_init)

  alpha_pred <- matrix(0, N, T); P_pred <- matrix(0, N, T)
  alpha_filt <- matrix(0, N, T); P_filt <- matrix(0, N, T); K <- matrix(0, N, T)
  alpha_smooth <- matrix(0, N, T); P_smooth <- matrix(0, N, T); J <- matrix(0, N, T)
  P_cov_smooth <- matrix(0, N, T)

  for (i in 1:N) {
    alpha_pred[i, 1] <- a1_init; P_pred[i, 1] <- P1_init
    K[i, 1] <- P_pred[i, 1] / (P_pred[i, 1] + sigma_eps2_init)
    alpha_filt[i, 1] <- alpha_pred[i, 1] + K[i, 1] * (y[i, 1] - (mu + alpha_pred[i, 1]))
    P_filt[i, 1] <- (1 - K[i, 1]) * P_pred[i, 1]
    for (t in 2:T) {
      alpha_pred[i, t] <- tanh(beta0_init + beta1_init * x[i]) * alpha_filt[i, t - 1]
      P_pred[i, t] <- (tanh(beta0_init + beta1_init * x[i]))^2 * P_filt[i, t - 1] + sigma_eta2_init
      K[i, t] <- P_pred[i, t] / (P_pred[i, t] + sigma_eps2_init)
      alpha_filt[i, t] <- alpha_pred[i, t] + K[i, t] * (y[i, t] - (mu + alpha_pred[i, t]))
      P_filt[i, t] <- (1 - K[i, t]) * P_pred[i, t]
    }
    alpha_smooth[i, T] <- alpha_filt[i, T]; P_smooth[i, T] <- P_filt[i, T]
    for (t in T:2) {
      J[i, t - 1] <- tanh(beta0_init + beta1_init * x[i]) * P_filt[i, t - 1] / P_pred[i, t]
      alpha_smooth[i, t - 1] <- alpha_filt[i, t - 1] + J[i, t - 1] * (alpha_smooth[i, t] - alpha_pred[i, t])
      P_smooth[i, t - 1] <- P_filt[i, t - 1] + J[i, t - 1]^2 * (P_smooth[i, t] - P_pred[i, t])
    }
    P_cov_smooth[i, T] <- (1 - K[i, T]) * tanh(beta0_init + beta1_init * x[i]) * P_filt[i, T - 1]
    for (t in T:3) {
      P_cov_smooth[i, t - 1] <- P_filt[i, t - 1] * J[i, t - 2] +
        J[i, t - 1] * (P_cov_smooth[i, t] - tanh(beta0_init + beta1_init * x[i]) * P_filt[i, t - 1]) * J[i, t - 2]
    }
  }
  list(alpha_pred = alpha_pred, P_pred = P_pred, alpha_filt = alpha_filt, P_filt = P_filt,
       alpha_smooth = alpha_smooth, P_smooth = P_smooth, P_cov_smooth = P_cov_smooth)
}

cat("=======================================================================\n")
cat("PART C: timing, original loop-over-(i,t) E-step vs vectorized kalman_run\n")
cat("=======================================================================\n")

timing_rows <- list()
for (NN in c(100, 500, 2000)) {
  Tt <- 35
  set.seed(42)
  dat <- generate_data_regression(NN, Tt, 0.4, -0.2, 0, 7, 15, 1, 1)
  init_v <- c(0.4, -0.2, 7, log(15), log(1), log(1))

  t_orig <- system.time(for (r in 1:5) em_e_step1_original(init_v, dat$x, dat$y, 0, NN, Tt))["elapsed"] / 5
  t_new  <- system.time(for (r in 1:5) kalman_run(dat$y, tanh(0.4 - 0.2 * dat$x), 0, 7, 15, 1, 1))["elapsed"] / 5

  cat(sprintf("N=%5d, T=%d:  original = %7.4f s   vectorized = %7.4f s   speed-up = %6.1fx\n",
              NN, Tt, t_orig, t_new, t_orig / t_new))
  timing_rows[[length(timing_rows) + 1]] <- data.frame(
    N = NN, T = Tt, seconds_original = t_orig, seconds_vectorized = t_new,
    speedup = t_orig / t_new
  )
}
timing_df <- do.call(rbind, timing_rows)

dir.create("results/simulation/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/simulation/figures", recursive = TRUE, showWarnings = FALSE)
write.csv(timing_df, "results/simulation/tables/part_c_timing.csv", row.names = FALSE)

png("results/simulation/figures/part_c_timing_speedup.png", width = 700, height = 500, res = 130)
plot(timing_df$N, timing_df$speedup, type = "b", pch = 19, col = "steelblue",
     xlab = "N (subjects)", ylab = "speed-up (original / vectorized)",
     main = "E-step vectorization speed-up, T = 35", log = "x")
abline(h = 1, lty = 2, col = "grey50")
dev.off()
cat("Saved: results/simulation/tables/part_c_timing.csv\n")
cat("Saved: results/simulation/figures/part_c_timing_speedup.png\n")

## ---------------------------------------------------------------------------
## PART D: coverage of 95% Wald CIs built from the analytic SEs
## ---------------------------------------------------------------------------
cat("\n=======================================================================\n")
cat("PART D: coverage simulation, homogeneous model\n")
cat("=======================================================================\n")

set.seed(100)
M <- 200; N <- 100; Tt <- 35; mu <- 0
truth_h <- c(phi = 0.6, a1 = 7, P1 = 15, sigma_eta2 = 1, sigma_eps2 = 1)

cover_h <- matrix(NA, M, 5, dimnames = list(NULL, names(truth_h)))
est_h <- matrix(NA, M, 5, dimnames = list(NULL, names(truth_h)))
se_store_h <- matrix(NA, M, 5, dimnames = list(NULL, names(truth_h)))

for (m in 1:M) {
  dat <- generate_data_homogeneous(N, Tt, truth_h["phi"], mu, truth_h["a1"], truth_h["P1"],
                                    truth_h["sigma_eta2"], truth_h["sigma_eps2"])
  fit <- em_homogeneous(dat$y, mu = mu,
                         init = list(phi = 0.3, a1 = mean(dat$y[, 1]), P1 = var(dat$y[, 1]),
                                     sigma_eta2 = 0.5, sigma_eps2 = 0.5),
                         thres = 1e-7, max_iter = 2000)
  se <- compute_SE_homogeneous(dat$y, mu, fit)
  est <- unlist(fit[c("phi", "a1", "P1", "sigma_eta2", "sigma_eps2")])
  est_h[m, ] <- est
  se_store_h[m, ] <- se$se[names(truth_h)]
  lo <- est - 1.96 * se$se[names(truth_h)]
  hi <- est + 1.96 * se$se[names(truth_h)]
  cover_h[m, ] <- (truth_h >= lo) & (truth_h <= hi)
}

cat(sprintf("M = %d replicates, N = %d, T = %d\n\n", M, N, Tt))
summary_h <- data.frame(
  parameter = names(truth_h), truth = truth_h,
  mean_est = colMeans(est_h), empirical_sd = apply(est_h, 2, sd),
  mean_analytic_se = colMeans(se_store_h),
  coverage_95 = colMeans(cover_h)
)
print(summary_h, row.names = FALSE, digits = 4)

## ---- save Part D (homogeneous): raw + summary tables, plots ----
M_h <- M
write.csv(summary_h, "results/simulation/tables/part_d_coverage_homogeneous_summary.csv", row.names = FALSE)
saveRDS(list(est = est_h, se = se_store_h, cover = cover_h, truth = truth_h, M = M_h, N = N, Tt = Tt),
        "results/simulation/tables/part_d_coverage_homogeneous_raw.rds")

plot_coverage <- function(summary_df, M, file, main) {
  mc_se <- sqrt(0.95 * 0.05 / M)
  png(file, width = 900, height = 550, res = 130)
  par(mar = c(6, 4, 3, 1))
  bp <- barplot(summary_df$coverage_95, ylim = c(0, 1), names.arg = rep("", nrow(summary_df)),
                col = "steelblue", ylab = "empirical coverage of 95% Wald CI", main = main)
  text(bp, par("usr")[3] - 0.03, labels = summary_df$parameter, srt = 30, adj = c(1, 1),
       xpd = TRUE, cex = 0.9)
  abline(h = 0.95, lty = 2, col = "red")
  rect(min(bp) - 1, 0.95 - 1.96 * mc_se, max(bp) + 1, 0.95 + 1.96 * mc_se,
       col = rgb(1, 0, 0, 0.1), border = NA)
  dev.off()
}
plot_coverage(summary_h, M_h, "results/simulation/figures/part_d_coverage_homogeneous.png",
              sprintf("Homogeneous model, M = %d", M_h))

plot_est_boxplot <- function(est, truth, file, main) {
  png(file, width = 900, height = 500, res = 130)
  par(mfrow = c(1, ncol(est)), mar = c(3, 3, 2, 1))
  for (j in seq_len(ncol(est))) {
    boxplot(est[, j], main = colnames(est)[j], ylab = "")
    abline(h = truth[j], col = "red", lty = 2)
  }
  mtext(main, outer = TRUE, line = -1.5)
  dev.off()
}
plot_est_boxplot(est_h, truth_h, "results/simulation/figures/part_d_estimates_homogeneous.png",
                  sprintf("Homogeneous model final estimates vs truth (dashed), M = %d", M_h))

cat("Saved: results/simulation/tables/part_d_coverage_homogeneous_summary.csv\n")
cat("Saved: results/simulation/tables/part_d_coverage_homogeneous_raw.rds\n")
cat("Saved: results/simulation/figures/part_d_coverage_homogeneous.png\n")
cat("Saved: results/simulation/figures/part_d_estimates_homogeneous.png\n")

cat("\n=======================================================================\n")
cat("PART D (cont.): coverage simulation, regression model (target: beta1)\n")
cat("=======================================================================\n")

set.seed(101)
M <- 100
truth_r <- c(beta0 = 0.4, beta1 = -0.2, a1 = 7, P1 = 15, sigma_eta2 = 1, sigma_eps2 = 1)

cover_r <- matrix(NA, M, 6, dimnames = list(NULL, names(truth_r)))
est_r <- matrix(NA, M, 6, dimnames = list(NULL, names(truth_r)))
se_store_r <- matrix(NA, M, 6, dimnames = list(NULL, names(truth_r)))

for (m in 1:M) {
  dat <- generate_data_regression(N, Tt, truth_r["beta0"], truth_r["beta1"], mu,
                                   truth_r["a1"], truth_r["P1"], truth_r["sigma_eta2"], truth_r["sigma_eps2"])
  fit <- em_regression(dat$y, dat$x, mu = mu,
                        init = list(beta0 = 0.1, beta1 = 0.1, a1 = mean(dat$y[, 1]),
                                    P1 = var(dat$y[, 1]), sigma_eta2 = 0.5, sigma_eps2 = 0.5),
                        thres = 1e-6, max_iter = 2000)
  se <- compute_SE_regression(dat$y, dat$x, mu, fit)
  est <- unlist(fit[c("beta0", "beta1", "a1", "P1", "sigma_eta2", "sigma_eps2")])
  est_r[m, ] <- est
  se_store_r[m, ] <- se$se[names(truth_r)]
  lo <- est - 1.96 * se$se[names(truth_r)]
  hi <- est + 1.96 * se$se[names(truth_r)]
  cover_r[m, ] <- (truth_r >= lo) & (truth_r <= hi)
}

cat(sprintf("M = %d replicates, N = %d, T = %d\n\n", M, N, Tt))
summary_r <- data.frame(
  parameter = names(truth_r), truth = truth_r,
  mean_est = colMeans(est_r), empirical_sd = apply(est_r, 2, sd),
  mean_analytic_se = colMeans(se_store_r),
  coverage_95 = colMeans(cover_r)
)
print(summary_r, row.names = FALSE, digits = 4)

## ---- save Part D (regression): raw + summary tables, plots ----
M_r <- M
write.csv(summary_r, "results/simulation/tables/part_d_coverage_regression_summary.csv", row.names = FALSE)
saveRDS(list(est = est_r, se = se_store_r, cover = cover_r, truth = truth_r, M = M_r, N = N, Tt = Tt),
        "results/simulation/tables/part_d_coverage_regression_raw.rds")

plot_coverage(summary_r, M_r, "results/simulation/figures/part_d_coverage_regression.png",
              sprintf("Regression model, M = %d", M_r))
plot_est_boxplot(est_r, truth_r, "results/simulation/figures/part_d_estimates_regression.png",
                  sprintf("Regression model final estimates vs truth (dashed), M = %d", M_r))

## ---- combined SE calibration plot (analytic mean SE vs empirical SD), both models ----
png("results/simulation/figures/part_d_se_calibration.png", width = 900, height = 500, res = 130)
par(mfrow = c(1, 2))
plot(summary_h$empirical_sd, summary_h$mean_analytic_se, pch = 19, col = "steelblue",
     xlab = "empirical SD across replicates", ylab = "mean analytic SE",
     main = "Homogeneous model")
abline(0, 1, lty = 2, col = "grey40")
text(summary_h$empirical_sd, summary_h$mean_analytic_se, summary_h$parameter, pos = 3, cex = 0.7)
plot(summary_r$empirical_sd, summary_r$mean_analytic_se, pch = 19, col = "darkorange",
     xlab = "empirical SD across replicates", ylab = "mean analytic SE",
     main = "Regression model")
abline(0, 1, lty = 2, col = "grey40")
text(summary_r$empirical_sd, summary_r$mean_analytic_se, summary_r$parameter, pos = 3, cex = 0.7)
dev.off()

cat("\nSaved: results/simulation/tables/part_d_coverage_regression_summary.csv\n")
cat("Saved: results/simulation/tables/part_d_coverage_regression_raw.rds\n")
cat("Saved: results/simulation/figures/part_d_coverage_regression.png\n")
cat("Saved: results/simulation/figures/part_d_estimates_regression.png\n")
cat("Saved: results/simulation/figures/part_d_se_calibration.png\n")
