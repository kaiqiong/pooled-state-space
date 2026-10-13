# Portable October 13 replication. Run from the repository root; see docs/SIMULATIONS.txt.
source("R/state_space_em.R")

## ---------------------------------------------------------------------------
## Data generators (vectorized)
## ---------------------------------------------------------------------------
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

cat("=======================================================================\n")
cat("PART A: homogeneous model -- analytic vs numeric I_obs = D + X\n")
cat("=======================================================================\n")

set.seed(1)
N <- 60; Tt <- 20; mu <- 0
truth_h <- list(phi = 0.6, a1 = 7, P1 = 15, sigma_eta2 = 1, sigma_eps2 = 1)
dat_h <- generate_data_homogeneous(N, Tt, truth_h$phi, mu, truth_h$a1, truth_h$P1,
                                    truth_h$sigma_eta2, truth_h$sigma_eps2)

fit_h <- em_homogeneous(dat_h$y, mu = mu,
                         init = list(phi = 0.3, a1 = mean(dat_h$y[, 1]), P1 = var(dat_h$y[, 1]),
                                     sigma_eta2 = 0.5, sigma_eps2 = 0.5))
cat(sprintf("converged: %s in %d iterations\n", fit_h$converged, fit_h$iter))
cat("estimates:\n"); print(unlist(fit_h[c("phi", "a1", "P1", "sigma_eta2", "sigma_eps2")]))

se_h <- compute_SE_homogeneous(dat_h$y, mu, fit_h)

theta_hat_h <- c(fit_h$phi, fit_h$a1, fit_h$P1, fit_h$sigma_eta2, fit_h$sigma_eps2)
names(theta_hat_h) <- c("phi", "a1", "P1", "sigma_eta2", "sigma_eps2")
num_h <- numeric_I_obs(function(th, th_star) Q_value_homogeneous(th, th_star, dat_h$y, mu), theta_hat_h)
dimnames(num_h$D) <- dimnames(se_h$D); dimnames(num_h$X) <- dimnames(se_h$X); dimnames(num_h$I_obs) <- dimnames(se_h$I_obs)

cat("\nmax|D_analytic - D_numeric| =", max(abs(se_h$D - num_h$D)), "\n")
cat("max|X_analytic - X_numeric| =", max(abs(se_h$X - num_h$X)), "\n")
cat("max|I_obs_analytic - I_obs_numeric| =", max(abs(se_h$I_obs - num_h$I_obs)), "\n")
cat("I_obs symmetry check, max|I_obs - t(I_obs)| =", max(abs(se_h$I_obs - t(se_h$I_obs))), "\n\n")

se_num_h <- sqrt(diag(solve(num_h$I_obs)))
cat("SE comparison (analytic vs numeric Hessian):\n")
print(round(rbind(analytic = se_h$se, numeric = se_num_h), 5))

dir.create("results/simulation/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/simulation/figures", recursive = TRUE, showWarnings = FALSE)

se_compare_h <- data.frame(
  model = "homogeneous", parameter = names(se_h$se),
  estimate = theta_hat_h, se_analytic = se_h$se, se_numeric = se_num_h,
  abs_diff = abs(se_h$se - se_num_h)
)
diag_h <- data.frame(
  model = "homogeneous",
  max_abs_diff_D = max(abs(se_h$D - num_h$D)),
  max_abs_diff_X = max(abs(se_h$X - num_h$X)),
  max_abs_diff_Iobs = max(abs(se_h$I_obs - num_h$I_obs)),
  Iobs_symmetry_max_abs = max(abs(se_h$I_obs - t(se_h$I_obs))),
  converged = fit_h$converged, n_iter = fit_h$iter
)

cat("\n=======================================================================\n")
cat("PART B: regression model (phi_i = tanh(beta0+beta1*x_i)) -- analytic vs numeric\n")
cat("=======================================================================\n")

set.seed(2)
truth_r <- list(beta0 = 0.4, beta1 = -0.2, a1 = 7, P1 = 15, sigma_eta2 = 1, sigma_eps2 = 1)
dat_r <- generate_data_regression(N, Tt, truth_r$beta0, truth_r$beta1, mu, truth_r$a1, truth_r$P1,
                                   truth_r$sigma_eta2, truth_r$sigma_eps2)

fit_r <- em_regression(dat_r$y, dat_r$x, mu = mu,
                        init = list(beta0 = 0.1, beta1 = 0.1, a1 = mean(dat_r$y[, 1]),
                                    P1 = var(dat_r$y[, 1]), sigma_eta2 = 0.5, sigma_eps2 = 0.5),
                        thres = 1e-7, max_iter = 5000)
cat(sprintf("converged: %s in %d iterations\n", fit_r$converged, fit_r$iter))
cat("estimates:\n"); print(unlist(fit_r[c("beta0", "beta1", "a1", "P1", "sigma_eta2", "sigma_eps2")]))

se_r <- compute_SE_regression(dat_r$y, dat_r$x, mu, fit_r)

theta_hat_r <- c(fit_r$beta0, fit_r$beta1, fit_r$a1, fit_r$P1, fit_r$sigma_eta2, fit_r$sigma_eps2)
names(theta_hat_r) <- c("beta0", "beta1", "a1", "P1", "sigma_eta2", "sigma_eps2")
num_r <- numeric_I_obs(function(th, th_star) Q_value_regression(th, th_star, dat_r$y, dat_r$x, mu), theta_hat_r)
dimnames(num_r$D) <- dimnames(se_r$D); dimnames(num_r$X) <- dimnames(se_r$X); dimnames(num_r$I_obs) <- dimnames(se_r$I_obs)

cat("\nmax|D_analytic - D_numeric| =", max(abs(se_r$D - num_r$D)), "\n")
cat("max|X_analytic - X_numeric| =", max(abs(se_r$X - num_r$X)), "\n")
cat("max|I_obs_analytic - I_obs_numeric| =", max(abs(se_r$I_obs - num_r$I_obs)), "\n")
cat("I_obs symmetry check, max|I_obs - t(I_obs)| =", max(abs(se_r$I_obs - t(se_r$I_obs))), "\n\n")

se_num_r <- sqrt(diag(solve(num_r$I_obs)))
cat("SE comparison (analytic vs numeric Hessian):\n")
print(round(rbind(analytic = se_r$se, numeric = se_num_r), 5))

se_compare_r <- data.frame(
  model = "regression", parameter = names(se_r$se),
  estimate = theta_hat_r, se_analytic = se_r$se, se_numeric = se_num_r,
  abs_diff = abs(se_r$se - se_num_r)
)
diag_r <- data.frame(
  model = "regression",
  max_abs_diff_D = max(abs(se_r$D - num_r$D)),
  max_abs_diff_X = max(abs(se_r$X - num_r$X)),
  max_abs_diff_Iobs = max(abs(se_r$I_obs - num_r$I_obs)),
  Iobs_symmetry_max_abs = max(abs(se_r$I_obs - t(se_r$I_obs))),
  converged = fit_r$converged, n_iter = fit_r$iter
)

## ---- save tables ----
se_compare_all <- rbind(se_compare_h, se_compare_r)
diag_all <- rbind(diag_h, diag_r)
write.csv(se_compare_all, "results/simulation/tables/part_a_se_analytic_vs_numeric.csv", row.names = FALSE)
write.csv(diag_all, "results/simulation/tables/part_a_hessian_diagnostics.csv", row.names = FALSE)
saveRDS(list(se_compare = se_compare_all, diagnostics = diag_all,
             homogeneous = list(fit = fit_h, se = se_h, numeric = num_h, data = dat_h),
             regression = list(fit = fit_r, se = se_r, numeric = num_r, data = dat_r)),
        "results/simulation/tables/part_a_full_objects.rds")

## ---- save a diagnostic plot: analytic vs numeric SE agreement ----
png("results/simulation/figures/part_a_analytic_vs_numeric_SE.png", width = 1000, height = 500, res = 130)
par(mfrow = c(1, 2))
plot(se_compare_h$se_numeric, se_compare_h$se_analytic, pch = 19, col = "steelblue",
     xlab = "numeric SE (numDeriv Hessian)", ylab = "analytic SE (Oakes)",
     main = "Homogeneous model")
abline(0, 1, lty = 2, col = "grey40")
text(se_compare_h$se_numeric, se_compare_h$se_analytic, se_compare_h$parameter, pos = 3, cex = 0.7)
plot(se_compare_r$se_numeric, se_compare_r$se_analytic, pch = 19, col = "darkorange",
     xlab = "numeric SE (numDeriv Hessian)", ylab = "analytic SE (Oakes)",
     main = "Regression model")
abline(0, 1, lty = 2, col = "grey40")
text(se_compare_r$se_numeric, se_compare_r$se_analytic, se_compare_r$parameter, pos = 3, cex = 0.7)
dev.off()

cat("\nSaved: results/simulation/tables/part_a_se_analytic_vs_numeric.csv\n")
cat("Saved: results/simulation/tables/part_a_hessian_diagnostics.csv\n")
cat("Saved: results/simulation/tables/part_a_full_objects.rds\n")
cat("Saved: results/simulation/figures/part_a_analytic_vs_numeric_SE.png\n")
