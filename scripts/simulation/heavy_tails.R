# Portable October 13 replication. Run from the repository root; see docs/SIMULATIONS.txt.
## Part L: robustness to heavy-tailed measurement error. Same design as Part E
## but epsilon_it ~ t_3 scaled to variance sigma_eps2 (t_3 has variance 3).
## Gaussian-likelihood EM is misspecified here; reports bias/RMSE/coverage
## for direct comparison with Part E's Gaussian results.
source("R/state_space_em.R")

generate_covariates <- function(N, sex_prob = 0.8, age_biomarker_cor = 0.3) {
  Sex <- rbinom(N, 1, sex_prob)
  Sigma <- matrix(c(1, age_biomarker_cor, age_biomarker_cor, 1), 2, 2)
  Z <- matrix(rnorm(N * 2), N, 2) %*% chol(Sigma)
  cbind(Intercept = 1, Age = Z[, 1], Sex = Sex, Biomarker = Z[, 2])
}
generate_data_t <- function(Xmat, beta, Tt, mu, a1, P1, sigma_eta2, sigma_eps2, df_t = 3) {
  N <- nrow(Xmat); phi <- tanh(as.vector(Xmat %*% beta))
  sc <- sqrt(sigma_eps2 / (df_t / (df_t - 2)))
  alpha <- matrix(0, N, Tt); y <- matrix(0, N, Tt)
  alpha[, 1] <- rnorm(N, a1, sqrt(P1)); y[, 1] <- mu + alpha[, 1] + sc * rt(N, df_t)
  for (t in 2:Tt) {
    alpha[, t] <- phi * alpha[, t - 1] + rnorm(N, 0, sqrt(sigma_eta2))
    y[, t] <- mu + alpha[, t] + sc * rt(N, df_t)
  }
  list(y = y, phi = phi)
}

dir.create("results/simulation/tables", recursive = TRUE, showWarnings = FALSE)
N <- 150; Tt <- 30; mu <- 0; M <- 150
truth <- c(Intercept = 0.3, Age = 0.15, Sex = -0.2, Biomarker = 0.25)
a1_true <- 7; P1_true <- 15; sigma_eta2_true <- 1; sigma_eps2_true <- 1
truth_full <- c(truth, a1 = a1_true, P1 = P1_true, sigma_eta2 = sigma_eta2_true, sigma_eps2 = sigma_eps2_true)
pn <- names(truth_full)
set.seed(3003)
est <- matrix(NA, M, length(pn), dimnames = list(NULL, pn))
se <- matrix(NA, M, length(pn), dimnames = list(NULL, pn))
cov <- matrix(NA, M, length(pn), dimnames = list(NULL, pn))
t0 <- Sys.time()
for (m in seq_len(M)) {
  Xm <- generate_covariates(N)
  dm <- generate_data_t(Xm, truth, Tt, mu, a1_true, P1_true, sigma_eta2_true, sigma_eps2_true)
  fm <- em_regression_mv(dm$y, Xm, mu = mu,
                          init = list(beta = c(0, 0, 0, 0), a1 = mean(dm$y[, 1]),
                                      P1 = max(var(dm$y[, 1]), 1), sigma_eta2 = 0.5, sigma_eps2 = 0.5),
                          thres = 1e-6, max_iter = 2000)
  sem <- tryCatch(compute_SE_regression_mv(dm$y, Xm, mu, fm), error = function(e) NULL)
  e <- c(fm$beta, fm$a1, fm$P1, fm$sigma_eta2, fm$sigma_eps2); names(e) <- pn
  est[m, ] <- e
  if (!is.null(sem)) { se[m, ] <- sem$se[pn]; cov[m, ] <- (truth_full >= e - 1.96 * se[m, ]) & (truth_full <= e + 1.96 * se[m, ]) }
}
cat(sprintf("done in %.1f s\n", as.numeric(Sys.time() - t0, units = "secs")))
df <- data.frame(parameter = pn, truth = truth_full, mean_est = colMeans(est),
                 bias = colMeans(est) - truth_full, rmse = sqrt(colMeans((est - rep(truth_full, each = M))^2)),
                 emp_sd = apply(est, 2, sd), mean_se = colMeans(se, na.rm = TRUE), coverage = colMeans(cov, na.rm = TRUE))
print(df, row.names = FALSE, digits = 4)
write.csv(df, "results/simulation/tables/part_l_heavytail_summary.csv", row.names = FALSE)
saveRDS(list(est = est, se = se, cov = cov, truth = truth_full, M = M), "results/simulation/tables/part_l_heavytail_raw.rds")
