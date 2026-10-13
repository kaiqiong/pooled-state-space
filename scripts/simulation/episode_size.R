# Portable October 13 replication. Run from the repository root; see docs/SIMULATIONS.txt.
## Part J: sample-size (N) and series-length (T) sweep for the joint estimator.
## Four-covariate design of Part E; reports RMSE, bias, SE calibration and
## 95% coverage for every parameter across N in {50,100,200,400} x T in
## {15,30,60}. The N=50, T=15 corner is roughly the per-subject regime of the
## CGMacros application.
source("R/state_space_em.R")

generate_covariates <- function(N, sex_prob = 0.8, age_biomarker_cor = 0.3) {
  Sex <- rbinom(N, 1, sex_prob)
  Sigma <- matrix(c(1, age_biomarker_cor, age_biomarker_cor, 1), 2, 2)
  Z <- matrix(rnorm(N * 2), N, 2) %*% chol(Sigma)
  cbind(Intercept = 1, Age = Z[, 1], Sex = Sex, Biomarker = Z[, 2])
}
generate_data_mv <- function(Xmat, beta, Tt, mu, a1, P1, sigma_eta2, sigma_eps2) {
  N <- nrow(Xmat); phi <- tanh(as.vector(Xmat %*% beta))
  alpha <- matrix(0, N, Tt); y <- matrix(0, N, Tt)
  alpha[, 1] <- rnorm(N, a1, sqrt(P1)); y[, 1] <- mu + alpha[, 1] + rnorm(N, 0, sqrt(sigma_eps2))
  for (t in 2:Tt) {
    alpha[, t] <- phi * alpha[, t - 1] + rnorm(N, 0, sqrt(sigma_eta2))
    y[, t] <- mu + alpha[, t] + rnorm(N, 0, sqrt(sigma_eps2))
  }
  list(y = y, phi = phi)
}

dir.create("results/simulation/tables", recursive = TRUE, showWarnings = FALSE)
truth <- c(Intercept = 0.3, Age = 0.15, Sex = -0.2, Biomarker = 0.25)
a1_true <- 7; P1_true <- 15; sigma_eta2_true <- 1; sigma_eps2_true <- 1; mu <- 0
truth_full <- c(truth, a1 = a1_true, P1 = P1_true, sigma_eta2 = sigma_eta2_true, sigma_eps2 = sigma_eps2_true)
pn <- names(truth_full)
M <- 100
N_grid <- c(50, 100, 200, 400); T_grid <- c(15, 30, 60)

rows <- list(); raw <- list()
t_start <- Sys.time()
for (N in N_grid) for (Tt in T_grid) {
  set.seed(1000 + N + Tt)
  est <- matrix(NA, M, length(pn), dimnames = list(NULL, pn))
  se <- matrix(NA, M, length(pn), dimnames = list(NULL, pn))
  cov <- matrix(NA, M, length(pn), dimnames = list(NULL, pn))
  iters <- integer(M)
  t0 <- Sys.time()
  for (m in seq_len(M)) {
    Xm <- generate_covariates(N)
    dm <- generate_data_mv(Xm, truth, Tt, mu, a1_true, P1_true, sigma_eta2_true, sigma_eps2_true)
    fm <- em_regression_mv(dm$y, Xm, mu = mu,
                            init = list(beta = c(0, 0, 0, 0), a1 = mean(dm$y[, 1]),
                                        P1 = max(var(dm$y[, 1]), 1), sigma_eta2 = 0.5, sigma_eps2 = 0.5),
                            thres = 1e-6, max_iter = 2000)
    sem <- tryCatch(compute_SE_regression_mv(dm$y, Xm, mu, fm), error = function(e) NULL)
    e <- c(fm$beta, fm$a1, fm$P1, fm$sigma_eta2, fm$sigma_eps2); names(e) <- pn
    est[m, ] <- e; iters[m] <- fm$iter
    if (!is.null(sem)) {
      se[m, ] <- sem$se[pn]
      cov[m, ] <- (truth_full >= e - 1.96 * se[m, ]) & (truth_full <= e + 1.96 * se[m, ])
    }
  }
  cat(sprintf("N=%3d T=%2d done in %6.1f s (median EM iters %.0f)\n", N, Tt,
              as.numeric(Sys.time() - t0, units = "secs"), median(iters)))
  raw[[paste(N, Tt)]] <- list(est = est, se = se, cov = cov, iters = iters)
  for (j in seq_along(pn)) rows[[length(rows) + 1]] <- data.frame(
    N = N, T = Tt, parameter = pn[j], truth = truth_full[j],
    mean_est = mean(est[, j]), bias = mean(est[, j]) - truth_full[j],
    rmse = sqrt(mean((est[, j] - truth_full[j])^2)),
    emp_sd = sd(est[, j]), mean_se = mean(se[, j], na.rm = TRUE),
    coverage = mean(cov[, j], na.rm = TRUE), n_se_ok = sum(!is.na(se[, j])))
}
df <- do.call(rbind, rows); rownames(df) <- NULL
cat(sprintf("\nTotal time %.1f min\n", as.numeric(Sys.time() - t_start, units = "mins")))
print(df, row.names = FALSE, digits = 4)
write.csv(df, "results/simulation/tables/part_j_NT_sweep_summary.csv", row.names = FALSE)
saveRDS(list(raw = raw, summary = df, truth = truth_full, M = M), "results/simulation/tables/part_j_NT_sweep_raw.rds")
