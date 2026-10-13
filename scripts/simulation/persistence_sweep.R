# Portable October 13 replication. Run from the repository root; see docs/SIMULATIONS.txt.
# Persistence-cap sweep; optional argument 'unguarded' disables the cap for
# the documented link-saturation diagnostic with otherwise identical fitting.
source("R/state_space_em.R")
args <- commandArgs(trailingOnly = TRUE)
unguarded <- length(args) && args[1] == "unguarded"
fit_cap <- if (unguarded) 1 else 0.995
prefix <- if (unguarded) "part_k_phi_sweep_unsafeguarded" else "part_k_v2_phi_sweep"

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
N <- 150; Tt <- 30; mu <- 0; M <- 100
a1_true <- 7; P1_true <- 15; sigma_eta2_true <- 1; sigma_eps2_true <- 1
cov_effects <- c(Age = 0.10, Sex = -0.10, Biomarker = 0.15)
phi_targets <- c(0.3, 0.6, 0.9, 0.97)

rows <- list(); raw <- list()
t_start <- Sys.time()
for (phi_t in phi_targets) {
  truth <- c(Intercept = atanh(phi_t), cov_effects)
  truth_full <- c(truth, a1 = a1_true, P1 = P1_true, sigma_eta2 = sigma_eta2_true, sigma_eps2 = sigma_eps2_true)
  pn <- names(truth_full)
  set.seed(2000 + round(100 * phi_t))
  est <- matrix(NA, M, length(pn), dimnames = list(NULL, pn))
  se <- matrix(NA, M, length(pn), dimnames = list(NULL, pn))
  cov <- matrix(NA, M, length(pn), dimnames = list(NULL, pn))
  pileup <- logical(M); maxphi <- numeric(M); iters <- integer(M); conv <- logical(M); nbt <- integer(M); status <- character(M)
  t0 <- Sys.time()
  for (m in seq_len(M)) {
    Xm <- generate_covariates(N)
    dm <- generate_data_mv(Xm, truth, Tt, mu, a1_true, P1_true, sigma_eta2_true, sigma_eps2_true)
    fm <- em_regression_mv(dm$y, Xm, mu = mu,
                            init = list(beta = c(atanh(0.5), 0, 0, 0), a1 = mean(dm$y[, 1]),
                                        P1 = max(var(dm$y[, 1]), 1), sigma_eta2 = 0.5, sigma_eps2 = 0.5),
                            thres = 1e-6, max_iter = 1500, phi_cap = fit_cap)
    phi_hat <- tanh(as.vector(Xm %*% fm$beta))
    maxphi[m] <- max(abs(phi_hat)); iters[m] <- fm$iter; conv[m] <- fm$converged
    ## "pile-up": sigma_eps2 driven to ~0, or phi_i stuck at the M-step cap (0.995);
    ## the safeguarded EM can no longer return phi_i = 1 exactly
    pileup[m] <- (fm$sigma_eps2 < 0.05 * sigma_eps2_true) || isTRUE(fm$at_cap)
    nbt[m] <- fm$n_backtrack; status[m] <- fm$status
    e <- c(fm$beta, fm$a1, fm$P1, fm$sigma_eta2, fm$sigma_eps2); names(e) <- pn
    est[m, ] <- e
    if (!pileup[m]) {
      sem <- tryCatch(compute_SE_regression_mv(dm$y, Xm, mu, fm), error = function(e) NULL)
      if (!is.null(sem)) {
        se[m, ] <- sem$se[pn]
        cov[m, ] <- (truth_full >= e - 1.96 * se[m, ]) & (truth_full <= e + 1.96 * se[m, ])
      }
    }
  }
  cat(sprintf("phi=%.2f done in %6.1f s | pile-up rate %.2f | converged %.2f | median iters %.0f | replicates needing M-step backtracking %.2f (mean %.1f backtracks)\n",
              phi_t, as.numeric(Sys.time() - t0, units = "secs"), mean(pileup), mean(conv), median(iters), mean(nbt > 0), mean(nbt)))
  raw[[as.character(phi_t)]] <- list(est = est, se = se, cov = cov, pileup = pileup, maxphi = maxphi,
                                     iters = iters, conv = conv, nbt = nbt, status = status, truth = truth_full)
  for (j in seq_along(pn)) rows[[length(rows) + 1]] <- data.frame(
    phi_target = phi_t, parameter = pn[j], truth = truth_full[j],
    pileup_rate = mean(pileup), conv_rate = mean(conv), backtrack_rate = mean(nbt > 0), constrained_rate = mean(status == "constrained"), stalled_rate = mean(status == "stalled"),
    mean_est_all = mean(est[, j]),
    mean_est_nonpileup = mean(est[!pileup, j]),
    bias_nonpileup = mean(est[!pileup, j]) - truth_full[j],
    rmse_nonpileup = sqrt(mean((est[!pileup, j] - truth_full[j])^2)),
    emp_sd_nonpileup = sd(est[!pileup, j]),
    mean_se = mean(se[, j], na.rm = TRUE), coverage = mean(cov[, j], na.rm = TRUE),
    n_nonpileup = sum(!pileup))
}
df <- do.call(rbind, rows); rownames(df) <- NULL
cat(sprintf("\nTotal time %.1f min\n", as.numeric(Sys.time() - t_start, units = "mins")))
print(df, row.names = FALSE, digits = 4)
write.csv(df, file.path("results/simulation/tables", paste0(prefix, "_summary.csv")), row.names = FALSE)
saveRDS(list(raw = raw, summary = df, M = M, N = N, Tt = Tt), file.path("results/simulation/tables", paste0(prefix, "_raw.rds")))
