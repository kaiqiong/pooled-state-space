source("R/simulation_design.R")
# Portable October 13 replication. Run from the repository root; see docs/SIMULATIONS.txt.
## Part P: comparison with a mixed-effects latent autoregression (mixed_latent_ar.R) on the application design
## (1174 windows, 44 participants, T = 14, no measurement error), same generating process and SEEDS as Part M
## (boundary): tau = 0 (working-independence model correct; mixed model has tau on its boundary) and tau = 0.4
## (mixed model correctly specified; working-independence model misspecified). Paired per replicate.
## Estimators: (a) proposed working-independence fit at sigma_eps2 = 0 with nuisance-aware cluster CR0 / t43;
## (b) mixed model with model-based (Hessian) SEs and normal quantiles. M = 500 per scenario.
## Rerun 2026-09-30 with the constrained M-step (EM) and the repaired mixed-model integration/Hessian; mixed-model
## bias/SD/coverage are reported over VALID fits only (optimizer converged, all cluster modes passed the score check,
## Hessian positive definite, finite SEs), with every denominator recorded explicitly.
source("R/state_space_em.R"); source("R/mixed_latent_ar.R"); suppressPackageStartupMessages(library(parallel))
f <- pss_simulation_design("positive"); pb <- pss_simulation_design("boundary"); Xmat <- f$Xmat; subj <- f$meta$subject; N <- nrow(Xmat); Tt <- ncol(f$Y); G <- length(unique(subj)); pn <- colnames(Xmat); p <- 6
beta0 <- setNames(pb$truth$beta, pn); a1_0 <- pb$truth$a1; P1_0 <- pb$truth$P1; se2_0 <- pb$truth$sigma_eta2
gen <- function(tau) { b_g <- setNames(rnorm(G, 0, tau), unique(subj)); phi <- tanh(as.vector(Xmat %*% beta0) + b_g[subj]); y <- matrix(0, N, Tt); y[, 1] <- rnorm(N, a1_0, sqrt(P1_0)); for (t in 2:Tt) y[, t] <- phi * y[, t - 1] + rnorm(N, 0, sqrt(se2_0)); y }
M <- 500; q_t <- qt(0.975, G - 1); q_n <- qnorm(0.975)
one <- function(m, tau) { set.seed(900000 + 10000 * round(10 * tau) + m); y <- gen(tau); nv <- naive_two_stage(y, Xmat); b_start <- c(atanh(min(0.98, median(pmax(pmin(nv$phi_hat, 0.99), -0.99), na.rm = TRUE))), rep(0, p - 1))
  t0 <- Sys.time(); fit <- em_regression_mv(y, Xmat, 0, list(beta = b_start, a1 = mean(y[, 1]), P1 = max(var(y[, 1]), 1), sigma_eta2 = 35, sigma_eps2 = 0), thres = 1e-8, max_iter = 3000, fix_sigma_eps2 = 0); wi_secs <- as.numeric(Sys.time() - t0, units = "secs")
  se_wi <- tryCatch(compute_SE_regression_mv_cluster_beta(y, Xmat, 0, fit, subj)$se[pn], error = function(e) rep(NA, p)); se_md <- tryCatch(compute_SE_regression_mv_cluster_MD(y, Xmat, 0, fit, subj)$se[pn], error = function(e) rep(NA, p))
  fm <- tryCatch(fit_mixed_ar(y, Xmat, subj, start_theta = b_start), error = function(e) NULL)
  if (is.null(fm)) return(list(m = m, wi = fit$beta, se_wi = se_wi, se_md = se_md, wi_conv = fit$converged, wi_status = fit$status, wi_secs = wi_secs, mixed_error = TRUE))
  list(m = m, wi = fit$beta, se_wi = se_wi, se_md = se_md, wi_conv = fit$converged, wi_status = fit$status, wi_secs = wi_secs, mixed_error = FALSE, mx = fm$theta, se_mx = fm$se_theta, tau = fm$tau, tau_floor = fm$tau_at_floor, mx_conv = fm$conv, mx_valid = fm$valid, mx_modes_ok = fm$modes_ok, mx_pd = fm$hessian_pd, mx_secs = fm$secs, mx_loglik = fm$loglik, wi_loglik = fit$loglik) }
rows <- list(); raw <- list(); rds <- "results/simulation/tables/part_p_mixed_comparator_v2_raw.rds"; if (file.exists(rds)) raw <- readRDS(rds)$raw
for (tau in c(0, 0.4)) { key <- as.character(tau); if (!(key %in% names(raw))) { t0 <- Sys.time(); raw[[key]] <- mclapply(seq_len(M), function(m) one(m, tau), mc.cores = as.integer(Sys.getenv("SIM_CORES", "1")), mc.preschedule = FALSE); cat(sprintf("tau=%.1f: %d replicates in %.1f min\n", tau, M, as.numeric(Sys.time() - t0, units = "mins"))); saveRDS(list(raw = raw, truth = beta0), rds) }
  r <- raw[[key]]; okm <- sapply(r, function(z) is.list(z) && !z$mixed_error && isTRUE(z$mx_valid)); retm <- sapply(r, function(z) is.list(z) && !z$mixed_error); W <- t(sapply(r, `[[`, "wi")); Sw <- t(sapply(r, `[[`, "se_wi")); Smd <- t(sapply(r, `[[`, "se_md")); Mx <- t(sapply(r[okm], `[[`, "mx")); Sm <- t(sapply(r[okm], `[[`, "se_mx"))
  for (j in seq_len(p)) rows[[length(rows) + 1]] <- data.frame(tau = tau, parameter = pn[j], truth = beta0[j], M = length(r), n_mixed = sum(okm),
    n_wi_se = sum(is.finite(Sw[, j])), wi_bias = mean(W[, j]) - beta0[j], wi_sd = sd(W[, j]), wi_se_cr0 = mean(Sw[, j], na.rm = TRUE), wi_cov_cr0 = mean(abs(W[, j] - beta0[j]) <= q_t * Sw[, j], na.rm = TRUE), wi_cov_md = mean(abs(W[, j] - beta0[j]) <= q_t * Smd[, j], na.rm = TRUE),
    mx_bias = mean(Mx[, j]) - beta0[j], mx_sd = sd(Mx[, j]), mx_se = mean(Sm[, j], na.rm = TRUE), mx_cov = mean(abs(Mx[, j] - beta0[j]) <= q_n * Sm[, j], na.rm = TRUE), mx_cov_t = mean(abs(Mx[, j] - beta0[j]) <= q_t * Sm[, j]), n_mx_se = sum(is.finite(Sm[, j])),
    rmse_ratio_mx_over_wi = sqrt(mean((Mx[, j] - beta0[j])^2)) / sqrt(mean((W[okm, j] - beta0[j])^2)),
    tau_mean = mean(sapply(r[okm], `[[`, "tau")), tau_at_floor = sum(sapply(r[okm], `[[`, "tau_floor")), mx_returned = sum(retm), mx_nonconv = sum(!sapply(r[retm], `[[`, "mx_conv")), mx_modes_fail = sum(!sapply(r[retm], `[[`, "mx_modes_ok")), mx_hessian_fail = sum(!sapply(r[retm], `[[`, "mx_pd")), mx_valid = sum(okm), mixed_errors = sum(!retm), wi_nonconv = sum(!sapply(r, `[[`, "wi_conv")), wi_constrained = sum(sapply(r, `[[`, "wi_status") == "constrained"), wi_stalled = sum(sapply(r, `[[`, "wi_status") == "stalled"),
    wi_median_secs = median(sapply(r, `[[`, "wi_secs")), mx_median_secs = median(sapply(r[okm], `[[`, "mx_secs")), lrt_tau_reject = mean(2 * (sapply(r[okm], `[[`, "mx_loglik") - sapply(r[okm], `[[`, "wi_loglik")) > qchisq(0.9, 1))) }
df <- do.call(rbind, rows); rownames(df) <- NULL; options(width = 250); print(df[, c("tau", "parameter", "wi_bias", "wi_sd", "wi_se_cr0", "wi_cov_cr0", "wi_cov_md", "mx_bias", "mx_sd", "mx_se", "mx_cov", "mx_cov_t", "rmse_ratio_mx_over_wi")], digits = 3, row.names = FALSE)
print(unique(df[, c("tau", "M", "n_mixed", "mx_returned", "mx_nonconv", "mx_modes_fail", "mx_hessian_fail", "mx_valid", "mixed_errors", "tau_mean", "tau_at_floor", "wi_nonconv", "wi_constrained", "wi_stalled", "wi_median_secs", "mx_median_secs", "lrt_tau_reject")]), digits = 3, row.names = FALSE)
write.csv(df, "results/simulation/tables/part_p_mixed_comparator_v2_summary.csv", row.names = FALSE); cat("saved\n")
