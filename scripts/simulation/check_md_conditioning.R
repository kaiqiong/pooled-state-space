source("R/simulation_design.R")
# Portable October 13 replication. Run from the repository root; see docs/SIMULATIONS.txt.
## Conditioning diagnostic for the Mancl-DeRouen-style adjustment (10 Oct 2026 audit item 2.3): for every heterogeneous
## boundary dataset, the largest real eigenvalue of H_g = I_g I^{-1} over participants (near 1 => Id - H_g nearly singular),
## the participant attaining it, its window count, and the MD SE; relates extreme MD SEs to near-singular corrections.
source("R/state_space_em.R"); suppressPackageStartupMessages(library(parallel))
f <- pss_simulation_design("positive"); Xmat <- f$Xmat; subj <- f$meta$subject; N <- nrow(Xmat); Tt <- ncol(f$Y); G <- 44; pn <- colnames(Xmat)
fb <- pss_simulation_design("boundary")$truth
beta0 <- fb$beta; a1_0 <- fb$a1; P1_0 <- fb$P1; se2_0 <- fb$sigma_eta2
gen <- function(tau) { b_g <- setNames(rnorm(G, 0, tau), unique(subj)); phi <- tanh(as.vector(Xmat %*% beta0) + b_g[subj]); alpha <- y <- matrix(0, N, Tt); alpha[, 1] <- rnorm(N, a1_0, sqrt(P1_0)); y[, 1] <- alpha[, 1] + rnorm(N, 0, 0); for (t in 2:Tt) { alpha[, t] <- phi * alpha[, t - 1] + rnorm(N, 0, sqrt(se2_0)); y[, t] <- alpha[, t] + rnorm(N, 0, 0) }; y }   # the zero-SD draws consume the RNG exactly as in the parent script
raw <- readRDS("results/simulation/tables/part_m_boundary_v2_raw.rds")$raw[["0.4"]]
one <- function(m) { set.seed(900000 + 10000 * round(10 * 0.4) + m); y <- gen(0.4); z <- raw[[m]]; fit <- em_regression_mv(y, Xmat, 0, list(beta = z$beta, a1 = mean(y[, 1]), P1 = max(var(y[, 1]), 1), sigma_eta2 = 35, sigma_eps2 = 0), thres = 1e-8, max_iter = 3000, fix_sigma_eps2 = 0)
  Ii <- compute_info_per_unit_fixed_eps2(y, Xmat, 0, fit); I_tot <- apply(Ii, c(2, 3), sum); I_inv <- solve(I_tot); ev <- sapply(unique(subj), function(g) { Ig <- apply(Ii[which(subj == g), , , drop = FALSE], c(2, 3), sum); max(Re(eigen(Ig %*% I_inv, only.values = TRUE)$values)) })
  md <- compute_SE_regression_mv_cluster_MD(y, Xmat, 0, fit, subj); mineig_I <- min(eigen(I_tot, symmetric = TRUE, only.values = TRUE)$values)
  data.frame(m = m, max_H_eig = max(ev), participant = unique(subj)[which.max(ev)], n_windows = sum(subj == unique(subj)[which.max(ev)]), n_H_eig_above_0.9 = sum(ev > 0.9), md_se_max = max(md$se[pn]), cr0_se_max = max(md$se_cr0[pn]), at_cap = isTRUE(fit$at_cap), min_eig_I = mineig_I, md_check = max(abs(md$se[pn] - z$se$md))) }
t0 <- Sys.time(); res <- do.call(rbind, mclapply(seq_along(raw), one, mc.cores = as.integer(Sys.getenv("SIM_CORES", "1")))); cat(sprintf("%d datasets in %.1f min; MD SEs reproduced within %.1e\n", nrow(res), as.numeric(Sys.time() - t0, units = "mins"), max(res$md_check)))
saveRDS(res, "results/simulation/tables/md_conditioning_boundary_tau0.4.rds"); write.csv(res, "results/simulation/tables/md_conditioning_boundary_tau0.4.csv", row.names = FALSE)
cat("Spearman corr(max H eigenvalue, max MD SE):", round(cor(res$max_H_eig, res$md_se_max, method = "spearman"), 3), "\n")
cat("max H eigenvalue: quantiles", round(quantile(res$max_H_eig, c(.5, .9, .99, 1)), 3), "\n")
big <- res[res$md_se_max > 1, ]; cat(sprintf("datasets with MD SE > 1: %d; their max H eigenvalue range %.3f-%.3f (others: median %.3f, max %.3f); cap-active among them: %d; datasets with MD SE > 5: max H eig %s\n", nrow(big), min(big$max_H_eig), max(big$max_H_eig), median(res$max_H_eig[res$md_se_max <= 1]), max(res$max_H_eig[res$md_se_max <= 1]), sum(big$at_cap), paste(round(res$max_H_eig[res$md_se_max > 5], 3), collapse = ",")))
print(res[order(-res$md_se_max), ][1:10, ], digits = 3)
print(res[res$m %in% c(724, 818), ], digits = 3)
