source("R/simulation_design.R")
# Portable October 13 replication. See docs/SIMULATIONS.txt.
## Rerun 2026-09-30 (v2) with the constrained M-step; a constrained fit has a nonzero raw nuisance score by construction, so the nuisance-score check is reported over unconstrained fits.
## Step 57: size and power of the one-sided boundary score test of Step 55 by simulation on the application
## design (1174 windows, 44 participants, T = 14), truth = boundary fit; white noise SD 0 (size), 1, 2, 3 mg/dL
## (power), M = 500 each; also SD 0 with a participant random effect tau = 0.4 on the link scale (size under
## within-participant dependence, which the clustered variance is meant to absorb). Analytic per-unit scores
## at r = 0; the nuisance projection uses the observed information for theta (analytic, from the library) and
## d(score_r)/d(theta) by central differences of the closed-form score.
source("R/state_space_em.R"); suppressPackageStartupMessages(library(parallel))
pb <- pss_simulation_design("boundary"); X <- pb$Xmat; cl <- pb$meta$subject; gi <- match(cl, unique(cl)); G <- length(unique(cl)); N <- nrow(X); Tt <- ncol(pb$Y); fit <- pb$truth
score_r_units <- function(th, Y) { phi <- as.vector(tanh(X %*% th[1:6])); a1 <- th[7]; P1 <- th[8]; q <- th[9]; m <- outer(phi, 0:(Tt - 1), function(p, k) a1 * p^k); d <- Y - m
  u <- cbind(d[, 1], d[, -1] - phi * d[, -Tt]); w <- u / rep(c(P1, rep(q, Tt - 1)), each = N); v <- cbind(w[, -Tt] - phi * w[, -1], w[, Tt]); (rowSums(v^2) - 1 / P1 - (Tt - 1) * (1 + phi^2) / q) / 2 }
score_th_units <- function(th, Y) { phi <- as.vector(tanh(X %*% th[1:6])); a1 <- th[7]; P1 <- th[8]; q <- th[9]; e <- Y[, -1] - phi * Y[, -Tt]; sb <- X * ((1 - phi^2) * rowSums(e * Y[, -Tt]) / q)
  cbind(sb, (Y[, 1] - a1) / P1, -1 / (2 * P1) + (Y[, 1] - a1)^2 / (2 * P1^2), -(Tt - 1) / (2 * q) + rowSums(e^2) / (2 * q^2)) }
stat <- function(Y) { f0 <- em_regression_mv(Y, X, 0, list(beta = fit$beta, a1 = mean(Y[, 1]), P1 = var(Y[, 1]), sigma_eta2 = 35, sigma_eps2 = 0), thres = 1e-9, max_iter = 2000, fix_sigma_eps2 = 0); th <- c(f0$beta, f0$a1, f0$P1, f0$sigma_eta2)
  sr <- score_r_units(th, Y); st <- score_th_units(th, Y); Ith <- compute_SE_regression_mv_fixed_eps2(Y, X, 0, f0)$I_obs
  Irt <- sapply(seq_along(th), function(j) { h <- 1e-5 * max(1, abs(th[j])); tp <- th; tp[j] <- tp[j] + h; tm <- th; tm[j] <- tm[j] - h; -(sum(score_r_units(tp, Y)) - sum(score_r_units(tm, Y))) / (2 * h) })
  adj <- as.vector(Irt %*% solve(Ith)); se <- sr - as.vector(st %*% adj); Sg <- rowsum(se, gi)[, 1]; c(score = sum(se), stat = sum(se) / sqrt(sum(Sg^2)), conv = f0$converged, constrained = identical(f0$status, "constrained"), max_score_th = max(abs(colSums(st)))) }
gen <- function(sd, tau) { b <- setNames(rnorm(G, 0, tau), unique(cl)); phi <- tanh(as.vector(X %*% fit$beta) + b[cl]); Y <- matrix(0, N, Tt); Y[, 1] <- rnorm(N, fit$a1, sqrt(fit$P1)); for (t in 2:Tt) Y[, t] <- phi * Y[, t - 1] + rnorm(N, 0, sqrt(fit$sigma_eta2)); Y + rnorm(N * Tt, 0, sd) }
M <- 500; scen <- list(c(0, 0), c(1, 0), c(2, 0), c(3, 0), c(0, 0.4)); out <- list()
for (sc in scen) { res <- do.call(rbind, mclapply(seq_len(M), function(m) { set.seed(570000 + 1000 * sc[1] + 100 * sc[2] + m); stat(gen(sc[1], sc[2])) }, mc.cores = as.integer(Sys.getenv("SIM_CORES", "1"))))
  out[[length(out) + 1]] <- data.frame(sd_eps = sc[1], tau = sc[2], M = M, reject_5pct = mean(res[, "stat"] > qnorm(0.95)), reject_1pct = mean(res[, "stat"] > qnorm(0.99)), mean_stat = mean(res[, "stat"]), sd_stat = sd(res[, "stat"]), q025 = quantile(res[, "stat"], .025), q975 = quantile(res[, "stat"], .975), nonconv = sum(res[, "conv"] == 0), constrained = sum(res[, "constrained"] == 1), max_nuisance_score_unconstrained = max(res[res[, "constrained"] == 0, "max_score_th"]))
  print(out[[length(out)]], digits = 3, row.names = FALSE) }
tab <- do.call(rbind, out); write.csv(tab, "results/simulation/tables/score_test_size_power_v2.csv", row.names = FALSE); cat("saved\n")
