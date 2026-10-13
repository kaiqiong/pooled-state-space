source("R/simulation_design.R")
# Portable October 13 replication. Run from the repository root; see docs/SIMULATIONS.txt.
source("R/state_space_em.R"); suppressPackageStartupMessages(library(parallel))
f <- pss_simulation_design("positive")
f <- pss_simulation_design("boundary"); f$sd_eps_primary <- 0
Xmat <- f$Xmat; subj <- f$meta$subject; N <- nrow(Xmat); Tt <- ncol(f$Y); G <- length(unique(subj)); pn <- colnames(Xmat); p <- length(pn)
beta0 <- setNames(f$truth$beta, pn); a1_0 <- f$truth$a1; P1_0 <- f$truth$P1; se2_0 <- f$truth$sigma_eta2; sp2_0 <- f$sd_eps_primary^2
cat(sprintf("N=%d G=%d T=%d | beta=(%s) a1=%.1f P1=%.0f sigma_eta2=%.1f sigma_eps2=%.0f | mean phi %.4f\n", N, G, Tt, paste(round(beta0, 3), collapse = ", "), a1_0, P1_0, se2_0, sp2_0, mean(tanh(Xmat %*% beta0))))
gen <- function(tau) { b_g <- setNames(rnorm(G, 0, tau), unique(subj)); phi <- tanh(as.vector(Xmat %*% beta0) + b_g[subj]); alpha <- y <- matrix(0, N, Tt)
  alpha[, 1] <- rnorm(N, a1_0, sqrt(P1_0)); y[, 1] <- alpha[, 1] + rnorm(N, 0, sqrt(sp2_0))
  for (t in 2:Tt) { alpha[, t] <- phi * alpha[, t - 1] + rnorm(N, 0, sqrt(se2_0)); y[, t] <- alpha[, t] + rnorm(N, 0, sqrt(sp2_0)) }; y }
try_se <- function(expr) { r <- tryCatch(expr, error = function(e) e); if (inherits(r, "error")) list(se = rep(NA_real_, p), err = conditionMessage(r)) else list(se = as.numeric(r[pn]), err = "") }
M <- 1000; M_boot <- 0; B_boot <- 199; q_t <- qt(0.975, G - 1); q_n <- qnorm(0.975)
one_rep <- function(m, tau) {
  set.seed(900000 + 10000 * round(10 * tau) + m); y <- gen(tau); nv <- naive_two_stage(y, Xmat)
  b_start <- c(atanh(min(0.98, median(pmax(pmin(nv$phi_hat, 0.99), -0.99), na.rm = TRUE))), rep(0, p - 1))
  fit <- tryCatch(em_regression_mv(y, Xmat, 0, list(beta = b_start, a1 = mean(y[, 1]), P1 = max(var(y[, 1]), 1), sigma_eta2 = 35, sigma_eps2 = sp2_0), thres = 1e-8, max_iter = 3000, fix_sigma_eps2 = sp2_0), error = function(e) e)
  if (inherits(fit, "error")) return(list(m = m, fit_error = conditionMessage(fit)))
  se <- list(model = try_se(compute_SE_regression_mv_fixed_eps2(y, Xmat, 0, fit)$se),
             cr0 = try_se(compute_SE_regression_mv_cluster_robust(y, Xmat, 0, fit, subj)$se),
             md = try_se(compute_SE_regression_mv_cluster_MD(y, Xmat, 0, fit, subj)$se))
  out <- list(m = m, fit_error = "", beta = fit$beta, naive = as.numeric(nv$beta_hat), converged = fit$converged, status = fit$status, iter = fit$iter, at_cap = isTRUE(fit$at_cap), n_backtrack = fit$n_backtrack,
              max_phi = max(abs(tanh(Xmat %*% fit$beta))), se = lapply(se, `[[`, "se"), se_err = sapply(se, `[[`, "err"))

  out }
dir.create("results/simulation/tables", showWarnings = FALSE, recursive = TRUE); rds <- "results/simulation/tables/part_m_boundary_v2_raw.rds"; raw <- if (file.exists(rds)) readRDS(rds)$raw else list()
for (tau in c(0, 0.4)) { key <- as.character(tau); if (key %in% names(raw)) { stopifnot(length(raw[[key]]) == M); next }; t0 <- Sys.time()
  raw[[key]] <- mclapply(seq_len(M), function(m) one_rep(m, tau), mc.cores = as.integer(Sys.getenv("SIM_CORES", "1")), mc.preschedule = FALSE)
  cat(sprintf("tau=%.1f: %d replicates in %.1f min\n", tau, M, as.numeric(Sys.time() - t0, units = "mins"))); saveRDS(list(raw = raw, truth = beta0, M = M, M_boot = M_boot, B_boot = B_boot, G = G, N = N), rds) }
