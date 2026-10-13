source("R/simulation_design.R")
# Portable October 13 replication. Run from the repository root; see docs/SIMULATIONS.txt.
## Matched cap sensitivity (10 Oct 2026 audit item 3.2): persistence caps 0.99, 0.995 (production) and 0.999 applied to
## IDENTICAL datasets and IDENTICAL participant-resampling indices, with the production optimization policy (EM with the
## constrained M-step, started as in the parent script; bootstrap refits started at the dataset's own cap-specific estimate).
## Sets: (i) the application windows (primary boundary analysis, sigma_eps2 = 0), B = 2000 with the primary analysis seed;
##       (ii) prespecified simulation subset: Part M v2 boundary datasets 1:100 in both scenarios, B = 999;
##       (iii) stress set: the 10 heterogeneous boundary datasets with the largest constrained-resample fractions at cap 0.995
##             (724, 746, 898, 558, 60, 255, 52, 601, 869, 207), B = 999.
## Recorded per dataset and cap: estimate, log-likelihood, status, max |phi|, model/CR0/MD SEs, bootstrap draws and statuses.
## Usage: Rscript sim_cap_sensitivity.R <cores>
cores <- as.integer(commandArgs(trailingOnly = TRUE)[1]); source("R/state_space_em.R"); suppressPackageStartupMessages(library(parallel)); caps <- c(0.99, 0.995, 0.999); B <- 999L
f <- pss_simulation_design("positive"); Xmat <- f$Xmat; subj <- f$meta$subject; N <- nrow(Xmat); Tt <- ncol(f$Y); G <- 44; pn <- colnames(Xmat); p <- 6
fb <- pss_simulation_design("boundary")$truth
beta0 <- setNames(fb$beta, pn); a1_0 <- fb$a1; P1_0 <- fb$P1; se2_0 <- fb$sigma_eta2
gen <- function(tau) { b_g <- setNames(rnorm(G, 0, tau), unique(subj)); phi <- tanh(as.vector(Xmat %*% beta0) + b_g[subj]); alpha <- y <- matrix(0, N, Tt); alpha[, 1] <- rnorm(N, a1_0, sqrt(P1_0)); y[, 1] <- alpha[, 1] + rnorm(N, 0, 0); for (t in 2:Tt) { alpha[, t] <- phi * alpha[, t - 1] + rnorm(N, 0, sqrt(se2_0)); y[, t] <- alpha[, t] + rnorm(N, 0, 0) }; y }
fit_cap <- function(y, cap, start) em_regression_mv(y, Xmat, 0, start, thres = 1e-8, max_iter = 3000, fix_sigma_eps2 = 0, phi_cap = cap)
se_all <- function(y, fit) list(model = tryCatch(compute_SE_regression_mv_fixed_eps2(y, Xmat, 0, fit)$se[pn], error = function(e) rep(NA, p)), cr0 = tryCatch(compute_SE_regression_mv_cluster_robust(y, Xmat, 0, fit, subj)$se[pn], error = function(e) rep(NA, p)), md = tryCatch(compute_SE_regression_mv_cluster_MD(y, Xmat, 0, fit, subj)$se[pn], error = function(e) rep(NA, p)))
run_set <- function(y, start, idx_all, label) { B <- length(idx_all); lapply(caps, function(cap) { fit <- fit_cap(y, cap, start); se <- se_all(y, fit); draws <- matrix(NA_real_, B, p, dimnames = list(NULL, pn)); st <- character(B); ll <- numeric(B); it <- integer(B)
    s1 <- list(beta = fit$beta, a1 = fit$a1, P1 = fit$P1, sigma_eta2 = fit$sigma_eta2, sigma_eps2 = 0)
    for (b in seq_len(B)) { idx <- idx_all[[b]]; fbb <- tryCatch(em_regression_mv(y[idx, ], Xmat[idx, ], 0, s1, thres = 1e-6, max_iter = 3000, fix_sigma_eps2 = 0, phi_cap = cap), error = function(e) NULL)
      if (is.null(fbb)) st[b] <- "error" else { st[b] <- if (fbb$converged) fbb$status else "not_converged"; ll[b] <- fbb$loglik; it[b] <- fbb$iter; if (fbb$converged) draws[b, ] <- fbb$beta } }
    list(cap = cap, beta = fit$beta, loglik = fit$loglik, status = fit$status, max_phi = max(abs(tanh(Xmat %*% fit$beta))), iter = fit$iter, se = se, draws = draws, boot_status = st, boot_loglik = ll, boot_iter = it) }) }
out_dir <- "results/simulation/cap_sensitivity"; dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
## (i) application
app <- function() { fn <- file.path(out_dir, "application.rds"); if (file.exists(fn)) return(invisible()); y <- f$Y; set.seed(20260929); idx_all <- lapply(seq_len(2000), function(b) unlist(lapply(sample(unique(subj), G, TRUE), function(s) which(subj == s))))
  start <- list(beta = f$truth$beta, a1 = mean(y[, 1]), P1 = var(y[, 1]), sigma_eta2 = 35, sigma_eps2 = 0); res <- run_set(y, start, idx_all, "application"); saveRDS(list(set = "application", res = res), fn) }
## (ii)+(iii) simulation datasets
sim_one <- function(m, tau, set) { fn <- file.path(out_dir, sprintf("%s_tau%s_m%03d.rds", set, tau, m)); if (file.exists(fn)) return(invisible())
  set.seed(900000 + 10000 * round(10 * tau) + m); y <- gen(tau); nv <- naive_two_stage(y, Xmat); b_start <- c(atanh(min(0.98, median(pmax(pmin(nv$phi_hat, 0.99), -0.99), na.rm = TRUE))), rep(0, p - 1)); start <- list(beta = b_start, a1 = mean(y[, 1]), P1 = max(var(y[, 1]), 1), sigma_eta2 = 35, sigma_eps2 = 0)
  set.seed(900000 + 5000000 + 10000 * round(10 * tau) + m); idx_all <- lapply(seq_len(B), function(b) unlist(lapply(sample(unique(subj), G, TRUE), function(s) which(subj == s))))   # the B = 999 study's resampling stream
  res <- run_set(y, start, idx_all, set); saveRDS(list(set = set, tau = tau, m = m, res = res), fn) }
stress <- c(724, 746, 898, 558, 60, 255, 52, 601, 869, 207); jobs <- rbind(data.frame(m = 0, tau = NA, set = "application"), data.frame(m = rep(1:100, 2), tau = rep(c(0, 0.4), each = 100), set = "subset"), data.frame(m = stress, tau = 0.4, set = "stress"))
t0 <- Sys.time(); invisible(mclapply(seq_len(nrow(jobs)), function(k) if (jobs$set[k] == "application") app() else sim_one(jobs$m[k], jobs$tau[k], jobs$set[k]), mc.cores = cores, mc.preschedule = FALSE)); cat(sprintf("cap sensitivity: %d jobs in %.1f min\n", nrow(jobs), as.numeric(Sys.time() - t0, units = "mins")))
