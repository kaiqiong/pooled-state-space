source("R/simulation_design.R")
# Portable October 13 replication. Run from the repository root; see docs/SIMULATIONS.txt.
## Part R (10 Oct 2026, audit item 3.3): participant-clustered inference with a FREELY ESTIMATED positive observation variance.
## Same datasets as Part Q (identical seeds: application design, truth = SD-20 fit, sigma_eps2 = 400 generated; tau in {0, 0.4})
## but sigma_eps2 is now ESTIMATED by EM (10 free parameters). Analytic intervals: model-based, nuisance-aware CR0 and MD-style
## from the per-episode free-variance information (free_variance_cluster_inference.R; validated against a numerical Hessian in
## pilot_free_variance.R). Target under tau = 0.4: the quadrature working target with r free (expected_likelihood_target_free_r.rds,
## interior: r* = 384, max phi 0.962). Bootstrap: B = 999 on the first 200 datasets with the SAME resampling stream as the Part Q
## B = 999 study, every refit estimating the same 10 parameters, started at the dataset's estimate; retry rule as in
## sim_bootstrap_b999.R; per-refit convergence, status, iterations and sigma_eps2 recorded. Output: results/simulation/freevar/tau<tau>/m<m>.rds.
## Usage: Rscript sim_clustered_freevar.R <M_analytic> <M_boot> <cores>
args <- commandArgs(trailingOnly = TRUE); M <- if (length(args) >= 1) as.integer(args[1]) else 500L; M_boot <- if (length(args) >= 2) as.integer(args[2]) else 200L; cores <- if (length(args) >= 3) as.integer(args[3]) else as.integer(Sys.getenv("SIM_CORES", "1")); B <- 999L
source("R/state_space_em.R"); source("R/free_variance_cluster_inference.R"); suppressPackageStartupMessages(library(parallel))
f <- pss_simulation_design("positive"); Xmat <- f$Xmat; subj <- f$meta$subject; N <- nrow(Xmat); Tt <- ncol(f$Y); G <- 44; pn <- colnames(Xmat); p <- 6
beta0 <- setNames(f$truth$beta, pn); a1_0 <- f$truth$a1; P1_0 <- f$truth$P1; se2_0 <- f$truth$sigma_eta2; r0 <- 400; seed_base <- 410000
gen <- function(tau, X = Xmat, cl = subj) { ug <- unique(cl); b <- setNames(rnorm(length(ug), 0, tau), ug); phi <- tanh(as.vector(X %*% beta0) + b[cl]); n <- nrow(X); al <- rnorm(n, a1_0, sqrt(P1_0)); y <- matrix(0, n, Tt); y[, 1] <- al + rnorm(n, 0, sqrt(r0)); for (t in 2:Tt) { al <- phi * al + rnorm(n, 0, sqrt(se2_0)); y[, t] <- al + rnorm(n, 0, sqrt(r0)) }; y }
fit_free <- function(y, X, start, thres, max_iter) tryCatch(em_regression_mv(y, X, 0, start, thres = thres, max_iter = max_iter), error = function(e) NULL)
start_generic <- function(y, b) list(beta = b, a1 = mean(y[, 1]), P1 = max(var(y[, 1]) - 300, 100), sigma_eta2 = 35, sigma_eps2 = 300)
one <- function(m, tau) { out_dir <- sprintf("results/simulation/freevar/tau%s", tau); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE); fn <- file.path(out_dir, sprintf("m%03d.rds", m)); if (file.exists(fn)) return(invisible())
  set.seed(seed_base + 10000 * round(10 * tau) + m); y <- gen(tau); t0 <- Sys.time(); fit <- fit_free(y, Xmat, start_generic(y, beta0), 1e-7, 20000); if (is.null(fit)) { saveRDS(list(m = m, tau = tau, err = TRUE), fn); return(invisible()) }
  ci <- tryCatch(cluster_inference_free(y, Xmat, 0, fit, subj), error = function(e) NULL)
  out <- list(m = m, tau = tau, err = FALSE, beta = fit$beta, a1 = fit$a1, P1 = fit$P1, sigma_eta2 = fit$sigma_eta2, sigma_eps2 = fit$sigma_eps2, conv = fit$converged, status = fit$status, iter = fit$iter, loglik = fit$loglik, max_phi = max(abs(tanh(Xmat %*% fit$beta))),
              se = if (is.null(ci)) NULL else list(model = ci$se_model, cr0 = ci$se_cr0, md = ci$se_md), min_eig = if (is.null(ci)) NA else ci$min_eig, score_max = if (is.null(ci)) NA else max(abs(ci$score_sum)), secs_fit = as.numeric(Sys.time() - t0, units = "secs"))
  if (m <= M_boot) { draws <- matrix(NA_real_, B, p, dimnames = list(NULL, pn)); diag <- data.frame(b = seq_len(B), converged = NA, status = NA_character_, iter = NA_integer_, attempt = NA_integer_, loglik = NA_real_, sigma_eps2 = NA_real_)
    set.seed(seed_base + 5000000 + 10000 * round(10 * tau) + m); s1 <- list(beta = fit$beta, a1 = fit$a1, P1 = fit$P1, sigma_eta2 = fit$sigma_eta2, sigma_eps2 = fit$sigma_eps2)
    for (b in seq_len(B)) { idx <- unlist(lapply(sample(unique(subj), G, TRUE), function(s) which(subj == s))); yb <- y[idx, ]; Xb <- Xmat[idx, ]; fb <- fit_free(yb, Xb, s1, 1e-6, 20000); att <- 1L
      if (is.null(fb) || !fb$converged) { fb <- fit_free(yb, Xb, start_generic(yb, fit$beta), 1e-6, 40000); att <- 2L }
      if (!is.null(fb)) { diag$converged[b] <- fb$converged; diag$status[b] <- fb$status; diag$iter[b] <- fb$iter; diag$loglik[b] <- fb$loglik; diag$sigma_eps2[b] <- fb$sigma_eps2; if (fb$converged) draws[b, ] <- fb$beta } else diag$status[b] <- "error"; diag$attempt[b] <- att }
    out$draws <- draws; out$boot_diag <- diag; out$n_unresolved <- sum(is.na(draws[, 1])) }
  out$minutes <- as.numeric(Sys.time() - t0, units = "mins"); saveRDS(out, fn); invisible() }
jobs <- expand.grid(tau = c(0, 0.4), m = seq_len(M)); jobs <- jobs[order(jobs$m > M_boot, jobs$m), ]   # bootstrap datasets first, both scenarios interleaved
if (length(args) >= 4 && args[4] == "reverse") jobs <- jobs[rev(seq_len(nrow(jobs))), ]   # 4th arg "reverse": second instance from the end (analytic-only datasets first, then bootstrap ones from m = 200 down)
t0 <- Sys.time(); invisible(mclapply(seq_len(nrow(jobs)), function(k) one(jobs$m[k], jobs$tau[k]), mc.cores = cores, mc.preschedule = FALSE)); cat(sprintf("Part R: %d jobs in %.1f min\n", nrow(jobs), as.numeric(Sys.time() - t0, units = "mins")))
