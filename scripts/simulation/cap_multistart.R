source("R/simulation_design.R")
# Portable October 13 replication. Run from the repository root; see docs/SIMULATIONS.txt.
## Targeted multistart for the cap-sensitivity study (11 Oct 2026): the feasible sets nest (0.99 < 0.995 < 0.999), so a fully
## maximized log-likelihood cannot decrease when the cap loosens; decreases in the stored results are optimizer failures.
## For every original fit and every bootstrap refit at each cap c, the solutions obtained at the OTHER caps that are feasible
## under c (max |phi| <= c) and differ from the current solution are used as additional EM starts at cap c; the best converged
## feasible solution is retained. Writes results/simulation/cap_sensitivity_ms/<same name>.rds with the updated solutions and a record of
## every replacement (which cap, improvement, coefficient change). Usage: Rscript cap_sensitivity_multistart.R <cores>
cores <- as.integer(commandArgs(trailingOnly = TRUE)[1]); source("R/state_space_em.R"); suppressPackageStartupMessages(library(parallel)); caps <- c(0.99, 0.995, 0.999)
f <- pss_simulation_design("positive"); Xmat <- f$Xmat; subj <- f$meta$subject; N <- nrow(Xmat); Tt <- ncol(f$Y); G <- 44; pn <- colnames(Xmat); p <- 6
fb <- pss_simulation_design("boundary")$truth; beta0 <- setNames(fb$beta, pn); a1_0 <- fb$a1; P1_0 <- fb$P1; se2_0 <- fb$sigma_eta2
gen <- function(tau) { b_g <- setNames(rnorm(G, 0, tau), unique(subj)); phi <- tanh(as.vector(Xmat %*% beta0) + b_g[subj]); alpha <- y <- matrix(0, N, Tt); alpha[, 1] <- rnorm(N, a1_0, sqrt(P1_0)); y[, 1] <- alpha[, 1] + rnorm(N, 0, 0); for (t in 2:Tt) { alpha[, t] <- phi * alpha[, t - 1] + rnorm(N, 0, sqrt(se2_0)); y[, t] <- alpha[, t] + rnorm(N, 0, 0) }; y }
se_all <- function(y, fit) list(model = tryCatch(compute_SE_regression_mv_fixed_eps2(y, Xmat, 0, fit)$se[pn], error = function(e) rep(NA, p)), cr0 = tryCatch(compute_SE_regression_mv_cluster_robust(y, Xmat, 0, fit, subj)$se[pn], error = function(e) rep(NA, p)), md = tryCatch(compute_SE_regression_mv_cluster_MD(y, Xmat, 0, fit, subj)$se[pn], error = function(e) rep(NA, p)))
dir.create("results/simulation/cap_sensitivity_ms", showWarnings = FALSE)
improve <- function(y, X, cur_beta, cur_ll, cap, cand, nuis, thres, max_iter) { best <- NULL; for (cb in cand) { if (max(abs(tanh(X %*% cb))) > cap + 1e-12 || max(abs(cb - cur_beta)) < 1e-8) next
    fit <- tryCatch(em_regression_mv(y, X, 0, list(beta = cb, a1 = nuis$a1, P1 = nuis$P1, sigma_eta2 = nuis$sigma_eta2, sigma_eps2 = 0), thres = thres, max_iter = max_iter, fix_sigma_eps2 = 0, phi_cap = cap), error = function(e) NULL)
    if (!is.null(fit) && fit$converged && fit$loglik > cur_ll + 1e-6 && (is.null(best) || fit$loglik > best$loglik)) best <- fit }; best }
one <- function(fn) { out <- file.path("results/simulation/cap_sensitivity_ms", basename(fn)); if (file.exists(out)) return(invisible()); z <- readRDS(fn); r <- z$res
  if (is.null(z$set) || z$set == "application") { y <- f$Y; set.seed(20260929); idx_all <- lapply(seq_len(2000), function(b) unlist(lapply(sample(unique(subj), G, TRUE), function(s) which(subj == s)))) } else { set.seed(900000 + 10000 * round(10 * z$tau) + z$m); y <- gen(z$tau); set.seed(900000 + 5000000 + 10000 * round(10 * z$tau) + z$m); idx_all <- lapply(seq_len(999), function(b) unlist(lapply(sample(unique(subj), G, TRUE), function(s) which(subj == s)))) }
  rec <- list(); nuis <- list(a1 = r[[2]]$beta * 0 + mean(y[, 1]), P1 = max(var(y[, 1]), 1), sigma_eta2 = 35)
  ## original fits
  for (ci in seq_along(caps)) { cand <- lapply(setdiff(seq_along(caps), ci), function(k) r[[k]]$beta); b <- improve(y, Xmat, r[[ci]]$beta, r[[ci]]$loglik, caps[ci], cand, list(a1 = r[[ci]]$beta[1] * 0 + mean(y[, 1]), P1 = max(var(y[, 1]), 1), sigma_eta2 = 35), 1e-8, 3000)
    if (!is.null(b)) { rec[[length(rec) + 1]] <- data.frame(kind = "original", cap = caps[ci], b = 0, improvement = b$loglik - r[[ci]]$loglik, max_beta_change = max(abs(b$beta - r[[ci]]$beta)), new_status = b$status); r[[ci]]$beta <- b$beta; r[[ci]]$loglik <- b$loglik; r[[ci]]$status <- b$status; r[[ci]]$max_phi <- max(abs(tanh(Xmat %*% b$beta))); r[[ci]]$se <- se_all(y, b) } }
  ## bootstrap refits
  B <- nrow(r[[1]]$draws); for (b in seq_len(B)) { idx <- idx_all[[b]]; yb <- y[idx, ]; Xb <- Xmat[idx, ]; nb <- list(a1 = mean(yb[, 1]), P1 = max(var(yb[, 1]), 1), sigma_eta2 = 35)
    for (ci in seq_along(caps)) { cur <- r[[ci]]$draws[b, ]; if (any(is.na(cur))) next; cand <- lapply(setdiff(seq_along(caps), ci), function(k) r[[k]]$draws[b, ]); cand <- cand[!sapply(cand, function(cb) any(is.na(cb)))]; if (!length(cand)) next
      bt <- improve(yb, Xb, cur, r[[ci]]$boot_loglik[b], caps[ci], cand, nb, 1e-6, 3000)
      if (!is.null(bt)) { rec[[length(rec) + 1]] <- data.frame(kind = "bootstrap", cap = caps[ci], b = b, improvement = bt$loglik - r[[ci]]$boot_loglik[b], max_beta_change = max(abs(bt$beta - cur)), new_status = bt$status); r[[ci]]$draws[b, ] <- bt$beta; r[[ci]]$boot_loglik[b] <- bt$loglik; r[[ci]]$boot_status[b] <- bt$status } } }
  z$res <- r; z$multistart <- if (length(rec)) do.call(rbind, rec) else NULL; saveRDS(z, out); invisible() }
fs <- c("results/simulation/cap_sensitivity/application.rds", list.files("results/simulation/cap_sensitivity", "^(subset|stress)_", full.names = TRUE)); t0 <- Sys.time(); invisible(mclapply(fs, one, mc.cores = cores, mc.preschedule = FALSE)); cat(sprintf("multistart: %d files in %.1f min\n", length(fs), as.numeric(Sys.time() - t0, units = "mins")))
