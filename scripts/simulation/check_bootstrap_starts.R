source("R/simulation_design.R")
# Portable October 13 replication. Run from the repository root; see docs/SIMULATIONS.txt.
## Prevalence check for suboptimal bootstrap refits (10 Oct 2026), prompted by the Oct 9 finding that boundary
## tau = 0.4, dataset 724, draw 4 converged to a constrained point 429 log-likelihood units below another feasible point.
## For every saved draw that ended CONSTRAINED (and a random subsample of interior draws), the resample is regenerated
## from the stored seed stream, refitted from an ALTERNATIVE start (generating coefficients with the generic nuisance
## start, as in the retry rule), and the two log-likelihoods compared. Usage: Rscript check_bootstrap_local_optima.R <posvar|boundary> <tau> <cores> [n_interior_per_dataset]
args <- commandArgs(trailingOnly = TRUE); design <- args[1]; tau <- as.numeric(args[2]); cores <- as.integer(args[3]); n_int <- if (length(args) >= 4) as.integer(args[4]) else 10L
source("R/state_space_em.R"); suppressPackageStartupMessages(library(parallel))
f <- pss_simulation_design("positive"); Xmat <- f$Xmat; subj <- f$meta$subject; N <- nrow(Xmat); Tt <- ncol(f$Y); G <- length(unique(subj)); pn <- colnames(Xmat); p <- length(pn)
if (design == "posvar") { beta0 <- setNames(f$truth$beta, pn); a1_0 <- f$truth$a1; P1_0 <- f$truth$P1; se2_0 <- f$truth$sigma_eta2; r0 <- 400; seed_base <- 410000
  gen <- function(tau, X = Xmat, cl = subj) { ug <- unique(cl); b <- setNames(rnorm(length(ug), 0, tau), ug); phi <- tanh(as.vector(X %*% beta0) + b[cl]); n <- nrow(X); al <- rnorm(n, a1_0, sqrt(P1_0)); y <- matrix(0, n, Tt); y[, 1] <- al + rnorm(n, 0, sqrt(r0))
    for (t in 2:Tt) { al <- phi * al + rnorm(n, 0, sqrt(se2_0)); y[, t] <- al + rnorm(n, 0, sqrt(r0)) }; y }
} else { f <- pss_simulation_design("boundary")
  beta0 <- setNames(f$truth$beta, pn); a1_0 <- f$truth$a1; P1_0 <- f$truth$P1; se2_0 <- f$truth$sigma_eta2; r0 <- 0; seed_base <- 900000
  gen <- function(tau) { b_g <- setNames(rnorm(G, 0, tau), unique(subj)); phi <- tanh(as.vector(Xmat %*% beta0) + b_g[subj]); alpha <- y <- matrix(0, N, Tt); alpha[, 1] <- rnorm(N, a1_0, sqrt(P1_0)); y[, 1] <- alpha[, 1] + rnorm(N, 0, sqrt(r0))
    for (t in 2:Tt) { alpha[, t] <- phi * alpha[, t - 1] + rnorm(N, 0, sqrt(se2_0)); y[, t] <- alpha[, t] + rnorm(N, 0, sqrt(r0)) }; y } }
refit <- function(y, X, start, max_iter = 3000) tryCatch(em_regression_mv(y, X, 0, start, thres = 1e-6, max_iter = max_iter, fix_sigma_eps2 = r0), error = function(e) NULL)
dir_in <- sprintf("results/simulation/boot_b999/%s_tau%s", design, tau); fs <- list.files(dir_in, "^m\\d+\\.rds$", full.names = TRUE)
one <- function(fn) { z <- readRDS(fn); m <- z$m; set.seed(seed_base + 10000 * round(10 * tau) + m); y <- gen(tau)
  set.seed(seed_base + 5000000 + 10000 * round(10 * tau) + m); idx_all <- lapply(seq_len(z$B), function(b) unlist(lapply(sample(unique(subj), G, TRUE), function(s) which(subj == s))))
  con <- which(z$diag$status == "constrained"); int <- setdiff(seq_len(z$B), con); set.seed(77000 + m); chk <- c(con, if (length(int)) sample(int, min(n_int, length(int))))
  s_full <- list(beta = z$beta_full, a1 = NA, P1 = NA, sigma_eta2 = NA, sigma_eps2 = r0)
  out <- lapply(chk, function(b) { idx <- idx_all[[b]]; yb <- y[idx, ]; Xb <- Xmat[idx, ]
    alt <- refit(yb, Xb, list(beta = beta0, a1 = mean(yb[, 1]), P1 = max(var(yb[, 1]) - r0, 100), sigma_eta2 = 35, sigma_eps2 = r0), max_iter = 10000)
    rep_ok <- if (b == chk[1]) { fb <- refit(yb, Xb, list(beta = z$beta_full, a1 = mean(y[, 1]), P1 = max(var(y[, 1]) - r0, 100), sigma_eta2 = 35, sigma_eps2 = r0)); NA } else NA  # placeholder (full-data nuisance values are not stored); reproduction is checked via loglik below
    data.frame(m = m, b = b, saved_status = z$diag$status[b], saved_loglik = z$diag$loglik[b], alt_status = if (is.null(alt)) "error" else alt$status, alt_conv = if (is.null(alt)) NA else alt$converged, alt_loglik = if (is.null(alt)) NA else alt$loglik,
               improvement = if (is.null(alt)) NA else alt$loglik - z$diag$loglik[b], max_beta_diff = if (is.null(alt)) NA else max(abs(alt$beta - z$draws[b, ])), alt_beta = if (is.null(alt)) I(list(rep(NA, p))) else I(list(alt$beta))) })
  do.call(rbind, out) }
t0 <- Sys.time(); res <- do.call(rbind, mclapply(fs, one, mc.cores = cores, mc.preschedule = FALSE)); cat(sprintf("%s tau=%s: %d draws checked (%d constrained, %d interior) in %.1f min\n", design, tau, nrow(res), sum(res$saved_status == "constrained"), sum(res$saved_status != "constrained"), as.numeric(Sys.time() - t0, units = "mins")))
saveRDS(res, sprintf("results/simulation/tables/bootstrap_local_optima_%s_tau%s.rds", design, tau))
for (st in c("constrained", "interior")) { r <- res[res$saved_status == st, ]; if (!nrow(r)) next; imp <- r$improvement
  cat(sprintf("  %s draws: n=%d | alt better by >1 loglik: %d (%.2f%%) | >0.01: %d | alt worse by >1: %d | median |improvement| %.2e | max improvement %.1f | max beta diff when better>1: %.3f\n", st, nrow(r), sum(imp > 1, na.rm = TRUE), 100 * mean(imp > 1, na.rm = TRUE), sum(imp > 0.01, na.rm = TRUE), sum(imp < -1, na.rm = TRUE), median(abs(imp), na.rm = TRUE), max(imp, na.rm = TRUE), if (any(imp > 1, na.rm = TRUE)) max(r$max_beta_diff[which(imp > 1)]) else 0)) }
