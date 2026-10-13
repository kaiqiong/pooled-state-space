source("R/simulation_design.R")
# Portable October 13 replication. Run from the repository root; see docs/SIMULATIONS.txt.
## Strengthened participant-bootstrap study (7 Oct 2026), run on the SAME simulated datasets as
##   * Part Q (positive observation variance, sigma_eps2 = 400 fixed; seeds 410000 + 10000*round(10*tau) + m), and
##   * Part M v2 (boundary, sigma_eps2 = 0; seeds 900000 + 10000*round(10*tau) + m),
## with B = 999 participant resamples per dataset and full per-refit diagnostics.
## Usage: Rscript sim_bootstrap_b999.R <design: posvar|boundary> <M_boot> <cores>
## Each dataset is regenerated from its seed, refitted exactly as in the parent script (the refit is checked
## against the saved estimate), and then bootstrapped. Every refit records converged / status / iterations / attempt,
## and the coefficient draws are saved. Retry procedure (documented in the supplement): a refit that errors or does
## not converge within 3000 iterations from the full-data estimate is retried once from a generic start
## (beta = full-data estimate, a1 = mean of first readings, P1 from their variance, sigma_eta2 = 35) with up to
## 10000 iterations; a refit that still fails is recorded as unresolved and excluded from the percentile interval.
## Output: results/simulation/boot_b999/<design>_tau<tau>/m<m>.rds (one file per dataset; the run is resumable).
args <- commandArgs(trailingOnly = TRUE); design <- if (length(args) >= 1) args[1] else "posvar"
M_boot <- if (length(args) >= 2) as.integer(args[2]) else 200L; cores <- if (length(args) >= 3) as.integer(args[3]) else 8L; B <- 999L
source("R/state_space_em.R"); suppressPackageStartupMessages(library(parallel))
f <- pss_simulation_design("positive"); Xmat <- f$Xmat; subj <- f$meta$subject; N <- nrow(Xmat); Tt <- ncol(f$Y); G <- length(unique(subj)); pn <- colnames(Xmat); p <- length(pn)
if (design == "posvar") {
  beta0 <- setNames(f$truth$beta, pn); a1_0 <- f$truth$a1; P1_0 <- f$truth$P1; se2_0 <- f$truth$sigma_eta2; r0 <- 400; seed_base <- 410000
  gen <- function(tau, X = Xmat, cl = subj) { ug <- unique(cl); b <- setNames(rnorm(length(ug), 0, tau), ug); phi <- tanh(as.vector(X %*% beta0) + b[cl]); n <- nrow(X); al <- rnorm(n, a1_0, sqrt(P1_0)); y <- matrix(0, n, Tt); y[, 1] <- al + rnorm(n, 0, sqrt(r0))
    for (t in 2:Tt) { al <- phi * al + rnorm(n, 0, sqrt(se2_0)); y[, t] <- al + rnorm(n, 0, sqrt(r0)) }; y }
  fit_full <- function(y) em_regression_mv(y, Xmat, 0, list(beta = beta0, a1 = mean(y[, 1]), P1 = max(var(y[, 1]) - r0, 100), sigma_eta2 = 35, sigma_eps2 = r0), thres = 1e-7, max_iter = 5000, fix_sigma_eps2 = r0)
  saved <- readRDS("results/simulation/tables/part_q_raw.rds")$raw; saved_beta <- function(tau, m) saved[[as.character(tau)]][[m]]$beta
} else if (design == "boundary") {
  f <- pss_simulation_design("boundary")
  beta0 <- setNames(f$truth$beta, pn); a1_0 <- f$truth$a1; P1_0 <- f$truth$P1; se2_0 <- f$truth$sigma_eta2; r0 <- 0; seed_base <- 900000
  gen <- function(tau) { b_g <- setNames(rnorm(G, 0, tau), unique(subj)); phi <- tanh(as.vector(Xmat %*% beta0) + b_g[subj]); alpha <- y <- matrix(0, N, Tt)
    alpha[, 1] <- rnorm(N, a1_0, sqrt(P1_0)); y[, 1] <- alpha[, 1] + rnorm(N, 0, sqrt(r0))
    for (t in 2:Tt) { alpha[, t] <- phi * alpha[, t - 1] + rnorm(N, 0, sqrt(se2_0)); y[, t] <- alpha[, t] + rnorm(N, 0, sqrt(r0)) }; y }
  fit_full <- function(y) { nv <- naive_two_stage(y, Xmat); b_start <- c(atanh(min(0.98, median(pmax(pmin(nv$phi_hat, 0.99), -0.99), na.rm = TRUE))), rep(0, p - 1))
    em_regression_mv(y, Xmat, 0, list(beta = b_start, a1 = mean(y[, 1]), P1 = max(var(y[, 1]), 1), sigma_eta2 = 35, sigma_eps2 = 0), thres = 1e-8, max_iter = 3000, fix_sigma_eps2 = 0) }
  saved <- readRDS("results/simulation/tables/part_m_boundary_v2_raw.rds")$raw; saved_beta <- function(tau, m) saved[[as.character(tau)]][[m]]$beta
} else stop("design must be posvar or boundary")
cat(sprintf("design %s: N=%d G=%d T=%d r0=%d | M_boot=%d B=%d cores=%d\n", design, N, G, Tt, r0, M_boot, B, cores))
refit <- function(y, X, start, thres = 1e-6, max_iter = 3000) tryCatch(em_regression_mv(y, X, 0, start, thres = thres, max_iter = max_iter, fix_sigma_eps2 = r0), error = function(e) NULL)
one <- function(m, tau) {
  out_dir <- sprintf("results/simulation/boot_b999/%s_tau%s", design, tau); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE); fn <- file.path(out_dir, sprintf("m%03d.rds", m)); if (file.exists(fn)) return(invisible(NULL))
  set.seed(seed_base + 10000 * round(10 * tau) + m); y <- gen(tau); t0 <- Sys.time(); fit <- fit_full(y)
  refit_dev <- max(abs(fit$beta - saved_beta(tau, m)))
  draws <- matrix(NA_real_, B, p, dimnames = list(NULL, pn)); diag <- data.frame(b = seq_len(B), converged = NA, status = NA_character_, iter = NA_integer_, attempt = NA_integer_, loglik = NA_real_)
  set.seed(seed_base + 5000000 + 10000 * round(10 * tau) + m)   # bootstrap resampling stream, separate from the data seed
  s1 <- list(beta = fit$beta, a1 = fit$a1, P1 = fit$P1, sigma_eta2 = fit$sigma_eta2, sigma_eps2 = r0)
  for (b in seq_len(B)) { idx <- unlist(lapply(sample(unique(subj), G, TRUE), function(s) which(subj == s))); yb <- y[idx, ]; Xb <- Xmat[idx, ]
    fb <- refit(yb, Xb, s1); att <- 1L
    if (is.null(fb) || !fb$converged) { fb <- refit(yb, Xb, list(beta = fit$beta, a1 = mean(yb[, 1]), P1 = max(var(yb[, 1]) - r0, 100), sigma_eta2 = 35, sigma_eps2 = r0), max_iter = 10000); att <- 2L }
    if (!is.null(fb)) { diag$converged[b] <- fb$converged; diag$status[b] <- fb$status; diag$iter[b] <- fb$iter; diag$loglik[b] <- fb$loglik; if (fb$converged) draws[b, ] <- fb$beta } else diag$status[b] <- "error"
    diag$attempt[b] <- att }
  saveRDS(list(design = design, tau = tau, m = m, B = B, beta_full = fit$beta, fit_status = fit$status, fit_iter = fit$iter, refit_dev_from_saved = refit_dev, draws = draws, diag = diag,
               n_unresolved = sum(is.na(draws[, 1])), minutes = as.numeric(Sys.time() - t0, units = "mins")), fn)
  invisible(NULL) }
## datasets are processed in the order m = 1, 2, ... with both scenarios interleaved, so partial results exist for both
jobs <- expand.grid(tau = c(0, 0.4), m = seq_len(M_boot)); if (length(args) >= 4 && args[4] == "reverse") jobs <- jobs[rev(seq_len(nrow(jobs))), ]; t0 <- Sys.time()   # 4th arg "reverse": a second instance working from the end
invisible(mclapply(seq_len(nrow(jobs)), function(k) one(jobs$m[k], jobs$tau[k]), mc.cores = cores, mc.preschedule = FALSE))
cat(sprintf("%s: %d datasets x 2 scenarios done in %.1f min\n", design, M_boot, as.numeric(Sys.time() - t0, units = "mins")))
cat("done\n")
