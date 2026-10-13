source("R/simulation_design.R")
# Portable October 13 replication. Run from the repository root; see docs/SIMULATIONS.txt.
source("R/state_space_em.R"); suppressPackageStartupMessages(library(parallel))
f <- pss_simulation_design("positive"); Xmat <- f$Xmat; subj <- f$meta$subject; N <- nrow(Xmat); Tt <- ncol(f$Y); G <- length(unique(subj)); pn <- colnames(Xmat); p <- 6
beta0 <- setNames(f$truth$beta, pn); a1_0 <- f$truth$a1; P1_0 <- f$truth$P1; se2_0 <- f$truth$sigma_eta2; r0 <- 400
cat(sprintf("truth: beta=(%s) a1=%.1f P1=%.0f sigma_eta2=%.1f sigma_eps2=%d\n", paste(round(beta0, 3), collapse = ", "), a1_0, P1_0, se2_0, r0))
gen <- function(tau, X = Xmat, cl = subj) { ug <- unique(cl); b <- setNames(rnorm(length(ug), 0, tau), ug); phi <- tanh(as.vector(X %*% beta0) + b[cl]); n <- nrow(X); al <- rnorm(n, a1_0, sqrt(P1_0)); y <- matrix(0, n, Tt); y[, 1] <- al + rnorm(n, 0, sqrt(r0))
  for (t in 2:Tt) { al <- phi * al + rnorm(n, 0, sqrt(se2_0)); y[, t] <- al + rnorm(n, 0, sqrt(r0)) }; y }
fitw <- function(y, X, start = beta0, thres = 1e-7) em_regression_mv(y, X, 0, list(beta = start, a1 = mean(y[, 1]), P1 = max(var(y[, 1]) - r0, 100), sigma_eta2 = 35, sigma_eps2 = r0), thres = thres, max_iter = 5000, fix_sigma_eps2 = r0)
## ---- replicates ----
M <- 500; M_boot <- 0; B_boot <- 99; q_t <- qt(0.975, G - 1); q_n <- qnorm(0.975)
one <- function(m, tau) { set.seed(410000 + 10000 * round(10 * tau) + m); y <- gen(tau); fit <- tryCatch(fitw(y, Xmat), error = function(e) NULL); if (is.null(fit)) return(list(m = m, err = TRUE))
  se <- list(model = tryCatch(compute_SE_regression_mv_fixed_eps2(y, Xmat, 0, fit)$se[pn], error = function(e) rep(NA, p)), cr0 = tryCatch(compute_SE_regression_mv_cluster_beta(y, Xmat, 0, fit, subj)$se[pn], error = function(e) rep(NA, p)), md = tryCatch(compute_SE_regression_mv_cluster_MD(y, Xmat, 0, fit, subj)$se[pn], error = function(e) rep(NA, p)))
  out <- list(m = m, err = FALSE, beta = fit$beta, conv = fit$converged, status = fit$status, iter = fit$iter, se = se)

  out }

# All final inference uses the separately seeded B=999 bootstrap.
raw <- list(); rds <- "results/simulation/tables/part_q_raw.rds"
if (file.exists(rds)) raw <- readRDS(rds)$raw
for (tau in c(0, 0.4)) {
 key <- as.character(tau)
 if (key %in% names(raw)) { stopifnot(length(raw[[key]]) == M); next }
 raw[[key]] <- mclapply(seq_len(M), function(m) one(m, tau), mc.cores = as.integer(Sys.getenv("SIM_CORES", "1")), mc.preschedule = FALSE)
 saveRDS(list(raw = raw, beta0 = beta0, M = M), rds)
}
