source("R/simulation_design.R")
# Portable October 13 replication. Run from the repository root; see docs/SIMULATIONS.txt.
## Part N (final): differential attenuation on the FINAL application design (1174 windows, 44 participants),
## with the validated nuisance-aware cluster CR0 sandwich (t_{G-1}) for the joint estimator and M = 500.
## Supersedes sim_differential_attenuation{,_v2,_v3}.R (earlier window sets; legacy CR2-type SE).
## Scenarios N1-N4 are HYPOTHETICAL white-noise scenarios (measurement SD 20 mg/dL on average; truth = the
## fit of the final windows with SD fixed at 20), fitted with SD fixed at 20:
##   N1: peak size a1 by group (45/55/72), group effects as fitted     N2: same, group effects zero
##   N3: noise SD by group (14/20/28), group effects as fitted         N4: same, group effects zero
## Scenarios N5-N6 mimic the PRIMARY analysis (no white measurement error; truth = boundary fit), fitted at
## sigma_eps2 = 0, with peak size by group:  N5: group effects as fitted   N6: group effects zero.
source("R/state_space_em.R"); suppressPackageStartupMessages(library(parallel))
f <- pss_simulation_design("positive")
Xmat <- f$Xmat; subj <- f$meta$subject; grp <- as.character(f$meta$group); N <- nrow(Xmat); Tt <- ncol(f$Y); G <- length(unique(subj)); pn <- colnames(Xmat)
fit0 <- pss_simulation_design("boundary")$truth
tr20 <- list(beta = setNames(f$truth$beta, pn), se2 = f$truth$sigma_eta2, P1 = f$truth$P1, a1 = f$truth$a1, fit_r = 400); tr0 <- list(beta = setNames(fit0$beta, pn), se2 = fit0$sigma_eta2, P1 = fit0$P1, a1 = fit0$a1, fit_r = 0)
a1_grp <- c(healthy = 45, prediabetic = 55, diabetic = 72); sp_grp <- c(healthy = 14, prediabetic = 20, diabetic = 28); cst <- function(v) c(healthy = v, prediabetic = v, diabetic = v)
scen <- list(N1 = list(tr = tr20, a1 = a1_grp, sp = cst(20), null = FALSE), N2 = list(tr = tr20, a1 = a1_grp, sp = cst(20), null = TRUE),
             N3 = list(tr = tr20, a1 = cst(tr20$a1), sp = sp_grp, null = FALSE), N4 = list(tr = tr20, a1 = cst(tr20$a1), sp = sp_grp, null = TRUE),
             N5 = list(tr = tr0, a1 = a1_grp, sp = cst(0), null = FALSE), N6 = list(tr = tr0, a1 = a1_grp, sp = cst(0), null = TRUE))
gen <- function(sc, beta) { phi <- tanh(as.vector(Xmat %*% beta)); a1 <- sc$a1[grp]; sp <- sc$sp[grp]; alpha <- y <- matrix(0, N, Tt)
  alpha[, 1] <- rnorm(N, a1, sqrt(sc$tr$P1)); y[, 1] <- alpha[, 1] + rnorm(N, 0, sp)
  for (t in 2:Tt) { alpha[, t] <- phi * alpha[, t - 1] + rnorm(N, 0, sqrt(sc$tr$se2)); y[, t] <- alpha[, t] + rnorm(N, 0, sp) }; y }
qt_ <- qt(0.975, G - 1); M <- 500
one <- function(m, s, sc, beta) { set.seed(8000000 + 10000 * match(s, names(scen)) + m); y <- gen(sc, beta); nv <- naive_two_stage(y, Xmat); z <- as.vector(nv$z); ok <- is.finite(z)
  Xn <- Xmat[ok, ]; XtXi <- solve(crossprod(Xn)); bh <- as.vector(XtXi %*% crossprod(Xn, z[ok])); res <- z[ok] - as.vector(Xn %*% bh); meat <- 0
  for (g in unique(subj[ok])) { i <- subj[ok] == g; sg <- crossprod(Xn[i, , drop = FALSE], res[i]); meat <- meat + sg %*% t(sg) }; nv_se <- sqrt(diag(XtXi %*% meat %*% XtXi))
  b_start <- c(atanh(min(0.98, median(pmax(pmin(nv$phi_hat, 0.99), -0.99), na.rm = TRUE))), rep(0, 5)); r <- sc$tr$fit_r
  fit <- tryCatch(em_regression_mv(y, Xmat, 0, list(beta = b_start, a1 = mean(y[, 1]), P1 = max(var(y[, 1]), 1), sigma_eta2 = 35, sigma_eps2 = r), thres = 1e-7, max_iter = 3000, fix_sigma_eps2 = r), error = function(e) e)
  if (inherits(fit, "error")) return(list(err = conditionMessage(fit)))
  se <- tryCatch(as.numeric(compute_SE_regression_mv_cluster_robust(y, Xmat, 0, fit, subj)$se[pn]), error = function(e) rep(NA_real_, 6))
  list(err = "", joint = as.numeric(fit$beta), naive = bh, joint_se = se, naive_se = nv_se, converged = fit$converged, at_cap = isTRUE(fit$at_cap), n_naive = sum(ok)) }
rows <- list(); raw <- list(); rds <- "results/simulation/tables/part_n_final_raw.rds"; if (file.exists(rds)) raw <- readRDS(rds)$raw
for (s in names(scen)) { sc <- scen[[s]]; beta <- sc$tr$beta; if (sc$null) beta[c("Prediabetic", "Diabetic")] <- 0
  if (!(s %in% names(raw))) { t0 <- Sys.time(); raw[[s]] <- mclapply(seq_len(M), function(m) one(m, s, sc, beta), mc.cores = as.integer(Sys.getenv("SIM_CORES", "1")), mc.preschedule = FALSE); cat(sprintf("%s done in %.1f min\n", s, as.numeric(Sys.time() - t0, units = "mins"))); saveRDS(list(raw = raw), rds) }
  res <- raw[[s]]; okr <- sapply(res, function(z) is.list(z) && identical(z$err, "")); rr <- res[okr]; J <- t(sapply(rr, `[[`, "joint")); Nv <- t(sapply(rr, `[[`, "naive")); Jse <- t(sapply(rr, `[[`, "joint_se")); Nse <- t(sapply(rr, `[[`, "naive_se"))
  for (j in seq_along(pn)) { u <- is.finite(Jse[, j]); rows[[length(rows) + 1]] <- data.frame(scenario = s, parameter = pn[j], truth = beta[j], M = length(res), n_fit = sum(okr), n_joint_se = sum(u), not_converged = sum(!sapply(rr, `[[`, "converged")), at_cap = sum(sapply(rr, `[[`, "at_cap")),
    joint_mean = mean(J[, j]), joint_bias = mean(J[, j]) - beta[j], joint_sd = sd(J[, j]), naive_mean = mean(Nv[, j]), naive_bias = mean(Nv[, j]) - beta[j], naive_sd = sd(Nv[, j]),
    joint_reject = mean(abs(J[u, j] / Jse[u, j]) > qt_), naive_reject = mean(abs(Nv[, j] / Nse[, j]) > qt_), joint_cover = mean(abs(J[u, j] - beta[j]) <= qt_ * Jse[u, j]), naive_cover = mean(abs(Nv[, j] - beta[j]) <= qt_ * Nse[, j])) } }
df <- do.call(rbind, rows); rownames(df) <- NULL; options(width = 220); print(df[df$parameter %in% c("Prediabetic", "Diabetic"), ], digits = 3, row.names = FALSE)
write.csv(df, "results/simulation/tables/part_n_final_summary.csv", row.names = FALSE); cat("saved\n")
