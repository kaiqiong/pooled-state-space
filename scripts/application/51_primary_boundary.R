## Step 51: PRIMARY ANALYSIS under option 1 (decided 2026-09-29): the measurement-error variance is
## estimated freely and goes to its boundary (0), so the primary fit is at sigma_eps2 = 0; fixed positive
## values are sensitivity analyses (profile, Steps 44/48). Inference: nuisance-aware cluster CR0 and a
## participant bootstrap (B = 2000) in which every resample is refitted at the boundary; each resample
## also records the log-likelihood at sigma_eps = 1 mg/dL to check that the boundary is preferred there too.
source(file.path(APP_ROOT, "R", "state_space_em.R"))
f <- readRDS("FINAL_fitted_models_final.rds"); Y <- f$Y; X <- f$Xmat; m <- f$meta; cl <- m$subject; bn <- colnames(X); G <- length(unique(cl)); qt_ <- qt(.975, G - 1); qn <- qnorm(.975)
fit0 <- function(Yb, Xb, r = 0, init = NULL) { if (is.null(init)) init <- list(beta = f$fit$beta, a1 = mean(Yb[, 1]), P1 = var(Yb[, 1]), sigma_eta2 = 35)
  em_regression_mv(Yb, Xb, 0, c(init[c("beta", "a1", "P1", "sigma_eta2")], list(sigma_eps2 = r)), thres = 1e-9, max_iter = 3000, fix_sigma_eps2 = r) }
fit <- fit0(Y, X); na <- compute_SE_regression_mv_cluster_beta(Y, X, mu = 0, fit = fit, cluster = cl, return_unit_info = TRUE); md <- compute_SE_regression_mv_cluster_MD(Y, X, 0, fit, cl)
ph <- as.vector(tanh(X %*% fit$beta)); hl <- 5 * log(2) / (-log(ph))
cat(sprintf("boundary fit: iter %d conv %s loglik %.2f | sigma_eta2 %.3f a1 %.2f P1 %.1f\n", fit$iter, fit$converged, fit$loglik, fit$sigma_eta2, fit$a1, fit$P1))
cat("phi by group:", round(tapply(ph, m$group, mean), 4), " | half-life (min) at group-mean phi:", round(5 * log(2) / (-log(tapply(ph, m$group, mean))), 1), " | mean of window half-lives:", round(tapply(hl, m$group, mean), 1), "\n")
## innovation diagnostics
kf <- fit$kf; e <- (Y - kf$alpha_pred) / sqrt(kf$P_pred); e <- e[, -1]; ec <- sweep(e, 2, colMeans(e)); cat("innovation lag-1 / lag-2 correlation (slot means removed):", round(cor(as.vector(ec[, -1]), as.vector(ec[, -ncol(ec)])), 3), round(cor(as.vector(ec[, -(1:2)]), as.vector(ec[, 1:(ncol(ec) - 2)])), 3), " | largest |slot mean| of standardized innovations:", round(max(abs(colMeans(e))), 3), "\n")
## naive two-stage
nv <- naive_two_stage(Y, X); z <- as.vector(nv$z); ok <- is.finite(z)
## leave-one-participant-out
loo <- t(sapply(unique(cl), function(g) { i <- cl != g; fit0(Y[i, ], X[i, ], init = list(beta = fit$beta, a1 = fit$a1, P1 = fit$P1, sigma_eta2 = fit$sigma_eta2))$beta })); colnames(loo) <- bn
cat("leave-one-participant-out Diabetic: range", round(range(loo[, "Diabetic"]), 3), " most influential:", names(which.max(abs(loo[, "Diabetic"] - fit$beta[6]))), round(loo[which.max(abs(loo[, "Diabetic"] - fit$beta[6])), "Diabetic"], 3), "\n")
## collinear positions treated as unobserved
Yna <- Y; Yna[f$flags != "kink"] <- NA; fna <- em_regression_mv_na(Yna, X, 0, list(beta = fit$beta, a1 = fit$a1, P1 = fit$P1, sigma_eta2 = fit$sigma_eta2, sigma_eps2 = 0), thres = 1e-9, max_iter = 3000, fix_sigma_eps2 = 0)
cat("collinear-as-missing: Diabetic", round(fna$beta[6], 4), " conv", fna$converged, "\n")
## bootstrap
set.seed(20260929); B <- pss_boot_count(2000); subs <- unique(cl); bj <- matrix(NA, B, 6, dimnames = list(NULL, bn)); bnv <- bj; bd <- matrix(NA, B, 6, dimnames = list(NULL, c("iter", "converged", "loglik0", "loglik1", "at_cap", "max_phi"))); fail <- character(B); t0 <- Sys.time()
for (b in seq_len(B)) { idx <- unlist(lapply(sample(subs, G, TRUE), function(s) which(cl == s))); Yb <- Y[idx, ]; Xb <- X[idx, ]
  fb <- tryCatch(fit0(Yb, Xb, init = list(beta = fit$beta, a1 = fit$a1, P1 = fit$P1, sigma_eta2 = fit$sigma_eta2)), error = function(e) e); if (inherits(fb, "error")) { fail[b] <- conditionMessage(fb); next }
  f1 <- tryCatch(fit0(Yb, Xb, r = 1, init = list(beta = fb$beta, a1 = fb$a1, P1 = fb$P1, sigma_eta2 = fb$sigma_eta2)), error = function(e) NULL)
  bj[b, ] <- fb$beta; bd[b, ] <- c(fb$iter, fb$converged, fb$loglik, if (is.null(f1)) NA else f1$loglik, isTRUE(fb$at_cap), max(abs(tanh(Xb %*% fb$beta)))); bnv[b, ] <- tryCatch(naive_two_stage(Yb, Xb)$beta_hat, error = function(e) rep(NA, 6))
  if (b %% 250 == 0) cat(sprintf("  boot %d/%d (%.1f min)\n", b, B, as.numeric(Sys.time() - t0, units = "mins"))) }
mcse <- function(b, pr) { q <- quantile(b, pr, na.rm = TRUE); d <- density(b[is.finite(b)]); sqrt(pr * (1 - pr) / sum(is.finite(b))) / approx(d$x, d$y, xout = q)$y }
est <- setNames(fit$beta, bn)
joint <- data.frame(parameter = bn, estimate = est, se_model = na$se_model_based[bn], se_cr0 = na$se[bn], se_md = md$se[bn], se_boot = apply(bj, 2, sd, na.rm = TRUE), cr0_lo_t = est - qt_ * na$se[bn], cr0_hi_t = est + qt_ * na$se[bn],
  boot_pct_lo = apply(bj, 2, quantile, .025, na.rm = TRUE), boot_pct_hi = apply(bj, 2, quantile, .975, na.rm = TRUE), boot_mcse_lo = apply(bj, 2, mcse, pr = .025), boot_frac_le0 = colMeans(bj <= 0, na.rm = TRUE), loo_min = apply(loo, 2, min), loo_max = apply(loo, 2, max), collinear_missing = fna$beta)
joint$boot_basic_lo <- 2 * est - joint$boot_pct_hi; joint$boot_basic_hi <- 2 * est - joint$boot_pct_lo
Xn <- X[ok, ]; XtXi <- solve(crossprod(Xn)); bh <- as.vector(XtXi %*% crossprod(Xn, z[ok])); res <- z[ok] - as.vector(Xn %*% bh); meat <- 0; for (g in unique(cl[ok])) { i <- cl[ok] == g; sg <- crossprod(Xn[i, , drop = FALSE], res[i]); meat <- meat + sg %*% t(sg) }; nse <- sqrt(diag(XtXi %*% meat %*% XtXi))
naive <- data.frame(parameter = bn, estimate = bh, se_model = as.vector(nv$se), se_cr0 = nse, se_boot = apply(bnv, 2, sd, na.rm = TRUE), cr0_lo_t = bh - qt_ * nse, cr0_hi_t = bh + qt_ * nse, boot_pct_lo = apply(bnv, 2, quantile, .025, na.rm = TRUE), boot_pct_hi = apply(bnv, 2, quantile, .975, na.rm = TRUE), n_windows_used = sum(ok))
diag <- data.frame(B = B, refit_errors = sum(fail != ""), not_converged = sum(bd[, "converged"] == 0, na.rm = TRUE), at_cap = sum(bd[, "at_cap"] == 1, na.rm = TRUE), iter_min = min(bd[, "iter"], na.rm = TRUE), iter_max = max(bd[, "iter"], na.rm = TRUE), max_phi = max(bd[, "max_phi"], na.rm = TRUE),
  resamples_preferring_boundary_to_sd1 = sum(bd[, "loglik0"] > bd[, "loglik1"], na.rm = TRUE), min_loglik_gain_boundary_vs_sd1 = min(bd[, "loglik0"] - bd[, "loglik1"], na.rm = TRUE), naive_failures = sum(!complete.cases(bnv)), G = G)
options(width = 230); print(joint, digits = 3, row.names = FALSE); print(naive, digits = 3, row.names = FALSE); print(diag, row.names = FALSE)
write.csv(joint, "primary_boundary_joint.csv", row.names = FALSE); write.csv(naive, "primary_boundary_naive.csv", row.names = FALSE); write.csv(diag, "primary_boundary_bootstrap_diagnostics.csv", row.names = FALSE); write.csv(data.frame(subject = rownames(loo), loo), "primary_boundary_leave_one_out.csv", row.names = FALSE)
saveRDS(list(fit = fit, nuisance_aware = na, md = md, naive_fit = nv, boot_joint = bj, boot_naive = bnv, boot_diag = bd, boot_fail = fail, loo = loo, fit_na = fna, phi = ph, halflife = hl, joint = joint, naive = naive, diagnostics = diag, Y = Y, Xmat = X, meta = m, flags = f$flags, session = sessionInfo()), "PRIMARY_boundary.rds"); cat("saved\n")
