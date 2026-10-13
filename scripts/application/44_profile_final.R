## Step 44 (final windows; derived from Steps 12/22/29). Keeps convergence metadata for EVERY start
## (sigma_eps2_profile_final_starts.csv) and continues the best unconstrained fit for up to 20000 further
## EM iterations, so that 'best attained likelihood' is not confused with a converged estimate.
start_log <- list()
## Step 22: Step 12 (safeguarded profile + unconstrained fit) on the corrected v2 windows.
## Step 12: redo the sigma_eps2 profile (Step 10) and the unconstrained fit
## with the SAFEGUARDED EM (phi_cap backtracking in em_regression_mv, added
## 2026-09-25) and multiple starting values, keeping the best log-likelihood.
##
## Why: a diagnostic (scratch check_profile_trap.R) showed that the Step-10
## fits with phi_hat = 1 at SD_eps <= 15 were optimizer artifacts -- from an
## interior start (b0 = atanh(0.94)) the same fixed-sigma_eps2 likelihoods have
## interior maxima with mean phi_hat ~0.94 and HIGHER log-likelihood
## (e.g. SD=10: -52913 vs -53792). The "phase transition" at 15-18 mg/dL was
## therefore where the poor start (b0 = atanh(0.3)) happened to escape the
## (1-phi^2) gradient trap, not a property of the likelihood.
source(file.path(APP_ROOT, "R", "state_space_em.R"))
f <- readRDS("FINAL_fitted_models_final.rds"); Y <- f$Y; Xmat <- f$Xmat; mu <- 0
N <- nrow(Y); Tt <- ncol(Y)
naive <- naive_two_stage(Y, Xmat)
b0_starts <- c(atanh(0.3), atanh(0.9), atanh(min(0.98, median(pmax(pmin(naive$phi_hat, 0.99), -0.99), na.rm = TRUE))))
cat("intercept starts:", round(b0_starts, 3), " (third = Fisher-z of median naive slope)\n")

fit_best <- function(fix) {
  best <- NULL
  for (b0 in b0_starts) for (se2_init in c(20, 60)) {
    init <- list(beta = c(b0, rep(0, ncol(Xmat) - 1)), a1 = mean(Y[, 1]), P1 = max(var(Y[, 1]), 1),
                 sigma_eta2 = se2_init, sigma_eps2 = if (is.null(fix)) 100 else fix)
    fit <- em_regression_mv(Y, Xmat, mu = mu, init = init, thres = 1e-7, max_iter = 4000, fix_sigma_eps2 = fix)
    fit$start <- c(b0 = b0, sigma_eta2 = se2_init)
    start_log[[length(start_log) + 1]] <<- data.frame(fixed_sigma_eps2 = if (is.null(fix)) NA else fix, b0 = b0, sigma_eta2_init = se2_init, iter = fit$iter, converged = fit$converged, loglik = fit$loglik, sigma_eps2 = fit$sigma_eps2, sigma_eta2 = fit$sigma_eta2, diabetic = fit$beta[6], n_backtrack = fit$n_backtrack, at_cap = fit$at_cap)
    if (is.null(best) || fit$loglik > best$loglik + 1e-6) best <- fit
    cat(sprintf("   start b0=%.2f se2=%2d: iter %4d conv %s loglik %.1f phi mean %.4f sig_eta2 %.2f sig_eps2 %.3f backtracks %d at_cap %s\n",
                b0, se2_init, fit$iter, fit$converged, fit$loglik, mean(tanh(Xmat %*% fit$beta)), fit$sigma_eta2, fit$sigma_eps2, fit$n_backtrack, fit$at_cap))
  }
  best
}

cat("\n=== Unconstrained fit (sigma_eps2 free) ===\n")
unc <- fit_best(NULL)
cat(sprintf("BEST: loglik %.1f | sigma_eps2 %.4f (SD %.3f mg/dL) | sigma_eta2 %.2f | phi mean %.4f range %.4f-%.4f\n",
            unc$loglik, unc$sigma_eps2, sqrt(unc$sigma_eps2), unc$sigma_eta2, mean(tanh(Xmat %*% unc$beta)), min(tanh(Xmat %*% unc$beta)), max(tanh(Xmat %*% unc$beta))))
cat("continuing the best unconstrained fit for up to 20000 further iterations...\n")
unc2 <- em_regression_mv(Y, Xmat, mu = mu, init = list(beta = unc$beta, a1 = unc$a1, P1 = unc$P1, sigma_eta2 = unc$sigma_eta2, sigma_eps2 = unc$sigma_eps2), thres = 1e-7, max_iter = 20000)
cat(sprintf("after continuation: iter %d converged %s sigma_eps2 %.6f sigma_eta2 %.3f loglik %.2f diabetic %.4f\n", unc2$iter, unc2$converged, unc2$sigma_eps2, unc2$sigma_eta2, unc2$loglik, unc2$beta[6]))
unc$continued <- unc2[c("beta", "a1", "P1", "sigma_eta2", "sigma_eps2", "iter", "converged", "loglik")]
saveRDS(unc, "unconstrained_fit_final.rds")

cat("\n=== Profile over fixed SD_eps, safeguarded EM, best of starts ===\n")
sd_grid <- c(2, 4, 6, 8, 10, 12, 15, 18, 20, 22, 26, 30)
rows <- list(); fits <- list()
for (sd in sd_grid) {
  cat(sprintf("-- SD = %d\n", sd))
  fit <- fit_best(sd^2); fits[[as.character(sd)]] <- fit
  ph <- tanh(as.vector(Xmat %*% fit$beta))
  rows[[length(rows) + 1]] <- data.frame(sd_eps = sd, sigma_eps2 = sd^2, converged = fit$converged, iter = fit$iter,
                                         sigma_eta2 = fit$sigma_eta2, a1 = fit$a1, P1 = fit$P1, loglik = fit$loglik,
                                         phi_mean = mean(ph), phi_min = min(ph), phi_max = max(ph),
                                         t(setNames(fit$beta, colnames(Xmat))), start_b0 = fit$start["b0"], start_se2 = fit$start["sigma_eta2"],
                                         n_backtrack = fit$n_backtrack, at_cap = fit$at_cap)
}
prof <- do.call(rbind, rows); rownames(prof) <- NULL
print(prof[, c("sd_eps", "iter", "loglik", "phi_mean", "phi_min", "phi_max", "sigma_eta2", "Intercept", "Age", "Sex", "BMI", "Prediabetic", "Diabetic", "n_backtrack", "at_cap")], digits = 4, row.names = FALSE)
write.csv(prof, "sigma_eps2_profile_final.csv", row.names = FALSE)
write.csv(do.call(rbind, start_log), "sigma_eps2_profile_final_starts.csv", row.names = FALSE)
saveRDS(list(profile = prof, fits = fits, unconstrained = unc), "sigma_eps2_profile_final_raw.rds")
cat("\nSaved sigma_eps2_profile_final.csv, sigma_eps2_profile_final_raw.rds, unconstrained_fit_final.rds\n")
