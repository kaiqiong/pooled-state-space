## Step 56: mixed-effects latent autoregression (participant random effect on the link scale) fitted to the
## primary windows at sigma_eps2 = 0 (implementation/mixed_latent_ar.R), as the comparator of Part P.
source(file.path(APP_ROOT, "R", "state_space_em.R")); source(file.path(APP_ROOT, "R", "mixed_latent_ar.R"))
pb <- readRDS("PRIMARY_boundary.rds"); Y <- pb$Y; X <- pb$Xmat; cl <- pb$meta$subject
fm <- fit_mixed_ar(Y, X, cl, start_theta = pb$fit$beta); fm2 <- fit_mixed_ar(Y, X, cl, start_theta = c(atanh(0.8), rep(0, 5)), tau_start = 0.05)
cat(sprintf("mixed model: loglik %.3f (second start %.3f) | tau %.4f (at floor %s) | sigma_eta2 %.3f | conv %s | %.1f s\n", fm$loglik, fm2$loglik, fm$tau, fm$tau_at_floor, fm$sigma_eta2, fm$conv, fm$secs))
tab <- data.frame(parameter = names(fm$theta), mixed_est = fm$theta, mixed_se = fm$se_theta, mixed_lo = fm$theta - qnorm(.975) * fm$se_theta, mixed_hi = fm$theta + qnorm(.975) * fm$se_theta, wi_est = pb$joint$estimate, wi_se_cr0 = pb$joint$se_cr0, wi_se_boot = pb$joint$se_boot, wi_boot_lo = pb$joint$boot_pct_lo, wi_boot_hi = pb$joint$boot_pct_hi)
print(tab, digits = 3, row.names = FALSE)
lrt <- 2 * (fm$loglik - pb$fit$loglik); cat(sprintf("LRT of tau = 0: %.2f (boundary null, p = %.4f from the 0.5 chi2_1 mixture) | SD of tau-hat: %.3f\n", lrt, 0.5 * pchisq(lrt, 1, lower.tail = FALSE), fm$se_all[8] * fm$tau))
cat("random-effect modes (link scale), range:", round(range(fm$modes), 3), "; mean phi by group with modes:", round(tapply(tanh(as.vector(X %*% fm$theta) + fm$modes[match(cl, unique(cl))]), pb$meta$group, mean), 4), "\n")
write.csv(tab, "mixed_model_application.csv", row.names = FALSE); saveRDS(list(fit = fm, fit2 = fm2, lrt = lrt), "mixed_model_application.rds"); cat("saved\n")
