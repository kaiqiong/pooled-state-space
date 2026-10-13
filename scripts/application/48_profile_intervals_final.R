## Step 48 (final windows): cluster CR0 (full joint, nuisance-aware) interval with t_{G-1} for every coefficient at each fixed SD_eps of the
## profile (v4), so that the dependence of the inference -- not only of the
## point estimate -- on the externally fixed sigma_eps is reported.
source(file.path(APP_ROOT, "R", "state_space_em.R"))
f <- readRDS("FINAL_fitted_models_final.rds"); r <- readRDS("sigma_eps2_profile_final_raw.rds"); subj <- f$meta$subject; q <- qt(0.975, length(unique(subj)) - 1)
out <- do.call(rbind, lapply(names(r$fits), function(k) { fit <- r$fits[[k]]; cr2 <- compute_SE_regression_mv_cluster_robust(f$Y, f$Xmat, 0, fit, subj)
  data.frame(sd_eps = as.numeric(k), parameter = colnames(f$Xmat), est = fit$beta, se_cr0 = cr2$se[colnames(f$Xmat)], lo = fit$beta - q * cr2$se[colnames(f$Xmat)], hi = fit$beta + q * cr2$se[colnames(f$Xmat)], loglik = fit$loglik, converged = fit$converged, iter = fit$iter) }))
print(out[out$parameter == "Diabetic", ], digits = 3, row.names = FALSE); write.csv(out, "sigma_eps2_profile_final_intervals.csv", row.names = FALSE)
u <- r$unconstrained; cat(sprintf("unconstrained: converged %s after %d iterations; sigma_eps2 %.5f; loglik %.1f\n", u$converged, u$iter, u$sigma_eps2, u$loglik))
