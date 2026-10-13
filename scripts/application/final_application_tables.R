## Assemble final Oct13 application tables and aggregate figure inputs.
## In particular, the primary 65-minute row reuses Step 51's bootstrap.
f <- readRDS("PRIMARY_boundary.rds")
joint <- read.csv("primary_boundary_joint.csv", check.names = FALSE)
sp <- list(`65` = readRDS("setpoint_boundary_65.rds"), `120` = readRDS("setpoint_boundary_120.rds"))
su <- read.csv("setpoint_bootstrap_uniform.csv", check.names = FALSE)
direct <- read.csv("setpoint_boundary_directML_summary.csv", check.names = FALSE)
prof <- read.csv("sigma_eps2_profile_final.csv", check.names = FALSE)
pi <- read.csv("sigma_eps2_profile_final_intervals.csv", check.names = FALSE)
diabetes <- joint[joint$parameter == "Diabetic", ]
stopifnot(nrow(diabetes) == 1L)
rows <- list()
for (dm in c(65, 120)) {
  order <- c("pre-meal", "pre-meal, windows of B", "180-min",
             if (dm == 65) "next-meal", "estimated set point, deviations", "estimated set point, full")
  for (v in order) {
    fixed <- sp[[as.character(dm)]]$res
    fixed <- fixed[fixed$parameter == "Diabetic" & fixed$variant == v, ]
    boot <- su[su$len == dm & su$variant == v, ]
    if (dm == 65 && v == "pre-meal") {
      row <- data.frame(len = dm, variant = v, n = nrow(f$Y), est = diabetes$estimate,
                         boot_lo = diabetes$boot_pct_lo, boot_hi = diabetes$boot_pct_hi,
                         cr0_lo = diabetes$cr0_lo_t, cr0_hi = diabetes$cr0_hi_t, B = nrow(f$boot_joint),
                         bootstrap_source = "Step 51 primary participant bootstrap")
    } else if (v == "pre-meal, windows of B") {
      stopifnot(nrow(fixed) == 1L)
      row <- data.frame(len = dm, variant = v, n = fixed$n, est = fixed$est,
                         boot_lo = NA_real_, boot_hi = NA_real_, cr0_lo = fixed$lo, cr0_hi = fixed$hi,
                         B = NA_integer_, bootstrap_source = "not computed")
    } else {
      stopifnot(nrow(boot) == 1L)
      row <- data.frame(len = dm, variant = v, n = boot$n, est = boot$est,
                         boot_lo = boot$boot_lo, boot_hi = boot$boot_hi,
                         cr0_lo = if (nrow(fixed)) fixed$lo else NA_real_,
                         cr0_hi = if (nrow(fixed)) fixed$hi else NA_real_, B = boot$B,
                         bootstrap_source = "Step 58 uniform participant bootstrap")
    }
    if (grepl("^estimated", v)) {
      detail <- direct[direct$len == dm & direct$variant == v & direct$parameter == "Diabetic", ]
      stopifnot(nrow(detail) == 1L)
      row$cr0_lo <- row$cr0_hi <- NA_real_
      row$phi_H <- detail$phi_H; row$phi_P <- detail$phi_P; row$phi_D <- detail$phi_D
      row$loglik <- detail$loglik
    } else if (dm == 65 && v == "pre-meal") {
      group_phi <- tapply(as.vector(tanh(f$Xmat %*% f$fit$beta)), f$meta$group, mean)
      row$phi_H <- unname(group_phi["healthy"]); row$phi_P <- unname(group_phi["prediabetic"])
      row$phi_D <- unname(group_phi["diabetic"]); row$loglik <- f$fit$loglik
    } else {
      stopifnot(nrow(fixed) == 1L)
      row$phi_H <- fixed$phi_H; row$phi_P <- fixed$phi_P; row$phi_D <- fixed$phi_D
      row$loglik <- fixed$loglik
    }
    rows[[length(rows) + 1L]] <- row
  }
}
references <- do.call(rbind, rows)
write.csv(references, "table_S14_reference_sensitivity.csv", row.names = FALSE)
write.csv(joint[, c("parameter", "estimate", "se_model", "se_cr0", "se_md", "se_boot",
                    "cr0_lo_t", "cr0_hi_t", "boot_pct_lo", "boot_pct_hi")],
          "table_S13_primary_regression.csv", row.names = FALSE)
profile <- data.frame(sd_eps = prof$sd_eps, sigma_eta2 = prof$sigma_eta2, loglik = prof$loglik,
                       phi_mean = prof$phi_mean, Diabetic = prof$Diabetic)
ci <- pi[pi$parameter == "Diabetic", ]
profile$cr0_lo <- ci$lo[match(profile$sd_eps, ci$sd_eps)]
profile$cr0_hi <- ci$hi[match(profile$sd_eps, ci$sd_eps)]
profile <- rbind(data.frame(sd_eps = 0, sigma_eta2 = f$fit$sigma_eta2, loglik = f$fit$loglik,
                            phi_mean = mean(tanh(f$Xmat %*% f$fit$beta)), Diabetic = diabetes$estimate,
                            cr0_lo = diabetes$cr0_lo_t, cr0_hi = diabetes$cr0_hi_t), profile)
write.csv(profile, "table_S12_variance_profile.csv", row.names = FALSE)
flow <- lapply(c(65, 120), function(dm) {
  w <- readRDS(sprintf("meal_windows_final_%dmin.rds", dm))
  data.frame(len = dm, step = c("logged_meals", "search_incomplete", "peak_at_end", "peak_at_onset",
                                "not_isolated", "grid_window_exclusion", "retained_episodes", "participants"),
              count = c(unname(w$flow[c("logged_meals", "search_incomplete", "peak_at_end", "peak_at_onset",
                                        "not_isolated")]), w$flow["candidates"] - w$flow["kept"],
                         nrow(w$Y), length(unique(w$meta$subject))))
})
write.csv(do.call(rbind, flow), "table_S11_cohort_flow.csv", row.names = FALSE)
labels <- c("Pre-meal (primary)", "Next meal", "180 min post-meal", "Estimated deviations*", "Estimated, intercept free")
variants <- c("pre-meal", "next-meal", "180-min", "estimated set point, deviations", "estimated set point, full")
sel <- references[references$len == 65, ]
sel <- sel[match(variants, sel$variant), ]
choices <- data.frame(choice = labels, est = sel$est, lo = sel$boot_lo, hi = sel$boot_hi,
                       n = sel$n, interval = "Bootstrap", B = sel$B,
                       display = paste0(labels, "\nN = ", sel$n))
write.csv(choices, "fig_reference_choices.csv", row.names = FALSE)
grp_lab <- c(healthy = "Healthy", prediabetic = "Prediabetic", diabetic = "Type 2 diabetes")
keys <- c("subject", "meal_time", "peak_time", "baseline_gl")
stopifnot(all(vapply(keys, function(k) isTRUE(all.equal(f$meta[[k]], sp[["65"]]$meta[[k]])), logical(1))))
meta <- sp[["65"]]$meta
keep <- is.finite(meta$gl_180) & is.finite(meta$baseline_gl)
stopifnot(sum(keep) == 1171L)
Y <- f$Y[keep, , drop = FALSE]
groups <- grp_lab[as.character(f$meta$group[keep])]
levels <- do.call(rbind, lapply(unname(grp_lab), function(group) {
  ix <- groups == group
  data.frame(group = group, n = sum(ix), reference = mean((meta$gl_180 - meta$baseline_gl)[keep][ix]))
}))
means <- do.call(rbind, lapply(unname(grp_lab), function(group) {
  data.frame(group = group, t = (0:13) * 5, y = colMeans(Y[groups == group, , drop = FALSE]))
}))
write.csv(means, "fig_reference_descriptive_means.csv", row.names = FALSE)
write.csv(levels, "fig_reference_descriptive_levels.csv", row.names = FALSE)
cat("Wrote Tables S11-S14 and the three aggregate recovery-reference figure inputs.\n")
