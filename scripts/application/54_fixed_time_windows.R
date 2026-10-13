## Step 54: fixed-time windows WITHOUT a peak criterion (diagnostic for the measurement-error variance; reproduces
## in this project the check reported in V2's smoother_peak_selection report). For each logged meal with a
## non-missing onset glucose, the window starts at the first grid position at or after t0 + ANCHOR minutes
## (within 5 minutes), has 14 slots at 5-minute spacing, all from one phase-identified block, with no other meal
## logged in (t0, window end]; centered on the onset glucose. Anchors 60, 90, 120. No coverage/peak rules.
## Then, for each anchor: EM with sigma_eps2 free (from two starts), fits with SD fixed at 0, 1, 2 mg/dL,
## innovation autocorrelation at the boundary, pooled AR(2) residual autocorrelation, and the boundary-fit
## coefficients with the nuisance-aware cluster CR0 SE.
source(file.path(APP_ROOT, "scripts", "application", "native_detect.R")); source(file.path(APP_ROOT, "R", "state_space_em.R"))
data_dir <- APP_DATA; N_SLOT <- 14; dirs <- list.dirs(data_dir, recursive = FALSE); dirs <- dirs[grepl("CGMacros-[0-9]+$", dirs)]
f <- readRDS("FINAL_fitted_models_final.rds"); mA <- mean(f$meta$Age); sA <- sd(f$meta$Age); mB <- mean(f$meta$BMI); sB <- sd(f$meta$BMI)   # standardization constants of the primary analysis
bio <- read.csv(file.path(APP_DATA, "bio.csv"), check.names = FALSE); bio$subject <- sprintf("CGMacros-%03d", bio$subject); bio$group <- cut(bio[["A1c PDL (Lab)"]], breaks = c(-Inf, 5.7, 6.5, Inf), labels = c("healthy", "prediabetic", "diabetic"), right = FALSE)
sub <- bio[, c("subject", "Age", "Gender", "BMI", "group")]; rownames(sub) <- sub$subject; chk <- unique(f$meta[, c("subject", "group")]); stopifnot(all(as.character(sub[chk$subject, "group"]) == as.character(chk$group)))   # all 45 participants, consistent with the primary analysis
build <- function(ANCHOR) { W <- list(); M <- list(); n_meals <- 0L; n_nobase <- 0L; n_meal_between <- 0L; n_grid <- 0L
  for (d in dirs) { sid <- basename(d); df <- read.csv(file.path(d, paste0(sid, ".csv")), stringsAsFactors = FALSE, check.names = FALSE); if (!all(c("Timestamp", "Dexcom GL", "Meal Type") %in% names(df))) next
    ts <- as.POSIXct(df$Timestamp, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"); gl <- suppressWarnings(as.numeric(df[["Dexcom GL"]])); meal <- df[["Meal Type"]]; o <- order(ts); ts <- ts[o]; gl <- gl[o]; meal <- meal[o]
    nat <- detect_native_grid(ts, gl); ok_nat <- !is.na(nat$phase_ok) & nat$phase_ok & !nat$boundary; meal_idx <- which(!is.na(meal) & meal != ""); meal_times <- ts[meal_idx]
    for (k in seq_along(meal_idx)) { n_meals <- n_meals + 1L; t0 <- ts[meal_idx[k]]; b <- which.min(abs(ts - t0)); if (is.na(gl[b]) || abs(as.numeric(ts[b] - t0, units = "mins")) > 1) { n_nobase <- n_nobase + 1L; next }
      a <- t0 + ANCHOR * 60; j0 <- which(ok_nat & nat$ts >= a & nat$ts <= a + 5 * 60); if (!length(j0)) { n_grid <- n_grid + 1L; next }; j0 <- j0[1]
      slot_t <- nat$ts[j0] + (0:(N_SLOT - 1)) * 5 * 60; js <- sapply(slot_t, function(s) { dm <- as.numeric(abs(nat$ts - s), units = "mins"); j <- which.min(dm); if (dm[j] <= 1e-6 && ok_nat[j]) j else NA })
      if (any(is.na(js)) || length(unique(nat$block[js])) > 1) { n_grid <- n_grid + 1L; next }
      if (any(meal_times > t0 & meal_times <= slot_t[N_SLOT])) { n_meal_between <- n_meal_between + 1L; next }
      W[[length(W) + 1]] <- nat$val[js] - gl[b]; M[[length(M) + 1]] <- data.frame(subject = sid, meal_time = t0, start_minute = as.numeric(nat$ts[j0] - t0, units = "mins"), baseline_gl = gl[b]) } }
  list(Y = do.call(rbind, W), meta = do.call(rbind, M), flow = c(logged = n_meals, no_onset_glucose = n_nobase, no_grid_window = n_grid, meal_between = n_meal_between, kept = length(W))) }
ac <- function(e) { e <- sweep(e, 2, colMeans(e)); c(lag1 = cor(as.vector(e[, -1]), as.vector(e[, -ncol(e)])), lag2 = cor(as.vector(e[, -(1:2)]), as.vector(e[, 1:(ncol(e) - 2)]))) }
res <- list(); out <- list()
for (A in c(60, 90, 120)) { w <- build(A); Y <- w$Y; m <- w$meta; s <- sub[m$subject, ]; X <- cbind(Intercept = 1, Age = (s$Age - mA) / sA, Sex = as.numeric(s$Gender == "F"), BMI = (s$BMI - mB) / sB, Prediabetic = as.numeric(s$group == "prediabetic"), Diabetic = as.numeric(s$group == "diabetic")); cl <- m$subject
  cat(sprintf("\n=== anchor %d min: %s | participants %d | mean excess at first/last slot by group: %s / %s\n", A, paste(names(w$flow), w$flow, collapse = ", "), length(unique(cl)), paste(round(tapply(Y[, 1], s$group, mean), 1), collapse = "/"), paste(round(tapply(Y[, N_SLOT], s$group, mean), 1), collapse = "/")))
  init <- function(b0, r) list(beta = c(atanh(b0), rep(0, 5)), a1 = mean(Y[, 1]), P1 = var(Y[, 1]), sigma_eta2 = 35, sigma_eps2 = r)
  free <- lapply(list(c(0.9, 100), c(0.6, 400)), function(z) em_regression_mv(Y, X, 0, init(z[1], z[2]), thres = 1e-9, max_iter = 2500))
  for (fr in free) cat(sprintf("  free: start r=%.0f -> sigma_eps2 %.4f (SD %.3f), sigma_eta2 %.2f, loglik %.2f, iter %d, converged %s, at_cap %s\n", NA, fr$sigma_eps2, sqrt(fr$sigma_eps2), fr$sigma_eta2, fr$loglik, fr$iter, fr$converged, isTRUE(fr$at_cap)))
  fx <- lapply(c(0, 1, 2), function(sd) em_regression_mv(Y, X, 0, init(0.9, sd^2), thres = 1e-9, max_iter = 4000, fix_sigma_eps2 = sd^2))
  ll <- sapply(fx, `[[`, "loglik"); cat(sprintf("  fixed SD 0/1/2: loglik %.2f / %.2f / %.2f (loss %.1f, %.1f); all converged %s\n", ll[1], ll[2], ll[3], ll[1] - ll[2], ll[1] - ll[3], all(sapply(fx, `[[`, "converged"))))
  f0 <- fx[[1]]; phi <- as.vector(tanh(X %*% f0$beta)); e <- Y[, -1] - phi * Y[, -N_SLOT]; a1c <- ac(e)
  y2 <- as.vector(Y[, -(1:2)]); y1 <- as.vector(Y[, 2:(N_SLOT - 1)]); y0 <- as.vector(Y[, 1:(N_SLOT - 2)]); f2 <- lm(y2 ~ y1 + y0 - 1); r2 <- matrix(residuals(f2), nrow(Y)); a2c <- ac(r2); rt <- polyroot(c(coef(f2)[2], coef(f2)[1], -1))   # characteristic roots z of z^2 - c1 z - c2 = 0 (an earlier version printed their reciprocals)
  se <- tryCatch(compute_SE_regression_mv_cluster_beta(Y, X, 0, f0, cl), error = function(e) list(se = c(Diabetic = NA, Prediabetic = NA))); cat(sprintf("  boundary fit: converged %s, at cap %s, max phi %.4f, n_backtrack %d\n", f0$converged, isTRUE(f0$at_cap), max(abs(phi)), f0$n_backtrack))
  cat(sprintf("  boundary fit: phi by group %s | innovation lag-1/2 corr %.3f / %.3f | pooled AR(2) coef %.3f %.3f, roots %s, residual lag-1/2 %.3f / %.3f\n  Diabetic %.3f (CR0 SE %.3f); Prediabetic %.3f (%.3f)\n", paste(round(tapply(phi, s$group, mean), 3), collapse = "/"), a1c[1], a1c[2], coef(f2)[1], coef(f2)[2], paste(round(Mod(rt), 3), collapse = "/"), a2c[1], a2c[2], f0$beta[6], se$se["Diabetic"], f0$beta[5], se$se["Prediabetic"]))
  out[[as.character(A)]] <- data.frame(anchor = A, n = nrow(Y), participants = length(unique(cl)), free_sigma_eps2 = sapply(free, `[[`, "sigma_eps2")[1], free_converged = free[[1]]$converged, free_iter = free[[1]]$iter, fit0_at_cap = isTRUE(f0$at_cap), max_phi = max(abs(phi)), loglik_sd0 = ll[1], loss_sd1 = ll[1] - ll[2], loss_sd2 = ll[1] - ll[3], innov_lag1 = a1c[1], innov_lag2 = a1c[2], ar2_resid_lag1 = a2c[1], ar2_root1 = max(Mod(rt)), ar2_root2 = min(Mod(rt)), phi_H = mean(phi[s$group == "healthy"]), phi_P = mean(phi[s$group == "prediabetic"]), phi_D = mean(phi[s$group == "diabetic"]), Diabetic = f0$beta[6], Diabetic_se_cr0 = se$se["Diabetic"], Prediabetic = f0$beta[5], Prediabetic_se_cr0 = se$se["Prediabetic"])
  res[[as.character(A)]] <- list(Y = Y, X = X, meta = m, flow = w$flow, free = free, fixed = fx, ar2 = coef(f2)) }
tab <- do.call(rbind, out); rownames(tab) <- NULL; options(width = 250); print(tab, digits = 4)
write.csv(tab, "fixed_time_windows_summary.csv", row.names = FALSE); saveRDS(res, "fixed_time_windows.rds"); cat("saved\n")
