## Step 52: final fixed-reference fits at sigma_eps2 = 0.
## Package adaptation omits the historical degenerate EM estimated-reference
## fits and pilot bootstrap. Step 53 estimates those references by direct ML;
## Step 58 supplies the reported B=2000 bootstrap intervals for every variant.
source(file.path(APP_ROOT, "R", "state_space_em.R"))
f <- readRDS("FINAL_fitted_models_final.rds"); sub <- unique(f$meta[, c("subject", "Age", "Gender", "BMI", "group")]); rownames(sub) <- sub$subject
mA <- mean(f$meta$Age); sA <- sd(f$meta$Age); mB <- mean(f$meta$BMI); sB <- sd(f$meta$BMI); se2 <- 0   # option 1: measurement variance at its boundary
later_baselines <- function(meta) { meta$gl_180 <- NA; meta$gl_next <- NA
  for (sid in unique(meta$subject)) { df <- read.csv(file.path(APP_DATA, sid, paste0(sid, ".csv")), stringsAsFactors = FALSE, check.names = FALSE)
    ts <- as.POSIXct(df$Timestamp, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"); gl <- suppressWarnings(as.numeric(df[["Dexcom GL"]])); ml <- df[["Meal Type"]]; o <- order(ts); ts <- ts[o]; gl <- gl[o]; ml <- ml[o]; mt <- ts[!is.na(ml) & ml != ""]
    for (i in which(meta$subject == sid)) { t0 <- meta$meal_time[i]; j <- which.min(abs(ts - (t0 + 180 * 60))); if (abs(as.numeric(ts[j] - (t0 + 180 * 60), units = "mins")) <= 5) meta$gl_180[i] <- gl[j]
      nx <- mt[mt > t0 + 60]; if (length(nx)) { k <- which.min(abs(ts - min(nx))); if (abs(as.numeric(ts[k] - min(nx), units = "mins")) <= 5) meta$gl_next[i] <- gl[k] } } }
  meta }
run_len <- function(dm) {
  w <- readRDS(sprintf("meal_windows_final_%dmin.rds", dm)); s <- sub[w$meta$subject, ]
  X <- cbind(Intercept = 1, Age = (s$Age - mA) / sA, Sex = as.numeric(s$Gender == "F"), BMI = (s$BMI - mB) / sB, Prediabetic = as.numeric(s$group == "prediabetic"), Diabetic = as.numeric(s$group == "diabetic"))
  subj <- w$meta$subject; G <- length(unique(subj)); q <- qt(0.975, G - 1); meta <- later_baselines(w$meta); Yraw <- w$Y; Y0 <- Yraw - meta$baseline_gl
  cat(sprintf("\n=== %d-min windows: N=%d | excess at last slot by group: %s | 180-min minus pre-meal: %s\n", dm, nrow(Y0), paste(round(tapply(Y0[, ncol(Y0)], s$group, mean), 1), collapse = "/"), paste(round(tapply(meta$gl_180 - meta$baseline_gl, s$group, function(z) mean(z, na.rm = TRUE)), 1), collapse = "/")))
  fitp <- function(Y, idx) em_regression_mv(Y[idx, ], X[idx, ], 0, list(beta = f$fit$beta, a1 = mean(Y[idx, 1]), P1 = max(var(Y[idx, 1]), 1), sigma_eta2 = 40, sigma_eps2 = se2), thres = 1e-7, max_iter = 4000, fix_sigma_eps2 = se2)
  summ <- function(fit, Y, idx, lab) { cr2 <- compute_SE_regression_mv_cluster_robust(Y[idx, ], X[idx, ], 0, fit, subj[idx]); ph <- tanh(X[idx, ] %*% fit$beta)
    data.frame(len = dm, variant = lab, parameter = colnames(X), est = fit$beta, se = cr2$se[colnames(X)], lo = fit$beta - q * cr2$se[colnames(X)], hi = fit$beta + q * cr2$se[colnames(X)], n = length(idx),
               phi_H = mean(ph[s$group[idx] == "healthy"]), phi_P = mean(ph[s$group[idx] == "prediabetic"]), phi_D = mean(ph[s$group[idx] == "diabetic"]), a1 = fit$a1, sigma_eta2 = fit$sigma_eta2, loglik = fit$loglik, conv = fit$converged) }
  res <- list(); all <- seq_len(nrow(Y0))
  fP <- fitp(Y0, all); res$P <- summ(fP, Y0, all, "pre-meal")
  iB <- which(!is.na(meta$gl_180)); YB <- Yraw - meta$gl_180; fPB <- fitp(Y0, iB); fB <- fitp(YB, iB); res$P_B <- summ(fPB, Y0, iB, "pre-meal, windows of B"); res$B <- summ(fB, YB, iB, "180-min"); fC <- NULL; iC <- NULL
  if (dm == 65) { iC <- which(!is.na(meta$gl_next)); YC <- Yraw - meta$gl_next; fC <- fitp(YC, iC); res$C <- summ(fC, YC, iC, "next-meal") }
  out <- do.call(rbind, res); rownames(out) <- NULL
  print(out[out$parameter == "Diabetic", c("variant", "n", "est", "lo", "hi", "phi_H", "phi_P", "phi_D", "a1", "loglik", "conv")], digits = 3, row.names = FALSE)
  list(res = out, fits = list(P = fP, P_B = fPB, B = fB, C = fC),
       idx = list(B = iB, C = iC), X = X, subj = subj, meta = meta)
}
r65 <- run_len(65); saveRDS(r65, "setpoint_boundary_65.rds"); r120 <- run_len(120); saveRDS(r120, "setpoint_boundary_120.rds")
write.csv(rbind(r65$res, r120$res), "setpoint_boundary_summary.csv", row.names = FALSE)
cat("saved fixed-reference fits; estimated references follow in Step 53\n")
