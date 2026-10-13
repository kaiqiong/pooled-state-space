## Step 50 (final windows; derived from Step 25): schematic example (Fig S2) on the corrected v2 windows: same event
## as Step 17 (CGMacros-006 lunch 2023-04-08 12:06) with the corrected detector.
source(file.path(APP_ROOT, "scripts", "application", "native_detect.R")); f <- readRDS("FINAL_fitted_models_final.rds"); w <- readRDS("meal_windows_final_65min.rds")
i <- which(w$meta$subject == "CGMacros-006" & format(w$meta$meal_time) == "2023-04-08 12:06:00"); stopifnot(length(i) == 1)
sid <- "CGMacros-006"; df <- read.csv(file.path(APP_DATA, sid, paste0(sid, ".csv")), stringsAsFactors = FALSE, check.names = FALSE)
ts <- as.POSIXct(df$Timestamp, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"); gl <- suppressWarnings(as.numeric(df[["Dexcom GL"]])); o <- order(ts); ts <- ts[o]; gl <- gl[o]
t0 <- w$meta$meal_time[i]; tp <- w$meta$peak_time[i]; sel <- ts >= t0 - 30 * 60 & ts <= t0 + 200 * 60
nat <- detect_native_grid(ts[sel], gl[sel]); dir <- ""
write.csv(data.frame(t_min = as.numeric(ts[sel] - t0, units = "mins"), gl = gl[sel]), paste0(dir, "schematic_example_series.csv"), row.names = FALSE)
write.csv(data.frame(t_min = as.numeric(nat$ts - t0, units = "mins"), gl = nat$val, kink = nat$kink), paste0(dir, "schematic_example_native.csv"), row.names = FALSE)
write.csv(data.frame(slot_min = (0:13) * 5, excess = f$Y[i, ], flag = w$flags[i, ]), paste0(dir, "schematic_example_window.csv"), row.names = FALSE)
mt <- ts[!is.na(df[["Meal Type"]][o]) & df[["Meal Type"]][o] != ""]; nx <- min(mt[mt > t0])
write.csv(data.frame(subject = sid, meal_type = w$meta$meal_type[i], meal_time = format(t0), peak_minute = w$meta$peak_minute[i], baseline_gl = w$meta$baseline_gl[i], next_meal_min = as.numeric(nx - t0, units = "mins"), group = "healthy", phi_hat = tanh(sum(f$Xmat[i, ] * f$fit$beta)), n_collinear = w$meta$n_collinear[i]), paste0(dir, "schematic_example_meta.csv"), row.names = FALSE)
cat("event flags:", w$flags[i, ], "\n")
