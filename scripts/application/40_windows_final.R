## Step 40: FINAL window construction (replaces Steps 8, 18, 20, 27). Adds to Step 27 the rule that the
## located peak must lie strictly inside the search window: peaks within one minute of meal onset are
## dropped (no post-meal rise was observed), symmetric with the existing rule for the end of the window.
## Step 27 rules: Changes from v2: (i) the
## peak-search window must have complete coverage (every minute 0-120 after
## onset present and non-missing), as in the original Step 8; (ii) slot
## timestamps and phase status are stored; (iii) a window is retained only if
## all slots come from ONE phase-resolved block, so that every interval is
## exactly 5 minutes; windows failing (iii) are counted by reason.
## Uses detect_native_grid() (native_detect.R): correct kink indexing plus the
## 5-min grid phase per contiguous block, so collinear native readings are no
## longer missed. Same rules as before: peak on the interpolated series within
## 0-120 min of onset (peaks at the search boundary dropped); no other meal
## between onset and window end; slot grid from the peak, slot filled by the
## grid reading within +/-2 min, else gap-filled from the interpolated column;
## at most MAX_GAPFILL gap-filled slots. Records per-slot flags: native-kink,
## native-collinear, gap-filled. Run with DECAY_MIN = 65 (T=14) and 120 (T=25).
DECAY_MIN <- APP_DECAY_MIN
source(file.path(APP_ROOT, "scripts", "application", "native_detect.R"))
data_dir <- APP_DATA; SEARCH_MIN <- 120; N_SLOT <- DECAY_MIN / 5 + 1; MAX_GAPFILL <- round(N_SLOT / 7); SLOT_TOL_MIN <- 2
cat(sprintf("DECAY_MIN=%d T=%d MAX_GAPFILL=%d\n", DECAY_MIN, N_SLOT, MAX_GAPFILL))
dirs <- list.dirs(data_dir, recursive = FALSE); dirs <- dirs[grepl("CGMacros-[0-9]+$", dirs)]
windows <- list(); flags <- list(); meta <- list(); ret <- list(); slot_ts <- list(); n_search_incomplete <- 0L; n_peak_at_end <- 0L; n_peak_at_onset <- 0L; n_meals <- 0L
for (d in dirs) {
  sid <- basename(d); df <- read.csv(file.path(d, paste0(sid, ".csv")), stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("Timestamp", "Dexcom GL", "Meal Type") %in% names(df))) next
  ts <- as.POSIXct(df$Timestamp, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"); gl <- suppressWarnings(as.numeric(df[["Dexcom GL"]])); meal <- df[["Meal Type"]]
  o <- order(ts); ts <- ts[o]; gl <- gl[o]; meal <- meal[o]
  nat <- detect_native_grid(ts, gl); meal_idx <- which(!is.na(meal) & meal != ""); meal_times <- ts[meal_idx]
  for (k in seq_along(meal_idx)) {
    n_meals <- n_meals + 1L; t0 <- ts[meal_idx[k]]; search_end <- t0 + SEARCH_MIN * 60; sw <- which(ts >= t0 & ts <= search_end)
    ## coverage requirement: all 121 per-minute rows of the search window present and non-missing
    if (length(sw) < SEARCH_MIN + 1 || any(is.na(gl[sw]))) { n_search_incomplete <- n_search_incomplete + 1L; next }
    pk <- sw[which.max(gl[sw])]; tp <- ts[pk]
    if (tp >= search_end - 60) { n_peak_at_end <- n_peak_at_end + 1L; next }     # peak at the end of the search window: ambiguous
    if (tp <= t0 + 60) { n_peak_at_onset <- n_peak_at_onset + 1L; next }         # peak at meal onset: no post-meal rise observed
    dend <- tp + DECAY_MIN * 60; if (any(meal_times > t0 & meal_times <= dend)) next
    slot_t <- tp + (0:(N_SLOT - 1)) * 5 * 60; vals <- numeric(N_SLOT); fl <- character(N_SLOT); ok <- TRUE; st <- rep(as.POSIXct(NA, tz = "UTC"), N_SLOT); blk <- rep(NA_integer_, N_SLOT)
    for (s in seq_len(N_SLOT)) {
      dm <- as.numeric(abs(nat$ts - slot_t[s]), units = "mins"); j <- which.min(dm)
      if (length(j) && dm[j] <= SLOT_TOL_MIN) { vals[s] <- nat$val[j]; st[s] <- nat$ts[j]; blk[s] <- nat$block[j]; fl[s] <- if (is.na(nat$phase_ok[j])) "unresolved" else if (nat$boundary[j]) "boundary" else if (nat$kink[j]) "kink" else "collinear" }
      else { r <- which.min(abs(ts - slot_t[s])); if (is.na(gl[r]) || abs(as.numeric(ts[r] - slot_t[s], units = "mins")) > 1) { ok <- FALSE; break }; vals[s] <- gl[r]; fl[s] <- "gapfill" }
    }
    n_gf <- sum(fl == "gapfill"); iv <- as.numeric(diff(st), units = "mins")
    unresolved <- any(fl %in% c("unresolved", "boundary")); irregular <- ok && n_gf == 0 && (any(abs(iv - 5) > 1e-6) || length(unique(blk)) > 1)
    reason <- if (!ok || n_gf > 0) "gapfill_needed" else if (unresolved) "unresolved_phase" else if (irregular) "irregular_interval" else "kept"
    complete <- reason == "kept"
    ret[[length(ret) + 1]] <- data.frame(subject = sid, meal_time = t0, peak_time = tp, n_gapfilled = n_gf, complete = complete, reason = reason)
    if (!complete) next
    b <- which.min(abs(ts - t0)); windows[[length(windows) + 1]] <- vals; flags[[length(flags) + 1]] <- fl; slot_ts[[length(slot_ts) + 1]] <- as.numeric(st)
    meta[[length(meta) + 1]] <- data.frame(subject = sid, meal_time = t0, meal_type = meal[meal_idx[k]], peak_time = tp, peak_minute = as.numeric(tp - t0, units = "mins"),
                                            baseline_gl = gl[b], n_gapfilled = n_gf, n_collinear = sum(fl == "collinear"))
  }
}
Y <- do.call(rbind, windows); FL <- do.call(rbind, flags); meta <- do.call(rbind, meta); ret <- do.call(rbind, ret)
by_subj <- aggregate(complete ~ subject, ret, function(z) c(n_candidate = length(z), n_kept = sum(z))); by_subj <- data.frame(subject = by_subj$subject, by_subj$complete); by_subj$retention_rate <- by_subj$n_kept / by_subj$n_candidate
n_not_isolated <- n_meals - n_search_incomplete - n_peak_at_end - n_peak_at_onset - nrow(ret)
cat(sprintf("logged meals %d | incomplete peak-search coverage %d | peak at end of search %d | peak at onset %d | another meal before window end %d\n", n_meals, n_search_incomplete, n_peak_at_end, n_peak_at_onset, n_not_isolated))
cat(sprintf("events dropped for incomplete peak-search coverage: %d\ncandidates %d, retained %d (%.1f%%), participants %d\n", n_search_incomplete, nrow(ret), nrow(Y), 100 * nrow(Y) / nrow(ret), length(unique(meta$subject)))); cat("exclusions by reason:"); print(table(ret$reason))
cat("slots: kink", round(mean(FL == "kink"), 3), " collinear", round(mean(FL == "collinear"), 3), " gap-filled", round(mean(FL == "gapfill"), 4), "\n")
cat("windows by n_gapfilled:"); print(table(meta$n_gapfilled)); cat("windows by n_collinear:"); print(table(meta$n_collinear))
out <- sprintf("meal_windows_final_%dmin.rds", DECAY_MIN)
saveRDS(list(Y = Y, flags = FL, slot_time = do.call(rbind, slot_ts), flow = c(logged_meals = n_meals, search_incomplete = n_search_incomplete, peak_at_end = n_peak_at_end, peak_at_onset = n_peak_at_onset, not_isolated = n_not_isolated, candidates = nrow(ret), kept = nrow(Y)), rules = list(search_min = SEARCH_MIN, decay_min = DECAY_MIN, slot_tol_min = SLOT_TOL_MIN, max_gapfill = MAX_GAPFILL, kink_tol = 0.05, max_gap_min = 1.5), n_search_incomplete = n_search_incomplete, meta = meta, retention = ret, by_subject = by_subj, decay_min = DECAY_MIN, n_slots = N_SLOT, max_gapfill = MAX_GAPFILL), out); cat("saved", out, "\n")
