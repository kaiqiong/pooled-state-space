## Step 42 (final windows; derived from Step 35): independent verification of the v3 windows from the raw files and
## the saved slot timestamps (does not reuse the builder's own bookkeeping).
source(file.path(APP_ROOT, "scripts", "application", "native_detect.R"))
for (dm in c(65, 120)) {
  w <- readRDS(sprintf("meal_windows_final_%dmin.rds", dm)); m <- w$meta; st <- w$slot_time; Tn <- ncol(w$Y)
  iv <- t(apply(st, 1, diff)) / 60
  cat(sprintf("\n=== %d-min windows: N = %d, T = %d, participants = %d\n", dm, nrow(w$Y), Tn, length(unique(m$subject))))
  cat("flags present:", paste(names(table(w$flags)), table(w$flags), collapse = "; "), "\n")
  cat("intervals != 5 min:", sum(abs(iv - 5) > 1e-9), " | first slot offset from peak time (s) range:", range(st[, 1] - as.numeric(m$peak_time)), "\n")
  chk <- do.call(rbind, lapply(unique(m$subject), function(sid) {
    df <- read.csv(file.path(APP_DATA, sid, paste0(sid, ".csv")), stringsAsFactors = FALSE, check.names = FALSE)
    ts <- as.POSIXct(df$Timestamp, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"); gl <- suppressWarnings(as.numeric(df[["Dexcom GL"]])); o <- order(ts); ts <- ts[o]; gl <- gl[o]
    ok <- !is.na(gl); tso <- ts[ok]; glo <- gl[ok]; blk <- cumsum(c(TRUE, as.numeric(diff(tso), units = "mins") > 1.5)); first <- !duplicated(blk); last <- !duplicated(blk, fromLast = TRUE)
    d <- diff(glo); kink <- c(FALSE, abs(diff(d)) > 0.05, FALSE)
    do.call(rbind, lapply(which(m$subject == sid), function(i) { j <- match(st[i, ], as.numeric(tso)); sw <- which(ts >= m$meal_time[i] & ts <= m$meal_time[i] + 7200)
      b <- blk[j]; ik <- which(blk == b[1] & kink & !first & !last); ph <- (as.integer(round(as.numeric(tso[ik]) / 60)) %% 5)
      data.frame(i = i, matched = !any(is.na(j)), value_ok = !any(is.na(j)) && max(abs(glo[j] - w$Y[i, ])) < 1e-9, one_block = length(unique(b)) == 1, n_interior_kinks = length(ik), phase_unique = length(unique(ph)) == 1,
                 on_phase = all((as.integer(round(st[i, ] / 60)) %% 5) == ph[1]), n_boundary = sum(first[j] | last[j]), n_true_kink = sum(kink[j] & !first[j] & !last[j]), flag_kink = sum(w$flags[i, ] == "kink"),
                 search_rows = length(sw), search_na = sum(is.na(gl[sw])), baseline_ok = abs(gl[which.min(abs(ts - m$meal_time[i]))] - m$baseline_gl[i]) < 1e-9) })) }))
  cat("slots matched to raw rows:", all(chk$matched), " values equal raw:", all(chk$value_ok), " single block:", all(chk$one_block), "\n")
  cat("min interior kinks in a window's block:", min(chk$n_interior_kinks), " phase unique:", all(chk$phase_unique), " all slots on phase:", all(chk$on_phase), "\n")
  cat("peak-search rows == 121:", all(chk$search_rows == 121), " NA in search:", sum(chk$search_na > 0), " baseline equals glucose at onset:", all(chk$baseline_ok), "\n")
  cat("windows containing block-boundary rows:", sum(chk$n_boundary > 0), " (slots:", sum(chk$n_boundary), ") | slots flagged kink but not an interior slope change:", sum(chk$flag_kink - chk$n_true_kink), "\n")
  cat("peak minute range:", range(m$peak_minute), " | peaks within 1 min of onset:", sum(m$peak_minute <= 1), " | within 1 min of search end:", sum(m$peak_minute >= 119), "\n")
  write.csv(chk, sprintf("window_verification_final_%dmin.csv", dm), row.names = FALSE)
}
