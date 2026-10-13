## Aggregate native-grid QC using the same final detector as Step 40.
source(file.path(APP_ROOT, "scripts", "application", "native_detect.R"))
dirs <- list.dirs(APP_DATA, recursive = FALSE)
dirs <- dirs[grepl("CGMacros-[0-9]+$", dirs)]
rows <- lapply(dirs, function(d) {
  sid <- basename(d)
  df <- read.csv(file.path(d, paste0(sid, ".csv")), stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("Timestamp", "Dexcom GL") %in% names(df))) return(NULL)
  ts <- as.POSIXct(df$Timestamp, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  gl <- suppressWarnings(as.numeric(df[["Dexcom GL"]]))
  nat <- detect_native_grid(ts, gl)
  data.frame(subject = sid, n_nonmissing = sum(!is.na(gl)), n_native = nrow(nat),
             frac_kink = mean(nat$kink, na.rm = TRUE))
})
write.csv(do.call(rbind, rows), "native_detection_qc_v2.csv", row.names = FALSE)
