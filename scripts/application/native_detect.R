## Corrected native-reading detection for the CGMacros per-minute Dexcom column.
## (1) Kinks: the interpolation knot is the row at which the incoming and
##     outgoing slopes differ, i.e. |d[j+1] - d[j]| > tol flags row j (the
##     original code flagged j+1 and a spurious row 2).
## (2) Grid phase: within each contiguous block (no gap > max_gap_min), the
##     device grid is every 5 minutes at a fixed phase; the phase is the modal
##     minute-mod-5 of the kinks. Every row on that phase is a native slot;
##     each is classified "kink" (confirmed native) or "collinear" (on-grid but
##     no slope change: a native reading collinear with its neighbours OR a
##     device dropout interpolated across -- indistinguishable from the file).
## Blocks with fewer than 3 interior kinks cannot establish a phase: their
## boundary/kink rows are returned with kink = NA and phase_ok = NA (status
## unresolved) and must not be used as confirmed readings.
## Returns a data frame of native slots with ts, val, kink flag, block id.
detect_native_grid <- function(ts, val, tol = 0.05, max_gap_min = 1.5, period = 5) {
  ok <- !is.na(val) & !is.na(ts); ts <- ts[ok]; val <- val[ok]; o <- order(ts); ts <- ts[o]; val <- val[o]
  n <- length(ts)
  ## fewer than 3 readings (including empty input): nothing can be established; return rows with unknown status
  if (n < 3) return(data.frame(ts = ts, val = val, kink = rep(NA, n), boundary = rep(TRUE, n), block = rep(1L, n), phase_ok = rep(NA, n)))
  gap <- as.numeric(diff(ts), units = "mins"); block <- cumsum(c(TRUE, gap > max_gap_min))
  d <- diff(val); kink <- c(FALSE, abs(diff(d)) > tol, FALSE)   # flags row j for |d[j+1]-d[j]| > tol, j = 2..n-1
  first <- !duplicated(block); last <- !duplicated(block, fromLast = TRUE); kink[first | last] <- TRUE
  mins <- as.numeric(ts, units = "secs") / 60
  out <- vector("list", max(block))
  for (b in seq_len(max(block))) {
    i <- which(block == b); k <- i[kink[i] & !first[i] & !last[i]]
    if (length(k) >= 3) {
      ph <- as.integer(round(mins[k])) %% period; tab <- tabulate(ph + 1, period); phase <- which.max(tab) - 1
      on_grid <- i[(as.integer(round(mins[i])) %% period) == phase]
      bd <- first[on_grid] | last[on_grid]
      ## a block-boundary row on the grid phase is NOT a confirmed knot (no slope change can be seen there): kink = NA, boundary = TRUE
      out[[b]] <- data.frame(ts = ts[on_grid], val = val[on_grid], kink = ifelse(bd, NA, kink[on_grid]), boundary = bd, block = b,
                             phase_ok = max(tab) / sum(tab) )
    } else out[[b]] <- data.frame(ts = ts[i[kink[i]]], val = val[i[kink[i]]], kink = NA, boundary = (first | last)[i[kink[i]]], block = b, phase_ok = NA)  # phase unresolved: status unknown, NOT a confirmed knot
  }
  do.call(rbind, out)
}
