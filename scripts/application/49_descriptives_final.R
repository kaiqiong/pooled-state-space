## Step 49: descriptive numbers and per-participant table for the final windows.
f <- readRDS("FINAL_fitted_models_final.rds"); w <- readRDS("meal_windows_final_65min.rds"); w2 <- readRDS("meal_windows_final_120min.rds"); m <- f$meta
cat("flow 65:\n"); print(w$flow); cat("flow 120:\n"); print(w2$flow)
cat("meal types:"); print(table(tolower(sub(" 1$", "", sub("s$", "", m$meal_type)))))
n <- table(m$subject); g <- tapply(as.character(m$group), m$subject, function(z) z[1])
print(do.call(rbind, lapply(split(as.integer(n), g[names(n)]), function(x) c(participants = length(x), windows = sum(x), mean = round(mean(x), 1), median = median(x), min = min(x), max = max(x)))))
cat(sprintf("all: mean %.1f median %.1f range %d-%d\n", mean(n), median(n), min(n), max(n)))
z <- as.vector(f$naive_fit$z); ok <- is.finite(z); sj <- factor(m$subject[ok]); zz <- z[ok]; nn <- table(sj); a <- anova(lm(zz ~ sj)); n0 <- (sum(nn) - sum(nn^2) / sum(nn)) / (length(nn) - 1)
icc <- (a[1, 3] - a[2, 3]) / (a[1, 3] + (n0 - 1) * a[2, 3]); r <- resid(lm(zz ~ Age + Sex + BMI + Prediabetic + Diabetic, data = as.data.frame(f$Xmat)[ok, ])); a2 <- anova(lm(r ~ sj)); icc2 <- (a2[1, 3] - a2[2, 3]) / (a2[1, 3] + (n0 - 1) * a2[2, 3])
cat(sprintf("ICC %.3f (adj %.3f) n0 %.1f DE %.2f; non-finite naive z %d\n", icc, icc2, n0, 1 + (n0 - 1) * icc, sum(!ok)))
cat("peak minute: median", median(m$peak_minute), " IQR", quantile(m$peak_minute, c(.25, .75)), " range", range(m$peak_minute), "\n")
cat("peak excess (first slot) by group:", round(tapply(f$Y[, 1], m$group, mean), 1), " | excess at 65 min:", round(tapply(f$Y[, 14], m$group, mean), 1), "\n")
cat("collinear share of slots:", round(mean(w$flags == "collinear"), 4), " | windows by number of collinear slots 0/1/2/3/4+:", table(pmin(m$n_collinear, 4)), "\n")
cat("120-min: by group", table(g[w2$meta$subject]), " collinear share", round(mean(w2$flags == "collinear"), 4), " excess at 120 min by group:", round(tapply((w2$Y - w2$meta$baseline_gl)[, 25], g[w2$meta$subject], mean), 1), "\n")
q <- read.csv("native_detection_qc_v2.csv"); bs <- merge(w$by_subject, q[, c("subject", "n_native", "frac_kink")], by = "subject", all.y = TRUE); bs$n_candidate[is.na(bs$n_candidate)] <- 0; bs$n_kept[is.na(bs$n_kept)] <- 0
gg <- c(healthy = "H", prediabetic = "P", diabetic = "D")[g[bs$subject]]; gg[is.na(gg)] <- "P"; bs$g <- gg
## logged meals per participant
lm_n <- sapply(bs$subject, function(sid) { df <- read.csv(file.path(APP_DATA, sid, paste0(sid, ".csv")), stringsAsFactors = FALSE, check.names = FALSE); sum(!is.na(df[["Meal Type"]]) & df[["Meal Type"]] != "") }); bs$n_meals <- lm_n
rows <- sprintf("%s & %s & %d & %d & %d & %.0f", sub("CGMacros-", "", bs$subject), bs$g, bs$n_meals, bs$n_kept, bs$n_native, 100 * bs$frac_kink)
out <- paste0(rows[1:23], " & ", c(rows[24:45], "& & & & &"), " \\\\"); writeLines(out, "retention_rows_final.tex"); bs$group <- g[bs$subject]; write.csv(bs, "per_subject_pipeline_table_final.csv", row.names = FALSE)
cat("grid readings total", sum(q$n_native), " nonmissing rows", sum(q$n_nonmissing), " share of grid readings that are knots (whole records):", round(weighted.mean(q$frac_kink, q$n_native), 4), "\n"); cat(out[7], "\n")
