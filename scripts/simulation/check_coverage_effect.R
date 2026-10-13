# Portable October 13 replication. Run from the repository root; see docs/SIMULATIONS.txt.
## Effect of substituting the better-likelihood solution for each affected bootstrap draw on the boundary tau = 0.4
## percentile intervals and their coverage. Usage: Rscript local_optima_coverage_effect.R <local_optima_rds>
args <- commandArgs(trailingOnly = TRUE); r <- readRDS(args[1]); design <- if (length(args) >= 2) args[2] else "boundary"; tg <- readRDS("results/simulation/tables/expected_likelihood_targets.rds")[[design]]; beta0 <- tg$generating$beta; bstar <- tg$beta_star; pn <- names(beta0)
aff <- r[!is.na(r$improvement) & r$improvement > 0.01, ]; cat("affected draws:", nrow(aff), "in", length(unique(aff$m)), "datasets; max per dataset", max(table(aff$m)), "| by saved status:", paste(names(table(aff$saved_status)), table(aff$saved_status)), "\n")
cat("improvement quantiles (50/90/100%):", round(quantile(aff$improvement, c(.5, .9, 1)), 1), "| max_beta_diff quantiles:", round(quantile(aff$max_beta_diff, c(.5, .9, 1)), 3), "\n")
fs <- list.files(paste0("results/simulation/boot_b999/", design, "_tau0.4"), "^m[0-9]+\\.rds$", full.names = TRUE); lo0 <- hi0 <- lo1 <- hi1 <- matrix(NA, length(fs), 6, dimnames = list(NULL, pn))
for (i in seq_along(fs)) { z <- readRDS(fs[i]); d0 <- z$draws; d1 <- d0; a <- aff[aff$m == z$m, ]; if (nrow(a)) for (k in seq_len(nrow(a))) d1[a$b[k], ] <- a$alt_beta[[k]]
  lo0[i, ] <- apply(d0, 2, quantile, .025); hi0[i, ] <- apply(d0, 2, quantile, .975); lo1[i, ] <- apply(d1, 2, quantile, .025); hi1[i, ] <- apply(d1, 2, quantile, .975) }
cv <- function(lo, hi, t) sweep(lo, 2, t, "<=") & sweep(hi, 2, t, ">=")
print(round(rbind(cov_working_saved = colMeans(cv(lo0, hi0, bstar)), cov_working_corrected = colMeans(cv(lo1, hi1, bstar)), cov_conditional_saved = colMeans(cv(lo0, hi0, beta0)), cov_conditional_corrected = colMeans(cv(lo1, hi1, beta0))), 3))
ch <- (lo0 != lo1) | (hi0 != hi1); dl <- c(lo0 - lo1, hi0 - hi1); chg <- dl != 0   # 10 Oct audit fix: index the 12000 endpoint changes by their own mask, not by the 6000-entry dataset mask
cat("datasets with any endpoint change:", sum(rowSums(ch) > 0), "of", length(fs), "; endpoints changed:", sum(chg), "of", length(dl), "; max |endpoint change|:", signif(max(abs(dl)), 3), "; median |change| over changed endpoints:", signif(median(abs(dl[chg])), 3), "; 90th pct:", signif(quantile(abs(dl[chg]), .9), 3), "\n")
w0 <- hi0 - lo0; w1 <- hi1 - lo1; cat("mean width change by coefficient:", signif(colMeans(w1 - w0), 2), "\n")
cat("dataset-coefficient coverage flips (working target):", sum(cv(lo0, hi0, bstar) != cv(lo1, hi1, bstar)), "; (conditional):", sum(cv(lo0, hi0, beta0) != cv(lo1, hi1, beta0)), "\n")
saveRDS(list(lo0 = lo0, hi0 = hi0, lo1 = lo1, hi1 = hi1, affected = aff), sub("\\.rds$", "_coverage_effect.rds", args[1]))
