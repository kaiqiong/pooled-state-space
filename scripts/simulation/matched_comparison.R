# Portable October 13 replication. Run from the repository root; see docs/SIMULATIONS.txt.
root <- getwd(); out <- file.path(root, "results/simulation/matched")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
options(width = 180)
impl <- root
target_free_path <- "results/simulation/tables/expected_likelihood_target_free_r.rds"
target_fixed_path <- "results/simulation/tables/expected_likelihood_targets.rds"
fixed_summary_path <- "results/simulation/tables/paired_posvar200_comparison.csv"
old_summary_path <- "results/simulation/tables/part_r_summary.csv"
free <- readRDS(target_free_path)$free
fixed <- readRDS(target_fixed_path)$posvar
pn <- names(fixed$generating$beta)
p <- length(pn)
stopifnot(p == 6L, identical(names(free$beta), pn))
all_rows <- list()
fit_diagnostics <- list()
bootstrap_diagnostics <- list()
for (tau in c(0, 0.4)) {
  target <- if (tau == 0) fixed$generating$beta else free$beta
  for (m in 1:500) {
    path <- file.path(impl, sprintf("results/simulation/freevar/tau%s/m%03d.rds", tau, m))
    z <- readRDS(path)
    stopifnot(z$m == m, z$tau == tau, !z$err, isTRUE(z$conv), z$status == "interior",
              identical(names(z$beta), pn), all(is.finite(z$beta)),
              is.finite(z$sigma_eps2), z$sigma_eps2 > 0,
              is.finite(z$min_eig), z$min_eig > 0)
    fit_diagnostics[[length(fit_diagnostics) + 1L]] <- data.frame(
      tau, m, converged = z$conv, status = z$status, observation_variance = z$sigma_eps2,
      min_information_eigenvalue = z$min_eig, in_matched200 = m <= 200)
    add_rows <- function(method, lo, hi, se) {
      stopifnot(length(lo) == p, length(hi) == p, length(se) == p,
                all(is.finite(c(lo, hi, se))), all(se > 0), all(hi > lo))
      data.frame(tau, m, parameter = pn, method, n_bootstrap_draws = if (method == "bootstrap") 999L else NA_integer_,
        in_matched200 = m <= 200, estimate = as.numeric(z$beta), target = as.numeric(target),
        se = as.numeric(se), lower = as.numeric(lo), upper = as.numeric(hi),
        length = as.numeric(hi - lo), covered = as.logical(target >= lo & target <= hi),
        below_target = as.logical(hi < target), above_target = as.logical(lo > target))
    }
    for (method in c("model", "cr0", "md")) {
      stopifnot(length(z$se[[method]]) == 10L,
                identical(names(z$se[[method]])[1:p], pn),
                all(is.finite(z$se[[method]]) & z$se[[method]] > 0))
      se <- z$se[[method]][pn]
      critical <- if (method == "model") qnorm(0.975) else qt(0.975, 43)
      all_rows[[length(all_rows) + 1L]] <- add_rows(method, z$beta - critical * se, z$beta + critical * se, se)
    }
    if (m <= 200) {
      stopifnot(identical(dim(z$draws), c(999L, 6L)), identical(colnames(z$draws), pn),
        all(is.finite(z$draws)), nrow(z$boot_diag) == 999L,
        identical(z$boot_diag$b, 1:999), all(z$boot_diag$converged),
        all(z$boot_diag$attempt == 1L), z$n_unresolved == 0L,
        all(is.finite(z$boot_diag$loglik)), all(is.finite(z$boot_diag$sigma_eps2)),
        all(z$boot_diag$sigma_eps2 > 0))
      limits <- apply(z$draws, 2, quantile, probs = c(0.025, 0.975), type = 7, na.rm = FALSE)
      boot_se <- apply(z$draws, 2, sd)
      all_rows[[length(all_rows) + 1L]] <- add_rows("bootstrap", limits[1, ], limits[2, ], boot_se)
      bootstrap_diagnostics[[length(bootstrap_diagnostics) + 1L]] <- data.frame(
        tau, m, B = 999L, nonconverged = sum(!z$boot_diag$converged),
        cap_active = sum(z$boot_diag$status == "constrained"),
        retries = sum(z$boot_diag$attempt > 1), unresolved = z$n_unresolved,
        min_observation_variance = min(z$boot_diag$sigma_eps2),
        max_observation_variance = max(z$boot_diag$sigma_eps2))
    } else stopifnot(is.null(z$draws))
  }
}
replicates <- do.call(rbind, all_rows)
row.names(replicates) <- NULL
matched <- replicates[replicates$in_matched200, ]
analytic <- replicates[replicates$method != "bootstrap", ]
stopifnot(nrow(replicates) == 20400L, nrow(matched) == 9600L, nrow(analytic) == 18000L)
write.csv(replicates, file.path(out, "replicate_intervals.csv"), row.names = FALSE)
write.csv(matched, file.path(out, "matched200_replicate_intervals.csv"), row.names = FALSE)
write.csv(do.call(rbind, fit_diagnostics), file.path(out, "fit_diagnostics.csv"), row.names = FALSE)
boot_diag <- do.call(rbind, bootstrap_diagnostics)
write.csv(boot_diag, file.path(out, "bootstrap_diagnostics.csv"), row.names = FALSE)

summarize <- function(records, n_required, include_bootstrap) {
  methods <- c("model", "cr0", "md", if (include_bootstrap) "bootstrap")
  rows <- list()
  for (tau in c(0, 0.4)) for (parameter in pn) for (method in methods) {
    z <- records[records$tau == tau & records$parameter == parameter & records$method == method, ]
    z <- z[order(z$m), ]
    stopifnot(nrow(z) == n_required, identical(z$m, seq_len(n_required)),
              length(unique(z$target)) == 1L, !anyNA(z$covered))
    n <- nrow(z)
    empirical_sd <- sd(z$estimate)
    coverage <- mean(z$covered)
    bias <- mean(z$estimate) - z$target[1]
    row <- data.frame(tau, parameter, method, n, coverage,
      coverage_mcse = sqrt(coverage * (1 - coverage) / n),
      mean_length = mean(z$length), length_mcse = sd(z$length) / sqrt(n),
      empirical_sd = empirical_sd, mean_se = mean(z$se),
      se_empirical_sd_ratio = mean(z$se) / empirical_sd,
      target = z$target[1], mean_estimate = mean(z$estimate), bias = bias,
      bias_mcse = empirical_sd / sqrt(n), median_length = median(z$length),
      p90_length = as.numeric(quantile(z$length, 0.90, type = 7)),
      below_target = mean(z$below_target), above_target = mean(z$above_target))
    if (include_bootstrap) {
      boot <- records[records$tau == tau & records$parameter == parameter & records$method == "bootstrap", ]
      boot <- boot[order(boot$m), ]
      stopifnot(identical(z$m, boot$m))
      difference <- as.numeric(boot$covered) - as.numeric(z$covered)
      row$empirical_sd_same200 <- empirical_sd
      row$bootstrap_minus_coverage <- mean(difference)
      row$paired_mcse <- sd(difference) / sqrt(n)
    }
    rows[[length(rows) + 1L]] <- row
  }
  result <- do.call(rbind, rows)
  row.names(result) <- NULL
  result
}
matched_summary <- summarize(matched, 200L, TRUE)
analytic_summary <- summarize(analytic, 500L, FALSE)
stopifnot(nrow(matched_summary) == 48L, nrow(analytic_summary) == 36L)
# Put the original figure input columns first for drop-in schema compatibility.
figure_columns <- c("tau", "parameter", "method", "n", "coverage", "coverage_mcse", "mean_length", "length_mcse",
  "empirical_sd_same200", "mean_se", "se_empirical_sd_ratio", "bootstrap_minus_coverage", "paired_mcse")
matched_summary <- matched_summary[, c(figure_columns, setdiff(names(matched_summary), figure_columns))]
write.csv(matched_summary, file.path(out, "paired_posvar200_free_variance_comparison.csv"), row.names = FALSE)
write.csv(analytic_summary, file.path(out, "free_variance_analytic500_summary.csv"), row.names = FALSE)

paired_rows <- list()
for (tau in c(0, 0.4)) for (parameter in pn) {
  z <- matched_summary[matched_summary$tau == tau & matched_summary$parameter == parameter, ]
  get <- function(method, column) z[z$method == method, column]
  paired_rows[[length(paired_rows) + 1L]] <- data.frame(tau, parameter, n = 200L,
    target = get("bootstrap", "target"), coverage_bootstrap = get("bootstrap", "coverage"),
    coverage_cr0 = get("cr0", "coverage"), coverage_md = get("md", "coverage"),
    bootstrap_minus_cr0 = get("cr0", "bootstrap_minus_coverage"),
    bootstrap_minus_md = get("md", "bootstrap_minus_coverage"),
    paired_mcse_bootstrap_cr0 = get("cr0", "paired_mcse"),
    paired_mcse_bootstrap_md = get("md", "paired_mcse"),
    mean_length_ratio_bootstrap_cr0 = get("bootstrap", "mean_length") / get("cr0", "mean_length"),
    mean_length_ratio_bootstrap_md = get("bootstrap", "mean_length") / get("md", "mean_length"),
    bootstrap_se_empirical_sd_ratio = get("bootstrap", "se_empirical_sd_ratio"))
}
paired <- do.call(rbind, paired_rows)
write.csv(paired, file.path(out, "bootstrap_paired_comparisons.csv"), row.names = FALSE)

old_fixed <- read.csv(fixed_summary_path)
stopifnot(nrow(old_fixed) == 48L, all(old_fixed$n == 200L))
stopifnot(file.copy(fixed_summary_path, file.path(out, "fixed_variance_reference_summary.csv"), overwrite = TRUE))
comparison <- merge(matched_summary, old_fixed, by = c("tau", "parameter", "method"),
                    suffixes = c("_free", "_fixed"), all = TRUE, sort = FALSE)
stopifnot(nrow(comparison) == 48L, !anyNA(comparison$coverage_free), !anyNA(comparison$coverage_fixed))
comparison$coverage_difference_free_minus_fixed <- comparison$coverage_free - comparison$coverage_fixed
comparison$mean_length_ratio_free_fixed <- comparison$mean_length_free / comparison$mean_length_fixed
comparison$mean_se_ratio_free_fixed <- comparison$mean_se_free / comparison$mean_se_fixed
comparison$se_empirical_sd_ratio_difference <- comparison$se_empirical_sd_ratio_free - comparison$se_empirical_sd_ratio_fixed
write.csv(comparison, file.path(out, "free_fixed_matched200_comparison.csv"), row.names = FALSE)

# The existing Part R summary is a valid reference for 500-fit analytic columns
# and 200-dataset bootstrap columns separately, never as a matched four-method table.
old_free <- read.csv(old_summary_path)
reproduction <- list()
for (i in seq_len(nrow(analytic_summary))) {
  z <- analytic_summary[i, ]
  ref <- old_free[old_free$tau == z$tau & old_free$parameter == z$parameter, ]
  stopifnot(nrow(ref) == 1L)
  errors <- c(coverage = z$coverage - ref[[paste0("cov_", z$method)]],
    mean_se = z$mean_se - ref[[paste0("se_", z$method)]], bias = z$bias - ref$bias_vs_target,
    empirical_sd = z$empirical_sd - ref$emp_sd, target = z$target - ref$target_free)
  stopifnot(max(abs(errors)) < 1e-12)
  reproduction[[length(reproduction) + 1L]] <- data.frame(tau = z$tau, parameter = z$parameter,
    method = z$method, n = 500L, max_difference_from_regenerated = max(abs(errors)))
}
for (i in which(matched_summary$method == "bootstrap")) {
  z <- matched_summary[i, ]
  ref <- old_free[old_free$tau == z$tau & old_free$parameter == z$parameter, ]
  errors <- c(coverage = z$coverage - ref$cov_boot, mean_se = z$mean_se - ref$se_boot,
              median_length = z$median_length - ref$medlen_boot)
  stopifnot(max(abs(errors)) < 1e-12)
  reproduction[[length(reproduction) + 1L]] <- data.frame(tau = z$tau, parameter = z$parameter,
    method = "bootstrap", n = 200L, max_difference_from_regenerated = max(abs(errors)))
}
write.csv(do.call(rbind, reproduction), file.path(out, "regenerated_summary_reproduction.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out, "R_session.txt"))

report <- capture.output({
  cat("Matched free-observation-variance inference: regenerated free-variance fit aggregation.\n")
  cat("First 200 original datasets per scenario; all six coefficients; B = 999.\n")
  cat("Model-based: normal quantile; CR0/MD: t43; bootstrap percentiles: quantile type 7.\n")
  cat("Every original fit and interval retained; no invalid results, refitting or new resampling.\n\n")
  for (tau in c(0, 0.4)) {
    cat(sprintf("tau = %.1f\n", tau))
    for (method in c("model", "cr0", "md", "bootstrap")) {
      z <- matched_summary[matched_summary$tau == tau & matched_summary$method == method, ]
      cat(sprintf("  %-9s coverage %.3f--%.3f; mean SE/empirical SD %.6f--%.6f\n",
                  method, min(z$coverage), max(z$coverage), min(z$se_empirical_sd_ratio), max(z$se_empirical_sd_ratio)))
    }
    z <- paired[paired$tau == tau, ]
    cat(sprintf("  Bootstrap minus CR0 coverage %.3f--%.3f; paired MCSE %.6f--%.6f\n",
      min(z$bootstrap_minus_cr0), max(z$bootstrap_minus_cr0), min(z$paired_mcse_bootstrap_cr0), max(z$paired_mcse_bootstrap_cr0)))
    cat(sprintf("  Bootstrap minus MD coverage %.3f to %.3f; paired MCSE %.6f--%.6f\n",
      min(z$bootstrap_minus_md), max(z$bootstrap_minus_md), min(z$paired_mcse_bootstrap_md), max(z$paired_mcse_bootstrap_md)))
    cat(sprintf("  Mean bootstrap/CR0 length %.6f--%.6f; bootstrap/MD length %.6f--%.6f\n",
      min(z$mean_length_ratio_bootstrap_cr0), max(z$mean_length_ratio_bootstrap_cr0),
      min(z$mean_length_ratio_bootstrap_md), max(z$mean_length_ratio_bootstrap_md)))
    bd <- boot_diag[boot_diag$tau == tau, ]
    cat(sprintf("  Bootstrap refits %d; nonconverged %d; unresolved %d; retries %d; cap-active %d; minimum observation variance %.6f\n\n",
      sum(bd$B), sum(bd$nonconverged), sum(bd$unresolved), sum(bd$retries), sum(bd$cap_active), min(bd$min_observation_variance)))
  }
  cat(sprintf("Maximum absolute matched free-minus-fixed coverage change: %.3f.\n",
              max(abs(comparison$coverage_difference_free_minus_fixed))))
  cat("All 1000 original records contribute to the separately saved 500-fit analytic summary; it contains no bootstrap columns.\n")
})
writeLines(report, file.path(out, "FINDINGS.txt"))
cat(report, sep = "\n")
