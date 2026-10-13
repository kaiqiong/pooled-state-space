source("R/simulation_design.R")
# Portable October 13 replication. Run from the repository root; see docs/SIMULATIONS.txt.
#!/usr/bin/env Rscript
# Independent algebra and summary audit. Reads saved information/scores/results;
# uses Cholesky solves and influence-vector norms, not parent covariance helpers.
project <- getwd()
work <- file.path(project, "results/simulation/free_variance_ablation")
records <- readRDS(file.path(work, "results/covariance_results.rds"))
replicate_df <- read.csv(file.path(work, "results/replicate_results.csv"))
summary_df <- read.csv(file.path(work, "results/summary.csv"))
design <- pss_simulation_design("positive")
free_target <- readRDS(file.path(project, "results/simulation/tables/expected_likelihood_target_free_r.rds"))$free$beta
pn <- colnames(design$Xmat)
p <- length(pn)
stopifnot(length(records) == 1000L, nrow(replicate_df) == 36000L,
          nrow(summary_df) == 72L, p == 6L)
for (tau in c(0, 0.4)) {
  ids <- vapply(records[vapply(records, function(x) x$tau == tau, logical(1))], function(x) x$m, numeric(1))
  stopifnot(length(ids) == 500L, identical(sort(as.integer(ids)), 1:500))
}

checks <- list()
max_expected_se_relative_error <- 0
max_rep_se_relative_error <- 0
max_nuisance_model_inflation <- -Inf
for (fit in records) {
  expected_target <- if (fit$tau == 0) design$truth$beta else free_target
  stopifnot(max(abs(fit$target - expected_target)) < 1e-13)
  variants <- list(full = list(J = fit$H, scores = fit$Sg),
    no_correction = list(J = fit$D, scores = fit$Sg),
    nuisance_fixed = list(J = fit$H[1:p, 1:p], scores = fit$Sg[, 1:p]))
  row_idx <- which(replicate_df$tau == fit$tau & replicate_df$m == fit$m)
  dataset_rows <- replicate_df[row_idx, ]
  stopifnot(length(row_idx) == 36L)
  for (name in names(variants)) {
    v <- variants[[name]]
    R <- chol(v$J)
    inverse_by_cholesky <- chol2inv(R)
    # J = R'R, so each column below is J^{-1} times a participant score.
    participant_influence <- backsolve(R, forwardsolve(t(R), t(v$scores)))
    model_se <- sqrt(diag(inverse_by_cholesky))
    cr0_se <- sqrt(rowSums(participant_influence^2))
    expected <- list(model = model_se[1:p], cr0 = cr0_se[1:p])
    recorded <- fit$variants[[name]]
    stopifnot(recorded$pd, recorded$finite,
      length(recorded$model_all) == ncol(v$J), length(recorded$cr0_all) == ncol(v$J))
    rel <- max(abs(model_se / recorded$model_all - 1), abs(cr0_se / recorded$cr0_all - 1))
    max_expected_se_relative_error <- max(max_expected_se_relative_error, rel)
    for (method in names(expected)) {
      z <- dataset_rows[dataset_rows$variant == name & dataset_rows$method == method, ]
      stopifnot(nrow(z) == p, identical(z$parameter, pn), all(z$valid))
      se_error <- max(abs(z$se / expected[[method]] - 1))
      max_rep_se_relative_error <- max(max_rep_se_relative_error, se_error)
      stopifnot(max(abs(z$estimate - fit$beta)) < 1e-13,
                max(abs(z$target - expected_target)) < 1e-13)
      critical_value <- if (method == "model") qnorm(0.975) else qt(0.975, 43)
      expected_covered <- abs(fit$beta - expected_target) <= critical_value * expected[[method]]
      stopifnot(identical(unname(expected_covered), z$covered))
      checks[[length(checks) + 1L]] <- data.frame(tau = fit$tau, m = fit$m,
        variant = name, method, cholesky_relative_error = rel, replicate_se_relative_error = se_error)
    }
  }
  # With positive definite H, ignoring nuisance estimation cannot increase the
  # model-based variance; there is no corresponding CR0 monotonicity assertion.
  max_nuisance_model_inflation <- max(max_nuisance_model_inflation,
    fit$variants$nuisance_fixed$model / fit$variants$full$model - 1)
}
stopifnot(max_expected_se_relative_error < 1e-8, max_rep_se_relative_error < 1e-8,
          max_nuisance_model_inflation < 1e-10)

numeric_checks <- list()
for (i in seq_len(nrow(summary_df))) {
  s <- summary_df[i, ]
  index <- which(replicate_df$tau == s$tau & replicate_df$variant == s$variant &
                 replicate_df$method == s$method & replicate_df$parameter == s$parameter)
  z <- replicate_df[index, ]
  stopifnot(nrow(z) == 500L, identical(sort(z$m), 1:500), all(z$valid))
  n <- nrow(z)
  empirical_sd <- sqrt(sum((z$estimate - mean(z$estimate))^2) / (n - 1L))
  critical_value <- if (s$method == "model") qnorm(0.975) else qt(0.975, 43)
  cover <- abs(z$estimate - z$target) <= critical_value * z$se
  cover_full <- abs(z$estimate - z$target) <= critical_value * z$full_se
  difference <- as.numeric(cover) - as.numeric(cover_full)
  coverage <- sum(cover) / n
  ratio <- z$se / z$full_se
  values <- c(n_total = n, n_valid = n, n_invalid = 0, target = z$target[1],
    mean_estimate = mean(z$estimate), empirical_sd = empirical_sd, mean_se = mean(z$se),
    ratio_to_emp_sd = mean(z$se) / empirical_sd, ratio_to_full = mean(z$se) / mean(z$full_se),
    coverage = coverage, valid_interval_coverage = coverage, valid_and_cover_fraction = coverage,
    coverage_mcse = sqrt(coverage * (1 - coverage) / n),
    coverage_difference_from_full = sum(difference) / n,
    paired_difference_mcse = sqrt(sum((difference - mean(difference))^2) / (n * (n - 1L))),
    lost_coverage = sum(cover_full & !cover), gained_coverage = sum(!cover_full & cover),
    individual_se_ratio_min = min(ratio), individual_se_ratio_median = median(ratio),
    individual_se_ratio_max = max(ratio))
  differences <- vapply(names(values), function(name) as.numeric(s[[name]]) - values[[name]], numeric(1))
  stopifnot(all(is.finite(differences)), max(abs(differences)) < 1e-10)
  numeric_checks[[length(numeric_checks) + 1L]] <- data.frame(tau = s$tau, variant = s$variant,
    method = s$method, parameter = s$parameter, max_summary_absolute_difference = max(abs(differences)))
}
write.csv(do.call(rbind, checks), file.path(work, "validation/independent_covariance_checks.csv"), row.names = FALSE)
write.csv(do.call(rbind, numeric_checks), file.path(work, "validation/independent_summary_checks.csv"), row.names = FALSE)
findings <- c(
  "PASS: all 1000 saved datasets and all 36000 coefficient/procedure rows independently audited.",
  "No optimizers, EM updates, simulated replacement datasets, or bootstrap draws were used.",
  "Model covariance was independently obtained by Cholesky inversion; CR0 SEs by norms of participant influence vectors.",
  sprintf("Maximum relative SE difference against helper results: %.12g.", max_expected_se_relative_error),
  sprintf("Maximum relative SE difference against coefficient-level CSV: %.12g.", max_rep_se_relative_error),
  sprintf("Maximum nuisance-fixed/full model-SE ratio minus one: %.12g.", max_nuisance_model_inflation),
  "All 72 summary rows reproduce from coefficient-level records, including empirical SD, both SE ratios, coverage and paired MCSE.",
  "All coefficient estimates and targets match saved covariance records and the generating/free working targets, respectively.",
  "All covariance evaluations are valid; every primary coverage denominator is 500.")
writeLines(findings, file.path(work, "validation/INDEPENDENT_AUDIT.txt"))
cat(findings, sep = "\n")
