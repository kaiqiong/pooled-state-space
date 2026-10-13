source("R/simulation_design.R")
# Portable October 13 replication. Run from the repository root; see docs/SIMULATIONS.txt.
root <- getwd(); out <- file.path(root, "results/simulation/free_variance_ablation")
source("R/state_space_em.R")
source("R/free_variance_cluster_inference.R")
source("scripts/simulation/helpers/ablation_covariance.R")
dir.create(file.path(out, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "validation"), recursive = TRUE, showWarnings = FALSE)
f <- pss_simulation_design("positive")
Xmat <- f$Xmat; subj <- f$meta$subject; N <- nrow(Xmat); Tt <- ncol(f$Y)
pn <- colnames(Xmat); p <- ncol(Xmat); G <- length(unique(subj))
stopifnot(N == 1174L, Tt == 14L, p == 6L, G == 44L)
beta0 <- setNames(f$truth$beta, pn); a1_0 <- f$truth$a1; P1_0 <- f$truth$P1
se2_0 <- f$truth$sigma_eta2; r0 <- 400
tg <- readRDS(file.path(root, "results/simulation/tables/expected_likelihood_target_free_r.rds"))$free
stopifnot(identical(names(tg$beta), pn))
# Extract only the generator definition, never execute the simulation runner.
exprs <- as.list(parse(file.path(root, "scripts/simulation/clustered_free_variance.R")))
gen_expr <- Filter(function(e) is.call(e) && identical(e[[1]], as.name("<-")) &&
                   identical(e[[2]], as.name("gen")), exprs)
stopifnot(length(gen_expr) == 1L); eval(gen_expr[[1]])

all_results <- list(); replicate_rows <- list(); diag_rows <- list()
all_start <- proc.time()[["elapsed"]]
for (tau in c(0, 0.4)) {
  cache <- file.path(root, "results/simulation/freevar", paste0("tau", tau))
  fs <- list.files(cache, "^m[0-9]+\\.rds$", full.names = TRUE)
  ids <- as.integer(sub("^m([0-9]+)\\.rds$", "\\1", basename(fs)))
  stopifnot(length(fs) == 500L, setequal(ids, seq_len(500L)))
  fs <- fs[order(ids)]
  target <- if (tau == 0) beta0 else tg$beta
  for (m in seq_len(500L)) {
    t0 <- proc.time()[["elapsed"]]
    fit <- readRDS(fs[m])
    stopifnot(fit$m == m, fit$tau == tau, isTRUE(fit$conv), !fit$err,
              identical(names(fit$beta), pn), !is.null(fit$se), fit$sigma_eps2 > 0)
    seed <- 410000 + 10000 * round(10 * tau) + m
    set.seed(seed); y <- gen(tau)
    fit$kf <- kalman_run(y, tanh(as.vector(Xmat %*% fit$beta)), 0,
                        fit$a1, fit$P1, fit$sigma_eta2, fit$sigma_eps2)
    Fv <- fit$kf$P_pred + fit$sigma_eps2; e <- y - fit$kf$alpha_pred
    loglik <- -0.5 * sum(log(2 * pi * Fv) + e^2 / Fv)
    ll_diff <- abs(loglik - fit$loglik)
    stopifnot(ll_diff < 1e-8)
    Ii <- info_per_unit_free(y, Xmat, 0, fit)
    Hraw <- apply(Ii, c(2, 3), sum)
    asym <- max(abs(Hraw - t(Hraw))) / max(abs(Hraw))
    H <- (Hraw + t(Hraw)) / 2
    D <- direct_information_free(y, Xmat, 0, fit)
    S <- scores_per_unit_free(y, Xmat, 0, fit)
    Sg <- rowsum(S, subj, reorder = FALSE)
    stopifnot(nrow(Sg) == G, identical(colnames(H), colnames(Sg)),
              identical(colnames(H), colnames(D)))
    vs <- evaluate_variants(H, D, Sg, p)
    stopifnot(vs$full$pd, vs$full$finite)
    saved_diff <- max(abs(vs$full$model_all - fit$se$model),
                      abs(vs$full$cr0_all - fit$se$cr0))
    saved_rel <- max(abs(vs$full$model_all / fit$se$model - 1),
                    abs(vs$full$cr0_all / fit$se$cr0 - 1))
    stopifnot(saved_diff < 1e-9, saved_rel < 1e-9)
    for (variant in names(vs)) for (method in c("model", "cr0")) {
      z <- vs[[variant]]; se <- z[[method]]; crit <- if (method == "model") qnorm(.975) else qt(.975, G - 1)
      valid <- is.finite(se) & se > 0
      covered <- ifelse(valid, abs(fit$beta - target) <= crit * se, NA)
      full_covered <- abs(fit$beta - target) <= crit * vs$full[[method]]
      if (variant == "full") {
        original_covered <- abs(fit$beta - target) <= crit * fit$se[[method]][pn]
        stopifnot(identical(as.logical(covered), as.logical(original_covered)))
      }
      replicate_rows[[length(replicate_rows) + 1L]] <- data.frame(
        tau = tau, m = m, parameter = pn, variant = variant, method = method,
        estimate = as.numeric(fit$beta), target = as.numeric(target),
        se = as.numeric(se), full_se = as.numeric(vs$full[[method]]),
        valid = valid, covered = as.logical(covered), full_covered = as.logical(full_covered))
      diag_rows[[length(diag_rows) + 1L]] <- data.frame(tau = tau, m = m,
        variant = variant, method = method, positive_definite = z$pd, finite_se = z$finite,
        min_eigenvalue = z$min_eigenvalue, condition_number = z$condition_number,
        loglik_abs_difference = ll_diff, full_se_max_abs_difference = saved_diff,
        full_se_max_relative_difference = saved_rel, information_asymmetry = asym,
        original_converged = fit$conv, original_status = fit$status,
        elapsed_seconds = proc.time()[["elapsed"]] - t0)
    }
    all_results[[length(all_results) + 1L]] <- list(tau = tau, m = m, seed = seed,
      beta = fit$beta, target = target, H = H, D = D, Sg = Sg, variants = vs,
      sigma_eps2 = fit$sigma_eps2, status = fit$status,
      loglik_abs_difference = ll_diff, full_se_max_abs_difference = saved_diff)
    if (m %% 100L == 0L) cat(sprintf("tau=%.1f: %d/500; elapsed %.1fs\n",
                    tau, m, proc.time()[["elapsed"]] - all_start))
  }
}
replicates <- do.call(rbind, replicate_rows); diagnostics <- do.call(rbind, diag_rows)
write.csv(replicates, file.path(out, "results/replicate_results.csv"), row.names = FALSE)
write.csv(diagnostics, file.path(out, "validation/dataset_diagnostics.csv"), row.names = FALSE)
saveRDS(all_results, file.path(out, "results/covariance_results.rds"))

keys <- unique(replicates[, c("tau", "variant", "method", "parameter")]); summaries <- list()
for (i in seq_len(nrow(keys))) {
  key <- keys[i, ]; z <- subset(replicates, tau == key$tau & variant == key$variant &
                               method == key$method & parameter == key$parameter)
  stopifnot(nrow(z) == 500L, length(unique(z$m)) == 500L)
  good <- z$valid; nv <- sum(good); nc <- sum(z$covered[good]); emp <- sd(z$estimate)
  # The primary summaries require all 500 intervals to be valid. Failures are
  # counted and retained, never silently removed from the coverage denominator.
  complete <- nv == nrow(z)
  delta <- if (complete) as.numeric(z$covered) - as.numeric(z$full_covered) else NA_real_
  ratio <- z$se[good] / z$full_se[good]
  summaries[[i]] <- cbind(key, data.frame(n_total = nrow(z), n_valid = nv, n_invalid = nrow(z) - nv,
    target = z$target[1], mean_estimate = mean(z$estimate), empirical_sd = emp,
    mean_se = if (complete) mean(z$se) else NA_real_,
    ratio_to_emp_sd = if (complete) mean(z$se) / emp else NA_real_,
    ratio_to_full = if (complete) mean(z$se) / mean(z$full_se) else NA_real_,
    coverage = if (complete) nc / nrow(z) else NA_real_,
    valid_interval_coverage = if (nv > 0) nc / nv else NA_real_,
    valid_and_cover_fraction = nc / nrow(z),
    coverage_mcse = if (complete) sqrt((nc / nv) * (1 - nc / nv) / nv) else NA_real_,
    coverage_difference_from_full = if (complete) mean(delta) else NA_real_,
    paired_difference_mcse = if (complete) sd(delta) / sqrt(nrow(z)) else NA_real_,
    lost_coverage = if (complete) sum(z$full_covered & !z$covered) else NA_integer_,
    gained_coverage = if (complete) sum(!z$full_covered & z$covered) else NA_integer_,
    individual_se_ratio_min = if (nv) min(ratio) else NA_real_,
    individual_se_ratio_median = if (nv) median(ratio) else NA_real_,
    individual_se_ratio_max = if (nv) max(ratio) else NA_real_))
}
summary <- do.call(rbind, summaries); rownames(summary) <- NULL
write.csv(summary, file.path(out, "results/summary.csv"), row.names = FALSE)

# Independent published-summary check in addition to per-fit saved-SE checks.
saved_summary <- read.csv(file.path(root, "results/simulation/tables/part_r_summary.csv"))
full_rows <- subset(summary, variant == "full")
summary_diffs <- list()
for (i in seq_len(nrow(full_rows))) {
  x <- full_rows[i, ]; old <- subset(saved_summary, tau == x$tau & parameter == x$parameter)
  stopifnot(nrow(old) == 1L)
  cov_diff <- x$coverage - old[[paste0("cov_", x$method)]]
  se_diff <- x$mean_se - old[[paste0("se_", x$method)]]
  stopifnot(abs(cov_diff) < 1e-12, abs(se_diff) < 1e-10,
            abs(x$target - old$target_free) < 1e-12)
  summary_diffs[[i]] <- data.frame(tau = x$tau, parameter = x$parameter,
                                  method = x$method, coverage_difference = cov_diff,
                                  mean_se_difference = se_diff)
}
write.csv(do.call(rbind, summary_diffs), file.path(out, "validation/full_summary_reproduction.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out, "validation/R_session.txt"))
elapsed <- proc.time()[["elapsed"]] - all_start
writeLines(c(sprintf("1000 saved fits; no EM or bootstrap; elapsed %.2f seconds.", elapsed),
  sprintf("All full SEs reproduced: max absolute difference %.12g.", max(diagnostics$full_se_max_abs_difference)),
  sprintf("Invalid covariance evaluations: %d of %d.", sum(!diagnostics$finite_se), nrow(diagnostics)),
  "Full-information coverage outcomes exactly match saved fits and the regenerated free-variance summary."),
  file.path(out, "validation/RUN_SUMMARY.txt"))
cat(readLines(file.path(out, "validation/RUN_SUMMARY.txt")), sep = "\n")
