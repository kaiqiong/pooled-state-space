#!/usr/bin/env Rscript
# Execute from any directory; all calculations write below results/simulation.
args <- commandArgs(trailingOnly = TRUE)
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
root <- normalizePath(file.path(dirname(script), "../.."), mustWork = TRUE)
setwd(root)
stage <- if (length(args)) args[1] else "help"
cores <- if (length(args) >= 2) as.integer(args[2]) else as.integer(Sys.getenv("SIM_CORES", "1"))
if (is.na(cores) || cores < 1) stop("cores must be a positive integer")
if (.Platform$OS.type == "windows" && cores > 1L) stop("Use cores=1 on Windows; multicore jobs use parallel::mclapply")
Sys.setenv(SIM_CORES = cores)
run <- function(name, ...) {
  extra <- as.character(c(...))
  cat("\nRunning", name, paste(extra, collapse = " "), "\n")
  status <- system2(file.path(R.home("bin"), "Rscript"),
                    c("--vanilla", shQuote(file.path("scripts/simulation", name)), shQuote(extra)))
  if (status != 0L) stop(name, " failed with status ", status)
}
app_required <- function() {
  required <- "results/application/application_design.rds"
  if (!file.exists(required)) stop("Run Rscript scripts/application/run.R preprocess first: missing ", required)
}
jobs <- list(
  smoke = function() { run("check_information.R"); run("check_generator_streams.R") },
  independent = function() {
    for (s in c("estimator_benchmark.R", "base_coverage.R", "interaction.R", "episode_size.R",
                "persistence_sweep.R", "heavy_tails.R", "simple_models.R")) run(s)
    run("persistence_sweep.R", "unguarded")
  },
  targets = function() {
    app_required()
    for (s in c("working_targets.R", "free_variance_target.R", "controlled_targets.R")) run(s)
  },
  clustered = function() {
    app_required()
    run("clustered_boundary.R"); run("clustered_fixed_variance.R")
    run("cluster_bootstrap.R", "boundary", 1000, cores)
    run("cluster_bootstrap.R", "posvar", 200, cores)
    run("clustered_free_variance.R", 500, 200, cores)
  },
  summaries = function() {
    for (s in c("summarize_fixed.R", "summarize_free_variance.R", "matched_comparison.R", "free_variance_ablation.R")) run(s)
  },
  robustness = function() {
    app_required()
    for (s in c("mixed_comparator.R", "differential_attenuation.R", "direct_optimization.R", "score_test_size_power.R")) run(s)
  },
  numerical = function() {
    app_required()
    if (!file.exists("results/application/FINAL_fitted_models_final.rds")) stop("Run the application prepare stage before actual-application numerical checks")
    for (s in c("check_information.R", "check_unit_information.R", "check_nuisance_inference.R",
                "check_mixed_model.R", "check_free_variance_derivatives.R", "check_ablation_hessian.R",
                "check_ablation_summary.R")) run(s)
  },
  diagnostics = function() {
    app_required()
    for (design in c("boundary", "posvar")) {
      run("check_original_starts.R", design, cores)
      # Manuscript: all 999 heterogeneous boundary draws; 5 interior draws
      # plus all cap-active heterogeneous fixed-positive-variance draws.
      run("check_bootstrap_starts.R", design, 0.4, cores, if (design == "boundary") 999 else 5)
    }
    run("check_coverage_effect.R", "results/simulation/tables/bootstrap_local_optima_boundary_tau0.4.rds", "boundary")
    run("check_coverage_effect.R", "results/simulation/tables/bootstrap_local_optima_posvar_tau0.4.rds", "posvar")
    run("check_md_conditioning.R")
  },
  caps = function() {
    app_required(); run("cap_sensitivity.R", cores); run("cap_multistart.R", cores)
    run("summarize_caps.R", "results/simulation/cap_sensitivity_ms", "_ms")
    run("check_caps.R")
  }
)
if (stage == "help") {
  cat("Usage: Rscript scripts/simulation/run.R STAGE [CORES]\n",
      "Stages: ", paste(c(names(jobs), "all"), collapse = ", "), "\n",
      "Full replication: all (substantial computation; no shipped fit caches).\n",
      "smoke checks numerical information and RNG reconstruction; it is not a Monte Carlo replication.\n",
      "See docs/SIMULATIONS.txt for inputs, output tables and precise stage order.\n", sep = "")
} else {
  dir.create("results/simulation/tables", recursive = TRUE, showWarnings = FALSE)
  dir.create("results/simulation/figures", recursive = TRUE, showWarnings = FALSE)
  dir.create("results/simulation/validation", recursive = TRUE, showWarnings = FALSE)
  if (stage == "all") {
    app_required()
    if (!file.exists("results/application/FINAL_fitted_models_final.rds")) stop("Run the application prepare stage before all, which includes actual-application numerical checks")
    for (key in c("independent", "targets", "clustered", "summaries", "robustness", "numerical", "diagnostics", "caps")) jobs[[key]]()
  } else if (stage %in% names(jobs)) jobs[[stage]]()
  else stop("Unknown stage: ", stage)
  capture.output(sessionInfo(), file = file.path("results/simulation/validation", paste0("session_", stage, ".txt")))
}
