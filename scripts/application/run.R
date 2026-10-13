#!/usr/bin/env Rscript
## Portable entry point. Run from any directory; all paths resolve at runtime.
args <- commandArgs(trailingOnly = TRUE)
stage <- if (length(args)) args[1] else "help"
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
APP_ROOT <- normalizePath(file.path(dirname(script), "..", ".."), mustWork = TRUE)
option <- function(name, fallback) {
  at <- match(name, args)
  if (is.na(at)) return(fallback)
  if (at == length(args)) stop("Missing value for ", name)
  args[at + 1L]
}
stages <- c("preprocess", "prepare", "primary", "profile", "reference", "diagnostics", "tables", "all")
if (!stage %in% stages) {
  cat("Usage: Rscript scripts/application/run.R STAGE [--data DIR] [--out DIR]\n",
      "Stages: preprocess, prepare, primary, profile, reference, diagnostics, tables, all\n",
      "Defaults: data/raw/CGMacros; results/application (relative to repository).\n",
      "Optional PSS_B >= 2 overrides bootstrap counts for development only.\n",
      "Optional PSS_CORES controls reference-bootstrap workers (default 8; Windows 1).\n")
  if (stage == "help") quit(status = 0)
  stop("Unknown stage: ", stage)
}
APP_DATA <- normalizePath(option("--data", Sys.getenv("PSS_RAW_DATA", file.path(APP_ROOT, "data", "raw", "CGMacros"))), mustWork = FALSE)
if (stage %in% c("preprocess", "prepare", "reference", "diagnostics", "all") &&
    !file.exists(file.path(APP_DATA, "bio.csv"))) stop("Expected bio.csv in data directory: ", APP_DATA)
out <- option("--out", Sys.getenv("PSS_APPLICATION_OUT", file.path(APP_ROOT, "results", "application")))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
APP_OUT <- normalizePath(out, mustWork = TRUE)
setwd(APP_OUT)
dir.create("logs", showWarnings = FALSE)
APP_CORES <- as.integer(Sys.getenv("PSS_CORES", if (.Platform$OS.type == "windows") "1" else "8"))
if (!is.finite(APP_CORES) || APP_CORES < 1L) stop("PSS_CORES must be a positive integer")
if (.Platform$OS.type == "windows") APP_CORES <- 1L
pss_boot_count <- function(default) {
  value <- Sys.getenv("PSS_B", "")
  if (!nzchar(value)) return(default)
  count <- suppressWarnings(as.integer(value))
  if (!is.finite(count) || count < 2L) stop("PSS_B must be an integer >= 2")
  count
}
if (nzchar(Sys.getenv("PSS_B"))) {
  pss_boot_count(2000)
  warning("PSS_B overrides published bootstrap counts. These are development results.")
}
run_step <- function(name, decay = NULL) {
  env <- new.env(parent = globalenv())
  env$APP_DECAY_MIN <- decay
  label <- paste0(sub("[.]R$", "", name), if (!is.null(decay)) paste0("_", decay))
  cat("Running ", label, "; log: ", file.path(APP_OUT, "logs", paste0(label, ".log")), "\n", sep = "")
  sink(file.path("logs", paste0(label, ".log")), split = TRUE)
  on.exit(sink(), add = TRUE)
  source(file.path(APP_ROOT, "scripts", "application", name), local = env, chdir = FALSE)
}
run_stage <- function(which) {
  switch(which,
    preprocess = {
      run_step("40_windows_final.R", 65)
      run_step("40_windows_final.R", 120)
      run_step("00_make_design.R")
      run_step("01_native_qc.R")
      run_step("42_verify_final.R")
    },
    prepare = run_step("41_refit_final.R"),
    primary = run_step("51_primary_boundary.R"),
    profile = {
      run_step("44_profile_final.R")
      run_step("48_profile_intervals_final.R")
    },
    reference = {
      run_step("52_setpoint_boundary.R")
      run_step("53_setpoint_boundary_directML.R")
      run_step("58_setpoint_bootstrap_uniform.R")
    },
    diagnostics = {
      run_step("49_descriptives_final.R")
      run_step("50_schematic_example_final.R")
      run_step("54_fixed_time_windows.R")
      run_step("55_boundary_score_test.R")
      run_step("56_mixed_model_application.R")
    },
    tables = run_step("final_application_tables.R"))
}
writeLines(c(paste("stage:", stage), paste("PSS_B:", Sys.getenv("PSS_B", "published defaults")),
             paste("PSS_CORES:", APP_CORES), capture.output(sessionInfo())),
           file.path("logs", paste0("session_", stage, ".txt")))
if (stage == "all") for (s in stages[stages != "all"]) run_stage(s) else run_stage(stage)
cat("Completed ", stage, ". Outputs: ", APP_OUT, "\n", sep = "")
