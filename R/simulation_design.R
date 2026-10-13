## Application-shaped simulation data and frozen aggregate generating parameters.
## Refitted application estimates remain in their own results/application RDS files.
pss_simulation_design <- function(design = c("positive", "boundary"), root = getwd(),
                                  application_path = NULL, include_application_fit = FALSE) {
  design <- match.arg(design)
  root <- normalizePath(root, mustWork = TRUE)
  if (is.null(application_path))
    application_path <- file.path(root, "results", "application", "application_design.rds")
  if (!file.exists(application_path))
    stop("Missing regenerated application design: ", application_path,
         ". Run Rscript scripts/application/run.R preprocess first.")
  config <- new.env(parent = baseenv())
  sys.source(file.path(root, "config", "generating_parameters.R"), envir = config)
  truth <- config$pss_generating_parameters[[design]]
  data <- readRDS(application_path)
  required <- c("Y", "Xmat", "meta", "flags", "slot_time")
  if (!all(required %in% names(data))) stop("Application design lacks required fields")
  stopifnot(is.matrix(data$Y), is.matrix(data$Xmat),
            identical(dim(data$Y), c(1174L, 14L)),
            identical(dim(data$Xmat), c(1174L, 6L)),
            length(data$meta$subject) == 1174L,
            length(unique(data$meta$subject)) == 44L,
            identical(colnames(data$Xmat), names(truth$beta)),
            all(is.finite(data$Y)), all(is.finite(data$Xmat)))
  out <- data[required]
  out$truth <- truth
  out$design <- design
  out$sd_eps_primary <- sqrt(truth$sigma_eps2)
  out$application_fit <- NULL
  if (isTRUE(include_application_fit)) {
    name <- if (design == "positive") "FINAL_fitted_models_final.rds" else "PRIMARY_boundary.rds"
    path <- file.path(dirname(application_path), name)
    if (!file.exists(path)) stop("Requested refitted application estimate is missing: ", path)
    out$application_fit <- readRDS(path)$fit
  }
  out
}
