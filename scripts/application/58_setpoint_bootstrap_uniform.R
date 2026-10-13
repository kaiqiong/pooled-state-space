## Step 58: participant-bootstrap percentile intervals (B = 2000, same procedure as the primary analysis) for EVERY
## set-point variant of Steps 52-53, at both window lengths, so that the reference-sensitivity figure and Table S22
## use one interval procedure. Fixed references: EM at sigma_eps2 = 0 (3 iterations). Estimated references: direct
## ML of the closed-form boundary likelihood (Step 53 code), with a1, P1, sigma_eta2 profiled; 2000 refits each.
source(file.path(APP_ROOT, "R", "state_space_em.R")); suppressPackageStartupMessages(library(parallel))
nll <- function(par, Y, X, Z) { p <- ncol(X); phi <- as.vector(tanh(X %*% par[1:p])); d <- Y - as.vector(Z %*% par[-(1:p)]); Tt <- ncol(Y); N <- nrow(Y)
  e <- d[, -1] - phi * d[, -Tt]; s2 <- mean(e^2); P1 <- mean((d[, 1] - mean(d[, 1]))^2); v <- N / 2 * (log(2 * pi * P1) + 1) + N * (Tt - 1) / 2 * (log(2 * pi * s2) + 1); if (is.finite(v)) v else 1e300 }
ml <- function(Y, X, Z, start) { o <- optim(start, nll, Y = Y, X = X, Z = Z, method = "BFGS", control = list(maxit = 2000, reltol = 1e-12)); o <- optim(o$par, nll, Y = Y, X = X, Z = Z, method = "BFGS", control = list(maxit = 2000, reltol = 1e-13)); list(par = o$par, conv = o$convergence == 0) }
B <- pss_boot_count(2000); out <- list(); diag <- list()
for (dm in c("65", "120")) { r52 <- readRDS(sprintf("setpoint_boundary_%s.rds", dm)); r53 <- readRDS("setpoint_boundary_directML.rds"); X <- r52$X; subj <- r52$subj; G <- length(unique(subj)); Yraw <- readRDS(sprintf("meal_windows_final_%smin.rds", dm))$Y; meta <- r52$meta
  variants <- list(`pre-meal` = list(Y = Yraw - meta$baseline_gl, idx = seq_len(nrow(Yraw))), `180-min` = list(Y = Yraw - meta$gl_180, idx = which(!is.na(meta$gl_180))))
  if (dm == "65") variants[["next-meal"]] <- list(Y = Yraw - meta$gl_next, idx = which(!is.na(meta$gl_next)))
  for (v in names(variants)) { Y <- variants[[v]]$Y[variants[[v]]$idx, ]; Xv <- X[variants[[v]]$idx, ]; cl <- subj[variants[[v]]$idx]; us <- unique(cl)   # covariates subset in step with the windows (a first run indexed the full X: wrong)
    fit <- r52$fits[[c(`pre-meal` = "P", `180-min` = "B", `next-meal` = "C")[v]]]
    set.seed(58000 + as.integer(dm) + nchar(v)); t0 <- Sys.time()
    bb <- do.call(rbind, mclapply(seq_len(B), function(b) { set.seed(58000000 + 10000 * as.integer(dm) + 100 * nchar(v) + b); idx <- unlist(lapply(sample(us, length(us), TRUE), function(s) which(cl == s)))
      f <- tryCatch(em_regression_mv(Y[idx, ], Xv[idx, ], 0, list(beta = fit$beta, a1 = fit$a1, P1 = fit$P1, sigma_eta2 = fit$sigma_eta2, sigma_eps2 = 0), thres = 1e-9, max_iter = 3000, fix_sigma_eps2 = 0), error = function(e) NULL)
      if (is.null(f)) rep(NA, 8) else c(f$beta, f$converged, isTRUE(f$at_cap)) }, mc.cores = APP_CORES))
    ok <- complete.cases(bb); est <- fit$beta[6]
    out[[paste(dm, v)]] <- data.frame(len = as.integer(dm), variant = v, n = nrow(Y), est = est, boot_se = sd(bb[ok, 6]), boot_lo = quantile(bb[ok, 6], .025), boot_hi = quantile(bb[ok, 6], .975), B = B, boot_fail = sum(!ok), boot_nonconv = sum(bb[ok, 7] == 0), boot_at_cap = sum(bb[ok, 8] == 1), method = "EM at sigma_eps2 = 0")
    cat(sprintf("%s-min %-10s: est %.3f, bootstrap (%.3f, %.3f), fail %d, nonconv %d, cap %d, %.1f min\n", dm, v, est, out[[paste(dm, v)]]$boot_lo, out[[paste(dm, v)]]$boot_hi, sum(!ok), sum(bb[ok, 7] == 0), sum(bb[ok, 8] == 1), as.numeric(Sys.time() - t0, units = "mins"))) }
  Y0 <- Yraw - meta$baseline_gl; cl <- subj; us <- unique(cl)
  for (v in c("deviations", "full")) { best <- r53$detail[[paste(dm, v)]]$best; Z <- if (v == "deviations") X[, -1] else X; start <- c(best$beta, best$gamma); t0 <- Sys.time()
    bb <- do.call(rbind, mclapply(seq_len(B), function(b) { set.seed(58100000 + 10000 * as.integer(dm) + 100 * nchar(v) + b); idx <- unlist(lapply(sample(us, length(us), TRUE), function(s) which(cl == s)))
      f <- tryCatch(ml(Y0[idx, ], X[idx, ], Z[idx, , drop = FALSE], start), error = function(e) NULL); if (is.null(f)) rep(NA, 7) else c(f$par[1:6], f$conv) }, mc.cores = APP_CORES))
    ok <- complete.cases(bb); est <- best$beta[6]
    out[[paste(dm, v)]] <- data.frame(len = as.integer(dm), variant = paste("estimated set point,", v), n = nrow(Y0), est = est, boot_se = sd(bb[ok, 6]), boot_lo = quantile(bb[ok, 6], .025), boot_hi = quantile(bb[ok, 6], .975), B = B, boot_fail = sum(!ok), boot_nonconv = sum(bb[ok, 7] == 0), boot_at_cap = NA, method = "direct ML, boundary likelihood")
    cat(sprintf("%s-min estimated %-10s: est %.3f, bootstrap (%.3f, %.3f), fail %d, nonconv %d, %.1f min\n", dm, v, est, out[[paste(dm, v)]]$boot_lo, out[[paste(dm, v)]]$boot_hi, sum(!ok), sum(bb[ok, 7] == 0), as.numeric(Sys.time() - t0, units = "mins"))) } }
res <- do.call(rbind, out); rownames(res) <- NULL; print(res, digits = 3); write.csv(res, "setpoint_bootstrap_uniform.csv", row.names = FALSE); cat("saved\n")
