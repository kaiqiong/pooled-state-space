## Step 53: estimated set point at the boundary (sigma_eps2 = 0) by DIRECT maximum likelihood.
## Package adaptation: pilot bootstrap omitted; Step 58 supplies final intervals.
## Reason: at sigma_eps2 = 0 the smoothed state equals y - mu exactly, so the EM update of the set-point
## coefficients gamma (Step 52, rows "estimated ...") returns its starting value: those rows of Step 52 are
## NOT estimates and are superseded by this script. At the boundary the likelihood is closed form:
##   d_it = y_it - Z_i'gamma;  d_i1 ~ N(a1, P1);  d_it = phi_i d_{i,t-1} + eta_it, eta ~ N(0, sigma_eta2),
## with a1, P1, sigma_eta2 profiled out. Maximised by BFGS from several starts.
r52 <- list(`65` = readRDS("setpoint_boundary_65.rds"), `120` = readRDS("setpoint_boundary_120.rds"))
nll <- function(par, Y, X, Z) { p <- ncol(X); phi <- as.vector(tanh(X %*% par[1:p])); d <- Y - as.vector(Z %*% par[-(1:p)]); Tt <- ncol(Y); N <- nrow(Y)
  e <- d[, -1] - phi * d[, -Tt]; s2 <- mean(e^2); P1 <- mean((d[, 1] - mean(d[, 1]))^2); v <- N / 2 * (log(2 * pi * P1) + 1) + N * (Tt - 1) / 2 * (log(2 * pi * s2) + 1); if (is.finite(v)) v else 1e300 }
nll0 <- function(b, Y, X) nll(c(b, 0), Y, X, matrix(0, nrow(Y), 1))
ml <- function(Y, X, Z, start) { o <- optim(start, nll, Y = Y, X = X, Z = Z, method = "BFGS", control = list(maxit = 2000, reltol = 1e-12)); o2 <- optim(o$par, nll, Y = Y, X = X, Z = Z, method = "BFGS", control = list(maxit = 2000, reltol = 1e-13))
  list(beta = setNames(o2$par[1:ncol(X)], colnames(X)), gamma = setNames(o2$par[-(1:ncol(X))], colnames(Z)), loglik = -o2$value, conv = o2$convergence == 0, grad = max(abs(numDeriv::grad(nll, o2$par, Y = Y, X = X, Z = Z)))) }
out <- list(); gam <- list(); keep <- list()
for (dm in c("65", "120")) { r <- r52[[dm]]; X <- r$X; subj <- r$subj; G <- length(unique(subj)); Y0 <- readRDS(sprintf("meal_windows_final_%smin.rds", dm))$Y - r$meta$baseline_gl; b0 <- r$fits$P$beta; grp <- ifelse(X[, "Diabetic"] == 1, "D", ifelse(X[, "Prediabetic"] == 1, "P", "H"))
  ## check: with gamma = 0 the direct likelihood must reproduce the EM boundary fit
  chk <- optim(b0, nll0, Y = Y0, X = X, method = "BFGS", control = list(reltol = 1e-13)); cat(sprintf("\n=== %s-min | check vs EM boundary fit: loglik %.2f vs %.2f, max |beta diff| %.1e\n", dm, -chk$value, r$fits$P$loglik, max(abs(chk$par - b0))))
  for (v in c("deviations", "full")) { Z <- if (v == "deviations") X[, -1] else X; k <- ncol(Z)
    set.seed(53)  # seed the random fourth start before constructing it
    starts <- list(c(b0, rep(0, k)), c(b0, if (v == "full") c(-30, rep(0, 5)) else rep(-5, k)), c(atanh(0.8), rep(0, 5), if (v == "full") c(20, rep(0, 5)) else rep(5, k)), c(b0 * 0.8, rnorm(k, 0, 10)))
    fits <- lapply(starts, function(s) tryCatch(ml(Y0, X, Z, s), error = function(e) NULL)); fits <- fits[!sapply(fits, is.null)]; ll <- sapply(fits, `[[`, "loglik"); best <- fits[[which.max(ll)]]
    cat(sprintf("%s: loglik by start %s | best %.2f (gain over gamma = 0: %.1f on %d df) | max |grad| %.1e\n   gamma: %s\n   beta: %s\n", v, paste(round(ll, 2), collapse = ", "), best$loglik, best$loglik - r$fits$P$loglik, k, best$grad, paste(names(best$gamma), round(best$gamma, 1), collapse = "  "), paste(names(best$beta), round(best$beta, 3), collapse = "  ")))
    # Pilot B=500/300 intervals are superseded by Step 58 (B=2000).
    ph <- as.vector(tanh(X %*% best$beta)); mu <- as.vector(Z %*% best$gamma)
    out[[paste(dm, v)]] <- data.frame(len = as.integer(dm), variant = paste("estimated set point,", v), parameter = colnames(X), est = best$beta, se = NA_real_, lo = NA_real_, hi = NA_real_, n = nrow(Y0),
      phi_H = mean(ph[grp == "H"]), phi_P = mean(ph[grp == "P"]), phi_D = mean(ph[grp == "D"]), mu_H = mean(mu[grp == "H"]), mu_P = mean(mu[grp == "P"]), mu_D = mean(mu[grp == "D"]), loglik = best$loglik, loglik_gamma0 = r$fits$P$loglik, B = 0L, boot_fail = NA_integer_)
    gam[[paste(dm, v)]] <- data.frame(len = as.integer(dm), variant = v, parameter = colnames(Z), gamma = best$gamma, se = NA_real_, lo = NA_real_, hi = NA_real_)
    keep[[paste(dm, v)]] <- list(best = best, fits = fits) } }
res <- do.call(rbind, out); rownames(res) <- NULL; g <- do.call(rbind, gam); rownames(g) <- NULL; options(width = 220)
print(res[res$parameter == "Diabetic", c("len", "variant", "est", "se", "lo", "hi", "phi_H", "phi_P", "phi_D", "mu_H", "mu_P", "mu_D", "loglik", "loglik_gamma0")], digits = 3, row.names = FALSE); print(g, digits = 3, row.names = FALSE)
write.csv(res, "setpoint_boundary_directML_summary.csv", row.names = FALSE); write.csv(g, "setpoint_boundary_directML_gamma.csv", row.names = FALSE); saveRDS(list(res = res, gamma = g, detail = keep), "setpoint_boundary_directML.rds"); cat("saved\n")
