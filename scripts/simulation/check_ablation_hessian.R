source("R/simulation_design.R")
# Portable October 13 replication. Run from the repository root; see docs/SIMULATIONS.txt.
#!/usr/bin/env Rscript
# Independent, bounded covariance-ablation validation. No EM or resampling.
# Numeric Hessian of negative Q holds the saved-fit smoother moments fixed.
options(warn = 1, width = 180)
project <- getwd()
work <- file.path(project, "results/simulation/free_variance_ablation")
output_dir <- file.path(work, "validation")
source(file.path(project, "R/state_space_em.R"))
source(file.path(project, "R/free_variance_cluster_inference.R"))
stopifnot(requireNamespace("numDeriv", quietly = TRUE))
source(file.path(project, "scripts/simulation/helpers/ablation_covariance.R"))
design <- pss_simulation_design("positive")
X <- design$Xmat
clusters <- design$meta$subject
Tobs <- ncol(design$Y)
p <- ncol(X)
parameter_names <- c(colnames(X), "a1", "P1", "sigma_eta2", "sigma_eps2")

# Mirrors the published data-generation ordering only, without calling its
# executable simulation runner or modifying any cache.
reconstruct <- function(tau, m = 1L) {
  set.seed(410000 + 10000 * round(10 * tau) + m)
  unique_clusters <- unique(clusters)
  b <- setNames(rnorm(length(unique_clusters), 0, tau), unique_clusters)
  phi <- tanh(as.vector(X %*% design$truth$beta) + b[clusters])
  state <- rnorm(nrow(X), design$truth$a1, sqrt(design$truth$P1))
  y <- matrix(0, nrow(X), Tobs)
  y[, 1] <- state + rnorm(nrow(X), 0, 20)
  for (tt in 2:Tobs) {
    state <- phi * state + rnorm(nrow(X), 0, sqrt(design$truth$sigma_eta2))
    y[, tt] <- state + rnorm(nrow(X), 0, 20)
  }
  y
}

pooled <- function(x) {
  if (length(dim(x)) == 3L) return(apply(x, c(2, 3), sum))
  stopifnot(is.matrix(x), all(dim(x) == p + 4L))
  x
}

# Independent construction: expected squared residuals are assembled at each
# time directly from frozen first/second moments, without compute_C.
make_negative_Q <- function(y, Xmat, mu, kf) {
  a <- kf$alpha_smooth
  v <- kf$P_smooth
  lag_covariance <- kf$P_cov_smooth
  n <- nrow(y)
  tt <- ncol(y)
  E_observation_square <- sum((y - mu - a)^2 + v)
  function(theta) {
    beta <- theta[seq_len(ncol(Xmat))]
    initial_mean <- theta[ncol(Xmat) + 1L]
    initial_variance <- theta[ncol(Xmat) + 2L]
    process_variance <- theta[ncol(Xmat) + 3L]
    observation_variance <- theta[ncol(Xmat) + 4L]
    stopifnot(min(initial_variance, process_variance, observation_variance) > 0)
    phi <- tanh(as.vector(Xmat %*% beta))
    E_initial_square <- sum(v[, 1] + (a[, 1] - initial_mean)^2)
    E_transition_square <- 0
    for (t in 2:tt) {
      E_transition_square <- E_transition_square + sum(
        v[, t] + a[, t]^2 -
          2 * phi * (lag_covariance[, t] + a[, t] * a[, t - 1L]) +
          phi^2 * (v[, t - 1L] + a[, t - 1L]^2)
      )
    }
    0.5 * (n * log(initial_variance) + E_initial_square / initial_variance +
      n * (tt - 1L) * log(process_variance) + E_transition_square / process_variance +
      n * tt * log(observation_variance) + E_observation_square / observation_variance)
  }
}

summary_rows <- list()
entry_rows <- list()
fit_rows <- list()
for (tau in c(0, 0.4)) {
  y <- reconstruct(tau)
  cache <- readRDS(file.path(project, "results/simulation/freevar", paste0("tau", tau), "m001.rds"))
  fit <- cache[c("beta", "a1", "P1", "sigma_eta2", "sigma_eps2")]
  fit$kf <- kalman_run(y, tanh(as.vector(X %*% fit$beta)), 0, fit$a1, fit$P1, fit$sigma_eta2, fit$sigma_eps2)
  Fv <- fit$kf$P_pred + fit$sigma_eps2
  observed_loglik <- -0.5 * sum(log(2 * pi * Fv) + (y - fit$kf$alpha_pred)^2 / Fv)
  stopifnot(abs(observed_loglik - cache$loglik) < 1e-8)
  theta0 <- setNames(c(unname(fit$beta), fit$a1, fit$P1, fit$sigma_eta2, fit$sigma_eps2), parameter_names)
  scales <- pmax(abs(theta0), 1)
  subsets <- list(pooled = seq_len(nrow(y)), subset17 = unique(round(seq(1, nrow(y), length.out = 17))))
  for (subset_name in names(subsets)) {
    idx <- subsets[[subset_name]]
    ys <- y[idx, , drop = FALSE]
    Xs <- X[idx, , drop = FALSE]
    fs <- fit
    # Episode factorization means subsetting these already computed moments is
    # identical to smoothing only those episodes, and explicitly freezes them.
    for (key in names(fit$kf)) {
      if (is.matrix(fit$kf[[key]])) fs$kf[[key]] <- fit$kf[[key]][idx, , drop = FALSE]
      else if (length(fit$kf[[key]]) == nrow(y)) fs$kf[[key]] <- fit$kf[[key]][idx]
    }
    negQ <- make_negative_Q(ys, Xs, 0, fs$kf)
    analytic <- pooled(direct_information_free(ys, Xs, 0, fs))
    analytic_scaled <- analytic * outer(scales, scales)
    for (step in c(0.01, 0.02)) {
      numeric_scaled <- numDeriv::hessian(function(u) negQ(theta0 + scales * (u - 1)),
        rep(1, length(theta0)), method = "Richardson", method.args = list(d = step, r = 6L))
      numerical <- numeric_scaled / outer(scales, scales)
      difference_scaled <- numeric_scaled - analytic_scaled
      # Entrywise normalization respects parameter units and handles theoretical
      # zero cross terms without misleading relative errors near zero.
      normalization <- outer(sqrt(pmax(abs(diag(analytic_scaled)), 1)),
                             sqrt(pmax(abs(diag(analytic_scaled)), 1)))
      scaled_entry_error <- abs(difference_scaled) / normalization
      relative_frobenius <- sqrt(sum(difference_scaled^2)) / sqrt(sum(analytic_scaled^2))
      max_entry_error <- max(scaled_entry_error)
      pass <- relative_frobenius < 2e-5 && max_entry_error < 1e-4
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(tau, subset = subset_name,
        episodes = length(idx), step, relative_frobenius, max_entry_error, pass)
      # Report every raw-unit entry, including cross derivatives that disappear
      # approximately at a pooled stationary fit but survive on the subset.
      for (j in seq_along(parameter_names)) for (k in j:length(parameter_names)) {
        entry_rows[[length(entry_rows) + 1L]] <- data.frame(tau, subset = subset_name, step,
          row = parameter_names[j], column = parameter_names[k], analytic = analytic[j, k],
          numeric = numerical[j, k], normalized_error = scaled_entry_error[j, k])
      }
    }
    if (subset_name == "subset17") {
      stopifnot(max(abs(analytic[seq_len(p), p + 3L])) > 1e-5,
                abs(analytic[p + 1L, p + 2L]) > 1e-7)
    }
  }
  fit_rows[[length(fit_rows) + 1L]] <- data.frame(tau, m = 1, loglik_error = observed_loglik - cache$loglik,
    observation_variance = cache$sigma_eps2, status = cache$status)
}
summary_df <- do.call(rbind, summary_rows)
entry_df <- do.call(rbind, entry_rows)
write.csv(summary_df, file.path(output_dir, "direct_hessian_summary.csv"), row.names = FALSE)
write.csv(entry_df, file.path(output_dir, "direct_hessian_entries.csv"), row.names = FALSE)
write.csv(do.call(rbind, fit_rows), file.path(output_dir, "direct_hessian_reconstruction.csv"), row.names = FALSE)
print(summary_df, row.names = FALSE)
cat("\nNonzero subset cross-parameter terms (first step):\n")
print(subset(entry_df, subset == "subset17" & step == 0.01 &
  ((row %in% colnames(X) & column == "sigma_eta2") | (row == "a1" & column == "P1"))), row.names = FALSE)
cat("\nVersions:\n")
print(sessionInfo())
stopifnot(all(summary_df$pass))
cat("\nPASS: independent fixed-moment negative-Q Hessians agree with the proposed direct information.\n")
