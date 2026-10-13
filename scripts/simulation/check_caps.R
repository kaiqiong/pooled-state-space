source("R/simulation_design.R")
# Verify the published feasible cross-cap likelihood comparisons on final fits.
f <- pss_simulation_design("positive")
X <- f$Xmat; cl <- f$meta$subject; G <- length(unique(cl))
paths <- c("results/simulation/cap_sensitivity_ms/application.rds",
           list.files("results/simulation/cap_sensitivity_ms", "^(subset|stress)_.*[.]rds$", full.names = TRUE))
stopifnot(length(paths) == 211L, all(file.exists(paths)))
rows <- list()
for (path in paths) {
  z <- readRDS(path); rr <- z$res; is_app <- identical(z$set, "application")
  B <- if (is_app) 2000L else 999L
  seed <- if (is_app) 20260929 else 900000 + 5000000 + 10000 * round(10 * z$tau) + z$m
  set.seed(seed)
  ids <- lapply(seq_len(B), function(b) unlist(lapply(sample(unique(cl), G, TRUE), function(g) which(cl == g))))
  for (b in 0:B) {
    XX <- if (b == 0L) X else X[ids[[b]], , drop = FALSE]
    betas <- lapply(rr, function(r) if (b == 0L) r$beta else r$draws[b, ])
    ll <- vapply(rr, function(r) if (b == 0L) r$loglik else r$boot_loglik[b], numeric(1))
    caps <- vapply(rr, `[[`, numeric(1), "cap")
    phi <- vapply(betas, function(beta) max(abs(tanh(XX %*% beta))), numeric(1))
    stopifnot(all(is.finite(ll)), all(is.finite(phi)), all(phi <= caps + 1e-7))
    for (i in seq_along(rr)) {
      feasible <- phi <= caps[i] + 1e-12
      deficit <- max(ll[feasible] - ll[i])
      rows[[length(rows) + 1L]] <- data.frame(file = basename(path), draw = b, cap = caps[i],
        max_phi = phi[i], feasible_candidate_deficit = deficit, pass = deficit <= 1e-3)
    }
  }
}
d <- do.call(rbind, rows)
write.csv(d, "results/simulation/validation/cap_feasibility_checks.csv", row.names = FALSE)
cat(nrow(d), "cap-specific fit checks; maximum feasible likelihood deficit", max(d$feasible_candidate_deficit), "\n")
stopifnot(all(d$pass))
