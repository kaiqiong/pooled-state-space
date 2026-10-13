# Portable October 13 replication. Run from the repository root; see docs/SIMULATIONS.txt.
source("scripts/simulation/helpers/target_functions.R")
beta0 <- setNames(f$truth$beta, pn); a1_0 <- f$truth$a1; P1_0 <- f$truth$P1; se2_0 <- f$truth$sigma_eta2; r0 <- 400
solve_target_free <- function(beta0, a0, P0, q0, r0, tau, nodes = 120) { TM <- true_moments(beta0, a0, P0, q0, r0, tau, nodes); st <- c(beta0, a0, log(P0), log(q0), log(r0))
  obj <- function(par) expected_nll(par[1:(p + 3)], TM, exp(par[p + 4])); o <- optim(st, obj, method = "BFGS", control = list(maxit = 3000, reltol = 1e-14)); o <- optim(o$par, obj, method = "BFGS", control = list(maxit = 3000, reltol = 1e-14))
  list(beta = setNames(o$par[1:p], pn), a1 = o$par[p + 1], P1 = exp(o$par[p + 2]), q = exp(o$par[p + 3]), r = exp(o$par[p + 4]), conv = o$convergence, max_phi = max(abs(tanh(W %*% o$par[1:p])))) }
chk <- solve_target_free(beta0, a1_0, P1_0, se2_0, r0, 0); cat(sprintf("\nfree-r target, tau=0: max |beta - beta0| %.1e, r %.2f (400), q %.3f (%.3f)\n", max(abs(chk$beta - beta0)), chk$r, chk$q, se2_0))
tf <- solve_target_free(beta0, a1_0, P1_0, se2_0, r0, 0.4); tfx <- readRDS("results/simulation/tables/expected_likelihood_targets.rds")$posvar
cat(sprintf("free-r target, tau=0.4: beta* = %s | r* = %.1f (400) q* = %.2f (38.5) a1* %.1f P1* %.0f | max phi %.3f | conv %d\nfixed-r target (r = 400):   %s\nratio free/generating slopes: %s\n", paste(round(tf$beta, 4), collapse = " "), tf$r, tf$q, tf$a1, tf$P1, tf$max_phi, tf$conv, paste(round(tfx$beta_star, 4), collapse = " "), paste(round(tf$beta[-1] / beta0[-1], 3), collapse = " ")))
saveRDS(list(free = tf, check_tau0 = chk), "results/simulation/tables/expected_likelihood_target_free_r.rds")
