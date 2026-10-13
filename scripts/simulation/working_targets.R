source("R/simulation_design.R")
# Portable October 13 replication. Run from the repository root; see docs/SIMULATIONS.txt.
## Deterministic working-likelihood target for Part Q (and the Part M v2 boundary design), 7 Oct 2026.
## The working model treats every episode as an independent Gaussian vector with mean m(theta) and covariance
## V(theta) (AR(1) state started at N(a1, P1), innovation variance q, observation variance r fixed). Its expected
## negative log-likelihood under the generating mixture (participant random effect b ~ N(0, tau^2) on the link scale)
## depends on the generating distribution only through the first two moments of each episode,
##   E_b[E(Y_t | b)] = a0 E_b[phi^t],   E_b[E(Y_t Y_s | b)] = (P0 + a0^2) E_b[phi^(t+s)] + q0 sum_{k=1}^{min(t,s)} E_b[phi^(t+s-2k)] + r0 1(t=s),
## (t, s = 0, ..., T-1), which we integrate over b by Gauss-Hermite quadrature. The pseudo-true value theta* minimizes
##   sum_g n_g { 0.5 log det V_g + 0.5 tr[V_g^{-1} (S_g - m_g mu_g' - mu_g m_g' + m_g m_g')] }
## over (beta, a1, log P1, log q), where g runs over the 44 participants (covariate rows), n_g is the number of windows,
## mu_g, S_g are the true first and second moments and m_g, V_g the working ones. Within-participant dependence does not
## enter: the expectation of a sum of episode log-likelihoods is the sum of expectations.
## Checks: tau = 0 returns the generating parameters; two starting values; 60 vs 120 nodes.
source("scripts/simulation/helpers/target_functions.R")
res <- list()
## ---- Part Q design (truth = SD-20 fit; r0 = 400 fixed) ----
beta0 <- setNames(f$truth$beta, pn); a0 <- f$truth$a1; P0 <- f$truth$P1; q0 <- f$truth$sigma_eta2; r0 <- 400
chk0 <- solve_target(beta0, a0, P0, q0, r0, tau = 0); cat(sprintf("Part Q tau=0 check: max |beta - beta0| = %.2e, a1 %.4f vs %.4f, P1 %.3f vs %.3f, q %.4f vs %.4f\n", max(abs(chk0$beta - beta0)), chk0$a1, a0, chk0$P1, P0, chk0$q, q0))
tq60 <- solve_target(beta0, a0, P0, q0, r0, tau = 0.4, nodes = 60); tq120 <- solve_target(beta0, a0, P0, q0, r0, tau = 0.4, nodes = 120)
tq_alt <- solve_target(beta0, a0, P0, q0, r0, tau = 0.4, nodes = 120, start = c(beta0 * 0.8, a0 * 1.1, log(P0 * 0.7), log(q0 * 1.5)))
cat(sprintf("Part Q tau=0.4: 60 vs 120 nodes max |diff| %.2e; second start max |diff| %.2e; conv %d/%d/%d; max phi %.4f\n", max(abs(tq60$beta - tq120$beta)), max(abs(tq_alt$beta - tq120$beta)), tq60$conv, tq120$conv, tq_alt$conv, tq120$max_phi))
res$posvar <- list(tau = 0.4, beta_star = tq120$beta, a1 = tq120$a1, P1 = tq120$P1, q = tq120$q, generating = list(beta = beta0, a1 = a0, P1 = P0, q = q0, r = r0), nodes = 120, check_tau0 = max(abs(chk0$beta - beta0)), node_sens = max(abs(tq60$beta - tq120$beta)), start_sens = max(abs(tq_alt$beta - tq120$beta)))
## ---- Part M v2 boundary design (truth = boundary fit; r0 = 0) ----
source("R/state_space_em.R"); fb <- pss_simulation_design("boundary")$truth
beta0b <- setNames(fb$beta, pn); chk0b <- solve_target(beta0b, fb$a1, fb$P1, fb$sigma_eta2, 0, tau = 0)
tb120 <- solve_target(beta0b, fb$a1, fb$P1, fb$sigma_eta2, 0, tau = 0.4, nodes = 120); tb60 <- solve_target(beta0b, fb$a1, fb$P1, fb$sigma_eta2, 0, tau = 0.4, nodes = 60)
cat(sprintf("Boundary tau=0 check: max |beta - beta0| = %.2e; tau=0.4: 60 vs 120 nodes %.2e; conv %d; max phi %.4f\n", max(abs(chk0b$beta - beta0b)), max(abs(tb60$beta - tb120$beta)), tb120$conv, tb120$max_phi))
print(round(rbind(generating = beta0b, quadrature = tb120$beta, rel_gap = (tb120$beta - beta0b) / beta0b), 4))
res$boundary <- list(tau = 0.4, beta_star = tb120$beta, a1 = tb120$a1, P1 = tb120$P1, q = tb120$q, generating = list(beta = beta0b, a1 = fb$a1, P1 = fb$P1, q = fb$sigma_eta2, r = 0), nodes = 120, check_tau0 = max(abs(chk0b$beta - beta0b)), node_sens = max(abs(tb60$beta - tb120$beta)))
saveRDS(res, "results/simulation/tables/expected_likelihood_targets.rds"); cat("saved results/simulation/tables/expected_likelihood_targets.rds\n")
