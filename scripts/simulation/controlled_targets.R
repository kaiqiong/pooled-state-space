# Portable October 13 replication. Run from the repository root; see docs/SIMULATIONS.txt.
## Controlled comparison: Part Q generating parameters (beta, a1, P1, q from the SD-20 fit), tau = 0.4, with the
## observation variance r0 set to 0, 100, 400 and 900 in BOTH generation and the working model (r fixed at r0).
## Isolates the effect of observation noise on the working-likelihood target from the other generating parameters.
source("scripts/simulation/helpers/target_functions.R")
beta0 <- setNames(f$truth$beta, pn); a0 <- f$truth$a1; P0 <- f$truth$P1; q0 <- f$truth$sigma_eta2
out <- t(sapply(c(0, 100, 400, 900), function(r0) { tg <- solve_target(beta0, a0, P0, q0, r0, tau = 0.4, nodes = 120); c(r0 = r0, tg$beta, ratio_diab = tg$beta["Diabetic"] / beta0["Diabetic"], ratio_age = tg$beta["Age"] / beta0["Age"], conv = tg$conv) }))
print(round(out, 4)); write.csv(out, "results/simulation/tables/expected_likelihood_target_by_r0.csv", row.names = FALSE)
