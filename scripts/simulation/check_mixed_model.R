source("R/simulation_design.R")
# Portable October 13 replication. Run from the repository root; see docs/SIMULATIONS.txt.
## Checks for mixed_latent_ar.R: (1) at tau -> 0 the marginal log-likelihood equals the working-independence
## AR(1) log-likelihood (closed form); (2) the adaptive quadrature agrees with a brute-force 1-D integral on a
## fine grid; (3) parameter recovery on one simulated clustered dataset.
source("R/state_space_em.R"); source("R/mixed_latent_ar.R")
f <- pss_simulation_design("positive"); X <- f$Xmat; cl <- f$meta$subject; pb <- pss_simulation_design("boundary"); Y <- pb$Y; S <- mixed_ar_suff(Y)
kf0 <- kalman_run(Y, tanh(as.vector(X %*% pb$truth$beta)), 0, pb$truth$a1, pb$truth$P1, pb$truth$sigma_eta2, 0)
ll0 <- -0.5 * sum(log(2 * pi * kf0$P_pred) + (Y - kf0$alpha_pred)^2 / kf0$P_pred)
th <- c(pb$truth$beta, log(pb$truth$sigma_eta2), log(1e-6)); cat(sprintf("tau -> 0: mixed loglik %.4f vs AR(1) boundary fit %.4f\n", mixed_ar_loglik(th, X, S, cl), ll0))
## brute force at tau = 0.3
th2 <- c(pb$truth$beta, log(pb$truth$sigma_eta2), log(0.3)); q <- pb$truth$sigma_eta2; eta <- as.vector(X %*% pb$truth$beta); gi <- match(cl, unique(cl)); grid <- seq(-2, 2, length.out = 4001); dg <- grid[2] - grid[1]
lint <- sapply(seq_along(unique(cl)), function(g) { i <- which(gi == g); lf <- sapply(grid, function(b) { ph <- tanh(eta[i] + b); sum(-(S$C0[i] - 2 * ph * S$C1[i] + ph^2 * S$C2[i]) / (2 * q)) + dnorm(b, 0, 0.3, log = TRUE) }); mx <- max(lf); mx + log(sum(exp(lf - mx)) * dg) })
a1 <- mean(S$y1); P1 <- mean((S$y1 - a1)^2); ll_bf <- sum(lint) - nrow(Y) * S$n / 2 * log(2 * pi * q) - nrow(Y) / 2 * (log(2 * pi * P1) + 1)
cat(sprintf("tau = 0.3: adaptive GH (K = 15) %.6f vs grid integration %.6f (K = 5: %.6f)\n", mixed_ar_loglik(th2, X, S, cl, K = 15), ll_bf, mixed_ar_loglik(th2, X, S, cl, K = 5)))
## recovery on simulated clustered data (Part M boundary design, tau = 0.4)
set.seed(1); subj <- cl; G <- length(unique(subj)); N <- nrow(X); Tt <- ncol(Y); b_g <- setNames(rnorm(G, 0, 0.4), unique(subj)); phi <- tanh(as.vector(X %*% pb$truth$beta) + b_g[subj]); ys <- matrix(0, N, Tt); ys[, 1] <- rnorm(N, pb$truth$a1, sqrt(pb$truth$P1)); for (t in 2:Tt) ys[, t] <- phi * ys[, t - 1] + rnorm(N, 0, sqrt(q))
fm <- fit_mixed_ar(ys, X, subj, start_theta = c(atanh(0.9), rep(0, 5))); cat("simulated tau = 0.4: theta_hat", round(fm$theta, 3), " (truth", round(pb$truth$beta, 3), ")\n  se", round(fm$se_theta, 3), " tau_hat", round(fm$tau, 3), " q_hat", round(fm$sigma_eta2, 2), " conv", fm$conv, " secs", round(fm$secs, 1), " evals", fm$evals, "\n")
cat("  correlation of estimated modes with true b_g:", round(cor(fm$modes, b_g[unique(subj)]), 3), "\n")
