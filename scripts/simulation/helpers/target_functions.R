source("R/simulation_design.R")
# Portable October 13 replication. Run from the repository root; see docs/SIMULATIONS.txt.
suppressPackageStartupMessages(library(statmod))
f <- pss_simulation_design("positive"); Xmat <- f$Xmat; subj <- f$meta$subject; Tt <- ncol(f$Y); pn <- colnames(Xmat); p <- length(pn)
rows <- !duplicated(subj); W <- Xmat[rows, , drop = FALSE]; n_g <- as.vector(table(subj)[subj[rows]]); stopifnot(sum(n_g) == nrow(Xmat))
ar1_moments <- function(phi, a, P, q, r) { tt <- 0:(Tt - 1); ph <- phi^tt; m <- a * ph; V <- P * outer(ph, ph)
  for (t in seq_len(Tt)) for (s in seq_len(Tt)) { k <- seq_len(min(t, s) - 1); if (length(k)) V[t, s] <- V[t, s] + q * sum(phi^((t - 1) + (s - 1) - 2 * k)) }
  diag(V) <- diag(V) + r; list(m = m, V = V) }
true_moments <- function(beta0, a0, P0, q0, r0, tau, nodes = 60) { gh <- gauss.quad.prob(nodes, "normal", mu = 0, sigma = 1); lapply(seq_len(nrow(W)), function(g) { eta <- sum(W[g, ] * beta0); mu <- numeric(Tt); S <- matrix(0, Tt, Tt)
  for (k in seq_len(nodes)) { b <- if (tau > 0) tau * gh$nodes[k] else 0; w <- gh$weights[k]; mo <- ar1_moments(tanh(eta + b), a0, P0, q0, r0); mu <- mu + w * mo$m; S <- S + w * (mo$V + outer(mo$m, mo$m)); if (tau == 0) break }
  if (tau == 0) { mo <- ar1_moments(tanh(eta), a0, P0, q0, r0); mu <- mo$m; S <- mo$V + outer(mo$m, mo$m) }; list(mu = mu, S = S) }) }
expected_nll <- function(par, TM, r, cap = 0.995) { beta <- par[1:p]; a <- par[p + 1]; P <- exp(par[p + 2]); q <- exp(par[p + 3]); eta <- as.vector(W %*% beta); if (max(abs(tanh(eta))) > cap) return(1e10); tot <- 0
  for (g in seq_len(nrow(W))) { mo <- ar1_moments(tanh(eta[g]), a, P, q, r); R <- chol(mo$V); d <- TM[[g]]$mu - mo$m; Sc <- TM[[g]]$S - outer(TM[[g]]$mu, TM[[g]]$mu) + outer(d, d)
    Vi <- chol2inv(R); tot <- tot + n_g[g] * (sum(log(diag(R))) + 0.5 * sum(Vi * Sc)) }; tot }
solve_target <- function(beta0, a0, P0, q0, r0, tau, nodes = 60, start = NULL) { TM <- true_moments(beta0, a0, P0, q0, r0, tau, nodes); st <- if (is.null(start)) c(beta0, a0, log(P0), log(q0)) else start
  o <- optim(st, expected_nll, TM = TM, r = r0, method = "BFGS", control = list(maxit = 2000, reltol = 1e-14)); o2 <- optim(o$par, expected_nll, TM = TM, r = r0, method = "BFGS", control = list(maxit = 2000, reltol = 1e-14))
  list(beta = setNames(o2$par[1:p], pn), a1 = o2$par[p + 1], P1 = exp(o2$par[p + 2]), q = exp(o2$par[p + 3]), value = o2$value, conv = o2$convergence, max_phi = max(abs(tanh(W %*% o2$par[1:p])))) }
