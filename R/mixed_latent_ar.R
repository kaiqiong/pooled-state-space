## Mixed-effects latent autoregression (comparator): phi_ig = tanh(W_i' theta + b_g), b_g ~ N(0, tau^2), units i
## nested in clusters g, with the measurement-error variance at its boundary (sigma_eps2 = 0), so that the
## per-unit likelihood is that of a first-order autoregression: y_i1 ~ N(a1, P1), y_it | y_{i,t-1} ~ N(phi_i y_{i,t-1}, q).
## The marginal likelihood integrates b_g out by ADAPTIVE Gauss-Hermite quadrature. a1 and P1 enter only through
## y_i1 and are profiled out. Parameters: (theta, log q, log tau). The working-independence model is tau = 0.
## For numerical stability, the per-cluster mode
## is now found by bracketed 1-D maximization followed by Newton polishing with a score check (|score| < 1e-8),
## the quadrature uses K = 20 nodes scaled by the observed curvature at the mode (expected curvature if that is
## not negative), the Hessian is computed with a Richardson scheme on a coarser step and checked for positive
## definiteness, and a fit is reported as valid only if the optimizer converged, every cluster mode passed the
## score check, and the Hessian is positive definite with finite standard errors.
gh_nodes <- function(K) { i <- seq_len(K - 1); a <- sqrt(i / 2); J <- matrix(0, K, K); J[cbind(i, i + 1)] <- a; J[cbind(i + 1, i)] <- a; e <- eigen(J, symmetric = TRUE); list(z = e$values, w = sqrt(pi) * e$vectors[1, ]^2) }
mixed_ar_suff <- function(y) { Tt <- ncol(y); list(C0 = rowSums(y[, -1]^2), C1 = rowSums(y[, -1] * y[, -Tt]), C2 = rowSums(y[, -Tt]^2), y1 = y[, 1], n = Tt - 1) }
mixed_ar_loglik <- function(par, X, S, cl, K = 20, return_modes = FALSE) {
  p <- ncol(X); theta <- par[1:p]; q <- exp(par[p + 1]); tau <- exp(par[p + 2]); eta <- as.vector(X %*% theta); ug <- unique(cl); G <- length(ug); gi <- match(cl, ug); gh <- gh_nodes(K)
  f_win <- function(phi) -(S$C0 - 2 * phi * S$C1 + phi^2 * S$C2) / (2 * q)
  ## per-cluster log-integrand and its first two derivatives in b
  fg <- function(b, idx) { ph <- tanh(eta[idx] + b); sum(f_win_i(ph, idx)) - b^2 / (2 * tau^2) }
  f_win_i <- function(ph, idx) -(S$C0[idx] - 2 * ph * S$C1[idx] + ph^2 * S$C2[idx]) / (2 * q)
  d12 <- function(b, idx) { ph <- tanh(eta[idx] + b); h <- 1 - ph^2; g1 <- (S$C1[idx] - ph * S$C2[idx]) / q; g2 <- -S$C2[idx] / q
    c(sum(g1 * h) - b / tau^2, sum(g2 * h^2 - 2 * ph * h * g1) - 1 / tau^2, -(sum(S$C2[idx] * h^2) / q + 1 / tau^2)) }   # score, observed curvature, expected curvature
  b <- numeric(G); sig <- numeric(G); ok <- logical(G); idx_list <- split(seq_along(cl), gi)
  for (g in seq_len(G)) { idx <- idx_list[[g]]; o <- optimize(function(bb) -fg(bb, idx), interval = c(-8 * tau - 1, 8 * tau + 1), tol = 1e-10); bm <- o$minimum
    for (it in 1:20) { d <- d12(bm, idx); if (abs(d[1]) < 1e-10) break; curv <- if (d[2] < 0) d[2] else d[3]; bm <- bm - d[1] / curv }
    d <- d12(bm, idx); ok[g] <- abs(d[1]) < 1e-8; b[g] <- bm; sig[g] <- 1 / sqrt(-(if (d[2] < 0) d[2] else d[3])) }
  nodes <- outer(sig, sqrt(2) * gh$z) + b
  lf <- matrix(0, G, K); for (k in seq_len(K)) { bk <- nodes[, k]; ph <- tanh(eta + bk[gi]); lf[, k] <- rowsum(f_win(ph), gi)[, 1] - bk^2 / (2 * tau^2) }
  mx <- apply(lf, 1, max); lint <- mx + log(rowSums(exp(lf - mx) * rep(exp(gh$z^2) * gh$w, each = G))) + log(sqrt(2) * sig) - 0.5 * log(2 * pi * tau^2)
  a1 <- mean(S$y1); P1 <- mean((S$y1 - a1)^2); N <- length(S$y1)
  ll <- sum(lint) - N * S$n / 2 * log(2 * pi * q) - N / 2 * (log(2 * pi * P1) + 1)
  if (return_modes) list(loglik = ll, b = b, sig = sig, modes_ok = ok, a1 = a1, P1 = P1) else ll }
fit_mixed_ar <- function(y, X, cl, start_theta, K = 20, tau_start = 0.3, tau_floor = 1e-3) {
  S <- mixed_ar_suff(y); p <- ncol(X); obj <- function(par) { v <- -mixed_ar_loglik(par, X, S, cl, K); if (is.finite(v)) v else 1e300 }
  q0 <- mean((y[, -1] - tanh(as.vector(X %*% start_theta)) * y[, -ncol(y)])^2); start <- c(start_theta, log(q0), log(tau_start)); t0 <- Sys.time()
  lower <- c(rep(-Inf, p + 1), log(tau_floor)); upper <- rep(Inf, p + 2)
  o <- optim(start, obj, method = "L-BFGS-B", lower = lower, upper = upper, control = list(maxit = 3000, factr = 1e5))
  o2 <- optim(o$par, obj, method = "L-BFGS-B", lower = lower, upper = upper, control = list(maxit = 3000, factr = 1e4))
  par <- o2$par; conv <- o2$convergence == 0
  H <- tryCatch(numDeriv::hessian(obj, par, method.args = list(eps = 1e-4, d = 0.01, r = 4)), error = function(e) NULL)
  pd <- !is.null(H) && all(is.finite(H)) && all(eigen(H, symmetric = TRUE, only.values = TRUE)$values > 0)
  se <- if (pd) sqrt(diag(solve(H))) else rep(NA_real_, p + 2)
  m <- mixed_ar_loglik(par, X, S, cl, K, return_modes = TRUE); valid <- conv && all(m$modes_ok) && pd && all(is.finite(se[1:p]))
  list(theta = setNames(par[1:p], colnames(X)), sigma_eta2 = exp(par[p + 1]), tau = exp(par[p + 2]), tau_at_floor = par[p + 2] <= log(tau_floor) + 1e-6, a1 = m$a1, P1 = m$P1, loglik = m$loglik,
       se_theta = setNames(se[1:p], colnames(X)), se_all = se, conv = conv, modes_ok = all(m$modes_ok), hessian_pd = pd, valid = valid, modes = m$b, mode_sd = m$sig, secs = as.numeric(Sys.time() - t0, units = "secs"), evals = o$counts[1] + o2$counts[1]) }
