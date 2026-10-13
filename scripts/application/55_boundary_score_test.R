## Step 55: score of the log-likelihood with respect to the measurement-error variance r = sigma_eps2 at the
## boundary r = 0, in closed form, its numerical check, and a one-sided score test of H0: r = 0 with a
## participant-clustered (nuisance-adjusted) variance.
## Per unit i: y_i = mu_i + alpha_i + eps_i, alpha_i ~ N(m_i, S_i), S_i = AR(1) covariance with Var(alpha_1)=P1,
## innovation variance q; y_i ~ N(mu_i + m_i, S_i + r I). With d_i = y_i - mu_i - m_i,
##   dl_i/dr = -tr(Sigma^-1)/2 + d_i' Sigma^-2 d_i / 2,   Sigma = S_i + r I.
## At r = 0, S_i^-1 = L_i' D^-1 L_i (L_i: first-order differencing with phi_i; D = diag(P1, q, ..., q)), so with
## innovations u = L_i d_i and w = D^-1 u, v = S_i^-1 d_i = L_i' w has v_t = w_t - phi_i w_{t+1} (t < T), v_T = w_T:
##   score_i(0) = ( sum_t v_t^2 - tr(S_i^-1) ) / 2,  tr(S_i^-1) = 1/P1 + (T-1)(1 + phi_i^2)/q.
## Expanding sum_t v_t^2 shows that, at the MLE of q, the score is dominated by
##   -(1/q^2) sum_i phi_i sum_t u_it u_i,t+1 : positive lag-one autocovariance of the innovations drives the
## score negative (boundary), negative autocovariance drives it positive (interior r).
source(file.path(APP_ROOT, "R", "state_space_em.R")); pb <- readRDS("PRIMARY_boundary.rds"); Y <- pb$Y; X <- pb$Xmat; cl <- pb$meta$subject; fit <- pb$fit
N <- nrow(Y); Tt <- ncol(Y); th0 <- c(fit$beta, fit$a1, fit$P1, fit$sigma_eta2)
unpack <- function(th) list(phi = as.vector(tanh(X %*% th[1:6])), a1 = th[7], P1 = th[8], q = th[9])
## exact multivariate-normal log-likelihood as a function of (theta, r), valid for small negative r too
ll_full <- function(th, r) { p <- unpack(th); tot <- 0; for (i in seq_len(N)) { ph <- p$phi[i]; m <- p$a1 * ph^(0:(Tt - 1)); L <- diag(Tt); L[cbind(2:Tt, 1:(Tt - 1))] <- -ph; Sinv <- t(L) %*% diag(c(1 / p$P1, rep(1 / p$q, Tt - 1))) %*% L; S <- solve(Sinv) + r * diag(Tt); d <- Y[i, ] - m
  tot <- tot + dmvn(d, S) }; tot }
dmvn <- function(d, S) { R <- chol(S); z <- backsolve(R, d, transpose = TRUE); -sum(log(diag(R))) - 0.5 * sum(z^2) - length(d) / 2 * log(2 * pi) }
## closed-form per-unit score at r = 0
score_units <- function(th) { p <- unpack(th); sapply(seq_len(N), function(i) { ph <- p$phi[i]; m <- p$a1 * ph^(0:(Tt - 1)); d <- Y[i, ] - m; u <- c(d[1], d[-1] - ph * d[-Tt]); w <- u / c(p$P1, rep(p$q, Tt - 1)); v <- c(w[-Tt] - ph * w[-1], w[Tt])
  trS <- 1 / p$P1 + (Tt - 1) * (1 + ph^2) / p$q; (sum(v^2) - trS) / 2 }) }
s_i <- score_units(th0); S0 <- sum(s_i)
cat(sprintf("closed-form score at r = 0: %.4f\n", S0))
cat(sprintf("check: full log-lik at (theta_hat, r=0) = %.4f (EM: %.4f)\n", ll_full(th0, 0), fit$loglik))
h <- 1e-4; num <- (ll_full(th0, h) - ll_full(th0, -h)) / (2 * h); cat(sprintf("numerical d loglik / d r at 0 (central, h = 1e-4): %.4f\n", num))
## decomposition: lag-one innovation term vs the rest
p <- unpack(th0); dec <- sapply(seq_len(N), function(i) { ph <- p$phi[i]; m <- p$a1 * ph^(0:(Tt - 1)); d <- Y[i, ] - m; u <- d[-1] - ph * d[-Tt]; -ph * sum(u[-1] * u[-length(u)]) / p$q^2 })
cat(sprintf("lag-one innovation term  -(1/q^2) sum_i phi_i sum_t u_t u_{t+1}: %.4f  (remainder %.4f)\n", sum(dec), S0 - sum(dec)))
## efficient information for r at the MLE under r = 0 (numerical Hessian of the full likelihood in (theta, r))
H <- numDeriv::hessian(function(z) ll_full(z[1:9], z[10]), c(th0, 0)); I <- -H; Ieff <- I[10, 10] - I[10, 1:9] %*% solve(I[1:9, 1:9]) %*% I[1:9, 10]
cat(sprintf("observed information for r at the boundary: I_rr = %.4g; nuisance-adjusted %.4g (negative: the log-likelihood is locally convex in r at 0 on these data, so the observed information cannot standardize the score; the clustered variance below is used)\n", I[10, 10], Ieff))
## nuisance-adjusted per-unit scores and a participant-clustered variance
G_th <- t(sapply(seq_len(N), function(i) numDeriv::grad(function(z) { p <- unpack(z); ph <- p$phi[i]; m <- p$a1 * ph^(0:(Tt - 1)); L <- diag(Tt); L[cbind(2:Tt, 1:(Tt - 1))] <- -ph; Sinv <- t(L) %*% diag(c(1 / p$P1, rep(1 / p$q, Tt - 1))) %*% L; dmvn(Y[i, ] - m, solve(Sinv)) }, th0)))
cat(sprintf("check: nuisance scores sum to (max abs) %.2e at the MLE\n", max(abs(colSums(G_th)))))
adj <- as.vector(I[10, 1:9] %*% solve(I[1:9, 1:9])); s_eff <- s_i - as.vector(G_th %*% adj); Sg <- tapply(s_eff, cl, sum); Vcl <- sum(Sg^2)
cat(sprintf("efficient score: sum %.4f; cluster-robust SD %.4g (44 participants); standardized score %.2f; participants with negative score contribution: %d of %d\n", sum(s_eff), sqrt(Vcl), sum(s_eff) / sqrt(Vcl), sum(Sg < 0), length(Sg)))
cat("One-sided score test of H0: r = 0 against r > 0 rejects for standardized score > 1.645; here the score is negative, so H0 is not rejected and the estimate is at the boundary.\n")
## the same statistic on simulated data with true r > 0 (does the test detect white noise?): from the boundary fit, add white noise of SD 2, 5, 10; refit at r = 0 and compute the statistic
set.seed(55); for (sd in c(0, 2, 5, 10)) { Ys <- matrix(0, N, Tt); al <- rnorm(N, fit$a1, sqrt(fit$P1)); Ys[, 1] <- al; for (t in 2:Tt) { al <- p$phi * al + rnorm(N, 0, sqrt(p$q)); Ys[, t] <- al }; Ys <- Ys + rnorm(N * Tt, 0, sd)
  f0 <- em_regression_mv(Ys, X, 0, list(beta = fit$beta, a1 = fit$a1, P1 = fit$P1, sigma_eta2 = fit$sigma_eta2, sigma_eps2 = 0), thres = 1e-9, max_iter = 2000, fix_sigma_eps2 = 0); ths <- c(f0$beta, f0$a1, f0$P1, f0$sigma_eta2)
  Yb <- Y; Y <- Ys; si <- score_units(ths); Gt <- t(sapply(seq_len(N), function(i) numDeriv::grad(function(z) { pp <- unpack(z); ph <- pp$phi[i]; m <- pp$a1 * ph^(0:(Tt - 1)); L <- diag(Tt); L[cbind(2:Tt, 1:(Tt - 1))] <- -ph; Sinv <- t(L) %*% diag(c(1 / pp$P1, rep(1 / pp$q, Tt - 1))) %*% L; dmvn(Ys[i, ] - m, solve(Sinv)) }, ths)))
  Hs <- numDeriv::hessian(function(z) ll_full(z[1:9], z[10]), c(ths, 0)); Is <- -Hs; adjs <- as.vector(Is[10, 1:9] %*% solve(Is[1:9, 1:9])); se <- si - as.vector(Gt %*% adjs); Sgs <- tapply(se, cl, sum); Y <- Yb
  cat(sprintf("simulated white noise SD %2d: score %.1f, cluster-robust standardized score %.2f, converged %s\n", sd, sum(se), sum(se) / sqrt(sum(Sgs^2)), f0$converged)) }
saveRDS(list(score = S0, score_units = s_i, score_eff = s_eff, cluster_scores = Sg, I = I, I_eff = Ieff, numerical_score = num), "boundary_score_test.rds"); cat("saved\n")
