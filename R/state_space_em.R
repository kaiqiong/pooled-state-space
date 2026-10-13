## Kalman filtering, smoothing, EM estimation and observed-information inference.
## Derivations: Supplementary Appendices A--C of the accompanying manuscript.
## All function bodies match the validated research implementation.
## Rows of y are episodes; cluster identifies participants.

## ---------------------------------------------------------------------------
## 1. Kalman filter + smoother + lag-one covariance smoother (vectorized over i)
## ---------------------------------------------------------------------------

#' @param y N x T matrix of observations
#' @param phi scalar (homogeneous) or length-N vector (phi_i, regression case)
kalman_run <- function(y, phi, mu, a1, P1, sigma_eta2, sigma_eps2) {
  N <- nrow(y); Tt <- ncol(y)
  if (length(phi) == 1L) phi <- rep(phi, N)

  alpha_pred <- matrix(0, N, Tt); P_pred <- matrix(0, N, Tt)
  alpha_filt <- matrix(0, N, Tt); P_filt <- matrix(0, N, Tt)
  K <- matrix(0, N, Tt)

  alpha_pred[, 1] <- a1
  P_pred[, 1] <- P1
  K[, 1] <- P_pred[, 1] / (P_pred[, 1] + sigma_eps2)
  alpha_filt[, 1] <- alpha_pred[, 1] + K[, 1] * (y[, 1] - (mu + alpha_pred[, 1]))
  P_filt[, 1] <- (1 - K[, 1]) * P_pred[, 1]

  for (t in 2:Tt) {
    alpha_pred[, t] <- phi * alpha_filt[, t - 1]
    P_pred[, t] <- phi^2 * P_filt[, t - 1] + sigma_eta2
    K[, t] <- P_pred[, t] / (P_pred[, t] + sigma_eps2)
    alpha_filt[, t] <- alpha_pred[, t] + K[, t] * (y[, t] - (mu + alpha_pred[, t]))
    P_filt[, t] <- (1 - K[, t]) * P_pred[, t]
  }

  alpha_smooth <- matrix(0, N, Tt); P_smooth <- matrix(0, N, Tt)
  J <- matrix(0, N, Tt)  # J[,s] = J_s, s = 1..T-1
  alpha_smooth[, Tt] <- alpha_filt[, Tt]
  P_smooth[, Tt] <- P_filt[, Tt]

  for (t in Tt:2) {
    J[, t - 1] <- phi * P_filt[, t - 1] / P_pred[, t]
    alpha_smooth[, t - 1] <- alpha_filt[, t - 1] + J[, t - 1] * (alpha_smooth[, t] - alpha_pred[, t])
    P_smooth[, t - 1] <- P_filt[, t - 1] + J[, t - 1]^2 * (P_smooth[, t] - P_pred[, t])
  }

  P_cov_smooth <- matrix(0, N, Tt)  # P_cov_smooth[,s] = P^T_{s,s-1}, s = 2..T
  P_cov_smooth[, Tt] <- (1 - K[, Tt]) * phi * P_filt[, Tt - 1]
  if (Tt >= 3) {
    for (t in Tt:3) {
      P_cov_smooth[, t - 1] <- P_filt[, t - 1] * J[, t - 2] +
        J[, t - 1] * (P_cov_smooth[, t] - phi * P_filt[, t - 1]) * J[, t - 2]
    }
  }

  list(alpha_pred = alpha_pred, P_pred = P_pred, K = K,
       alpha_filt = alpha_filt, P_filt = P_filt, J = J,
       alpha_smooth = alpha_smooth, P_smooth = P_smooth, P_cov_smooth = P_cov_smooth,
       phi = phi)
}

## Per-subject sufficient statistics C0_i, C1_i, C2_i (length-N vectors); sum over i
## for the homogeneous model's scalar C0, C1, C2.
compute_C <- function(alpha_smooth, P_smooth, P_cov_smooth) {
  Tt <- ncol(alpha_smooth)
  idx_t <- 2:Tt; idx_tm1 <- 1:(Tt - 1)
  a_t <- alpha_smooth[, idx_t, drop = FALSE]; a_tm1 <- alpha_smooth[, idx_tm1, drop = FALSE]
  P_t <- P_smooth[, idx_t, drop = FALSE]; P_tm1 <- P_smooth[, idx_tm1, drop = FALSE]
  Pcov_t <- P_cov_smooth[, idx_t, drop = FALSE]
  list(
    C0 = rowSums(a_t^2 + P_t),
    C1 = rowSums(a_t * a_tm1 + Pcov_t),
    C2 = rowSums(a_tm1^2 + P_tm1)
  )
}

## Product-rule ("dot") version of compute_C: given the base smoothed arrays and their
## sensitivities (d.../dtheta*), returns dC0_i, dC1_i, dC2_i.
compute_C_dot <- function(alpha_smooth, P_smooth, P_cov_smooth,
                           dalpha_smooth, dP_smooth, dP_cov_smooth) {
  Tt <- ncol(alpha_smooth)
  idx_t <- 2:Tt; idx_tm1 <- 1:(Tt - 1)
  a_t <- alpha_smooth[, idx_t, drop = FALSE]; a_tm1 <- alpha_smooth[, idx_tm1, drop = FALSE]
  da_t <- dalpha_smooth[, idx_t, drop = FALSE]; da_tm1 <- dalpha_smooth[, idx_tm1, drop = FALSE]
  dP_t <- dP_smooth[, idx_t, drop = FALSE]; dP_tm1 <- dP_smooth[, idx_tm1, drop = FALSE]
  dPcov_t <- dP_cov_smooth[, idx_t, drop = FALSE]
  list(
    dC0 = rowSums(2 * a_t * da_t + dP_t),
    dC1 = rowSums(da_t * a_tm1 + a_t * da_tm1 + dPcov_t),
    dC2 = rowSums(2 * a_tm1 * da_tm1 + dP_tm1)
  )
}

## ---------------------------------------------------------------------------
## 2. EM drivers
## ---------------------------------------------------------------------------

## Homogeneous model: theta = (phi, a1, P1, sigma_eta2, sigma_eps2), all shared.
## Every M-step update is closed-form -- no `optim` call anywhere.
em_homogeneous <- function(y, mu = 0, init, thres = 1e-8, max_iter = 500) {
  N <- nrow(y); Tt <- ncol(y)
  phi <- init$phi; a1 <- init$a1; P1 <- init$P1
  sigma_eta2 <- init$sigma_eta2; sigma_eps2 <- init$sigma_eps2

  for (iter in 1:max_iter) {
    kf <- kalman_run(y, phi, mu, a1, P1, sigma_eta2, sigma_eps2)
    Cs <- compute_C(kf$alpha_smooth, kf$P_smooth, kf$P_cov_smooth)
    C0 <- sum(Cs$C0); C1 <- sum(Cs$C1); C2 <- sum(Cs$C2)

    a1_new <- mean(kf$alpha_smooth[, 1])
    P1_new <- mean(kf$P_smooth[, 1] + (kf$alpha_smooth[, 1] - a1_new)^2)
    phi_new <- C1 / C2
    sigma_eta2_new <- (C0 - 2 * phi_new * C1 + phi_new^2 * C2) / (N * (Tt - 1))
    R <- sum(kf$P_smooth + (y - mu - kf$alpha_smooth)^2)
    sigma_eps2_new <- R / (N * Tt)

    delta <- max(abs(c(phi_new - phi, a1_new - a1, P1_new - P1,
                        sigma_eta2_new - sigma_eta2, sigma_eps2_new - sigma_eps2)))
    phi <- phi_new; a1 <- a1_new; P1 <- P1_new
    sigma_eta2 <- sigma_eta2_new; sigma_eps2 <- sigma_eps2_new
    if (delta < thres) break
  }

  list(phi = phi, a1 = a1, P1 = P1, sigma_eta2 = sigma_eta2, sigma_eps2 = sigma_eps2,
       iter = iter, converged = delta < thres,
       kf = kalman_run(y, phi, mu, a1, P1, sigma_eta2, sigma_eps2))
}

## Regression model: phi_i = tanh(beta0 + beta1*x_i). a1, P1, sigma_eps2 closed-form;
## (beta0,beta1) via `optim` with sigma_eta2 profiled out analytically (as in the
## original code), then sigma_eta2 recomputed at the optimum.
em_regression <- function(y, x, mu = 0, init, thres = 1e-8, max_iter = 500) {
  N <- nrow(y); Tt <- ncol(y)
  beta0 <- init$beta0; beta1 <- init$beta1; a1 <- init$a1; P1 <- init$P1
  sigma_eta2 <- init$sigma_eta2; sigma_eps2 <- init$sigma_eps2

  neg_Q_beta_profiled <- function(beta, x, alpha_smooth, P_smooth, P_cov_smooth) {
    phi_i <- tanh(beta[1] + beta[2] * x)
    Cs <- compute_C(alpha_smooth, P_smooth, P_cov_smooth)
    total <- sum(Cs$C0 - 2 * phi_i * Cs$C1 + phi_i^2 * Cs$C2)
    sigma_eta2_hat <- total / (N * (Tt - 1))
    val <- (N * (Tt - 1) / 2) * log(sigma_eta2_hat) + total / (2 * sigma_eta2_hat)
    if (!is.finite(val)) 1e300 else val
  }

  for (iter in 1:max_iter) {
    phi <- tanh(beta0 + beta1 * x)
    kf <- kalman_run(y, phi, mu, a1, P1, sigma_eta2, sigma_eps2)

    a1_new <- mean(kf$alpha_smooth[, 1])
    P1_new <- mean(kf$P_smooth[, 1] + (kf$alpha_smooth[, 1] - a1_new)^2)
    R <- sum(kf$P_smooth + (y - mu - kf$alpha_smooth)^2)
    sigma_eps2_new <- R / (N * Tt)

    opt <- optim(c(beta0, beta1), neg_Q_beta_profiled, x = x,
                 alpha_smooth = kf$alpha_smooth, P_smooth = kf$P_smooth,
                 P_cov_smooth = kf$P_cov_smooth,
                 method = "BFGS", control = list(reltol = 1e-12))
    beta0_new <- opt$par[1]; beta1_new <- opt$par[2]
    phi_new <- tanh(beta0_new + beta1_new * x)
    Cs <- compute_C(kf$alpha_smooth, kf$P_smooth, kf$P_cov_smooth)
    total <- sum(Cs$C0 - 2 * phi_new * Cs$C1 + phi_new^2 * Cs$C2)
    sigma_eta2_new <- total / (N * (Tt - 1))

    delta <- max(abs(c(beta0_new - beta0, beta1_new - beta1, a1_new - a1, P1_new - P1,
                        sigma_eta2_new - sigma_eta2, sigma_eps2_new - sigma_eps2)))
    beta0 <- beta0_new; beta1 <- beta1_new; a1 <- a1_new; P1 <- P1_new
    sigma_eta2 <- sigma_eta2_new; sigma_eps2 <- sigma_eps2_new
    if (delta < thres) break
  }

  phi <- tanh(beta0 + beta1 * x)
  list(beta0 = beta0, beta1 = beta1, a1 = a1, P1 = P1,
       sigma_eta2 = sigma_eta2, sigma_eps2 = sigma_eps2,
       iter = iter, converged = delta < thres,
       kf = kalman_run(y, phi, mu, a1, P1, sigma_eta2, sigma_eps2))
}

## ---------------------------------------------------------------------------
## 3. Sensitivity recursions (SE_derivation/louis_oakes_SE.tex, Section 4.1)
## ---------------------------------------------------------------------------

## Given the base filter/smoother output `kf` (already evaluated at theta* = theta_hat)
## and a seed direction (dphi, dsigma_eta2, dsigma_eps2, da1, dP1) specifying which
## theta*_l we differentiate with respect to, returns d(smoothed quantities)/d theta*_l.
kalman_sensitivity <- function(kf, mu, sigma_eta2, sigma_eps2, y,
                                dphi, dsigma_eta2, dsigma_eps2, da1, dP1) {
  N <- nrow(y); Tt <- ncol(y)
  phi <- kf$phi
  if (length(dphi) == 1L) dphi <- rep(dphi, N)

  alpha_pred <- kf$alpha_pred; P_pred <- kf$P_pred; K <- kf$K
  alpha_filt <- kf$alpha_filt; P_filt <- kf$P_filt; J <- kf$J
  alpha_smooth <- kf$alpha_smooth; P_smooth <- kf$P_smooth; P_cov_smooth <- kf$P_cov_smooth

  dalpha_pred <- matrix(0, N, Tt); dP_pred <- matrix(0, N, Tt)
  dalpha_filt <- matrix(0, N, Tt); dP_filt <- matrix(0, N, Tt); dK <- matrix(0, N, Tt)

  dalpha_pred[, 1] <- da1
  dP_pred[, 1] <- dP1
  dK[, 1] <- (dP_pred[, 1] * sigma_eps2 - P_pred[, 1] * dsigma_eps2) / (P_pred[, 1] + sigma_eps2)^2
  e1 <- y[, 1] - mu - alpha_pred[, 1]
  dalpha_filt[, 1] <- (1 - K[, 1]) * dalpha_pred[, 1] + dK[, 1] * e1
  dP_filt[, 1] <- (1 - K[, 1]) * dP_pred[, 1] - dK[, 1] * P_pred[, 1]

  for (t in 2:Tt) {
    dalpha_pred[, t] <- dphi * alpha_filt[, t - 1] + phi * dalpha_filt[, t - 1]
    dP_pred[, t] <- 2 * phi * dphi * P_filt[, t - 1] + phi^2 * dP_filt[, t - 1] + dsigma_eta2
    dK[, t] <- (dP_pred[, t] * sigma_eps2 - P_pred[, t] * dsigma_eps2) / (P_pred[, t] + sigma_eps2)^2
    et <- y[, t] - mu - alpha_pred[, t]
    dalpha_filt[, t] <- (1 - K[, t]) * dalpha_pred[, t] + dK[, t] * et
    dP_filt[, t] <- (1 - K[, t]) * dP_pred[, t] - dK[, t] * P_pred[, t]
  }

  dalpha_smooth <- matrix(0, N, Tt); dP_smooth <- matrix(0, N, Tt); dJ <- matrix(0, N, Tt)
  dalpha_smooth[, Tt] <- dalpha_filt[, Tt]
  dP_smooth[, Tt] <- dP_filt[, Tt]

  for (t in Tt:2) {
    dJ[, t - 1] <- ((dphi * P_filt[, t - 1] + phi * dP_filt[, t - 1]) * P_pred[, t] -
                       phi * P_filt[, t - 1] * dP_pred[, t]) / P_pred[, t]^2
    dalpha_smooth[, t - 1] <- dalpha_filt[, t - 1] +
      dJ[, t - 1] * (alpha_smooth[, t] - alpha_pred[, t]) +
      J[, t - 1] * (dalpha_smooth[, t] - dalpha_pred[, t])
    dP_smooth[, t - 1] <- dP_filt[, t - 1] +
      2 * J[, t - 1] * dJ[, t - 1] * (P_smooth[, t] - P_pred[, t]) +
      J[, t - 1]^2 * (dP_smooth[, t] - dP_pred[, t])
  }

  dP_cov_smooth <- matrix(0, N, Tt)
  dP_cov_smooth[, Tt] <- -dK[, Tt] * phi * P_filt[, Tt - 1] +
    (1 - K[, Tt]) * (dphi * P_filt[, Tt - 1] + phi * dP_filt[, Tt - 1])
  if (Tt >= 3) {
    for (t in Tt:3) {
      M <- P_cov_smooth[, t] - phi * P_filt[, t - 1]
      dM <- dP_cov_smooth[, t] - dphi * P_filt[, t - 1] - phi * dP_filt[, t - 1]
      dP_cov_smooth[, t - 1] <- dP_filt[, t - 1] * J[, t - 2] + P_filt[, t - 1] * dJ[, t - 2] +
        dJ[, t - 1] * M * J[, t - 2] + J[, t - 1] * dM * J[, t - 2] + J[, t - 1] * M * dJ[, t - 2]
    }
  }

  list(dalpha_smooth = dalpha_smooth, dP_smooth = dP_smooth, dP_cov_smooth = dP_cov_smooth)
}

## ---------------------------------------------------------------------------
## 4. Observed information / SEs via Oakes' identity: I_obs = D + X
## ---------------------------------------------------------------------------

compute_SE_homogeneous <- function(y, mu, fit) {
  N <- nrow(y); Tt <- ncol(y)
  phi <- fit$phi; a1 <- fit$a1; P1 <- fit$P1
  sigma_eta2 <- fit$sigma_eta2; sigma_eps2 <- fit$sigma_eps2
  kf <- fit$kf

  Cs <- compute_C(kf$alpha_smooth, kf$P_smooth, kf$P_cov_smooth)
  C2 <- sum(Cs$C2)

  pn <- c("phi", "a1", "P1", "sigma_eta2", "sigma_eps2")
  D <- matrix(0, 5, 5, dimnames = list(pn, pn))
  D["a1", "a1"] <- N / P1
  D["P1", "P1"] <- N / (2 * P1^2)
  D["phi", "phi"] <- C2 / sigma_eta2
  D["sigma_eta2", "sigma_eta2"] <- N * (Tt - 1) / (2 * sigma_eta2^2)
  D["sigma_eps2", "sigma_eps2"] <- N * Tt / (2 * sigma_eps2^2)

  seeds <- list(
    phi        = list(dphi = 1, dsigma_eta2 = 0, dsigma_eps2 = 0, da1 = 0, dP1 = 0),
    a1         = list(dphi = 0, dsigma_eta2 = 0, dsigma_eps2 = 0, da1 = 1, dP1 = 0),
    P1         = list(dphi = 0, dsigma_eta2 = 0, dsigma_eps2 = 0, da1 = 0, dP1 = 1),
    sigma_eta2 = list(dphi = 0, dsigma_eta2 = 1, dsigma_eps2 = 0, da1 = 0, dP1 = 0),
    sigma_eps2 = list(dphi = 0, dsigma_eta2 = 0, dsigma_eps2 = 1, da1 = 0, dP1 = 0)
  )

  X <- matrix(0, 5, 5, dimnames = list(pn, pn))
  for (l in pn) {
    s <- seeds[[l]]
    sens <- kalman_sensitivity(kf, mu, sigma_eta2, sigma_eps2, y,
                                dphi = s$dphi, dsigma_eta2 = s$dsigma_eta2,
                                dsigma_eps2 = s$dsigma_eps2, da1 = s$da1, dP1 = s$dP1)
    dCs <- compute_C_dot(kf$alpha_smooth, kf$P_smooth, kf$P_cov_smooth,
                          sens$dalpha_smooth, sens$dP_smooth, sens$dP_cov_smooth)
    dC0 <- sum(dCs$dC0); dC1 <- sum(dCs$dC1); dC2 <- sum(dCs$dC2)
    dR <- sum(sens$dP_smooth - 2 * (y - mu - kf$alpha_smooth) * sens$dalpha_smooth)
    dS <- sum(sens$dP_smooth[, 1] + 2 * (kf$alpha_smooth[, 1] - a1) * sens$dalpha_smooth[, 1])
    dSa1 <- sum(sens$dalpha_smooth[, 1])

    X["phi", l]        <- -(1 / sigma_eta2) * (dC1 - phi * dC2)
    X["a1", l]         <- -(1 / P1) * dSa1
    X["P1", l]         <- -(1 / (2 * P1^2)) * dS
    X["sigma_eta2", l] <- -(1 / (2 * sigma_eta2^2)) * (dC0 - 2 * phi * dC1 + phi^2 * dC2)
    X["sigma_eps2", l] <- -(1 / (2 * sigma_eps2^2)) * dR
  }

  I_obs <- D + X
  vcov <- solve(I_obs)
  list(I_obs = I_obs, D = D, X = X, vcov = vcov, se = setNames(sqrt(diag(vcov)), pn))
}

compute_SE_regression <- function(y, x, mu, fit) {
  N <- nrow(y); Tt <- ncol(y)
  beta0 <- fit$beta0; beta1 <- fit$beta1; a1 <- fit$a1; P1 <- fit$P1
  sigma_eta2 <- fit$sigma_eta2; sigma_eps2 <- fit$sigma_eps2
  phi <- tanh(beta0 + beta1 * x)
  kf <- fit$kf

  Cs <- compute_C(kf$alpha_smooth, kf$P_smooth, kf$P_cov_smooth)  # per-subject
  rho <- Cs$C1 - phi * Cs$C2
  w <- (1 - phi^2) * ((1 - phi^2) * Cs$C2 + 2 * phi * rho)

  pn <- c("beta0", "beta1", "a1", "P1", "sigma_eta2", "sigma_eps2")
  D <- matrix(0, 6, 6, dimnames = list(pn, pn))
  D["beta0", "beta0"] <- sum(w) / sigma_eta2
  D["beta0", "beta1"] <- D["beta1", "beta0"] <- sum(x * w) / sigma_eta2
  D["beta1", "beta1"] <- sum(x^2 * w) / sigma_eta2
  D["a1", "a1"] <- N / P1
  D["P1", "P1"] <- N / (2 * P1^2)
  D["sigma_eta2", "sigma_eta2"] <- N * (Tt - 1) / (2 * sigma_eta2^2)
  D["sigma_eps2", "sigma_eps2"] <- N * Tt / (2 * sigma_eps2^2)

  seeds <- list(
    beta0      = list(dphi = (1 - phi^2),       dsigma_eta2 = 0, dsigma_eps2 = 0, da1 = 0, dP1 = 0),
    beta1      = list(dphi = (1 - phi^2) * x,   dsigma_eta2 = 0, dsigma_eps2 = 0, da1 = 0, dP1 = 0),
    a1         = list(dphi = rep(0, N),         dsigma_eta2 = 0, dsigma_eps2 = 0, da1 = 1, dP1 = 0),
    P1         = list(dphi = rep(0, N),         dsigma_eta2 = 0, dsigma_eps2 = 0, da1 = 0, dP1 = 1),
    sigma_eta2 = list(dphi = rep(0, N),         dsigma_eta2 = 1, dsigma_eps2 = 0, da1 = 0, dP1 = 0),
    sigma_eps2 = list(dphi = rep(0, N),         dsigma_eta2 = 0, dsigma_eps2 = 1, da1 = 0, dP1 = 0)
  )

  X <- matrix(0, 6, 6, dimnames = list(pn, pn))
  for (l in pn) {
    s <- seeds[[l]]
    sens <- kalman_sensitivity(kf, mu, sigma_eta2, sigma_eps2, y,
                                dphi = s$dphi, dsigma_eta2 = s$dsigma_eta2,
                                dsigma_eps2 = s$dsigma_eps2, da1 = s$da1, dP1 = s$dP1)
    dCs <- compute_C_dot(kf$alpha_smooth, kf$P_smooth, kf$P_cov_smooth,
                          sens$dalpha_smooth, sens$dP_smooth, sens$dP_cov_smooth)  # per-subject
    drho <- dCs$dC1 - phi * dCs$dC2
    dR <- sum(sens$dP_smooth - 2 * (y - mu - kf$alpha_smooth) * sens$dalpha_smooth)
    dSa1 <- sum(sens$dalpha_smooth[, 1])
    dS <- sum(sens$dP_smooth[, 1] + 2 * (kf$alpha_smooth[, 1] - a1) * sens$dalpha_smooth[, 1])

    X["beta0", l] <- -(1 / sigma_eta2) * sum((1 - phi^2) * drho)
    X["beta1", l] <- -(1 / sigma_eta2) * sum(x * (1 - phi^2) * drho)
    X["a1", l]    <- -(1 / P1) * dSa1
    X["P1", l]    <- -(1 / (2 * P1^2)) * dS
    X["sigma_eta2", l] <- -(1 / (2 * sigma_eta2^2)) * sum(dCs$dC0 - 2 * phi * dCs$dC1 + phi^2 * dCs$dC2)
    X["sigma_eps2", l] <- -(1 / (2 * sigma_eps2^2)) * dR
  }

  I_obs <- D + X
  vcov <- solve(I_obs)
  list(I_obs = I_obs, D = D, X = X, vcov = vcov, se = setNames(sqrt(diag(vcov)), pn))
}

## ---------------------------------------------------------------------------
## 5. Q(theta | theta*) as an explicit bivariate function, for numerical validation
##    of D and X against numDeriv (see validate_and_simulate.R).
## ---------------------------------------------------------------------------

Q_value_homogeneous <- function(theta, theta_star, y, mu) {
  N <- nrow(y); Tt <- ncol(y)
  kf <- kalman_run(y, theta_star[1], mu, theta_star[2], theta_star[3], theta_star[4], theta_star[5])
  phi <- theta[1]; a1 <- theta[2]; P1 <- theta[3]; sigma_eta2 <- theta[4]; sigma_eps2 <- theta[5]

  Q1 <- -N / 2 * log(2 * P1) - (1 / (2 * P1)) * sum(kf$P_smooth[, 1] + (kf$alpha_smooth[, 1] - a1)^2)
  Cs <- compute_C(kf$alpha_smooth, kf$P_smooth, kf$P_cov_smooth)
  C0 <- sum(Cs$C0); C1 <- sum(Cs$C1); C2 <- sum(Cs$C2)
  Q2 <- -N * (Tt - 1) / 2 * log(sigma_eta2) - (1 / (2 * sigma_eta2)) * (C0 - 2 * phi * C1 + phi^2 * C2)
  R <- sum(kf$P_smooth + (y - mu - kf$alpha_smooth)^2)
  Q3 <- -N * Tt / 2 * log(sigma_eps2) - (1 / (2 * sigma_eps2)) * R
  Q1 + Q2 + Q3
}

Q_value_regression <- function(theta, theta_star, y, x, mu) {
  N <- nrow(y); Tt <- ncol(y)
  phi_star <- tanh(theta_star[1] + theta_star[2] * x)
  kf <- kalman_run(y, phi_star, mu, theta_star[3], theta_star[4], theta_star[5], theta_star[6])

  beta0 <- theta[1]; beta1 <- theta[2]; a1 <- theta[3]; P1 <- theta[4]
  sigma_eta2 <- theta[5]; sigma_eps2 <- theta[6]
  phi <- tanh(beta0 + beta1 * x)

  Q1 <- -N / 2 * log(2 * P1) - (1 / (2 * P1)) * sum(kf$P_smooth[, 1] + (kf$alpha_smooth[, 1] - a1)^2)
  Cs <- compute_C(kf$alpha_smooth, kf$P_smooth, kf$P_cov_smooth)
  Q2 <- -N * (Tt - 1) / 2 * log(sigma_eta2) -
    (1 / (2 * sigma_eta2)) * sum(Cs$C0 - 2 * phi * Cs$C1 + phi^2 * Cs$C2)
  R <- sum(kf$P_smooth + (y - mu - kf$alpha_smooth)^2)
  Q3 <- -N * Tt / 2 * log(sigma_eps2) - (1 / (2 * sigma_eps2)) * R
  Q1 + Q2 + Q3
}

## ---------------------------------------------------------------------------
## 6. Multivariate regression model: phi_i = tanh(Xmat[i,] %*% beta), Xmat is an
##    N x p design matrix (first column conventionally an intercept of 1s; any
##    mix of continuous/binary covariates in the remaining columns). Generalizes
##    em_regression / compute_SE_regression / Q_value_regression (p = 2, fixed
##    intercept + one binary covariate) to arbitrary p. Named columns of Xmat
##    (colnames(Xmat)) become the parameter names in the returned SE/vcov.
## ---------------------------------------------------------------------------

## fix_sigma_eps2: if non-NULL, sigma_eps2 is held at this fixed value for the
## entire EM run (not re-estimated in the M-step) -- used for profiling over
## externally-fixed measurement-error variances when the plain MLE for
## sigma_eps2 piles up at (or very near) the zero boundary, a well-documented
## identifiability weakness of local-level/AR(1)+noise models in short,
## smooth series (Harvey 1989). Everything else (beta, a1, P1, sigma_eta2) is
## still estimated by EM conditional on the fixed sigma_eps2.
## phi_cap: safeguard on the M-step for theta_phi. The profiled Q-function is a
## quadratic in phi_i whose unconstrained maximizer C1/C2 can exceed 1 early in
## EM (poor start, large initial transient); BFGS on the tanh scale then runs
## beta off to +/-Inf, phi_i = 1 exactly, and the (1 - phi_i^2) factor in the
## score makes the gradient vanish identically, freezing EM at a point whose
## log-likelihood can be far below the MLE (diagnosed 2026-09-25, Part K). The
## safeguard backtracks the M-step along beta_old -> beta_new until
## max|phi_i| <= phi_cap and Q has not decreased (a generalized EM step).

## ---- Constrained M-step for theta_phi (added 2026-09-30) ----------------------------------------------
## The cap |phi_i| <= phi_cap is the set of LINEAR constraints |W_i' beta| <= atanh(phi_cap). When the
## unconstrained BFGS minimizer of the profiled Q violates it, the M-step is solved exactly under the
## constraints by an adaptive log-barrier method (constrOptim, BFGS inner iterations, analytic gradient).
## Earlier versions halved the proposed step until feasible and, when the step became tiny, kept the old
## value; the small parameter change was then mistaken for convergence ("stall"). The M-step status is
## recorded: "interior", "constrained" (solution on the cap), "stalled" (barrier solve failed and the
## fallback did not improve Q), or "iteration_limit".
constrained_mstep_beta <- function(beta, negQ, gradQ, Xmat, phi_cap, ...) {
  bound <- atanh(phi_cap); ui <- rbind(-Xmat, Xmat); ci <- rep(-bound, 2 * nrow(Xmat))
  start <- beta; eta <- as.vector(Xmat %*% start); m <- max(abs(eta))
  if (m >= bound * (1 - 1e-7)) start <- start * (bound * (1 - 1e-6) / m)   # strictly feasible start
  res <- tryCatch(constrOptim(start, function(b) negQ(b, Xmat, ...), function(b) gradQ(b, Xmat, ...), ui = ui, ci = ci,
                              method = "BFGS", outer.iterations = 200, outer.eps = 1e-10, control = list(reltol = 1e-12, maxit = 500)),
                  error = function(e) NULL)
  if (is.null(res) || !all(is.finite(res$par))) return(list(beta = beta, status = "stalled"))
  bnew <- res$par; eta <- as.vector(Xmat %*% bnew); m <- max(abs(eta)); if (m > bound) bnew <- bnew * (bound / m)
  q_new <- negQ(bnew, Xmat, ...); q_old <- negQ(beta, Xmat, ...)
  if (q_new > q_old - 1e-10 * abs(q_old)) bnew <- beta     # improvement below tolerance (or none): keep the current feasible value so EM can declare convergence
  list(beta = bnew, status = if (max(abs(as.vector(Xmat %*% bnew))) >= bound * (1 - 1e-6)) "constrained" else "interior")
}
gradQ_beta_profiled <- function(beta, Xmat, alpha_smooth, P_smooth, P_cov_smooth) {
  phi_i <- tanh(as.vector(Xmat %*% beta)); Cs <- compute_C(alpha_smooth, P_smooth, P_cov_smooth)
  total <- sum(Cs$C0 - 2 * phi_i * Cs$C1 + phi_i^2 * Cs$C2); n <- length(phi_i) * (ncol(alpha_smooth) - 1)
  dtot <- colSums(Xmat * ((-2 * Cs$C1 + 2 * phi_i * Cs$C2) * (1 - phi_i^2)))
  (n / 2) * dtot / total
}

em_regression_mv <- function(y, Xmat, mu = 0, init, thres = 1e-8, max_iter = 500,
                              fix_sigma_eps2 = NULL, phi_cap = 0.995) {
  N <- nrow(y); Tt <- ncol(y)
  beta <- init$beta; a1 <- init$a1; P1 <- init$P1
  n_backtrack <- 0L; delta <- Inf; mstep_status <- "interior"
  sigma_eta2 <- init$sigma_eta2
  sigma_eps2 <- if (!is.null(fix_sigma_eps2)) fix_sigma_eps2 else init$sigma_eps2

  neg_Q_beta_profiled <- function(beta, Xmat, alpha_smooth, P_smooth, P_cov_smooth) {
    phi_i <- tanh(as.vector(Xmat %*% beta))
    Cs <- compute_C(alpha_smooth, P_smooth, P_cov_smooth)
    total <- sum(Cs$C0 - 2 * phi_i * Cs$C1 + phi_i^2 * Cs$C2)
    sigma_eta2_hat <- total / (N * (Tt - 1))
    val <- (N * (Tt - 1) / 2) * log(sigma_eta2_hat) + total / (2 * sigma_eta2_hat)
    if (!is.finite(val)) 1e300 else val
  }

  for (iter in 1:max_iter) {
    phi <- tanh(as.vector(Xmat %*% beta))
    kf <- kalman_run(y, phi, mu, a1, P1, sigma_eta2, sigma_eps2)

    a1_new <- mean(kf$alpha_smooth[, 1])
    P1_new <- mean(kf$P_smooth[, 1] + (kf$alpha_smooth[, 1] - a1_new)^2)
    if (is.null(fix_sigma_eps2)) {
      R <- sum(kf$P_smooth + (y - mu - kf$alpha_smooth)^2)
      sigma_eps2_new <- max(R / (N * Tt), 1e-6)   # floor: variance pile-up can drive this to 0
    } else {
      sigma_eps2_new <- fix_sigma_eps2
    }
    if (!all(is.finite(c(kf$alpha_smooth, kf$P_smooth, kf$P_cov_smooth)))) { warning("non-finite smoother output; stopping EM"); break }

    opt <- optim(beta, neg_Q_beta_profiled, Xmat = Xmat,
                 alpha_smooth = kf$alpha_smooth, P_smooth = kf$P_smooth,
                 P_cov_smooth = kf$P_cov_smooth, method = "BFGS",
                 control = list(reltol = 1e-12, maxit = 500))
    beta_new <- opt$par; if (!all(is.finite(beta_new))) beta_new <- beta
    phi_new <- tanh(as.vector(Xmat %*% beta_new))
    mstep_status <- "interior"
    if (max(abs(phi_new)) > phi_cap) {
      n_backtrack <- n_backtrack + 1L
      cm <- constrained_mstep_beta(beta, neg_Q_beta_profiled, gradQ_beta_profiled, Xmat, phi_cap,
                                   alpha_smooth = kf$alpha_smooth, P_smooth = kf$P_smooth, P_cov_smooth = kf$P_cov_smooth)
      beta_new <- cm$beta; phi_new <- tanh(as.vector(Xmat %*% beta_new)); mstep_status <- cm$status
    }
    Cs <- compute_C(kf$alpha_smooth, kf$P_smooth, kf$P_cov_smooth)
    total <- sum(Cs$C0 - 2 * phi_new * Cs$C1 + phi_new^2 * Cs$C2)
    sigma_eta2_new <- total / (N * (Tt - 1))

    delta <- max(abs(c(beta_new - beta, a1_new - a1, P1_new - P1,
                        sigma_eta2_new - sigma_eta2, sigma_eps2_new - sigma_eps2)))
    beta <- beta_new; a1 <- a1_new; P1 <- P1_new
    sigma_eta2 <- sigma_eta2_new; sigma_eps2 <- sigma_eps2_new
    if (delta < thres) break
  }

  phi <- tanh(as.vector(Xmat %*% beta))
  kf_final <- kalman_run(y, phi, mu, a1, P1, sigma_eta2, sigma_eps2)
  Fv <- kf_final$P_pred + sigma_eps2
  e <- y - mu - kf_final$alpha_pred
  loglik <- -0.5 * sum(log(2 * pi * Fv) + e^2 / Fv)

  converged <- delta < thres
  status <- if (!converged) "iteration_limit" else if (mstep_status == "stalled") "stalled" else if (max(abs(phi)) >= phi_cap - 1e-6) "constrained" else "interior"
  list(beta = beta, a1 = a1, P1 = P1, sigma_eta2 = sigma_eta2, sigma_eps2 = sigma_eps2,
       iter = iter, converged = converged && mstep_status != "stalled", loglik = loglik,
       sigma_eps2_fixed = !is.null(fix_sigma_eps2), kf = kf_final,
       n_backtrack = n_backtrack, at_cap = max(abs(phi)) >= phi_cap - 1e-6, status = status)
}

compute_SE_regression_mv <- function(y, Xmat, mu, fit) {
  N <- nrow(y); Tt <- ncol(y); p <- ncol(Xmat)
  beta <- fit$beta; a1 <- fit$a1; P1 <- fit$P1
  sigma_eta2 <- fit$sigma_eta2; sigma_eps2 <- fit$sigma_eps2
  phi <- tanh(as.vector(Xmat %*% beta))
  kf <- fit$kf

  Cs <- compute_C(kf$alpha_smooth, kf$P_smooth, kf$P_cov_smooth)  # per-subject
  rho <- Cs$C1 - phi * Cs$C2
  w <- (1 - phi^2) * ((1 - phi^2) * Cs$C2 + 2 * phi * rho)

  beta_names <- if (!is.null(colnames(Xmat))) colnames(Xmat) else paste0("beta", seq_len(p) - 1)
  pn <- c(beta_names, "a1", "P1", "sigma_eta2", "sigma_eps2")
  np <- length(pn)
  D <- matrix(0, np, np, dimnames = list(pn, pn))
  D[beta_names, beta_names] <- t(Xmat) %*% (Xmat * w) / sigma_eta2
  D["a1", "a1"] <- N / P1
  D["P1", "P1"] <- N / (2 * P1^2)
  D["sigma_eta2", "sigma_eta2"] <- N * (Tt - 1) / (2 * sigma_eta2^2)
  D["sigma_eps2", "sigma_eps2"] <- N * Tt / (2 * sigma_eps2^2)

  seed_list <- list()
  for (j in seq_len(p)) {
    seed_list[[beta_names[j]]] <- list(dphi = (1 - phi^2) * Xmat[, j],
                                        dsigma_eta2 = 0, dsigma_eps2 = 0, da1 = 0, dP1 = 0)
  }
  seed_list[["a1"]] <- list(dphi = rep(0, N), dsigma_eta2 = 0, dsigma_eps2 = 0, da1 = 1, dP1 = 0)
  seed_list[["P1"]] <- list(dphi = rep(0, N), dsigma_eta2 = 0, dsigma_eps2 = 0, da1 = 0, dP1 = 1)
  seed_list[["sigma_eta2"]] <- list(dphi = rep(0, N), dsigma_eta2 = 1, dsigma_eps2 = 0, da1 = 0, dP1 = 0)
  seed_list[["sigma_eps2"]] <- list(dphi = rep(0, N), dsigma_eta2 = 0, dsigma_eps2 = 1, da1 = 0, dP1 = 0)

  X <- matrix(0, np, np, dimnames = list(pn, pn))
  for (l in pn) {
    s <- seed_list[[l]]
    sens <- kalman_sensitivity(kf, mu, sigma_eta2, sigma_eps2, y,
                                dphi = s$dphi, dsigma_eta2 = s$dsigma_eta2,
                                dsigma_eps2 = s$dsigma_eps2, da1 = s$da1, dP1 = s$dP1)
    dCs <- compute_C_dot(kf$alpha_smooth, kf$P_smooth, kf$P_cov_smooth,
                          sens$dalpha_smooth, sens$dP_smooth, sens$dP_cov_smooth)  # per-subject
    drho <- dCs$dC1 - phi * dCs$dC2
    dR <- sum(sens$dP_smooth - 2 * (y - mu - kf$alpha_smooth) * sens$dalpha_smooth)
    dSa1 <- sum(sens$dalpha_smooth[, 1])
    dS <- sum(sens$dP_smooth[, 1] + 2 * (kf$alpha_smooth[, 1] - a1) * sens$dalpha_smooth[, 1])

    for (j in seq_len(p)) {
      X[beta_names[j], l] <- -(1 / sigma_eta2) * sum(Xmat[, j] * (1 - phi^2) * drho)
    }
    X["a1", l]    <- -(1 / P1) * dSa1
    X["P1", l]    <- -(1 / (2 * P1^2)) * dS
    X["sigma_eta2", l] <- -(1 / (2 * sigma_eta2^2)) * sum(dCs$dC0 - 2 * phi * dCs$dC1 + phi^2 * dCs$dC2)
    X["sigma_eps2", l] <- -(1 / (2 * sigma_eps2^2)) * dR
  }

  I_obs <- D + X
  vcov <- solve(I_obs)
  list(I_obs = I_obs, D = D, X = X, vcov = vcov, se = setNames(sqrt(diag(vcov)), pn))
}

## Oakes information for fixed sigma_eps2: theta = (beta, a1, P1, sigma_eta2).
## The direct block below is unsimplified, so it is valid away from a fixed
## point and is the exact sum of compute_info_per_unit_fixed_eps2(). A fixed
## measurement variance has no parameter row/column; varying that assumption
## requires a separate sensitivity analysis. See SE_derivation/PER_UNIT_INFORMATION.md.
compute_SE_regression_mv_fixed_eps2 <- function(y, Xmat, mu, fit) {
  N <- nrow(y); Tt <- ncol(y); p <- ncol(Xmat)
  beta <- fit$beta; a1 <- fit$a1; P1 <- fit$P1
  sigma_eta2 <- fit$sigma_eta2; sigma_eps2 <- fit$sigma_eps2
  phi <- tanh(as.vector(Xmat %*% beta))
  kf <- fit$kf

  Cs <- compute_C(kf$alpha_smooth, kf$P_smooth, kf$P_cov_smooth)
  rho <- Cs$C1 - phi * Cs$C2
  w <- (1 - phi^2) * ((1 - phi^2) * Cs$C2 + 2 * phi * rho)

  beta_names <- if (!is.null(colnames(Xmat))) colnames(Xmat) else paste0("beta", seq_len(p) - 1)
  pn <- c(beta_names, "a1", "P1", "sigma_eta2")  # sigma_eps2 excluded: fixed, not estimated
  np <- length(pn)
  D <- matrix(0, np, np, dimnames = list(pn, pn))
  D[beta_names, beta_names] <- t(Xmat) %*% (Xmat * w) / sigma_eta2
  D["a1", "a1"] <- N / P1
  S_initial <- kf$P_smooth[, 1] + (kf$alpha_smooth[, 1] - a1)^2
  R_transition <- Cs$C0 - 2 * phi * Cs$C1 + phi^2 * Cs$C2
  D["P1", "P1"] <- -N / (2 * P1^2) + sum(S_initial) / P1^3
  D["sigma_eta2", "sigma_eta2"] <- -N * (Tt - 1) / (2 * sigma_eta2^2) +
    sum(R_transition) / sigma_eta2^3
  D["a1", "P1"] <- D["P1", "a1"] <- sum(kf$alpha_smooth[, 1] - a1) / P1^2
  beta_q <- as.vector(crossprod(Xmat, (1 - phi^2) * rho)) / sigma_eta2^2
  D[beta_names, "sigma_eta2"] <- beta_q
  D["sigma_eta2", beta_names] <- beta_q

  seed_list <- list()
  for (j in seq_len(p)) {
    seed_list[[beta_names[j]]] <- list(dphi = (1 - phi^2) * Xmat[, j],
                                        dsigma_eta2 = 0, dsigma_eps2 = 0, da1 = 0, dP1 = 0)
  }
  seed_list[["a1"]] <- list(dphi = rep(0, N), dsigma_eta2 = 0, dsigma_eps2 = 0, da1 = 1, dP1 = 0)
  seed_list[["P1"]] <- list(dphi = rep(0, N), dsigma_eta2 = 0, dsigma_eps2 = 0, da1 = 0, dP1 = 1)
  seed_list[["sigma_eta2"]] <- list(dphi = rep(0, N), dsigma_eta2 = 1, dsigma_eps2 = 0, da1 = 0, dP1 = 0)
  ## no "sigma_eps2" seed: it is a fixed constant, not a direction in theta*

  X <- matrix(0, np, np, dimnames = list(pn, pn))
  for (l in pn) {
    s <- seed_list[[l]]
    sens <- kalman_sensitivity(kf, mu, sigma_eta2, sigma_eps2, y,
                                dphi = s$dphi, dsigma_eta2 = s$dsigma_eta2,
                                dsigma_eps2 = s$dsigma_eps2, da1 = s$da1, dP1 = s$dP1)
    dCs <- compute_C_dot(kf$alpha_smooth, kf$P_smooth, kf$P_cov_smooth,
                          sens$dalpha_smooth, sens$dP_smooth, sens$dP_cov_smooth)
    drho <- dCs$dC1 - phi * dCs$dC2
    dSa1 <- sum(sens$dalpha_smooth[, 1])
    dS <- sum(sens$dP_smooth[, 1] + 2 * (kf$alpha_smooth[, 1] - a1) * sens$dalpha_smooth[, 1])

    for (j in seq_len(p)) {
      X[beta_names[j], l] <- -(1 / sigma_eta2) * sum(Xmat[, j] * (1 - phi^2) * drho)
    }
    X["a1", l]    <- -(1 / P1) * dSa1
    X["P1", l]    <- -(1 / (2 * P1^2)) * dS
    X["sigma_eta2", l] <- -(1 / (2 * sigma_eta2^2)) * sum(dCs$dC0 - 2 * phi * dCs$dC1 + phi^2 * dCs$dC2)
  }

  I_obs <- D + X
  vcov <- solve(I_obs)
  list(I_obs = I_obs, D = D, X = X, vcov = vcov, se = setNames(sqrt(diag(vcov)), pn))
}

## Cluster-robust ("sandwich") correction on top of compute_SE_regression_mv_
## fixed_eps2's Oakes-identity information matrix, for when observation units
## (e.g. meal-response windows) are grouped into clusters (e.g. subjects) that
## share covariates and plausibly correlated residual variation the model
## doesn't explicitly represent. Uses the standard CR0 sandwich: bread =
## I_obs^{-1} (already-validated observed information), meat = sum over
## clusters of (summed within-cluster score) outer-product. Per-unit scores
## are the per-i contributions to the same Fisher's-identity gradient that
## feeds D's diagonal entries in compute_SE_regression_mv_fixed_eps2 -- each
## is exactly that function's per-i term, simply left unsummed over i.
## `cluster` must be a length-N vector/factor of cluster ids aligned to y's
## rows. Singleton clusters give the ordinary unit-level sandwich. Under a
## correctly specified independent-unit model it agrees asymptotically with
## model-based uncertainty, not exactly in a finite sample.
compute_SE_regression_mv_cluster_robust <- function(y, Xmat, mu, fit, cluster) {
  N <- nrow(y); Tt <- ncol(y); p <- ncol(Xmat)
  beta <- fit$beta; a1 <- fit$a1; P1 <- fit$P1
  sigma_eta2 <- fit$sigma_eta2; sigma_eps2 <- fit$sigma_eps2
  phi <- tanh(as.vector(Xmat %*% beta))
  kf <- fit$kf
  stopifnot(length(cluster) == N)

  base <- compute_SE_regression_mv_fixed_eps2(y, Xmat, mu, fit)
  pn <- rownames(base$I_obs)
  beta_names <- pn[seq_len(p)]

  Cs <- compute_C(kf$alpha_smooth, kf$P_smooth, kf$P_cov_smooth)  # per-i: C0,C1,C2
  rho <- Cs$C1 - phi * Cs$C2

  ## per-window (per-i) score contributions, each a length-N vector
  score <- matrix(0, N, length(pn), dimnames = list(NULL, pn))
  for (j in seq_len(p)) {
    score[, beta_names[j]] <- (1 / sigma_eta2) * Xmat[, j] * (1 - phi^2) * rho
  }
  score[, "a1"] <- (1 / P1) * (kf$alpha_smooth[, 1] - a1)
  score[, "P1"] <- -1 / (2 * P1) + (1 / (2 * P1^2)) * (kf$P_smooth[, 1] + (kf$alpha_smooth[, 1] - a1)^2)
  score[, "sigma_eta2"] <- -(Tt - 1) / (2 * sigma_eta2) +
    (1 / (2 * sigma_eta2^2)) * (Cs$C0 - 2 * phi * Cs$C1 + phi^2 * Cs$C2)

  ## sanity check: the per-i scores should sum (over ALL i) to ~0 at the MLE
  total_score <- colSums(score)

  clusters <- unique(cluster)
  meat <- matrix(0, length(pn), length(pn), dimnames = list(pn, pn))
  for (g in clusters) {
    Sg <- colSums(score[cluster == g, , drop = FALSE])
    meat <- meat + outer(Sg, Sg)
  }

  bread <- base$vcov  # = I_obs^{-1}
  vcov_robust <- bread %*% meat %*% bread
  se_robust <- sqrt(diag(vcov_robust))

  list(vcov = vcov_robust, se = setNames(se_robust, pn), meat = meat, bread = bread,
       total_score_check = total_score, n_clusters = length(clusters),
       se_model_based = base$se, score = score)
}

## Exact nuisance reduction of a JOINT cluster sandwich, in any parameter order.
## I_obs and score must use the same named parameter coordinates. No small-sample
## adjustment or reference distribution is implied. See NUISANCE_AWARE_INFERENCE.md.
reduce_regression_cluster_scores <- function(I_obs, score, cluster, beta_names,
                                             unit_info = NULL, tolerance = 1e-8) {
  pn <- colnames(I_obs)
  if (!is.matrix(I_obs) || nrow(I_obs) != ncol(I_obs) || is.null(pn) ||
      anyDuplicated(pn) || !identical(rownames(I_obs), pn) ||
      !is.matrix(score) || !identical(colnames(score), pn) ||
      any(!is.finite(I_obs)) || any(!is.finite(score)))
    stop("Information and scores must be finite matrices in identical named coordinates.")
  if (!length(beta_names) || anyDuplicated(beta_names) || !all(beta_names %in% pn))
    stop("beta_names must select distinct information parameters.")
  if (length(cluster) != nrow(score) || anyNA(cluster) || length(unique(cluster)) < 2L)
    stop("At least two nonmissing clusters aligned with score rows are required.")
  if (length(tolerance) != 1L || !is.finite(tolerance) || tolerance <= 0)
    stop("tolerance must be a positive finite scalar.")
  # Measure symmetry in coordinates scaled to unit diagonal before rounding.
  if (any(diag(I_obs) <= 0)) stop("Joint information is not positive definite.")
  information_scale <- sqrt(diag(I_obs))
  standardized <- I_obs / outer(information_scale, information_scale)
  asymmetry <- max(abs(standardized - t(standardized)))
  if (asymmetry > tolerance) stop("Joint information is not symmetric to tolerance.")
  H <- (I_obs + t(I_obs)) / 2
  H_scaled <- H / outer(information_scale, information_scale)
  if (inherits(try(chol(H_scaled), silent = TRUE), "try-error"))
    stop("Joint information is not positive definite; do not repair it by clipping eigenvalues.")
  nuisance <- setdiff(pn, beta_names); p <- length(beta_names)
  projection <- matrix(0, p, length(pn), dimnames = list(beta_names, pn))
  projection[, beta_names] <- diag(p)
  if (length(nuisance)) {
    C <- t(solve(H[nuisance, nuisance, drop = FALSE], H[nuisance, beta_names, drop = FALSE]))
    projection[, nuisance] <- -C
    K <- H[beta_names, beta_names, drop = FALSE] - C %*% H[nuisance, beta_names, drop = FALSE]
  } else {
    C <- matrix(0, p, 0, dimnames = list(beta_names, character()))
    K <- H[beta_names, beta_names, drop = FALSE]
  }
  K <- (K + t(K)) / 2
  K_inv <- solve(K)
  adjusted_score <- score %*% t(projection)
  cluster_score <- rowsum(adjusted_score, group = as.character(cluster), reorder = FALSE)
  influence <- cluster_score %*% t(K_inv)
  vc <- crossprod(influence)
  # Independent full-block calculation is retained as an exact algebraic check.
  joint_cluster_score <- rowsum(score, group = as.character(cluster), reorder = FALSE)
  joint_bread <- solve(H)
  joint_influence <- joint_cluster_score %*% t(joint_bread)
  full_beta_vcov <- crossprod(joint_influence[, beta_names, drop = FALSE])
  projected_unit <- projected_cluster <- NULL
  projection_sum_error <- NA_real_
  if (!is.null(unit_info)) {
    if (!is.array(unit_info) || !identical(dim(unit_info), c(nrow(score), length(pn), length(pn))) ||
        !identical(dimnames(unit_info)[[2]], pn) || !identical(dimnames(unit_info)[[3]], pn) ||
        any(!is.finite(unit_info))) stop("unit_info must align with score rows and joint parameters.")
    if (max(abs(apply(unit_info, c(2, 3), sum) - H) / (1 + abs(H))) > tolerance)
      stop("Per-unit information does not sum to joint information.")
    projected_unit <- array(0, c(nrow(score), p, p), dimnames = list(NULL, beta_names, beta_names))
    for (i in seq_len(nrow(score))) {
      Hi <- matrix(unit_info[i, , ], length(pn), length(pn))
      projected_unit[i, , ] <- projection %*% Hi %*% t(projection)
    }
    groups <- rownames(cluster_score)
    projected_cluster <- array(0, c(length(groups), p, p), dimnames = list(groups, beta_names, beta_names))
    for (g in seq_along(groups)) projected_cluster[g, , ] <-
      apply(projected_unit[as.character(cluster) == groups[g], , , drop = FALSE], c(2, 3), sum)
    projection_sum_error <- max(abs(apply(projected_unit, c(2, 3), sum) - K))
  }
  list(vcov = vc, se = setNames(sqrt(diag(vc)), beta_names),
       vcov_model_based = K_inv, se_model_based = setNames(sqrt(diag(K_inv)), beta_names),
       I_eff = K, projection = projection, nuisance_coefficient = C,
       beta_names = beta_names, nuisance_names = nuisance,
       score = adjusted_score, cluster_score = cluster_score, influence = influence,
       joint_cluster_score = joint_cluster_score,
       projected_info_per_unit = projected_unit, projected_info_per_cluster = projected_cluster,
       n_clusters = nrow(cluster_score),
       diagnostics = list(full_sandwich_error = max(abs(vc - full_beta_vcov)),
         full_influence_error = max(abs(influence - joint_influence[, beta_names, drop = FALSE])),
         full_model_block_error = max(abs(K_inv - joint_bread[beta_names, beta_names, drop = FALSE])),
         information_asymmetry_scaled = asymmetry,
         information_rcond_scaled = rcond(H_scaled),
         projected_information_sum_error = projection_sum_error,
         total_adjusted_score = colSums(adjusted_score),
         joint_newton_step = as.vector(joint_bread %*% colSums(score))),
       method = "Nuisance-adjusted cluster CR0 (exact full-sandwich beta block)",
       small_sample_adjustment = "none")
}

## Primary fixed-measurement-variance regression wrapper. All estimated a1, P1,
## sigma_eta2 parameters enter the nuisance adjustment. mu is known/fixed here.
compute_SE_regression_mv_cluster_beta <- function(y, Xmat, mu, fit, cluster,
                                                  return_unit_info = FALSE) {
  if (anyNA(y) || any(!is.finite(y)) || ncol(y) < 2L || nrow(Xmat) != nrow(y))
    stop("Complete finite observations (T >= 2) and aligned covariates are required.")
  if (!isTRUE(fit$sigma_eps2_fixed))
    stop("fit must be marked sigma_eps2_fixed=TRUE; a free-variance fit needs its full nuisance block.")
  if (!isTRUE(fit$converged) || isTRUE(fit$at_cap))
    stop("Inference requires a converged interior fit; inspect nonconvergence or an active cap first.")
  cr0 <- compute_SE_regression_mv_cluster_robust(y, Xmat, mu, fit, cluster)
  bn <- colnames(cr0$score)[seq_len(ncol(Xmat))]
  H <- solve(cr0$bread)
  Ii <- if (return_unit_info) compute_info_per_unit_fixed_eps2(y, Xmat, mu, fit) else NULL
  ans <- reduce_regression_cluster_scores(H, cr0$score, cluster, bn, Ii)
  ans$beta <- setNames(as.vector(fit$beta), bn)
  ans$diagnostics$existing_full_CR0_error <- max(abs(ans$vcov - cr0$vcov[bn, bn, drop = FALSE]))
  ans$diagnostics$total_joint_score <- cr0$total_score_check
  ans
}

## Explicit reference distribution for a single linear contrast. No automatic
## t_(G-1) or CR2/Satterthwaite claim is made. t requires a user-supplied df.
wald_regression_cluster_contrast <- function(object, contrast, null = 0,
                                             level = .95, reference = c("normal", "t"), df = NULL) {
  reference <- match.arg(reference)
  bn <- names(object$beta)
  if (!is.numeric(contrast) || !length(contrast) || any(!is.finite(contrast)) ||
      !any(contrast != 0)) stop("contrast must be a nonzero finite numeric vector.")
  if (!is.null(names(contrast))) {
    if (anyDuplicated(names(contrast)) || !all(names(contrast) %in% bn)) stop("Unknown or duplicated contrast names.")
    L <- setNames(numeric(length(bn)), bn); L[names(contrast)] <- contrast
  } else {
    if (length(contrast) != length(bn)) stop("Unnamed contrast must match all regression parameters.")
    L <- contrast
  }
  if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1 ||
      length(null) != 1L || !is.finite(null)) stop("Invalid level or null value.")
  estimate <- sum(L * object$beta); variance <- as.numeric(t(L) %*% object$vcov %*% L)
  if (!is.finite(variance) || variance <= 0) stop("Contrast variance must be positive.")
  se <- sqrt(variance); statistic <- (estimate - null) / se
  if (reference == "t") {
    if (length(df) != 1L || !is.finite(df) || df <= 0) stop("A positive df must be explicitly supplied for a t reference.")
    critical <- qt((1 + level) / 2, df); p_value <- 2 * pt(-abs(statistic), df)
  } else {
    if (!is.null(df)) stop("df only applies to a t reference.")
    critical <- qnorm((1 + level) / 2); p_value <- 2 * pnorm(-abs(statistic)); df <- Inf
  }
  data.frame(estimate = estimate, se = se, statistic = statistic, null = null,
             lower = estimate - critical * se, upper = estimate + critical * se,
             level = level, reference = reference, df = df, p_value = p_value)
}

## Per-unit contributions to the Oakes-identity observed information (fixed
## sigma_eps2 case). Because each unit's filter/smoother recursion is
## independent of the others given theta, both the direct block D and the
## cross block X are exact sums over units of per-unit terms; this function
## returns them as an N x np x np array (used for the bias-corrected sandwich).
compute_info_per_unit_fixed_eps2 <- function(y, Xmat, mu, fit) {
  N <- nrow(y); Tt <- ncol(y); p <- ncol(Xmat)
  beta <- fit$beta; a1 <- fit$a1; P1 <- fit$P1
  sigma_eta2 <- fit$sigma_eta2; sigma_eps2 <- fit$sigma_eps2
  phi <- tanh(as.vector(Xmat %*% beta)); kf <- fit$kf
  Cs <- compute_C(kf$alpha_smooth, kf$P_smooth, kf$P_cov_smooth)
  rho <- Cs$C1 - phi * Cs$C2
  w <- (1 - phi^2) * ((1 - phi^2) * Cs$C2 + 2 * phi * rho)
  beta_names <- if (!is.null(colnames(Xmat))) colnames(Xmat) else paste0("beta", seq_len(p) - 1)
  pn <- c(beta_names, "a1", "P1", "sigma_eta2"); np <- length(pn)
  Ii <- array(0, c(N, np, np), dimnames = list(NULL, pn, pn))
  for (j in seq_len(p)) for (k in seq_len(p))
    Ii[, beta_names[j], beta_names[k]] <- Xmat[, j] * Xmat[, k] * w / sigma_eta2
  Ii[, "a1", "a1"] <- 1 / P1
  ## Keep per-unit residual terms: only their FULL-DATA sums vanish at an
  ## interior optimum. Never substitute the global score equations per unit.
  S_initial <- kf$P_smooth[, 1] + (kf$alpha_smooth[, 1] - a1)^2
  R_transition <- Cs$C0 - 2 * phi * Cs$C1 + phi^2 * Cs$C2
  Ii[, "P1", "P1"] <- -1 / (2 * P1^2) + S_initial / P1^3
  Ii[, "sigma_eta2", "sigma_eta2"] <- -(Tt - 1) / (2 * sigma_eta2^2) +
    R_transition / sigma_eta2^3
  Ii[, "a1", "P1"] <- Ii[, "P1", "a1"] <- (kf$alpha_smooth[, 1] - a1) / P1^2
  for (j in seq_len(p)) {
    beta_q <- Xmat[, j] * (1 - phi^2) * rho / sigma_eta2^2
    Ii[, beta_names[j], "sigma_eta2"] <- beta_q
    Ii[, "sigma_eta2", beta_names[j]] <- beta_q
  }
  seed_list <- list()
  for (j in seq_len(p)) seed_list[[beta_names[j]]] <- list(dphi = (1 - phi^2) * Xmat[, j], dsigma_eta2 = 0, dsigma_eps2 = 0, da1 = 0, dP1 = 0)
  seed_list[["a1"]] <- list(dphi = rep(0, N), dsigma_eta2 = 0, dsigma_eps2 = 0, da1 = 1, dP1 = 0)
  seed_list[["P1"]] <- list(dphi = rep(0, N), dsigma_eta2 = 0, dsigma_eps2 = 0, da1 = 0, dP1 = 1)
  seed_list[["sigma_eta2"]] <- list(dphi = rep(0, N), dsigma_eta2 = 1, dsigma_eps2 = 0, da1 = 0, dP1 = 0)
  for (l in pn) {
    s <- seed_list[[l]]
    sens <- kalman_sensitivity(kf, mu, sigma_eta2, sigma_eps2, y, dphi = s$dphi, dsigma_eta2 = s$dsigma_eta2,
                                dsigma_eps2 = s$dsigma_eps2, da1 = s$da1, dP1 = s$dP1)
    dCs <- compute_C_dot(kf$alpha_smooth, kf$P_smooth, kf$P_cov_smooth, sens$dalpha_smooth, sens$dP_smooth, sens$dP_cov_smooth)
    drho <- dCs$dC1 - phi * dCs$dC2
    for (j in seq_len(p)) Ii[, beta_names[j], l] <- Ii[, beta_names[j], l] - (1 / sigma_eta2) * Xmat[, j] * (1 - phi^2) * drho
    Ii[, "a1", l] <- Ii[, "a1", l] - (1 / P1) * sens$dalpha_smooth[, 1]
    Ii[, "P1", l] <- Ii[, "P1", l] - (1 / (2 * P1^2)) * (sens$dP_smooth[, 1] + 2 * (kf$alpha_smooth[, 1] - a1) * sens$dalpha_smooth[, 1])
    Ii[, "sigma_eta2", l] <- Ii[, "sigma_eta2", l] - (1 / (2 * sigma_eta2^2)) * (dCs$dC0 - 2 * phi * dCs$dC1 + phi^2 * dCs$dC2)
  }
  Ii
}

## Bias-corrected cluster-robust sandwich in the spirit of Mancl & DeRouen
## (2001, Biometrics 57:126-134): each cluster's summed score S_g is inflated
## by (I - H_g)^{-1}, H_g = I_g I^{-1}, where I_g is the cluster's contribution
## to the observed information. Reduces the downward bias of CR0 when the
## number of clusters is small. Returns CR0 and model-based SEs alongside;
## use with t_{G-1} critical values.
compute_SE_regression_mv_cluster_MD <- function(y, Xmat, mu, fit, cluster) {
  cr0 <- compute_SE_regression_mv_cluster_robust(y, Xmat, mu, fit, cluster)
  Ii <- compute_info_per_unit_fixed_eps2(y, Xmat, mu, fit)
  pn <- colnames(cr0$score); np <- length(pn)
  I_tot <- apply(Ii, c(2, 3), sum)
  I_inv <- solve(I_tot)
  clusters <- unique(cluster); G <- length(clusters)
  meat <- matrix(0, np, np, dimnames = list(pn, pn))
  for (g in clusters) {
    idx <- which(cluster == g)
    Sg <- colSums(cr0$score[idx, , drop = FALSE])
    Ig <- apply(Ii[idx, , , drop = FALSE], c(2, 3), sum)
    Ag <- solve(diag(np) - Ig %*% I_inv)
    v <- Ag %*% Sg
    meat <- meat + v %*% t(v)
  }
  vcov_md <- I_inv %*% meat %*% I_inv
  list(vcov = vcov_md, se = setNames(sqrt(diag(vcov_md)), pn), se_cr0 = cr0$se, se_model_based = cr0$se_model_based,
       n_clusters = G, I_check = max(abs(I_tot - solve(cr0$bread))))
}

## LEGACY: this beta-only adjustment omits nuisance estimation and is not a
## validated CR2 implementation. Retained for reproducing previous results.
## Use compute_SE_regression_mv_cluster_beta() for exact nuisance-aware CR0.
## A small-sample replacement requires a working covariance derivation.
## CR2-type (Bell-McCaffrey-style) correction restricted to the regression
## block: A_g = I^{1/2} (I_p - B_g)^{-1/2} I^{-1/2} with B_g = I^{-1/2} I_g I^{-1/2}
## (symmetrized per-cluster information), applied to each cluster's summed
## score. Milder than Mancl-DeRouen ((I - H_g)^{-1}); both are returned by the
## real-data scripts so that the sensitivity of the conclusions to the
## small-sample correction can be reported. `leverage` gives max eigenvalue of
## B_g per cluster (the share of the information a cluster carries in its
## most influential direction).
compute_SE_regression_mv_cluster_CR2 <- function(y, Xmat, mu, fit, cluster) {
  cr0 <- compute_SE_regression_mv_cluster_robust(y, Xmat, mu, fit, cluster)
  Ii <- compute_info_per_unit_fixed_eps2(y, Xmat, mu, fit)
  bn <- colnames(Xmat)
  I_tot <- apply(Ii, c(2, 3), sum); Ib <- (I_tot[bn, bn] + t(I_tot[bn, bn])) / 2; Iinv <- solve(Ib)
  es <- eigen(Ib, symmetric = TRUE); Ihalf <- es$vectors %*% diag(sqrt(es$values)) %*% t(es$vectors); Imhalf <- solve(Ihalf)
  clusters <- unique(cluster); meat <- 0; lev <- setNames(numeric(length(clusters)), clusters)
  for (g in clusters) {
    idx <- which(cluster == g)
    Ig <- apply(Ii[idx, bn, bn, drop = FALSE], c(2, 3), sum); Ig <- (Ig + t(Ig)) / 2
    Bg <- Imhalf %*% Ig %*% Imhalf; Bg <- (Bg + t(Bg)) / 2; eb <- eigen(Bg, symmetric = TRUE); lev[g] <- max(eb$values)
    Ag <- Ihalf %*% (eb$vectors %*% diag(1 / sqrt(pmax(1 - eb$values, 1e-8)), length(bn)) %*% t(eb$vectors)) %*% Imhalf
    v <- Ag %*% colSums(cr0$score[idx, bn, drop = FALSE]); meat <- meat + v %*% t(v)
  }
  vc <- Iinv %*% meat %*% Iinv
  list(vcov = vc, se = setNames(sqrt(diag(vc)), bn), se_cr0 = cr0$se[bn], se_model_based = cr0$se_model_based[bn],
       leverage = lev, n_clusters = length(clusters),
       method = "Legacy beta-only information adjustment", nuisance_adjusted = FALSE, validated_CR2 = FALSE)
}

Q_value_regression_mv <- function(theta, theta_star, y, Xmat, mu) {
  N <- nrow(y); Tt <- ncol(y); p <- ncol(Xmat)
  beta_star <- theta_star[1:p]
  phi_star <- tanh(as.vector(Xmat %*% beta_star))
  kf <- kalman_run(y, phi_star, mu, theta_star[p + 1], theta_star[p + 2],
                    theta_star[p + 3], theta_star[p + 4])

  beta <- theta[1:p]; a1 <- theta[p + 1]; P1 <- theta[p + 2]
  sigma_eta2 <- theta[p + 3]; sigma_eps2 <- theta[p + 4]
  phi <- tanh(as.vector(Xmat %*% beta))

  Q1 <- -N / 2 * log(2 * P1) - (1 / (2 * P1)) * sum(kf$P_smooth[, 1] + (kf$alpha_smooth[, 1] - a1)^2)
  Cs <- compute_C(kf$alpha_smooth, kf$P_smooth, kf$P_cov_smooth)
  Q2 <- -N * (Tt - 1) / 2 * log(sigma_eta2) -
    (1 / (2 * sigma_eta2)) * sum(Cs$C0 - 2 * phi * Cs$C1 + phi^2 * Cs$C2)
  R <- sum(kf$P_smooth + (y - mu - kf$alpha_smooth)^2)
  Q3 <- -N * Tt / 2 * log(sigma_eps2) - (1 / (2 * sigma_eps2)) * R
  Q1 + Q2 + Q3
}

## Naive two-stage competitor: per-subject OLS AR(1) slope (ignores measurement
## error entirely), Fisher-z transformed, regressed on Xmat via OLS. Used as the
## "what would a non-specialist do" benchmark in sim_benchmark_naive.R.
naive_two_stage <- function(y, Xmat) {
  N <- nrow(y); Tt <- ncol(y)
  phi_hat <- vapply(seq_len(N), function(i) {
    coef(lm(y[i, 2:Tt] ~ y[i, 1:(Tt - 1)]))[2]
  }, numeric(1))
  phi_hat_clipped <- pmin(pmax(phi_hat, -0.999), 0.999)
  z <- atanh(phi_hat_clipped)
  fit <- lm(z ~ Xmat - 1)
  list(phi_hat = phi_hat, z = z, beta_hat = coef(fit), se = summary(fit)$coefficients[, "Std. Error"])
}

## Second competitor: fit each subject's OWN state-space model separately by
## direct numerical MLE (Kalman-filter prediction-error-decomposition
## likelihood, optimized over an unconstrained reparameterization), then
## Fisher-z transform and regress on Xmat via OLS -- same second stage as
## naive_two_stage(), but with correctly measurement-error-corrected per-subject
## phi_hat instead of the naive AR(1) OLS slope. Isolates how much of the joint
## EM's advantage is "correcting for measurement error" (this competitor should
## fix most of naive_two_stage's attenuation) vs. "pooling variance components
## and beta across subjects" (only the joint EM does that).
single_subject_negloglik <- function(par, y, mu) {
  phi <- tanh(par[1]); a1 <- par[2]; P1 <- exp(par[3])
  sigma_eta2 <- exp(par[4]); sigma_eps2 <- exp(par[5])
  Tt <- length(y); ymat <- matrix(y, 1, Tt)
  kf <- kalman_run(ymat, phi, mu, a1, P1, sigma_eta2, sigma_eps2)
  Fv <- pmax(kf$P_pred[1, ] + sigma_eps2, 1e-8)
  e <- ymat[1, ] - mu - kf$alpha_pred[1, ]
  val <- 0.5 * sum(log(2 * pi * Fv) + e^2 / Fv)
  if (!is.finite(val)) val <- 1e10
  val
}

fit_single_subject_MLE <- function(y_row, mu = 0) {
  Tt <- length(y_row)
  phi0 <- tryCatch({
    p <- coef(lm(y_row[2:Tt] ~ y_row[1:(Tt - 1)]))[2]
    min(max(p, -0.9), 0.9)
  }, error = function(e) 0.2)
  init <- c(atanh(phi0), mean(y_row), log(var(y_row) / 2 + 0.1),
            log(var(y_row) / 2 + 0.1), log(var(y_row) / 2 + 0.1))
  opt <- optim(init, single_subject_negloglik, y = y_row, mu = mu,
               method = "Nelder-Mead", control = list(maxit = 1000, reltol = 1e-12))
  list(phi = tanh(opt$par[1]), a1 = opt$par[2], P1 = exp(opt$par[3]),
       sigma_eta2 = exp(opt$par[4]), sigma_eps2 = exp(opt$par[5]),
       convergence = opt$convergence)
}

per_subject_kalman_two_stage <- function(y, Xmat, mu = 0) {
  N <- nrow(y)
  fits <- lapply(seq_len(N), function(i) fit_single_subject_MLE(y[i, ], mu))
  phi_hat <- vapply(fits, function(f) f$phi, numeric(1))
  phi_hat_clipped <- pmin(pmax(phi_hat, -0.999), 0.999)
  z <- atanh(phi_hat_clipped)
  fit <- lm(z ~ Xmat - 1)
  list(phi_hat = phi_hat, z = z, beta_hat = coef(fit),
       se = summary(fit)$coefficients[, "Std. Error"],
       n_nonconverged = sum(vapply(fits, function(f) f$convergence, numeric(1)) != 0))
}

## Third competitor: a "reasonable applied statistician's" marginal-model
## baseline -- nlme::gls with an AR(1) residual correlation structure, fit
## SEPARATELY within each level of a categorical covariate (here: Sex). This
## has no latent-state/measurement-error decomposition at all (it models the
## OBSERVED series' own autocorrelation) and, more importantly, has no way to
## use a *continuous* covariate to explain the AR parameter -- it can only ever
## report one rho per stratum. Returns rho per stratum, for comparison against
## each other method's stratum-average phi (see sim_benchmark_extended.R).
stratified_gls_corAR1 <- function(y, strata) {
  N <- nrow(y); Tt <- ncol(y)
  long_df <- data.frame(subject = rep(seq_len(N), times = Tt),
                         time = rep(seq_len(Tt), each = N),
                         y = as.vector(y),
                         stratum = rep(strata, times = Tt))
  levels_s <- sort(unique(strata))
  rho <- vapply(levels_s, function(s) {
    sub <- long_df[long_df$stratum == s, ]
    m <- nlme::gls(y ~ 1, data = sub, correlation = nlme::corAR1(form = ~ time | subject))
    as.numeric(coef(m$modelStruct$corStruct, unconstrained = FALSE))
  }, numeric(1))
  setNames(rho, levels_s)
}

## Numerically compute I_obs = D + X via numDeriv, for cross-checking the analytic
## formulas above. `theta_hat` is the fixed point; `Qfun(theta, theta_star)` must be
## one of Q_value_homogeneous / Q_value_regression with data arguments pre-bound.
numeric_I_obs <- function(Qfun, theta_hat) {
  D_num <- -numDeriv::hessian(function(th) Qfun(th, theta_hat), theta_hat)
  grad_theta <- function(theta_star_val) numDeriv::grad(function(th) Qfun(th, theta_star_val), theta_hat)
  X_num <- -numDeriv::jacobian(grad_theta, theta_hat)
  list(D = D_num, X = X_num, I_obs = D_num + X_num)
}

## Sensitivity variant: the set-point mu_i is not fixed but modeled as
## mu_i = Zmu_i' gamma and estimated inside EM. The M-step for gamma is the
## closed-form least-squares update of rowMeans(y - alpha_smooth) on Zmu
## (the Q-function is quadratic in gamma with equal weights across t).
## Everything else is identical to em_regression_mv (including the phi cap).
## Intended for sensitivity analysis: no Oakes SE is derived for gamma, use
## a cluster bootstrap for uncertainty.
em_regression_mv_setpoint <- function(y, Xmat, Zmu, init, thres = 1e-8, max_iter = 500,
                                      fix_sigma_eps2 = NULL, phi_cap = 0.995) {
  N <- nrow(y); Tt <- ncol(y)
  beta <- init$beta; a1 <- init$a1; P1 <- init$P1
  sigma_eta2 <- init$sigma_eta2
  sigma_eps2 <- if (!is.null(fix_sigma_eps2)) fix_sigma_eps2 else init$sigma_eps2
  gamma <- if (!is.null(init$gamma)) init$gamma else rep(0, ncol(Zmu))
  ZtZinv <- solve(crossprod(Zmu))
  n_backtrack <- 0L
  neg_Q_beta_profiled <- function(beta, Xmat, alpha_smooth, P_smooth, P_cov_smooth) {
    phi_i <- tanh(as.vector(Xmat %*% beta))
    Cs <- compute_C(alpha_smooth, P_smooth, P_cov_smooth)
    total <- sum(Cs$C0 - 2 * phi_i * Cs$C1 + phi_i^2 * Cs$C2)
    sigma_eta2_hat <- total / (N * (Tt - 1))
    val <- (N * (Tt - 1) / 2) * log(sigma_eta2_hat) + total / (2 * sigma_eta2_hat)
    if (!is.finite(val)) 1e300 else val
  }
  for (iter in 1:max_iter) {
    phi <- tanh(as.vector(Xmat %*% beta)); mu <- as.vector(Zmu %*% gamma)
    kf <- kalman_run(y, phi, mu, a1, P1, sigma_eta2, sigma_eps2)
    a1_new <- mean(kf$alpha_smooth[, 1])
    P1_new <- mean(kf$P_smooth[, 1] + (kf$alpha_smooth[, 1] - a1_new)^2)
    gamma_new <- as.vector(ZtZinv %*% crossprod(Zmu, rowMeans(y - kf$alpha_smooth)))
    mu_new <- as.vector(Zmu %*% gamma_new)
    if (is.null(fix_sigma_eps2)) {
      sigma_eps2_new <- sum(kf$P_smooth + (y - mu_new - kf$alpha_smooth)^2) / (N * Tt)
    } else sigma_eps2_new <- fix_sigma_eps2
    opt <- optim(beta, neg_Q_beta_profiled, Xmat = Xmat, alpha_smooth = kf$alpha_smooth,
                 P_smooth = kf$P_smooth, P_cov_smooth = kf$P_cov_smooth, method = "BFGS",
                 control = list(reltol = 1e-12, maxit = 500))
    beta_new <- opt$par; phi_new <- tanh(as.vector(Xmat %*% beta_new))
    if (max(abs(phi_new)) > phi_cap) {
      n_backtrack <- n_backtrack + 1L
      cm <- constrained_mstep_beta(beta, neg_Q_beta_profiled, gradQ_beta_profiled, Xmat, phi_cap,
                                   alpha_smooth = kf$alpha_smooth, P_smooth = kf$P_smooth, P_cov_smooth = kf$P_cov_smooth)
      beta_new <- cm$beta; phi_new <- tanh(as.vector(Xmat %*% beta_new))
    }
    Cs <- compute_C(kf$alpha_smooth, kf$P_smooth, kf$P_cov_smooth)
    sigma_eta2_new <- sum(Cs$C0 - 2 * phi_new * Cs$C1 + phi_new^2 * Cs$C2) / (N * (Tt - 1))
    delta <- max(abs(c(beta_new - beta, gamma_new - gamma, a1_new - a1, P1_new - P1,
                        sigma_eta2_new - sigma_eta2, sigma_eps2_new - sigma_eps2)))
    beta <- beta_new; gamma <- gamma_new; a1 <- a1_new; P1 <- P1_new
    sigma_eta2 <- sigma_eta2_new; sigma_eps2 <- sigma_eps2_new
    if (delta < thres) break
  }
  phi <- tanh(as.vector(Xmat %*% beta)); mu <- as.vector(Zmu %*% gamma)
  kf_final <- kalman_run(y, phi, mu, a1, P1, sigma_eta2, sigma_eps2)
  Fv <- kf_final$P_pred + sigma_eps2; e <- y - mu - kf_final$alpha_pred
  list(beta = beta, gamma = gamma, mu = mu, a1 = a1, P1 = P1, sigma_eta2 = sigma_eta2, sigma_eps2 = sigma_eps2,
       iter = iter, converged = delta < thres, loglik = -0.5 * sum(log(2 * pi * Fv) + e^2 / Fv),
       sigma_eps2_fixed = !is.null(fix_sigma_eps2), kf = kf_final, n_backtrack = n_backtrack,
       at_cap = max(abs(phi)) >= phi_cap - 1e-6)
}

## ---------------------------------------------------------------------------
## Missing observations: kalman_run_na() is kalman_run() with NA entries in y
## handled by skipping the measurement update (K_it = 0, filtered = predicted)
## at those (i,t); the smoother and lag-one covariance are unchanged. Used for
## the sensitivity analysis in which grid slots that are not confirmed native
## readings (collinear points / gap-filled slots) are treated as missing.
## em_regression_mv_na() is the EM for fixed sigma_eps2 with this filter; the
## sigma_eps2 M-step (if not fixed) sums only over observed entries.
## ---------------------------------------------------------------------------
kalman_run_na <- function(y, phi, mu, a1, P1, sigma_eta2, sigma_eps2) {
  stopifnot(is.matrix(y), ncol(y) >= 2)   # supported domain: at least two time points per unit
  N <- nrow(y); Tt <- ncol(y); if (length(phi) == 1L) phi <- rep(phi, N); if (length(mu) == 1L) mu <- rep(mu, N)
  obs <- !is.na(y); y0 <- y; y0[!obs] <- 0
  alpha_pred <- P_pred <- alpha_filt <- P_filt <- K <- matrix(0, N, Tt)
  alpha_pred[, 1] <- a1; P_pred[, 1] <- P1
  for (t in 1:Tt) {
    if (t > 1) { alpha_pred[, t] <- phi * alpha_filt[, t - 1]; P_pred[, t] <- phi^2 * P_filt[, t - 1] + sigma_eta2 }
    Kt <- P_pred[, t] / (P_pred[, t] + sigma_eps2); Kt[!obs[, t]] <- 0; K[, t] <- Kt
    alpha_filt[, t] <- alpha_pred[, t] + Kt * (y0[, t] - mu - alpha_pred[, t]); P_filt[, t] <- (1 - Kt) * P_pred[, t]
  }
  alpha_smooth <- alpha_filt; P_smooth <- P_filt; J <- matrix(0, N, Tt)
  for (t in (Tt - 1):1) {
    J[, t] <- P_filt[, t] * phi / P_pred[, t + 1]
    alpha_smooth[, t] <- alpha_filt[, t] + J[, t] * (alpha_smooth[, t + 1] - alpha_pred[, t + 1])
    P_smooth[, t] <- P_filt[, t] + J[, t]^2 * (P_smooth[, t + 1] - P_pred[, t + 1])
  }
  P_cov_smooth <- matrix(0, N, Tt)   # P_cov_smooth[, t] = Cov(alpha_t, alpha_{t-1} | Y), t >= 2
  P_cov_smooth[, Tt] <- (1 - K[, Tt]) * phi * P_filt[, Tt - 1]
  if (Tt > 2) for (t in (Tt - 1):2) P_cov_smooth[, t] <- P_filt[, t] * J[, t - 1] + J[, t] * (P_cov_smooth[, t + 1] - phi * P_filt[, t]) * J[, t - 1]
  list(alpha_pred = alpha_pred, P_pred = P_pred, K = K, alpha_filt = alpha_filt, P_filt = P_filt, J = J,
       alpha_smooth = alpha_smooth, P_smooth = P_smooth, P_cov_smooth = P_cov_smooth, phi = phi, obs = obs)
}
em_regression_mv_na <- function(y, Xmat, mu = 0, init, thres = 1e-8, max_iter = 500, fix_sigma_eps2 = NULL, phi_cap = 0.995) {
  N <- nrow(y); Tt <- ncol(y); obs <- !is.na(y); n_obs <- sum(obs)
  beta <- init$beta; a1 <- init$a1; P1 <- init$P1; sigma_eta2 <- init$sigma_eta2
  sigma_eps2 <- if (!is.null(fix_sigma_eps2)) fix_sigma_eps2 else init$sigma_eps2; n_backtrack <- 0L
  nQ <- function(beta, alpha_smooth, P_smooth, P_cov_smooth) { phi_i <- tanh(as.vector(Xmat %*% beta)); Cs <- compute_C(alpha_smooth, P_smooth, P_cov_smooth)
    tot <- sum(Cs$C0 - 2 * phi_i * Cs$C1 + phi_i^2 * Cs$C2); s2 <- tot / (N * (Tt - 1)); v <- (N * (Tt - 1) / 2) * log(s2) + tot / (2 * s2); if (is.finite(v)) v else 1e300 }
  for (iter in 1:max_iter) {
    phi <- tanh(as.vector(Xmat %*% beta)); kf <- kalman_run_na(y, phi, mu, a1, P1, sigma_eta2, sigma_eps2)
    a1_new <- mean(kf$alpha_smooth[, 1]); P1_new <- mean(kf$P_smooth[, 1] + (kf$alpha_smooth[, 1] - a1_new)^2)
    sigma_eps2_new <- if (is.null(fix_sigma_eps2)) sum((kf$P_smooth + (ifelse(obs, y, 0) - mu - kf$alpha_smooth)^2)[obs]) / n_obs else fix_sigma_eps2
    opt <- optim(beta, nQ, alpha_smooth = kf$alpha_smooth, P_smooth = kf$P_smooth, P_cov_smooth = kf$P_cov_smooth, method = "BFGS", control = list(reltol = 1e-12, maxit = 500))
    beta_new <- opt$par; phi_new <- tanh(as.vector(Xmat %*% beta_new))
    if (max(abs(phi_new)) > phi_cap) { n_backtrack <- n_backtrack + 1L
      cm <- constrained_mstep_beta(beta, function(b, Xmat, alpha_smooth, P_smooth, P_cov_smooth) nQ(b, alpha_smooth, P_smooth, P_cov_smooth), gradQ_beta_profiled, Xmat, phi_cap,
                                   alpha_smooth = kf$alpha_smooth, P_smooth = kf$P_smooth, P_cov_smooth = kf$P_cov_smooth)
      beta_new <- cm$beta; phi_new <- tanh(as.vector(Xmat %*% beta_new)) }
    Cs <- compute_C(kf$alpha_smooth, kf$P_smooth, kf$P_cov_smooth); sigma_eta2_new <- sum(Cs$C0 - 2 * phi_new * Cs$C1 + phi_new^2 * Cs$C2) / (N * (Tt - 1))
    delta <- max(abs(c(beta_new - beta, a1_new - a1, P1_new - P1, sigma_eta2_new - sigma_eta2, sigma_eps2_new - sigma_eps2)))
    beta <- beta_new; a1 <- a1_new; P1 <- P1_new; sigma_eta2 <- sigma_eta2_new; sigma_eps2 <- sigma_eps2_new
    if (delta < thres) break
  }
  phi <- tanh(as.vector(Xmat %*% beta)); kf <- kalman_run_na(y, phi, mu, a1, P1, sigma_eta2, sigma_eps2)
  Fv <- kf$P_pred + sigma_eps2; e <- ifelse(obs, y, 0) - mu - kf$alpha_pred
  list(beta = beta, a1 = a1, P1 = P1, sigma_eta2 = sigma_eta2, sigma_eps2 = sigma_eps2, iter = iter, converged = delta < thres,
       loglik = -0.5 * sum((log(2 * pi * Fv) + e^2 / Fv)[obs]), sigma_eps2_fixed = !is.null(fix_sigma_eps2), kf = kf, n_backtrack = n_backtrack, at_cap = max(abs(phi)) >= phi_cap - 1e-6)
}
