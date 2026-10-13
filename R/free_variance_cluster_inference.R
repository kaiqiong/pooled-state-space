## Per-episode observed information and scores with a FREELY ESTIMATED observation variance.
## Parameter vector (beta, a1, P1, sigma_eta2, sigma_eps2). The entries are the unsimplified per-episode terms (no first-order
## conditions substituted), as in compute_info_per_unit_fixed_eps2, plus the sigma_eps2 row and column: with
## R_obs,i = sum_t [(y_it - alpha_hat_it)^2 + P_hat_it], the complete-data terms are
##   -d2Q/d(sigma_eps2)^2 = -T/(2 r^2) + R_obs/r^3,  cross terms zero,
## and the smoothed-moment correction of row sigma_eps2 against parameter l is -(1/(2 r^2)) dR_obs/dl, while the sigma_eps2
## column of the other rows uses the sensitivity seed dsigma_eps2 = 1. Per-episode scores: dQ/dr = -T/(2r) + R_obs/(2r^2).
## Functions: info_per_unit_free(y, Xmat, mu, fit), scores_per_unit_free(y, Xmat, mu, fit), cluster_inference_free(y, Xmat, mu, fit, cluster).
info_per_unit_free <- function(y, Xmat, mu, fit) {
  N <- nrow(y); Tt <- ncol(y); p <- ncol(Xmat); beta <- fit$beta; a1 <- fit$a1; P1 <- fit$P1; q <- fit$sigma_eta2; r <- fit$sigma_eps2
  phi <- tanh(as.vector(Xmat %*% beta)); kf <- fit$kf; Cs <- compute_C(kf$alpha_smooth, kf$P_smooth, kf$P_cov_smooth); rho <- Cs$C1 - phi * Cs$C2; w <- (1 - phi^2) * ((1 - phi^2) * Cs$C2 + 2 * phi * rho)
  bn <- colnames(Xmat); pn <- c(bn, "a1", "P1", "sigma_eta2", "sigma_eps2"); np <- length(pn); Ii <- array(0, c(N, np, np), dimnames = list(NULL, pn, pn))
  for (j in seq_len(p)) for (k in seq_len(p)) Ii[, bn[j], bn[k]] <- Xmat[, j] * Xmat[, k] * w / q
  Ii[, "a1", "a1"] <- 1 / P1; S_initial <- kf$P_smooth[, 1] + (kf$alpha_smooth[, 1] - a1)^2; R_tr <- Cs$C0 - 2 * phi * Cs$C1 + phi^2 * Cs$C2; resid <- y - mu - kf$alpha_smooth; R_obs <- rowSums(resid^2 + kf$P_smooth)
  Ii[, "P1", "P1"] <- -1 / (2 * P1^2) + S_initial / P1^3; Ii[, "sigma_eta2", "sigma_eta2"] <- -(Tt - 1) / (2 * q^2) + R_tr / q^3; Ii[, "sigma_eps2", "sigma_eps2"] <- -Tt / (2 * r^2) + R_obs / r^3
  Ii[, "a1", "P1"] <- Ii[, "P1", "a1"] <- (kf$alpha_smooth[, 1] - a1) / P1^2
  for (j in seq_len(p)) { bq <- Xmat[, j] * (1 - phi^2) * rho / q^2; Ii[, bn[j], "sigma_eta2"] <- bq; Ii[, "sigma_eta2", bn[j]] <- bq }
  seed_list <- list(); for (j in seq_len(p)) seed_list[[bn[j]]] <- list(dphi = (1 - phi^2) * Xmat[, j], dsigma_eta2 = 0, dsigma_eps2 = 0, da1 = 0, dP1 = 0)
  seed_list[["a1"]] <- list(dphi = rep(0, N), dsigma_eta2 = 0, dsigma_eps2 = 0, da1 = 1, dP1 = 0); seed_list[["P1"]] <- list(dphi = rep(0, N), dsigma_eta2 = 0, dsigma_eps2 = 0, da1 = 0, dP1 = 1)
  seed_list[["sigma_eta2"]] <- list(dphi = rep(0, N), dsigma_eta2 = 1, dsigma_eps2 = 0, da1 = 0, dP1 = 0); seed_list[["sigma_eps2"]] <- list(dphi = rep(0, N), dsigma_eta2 = 0, dsigma_eps2 = 1, da1 = 0, dP1 = 0)
  for (l in pn) { s <- seed_list[[l]]; sens <- kalman_sensitivity(kf, mu, q, r, y, dphi = s$dphi, dsigma_eta2 = s$dsigma_eta2, dsigma_eps2 = s$dsigma_eps2, da1 = s$da1, dP1 = s$dP1)
    dCs <- compute_C_dot(kf$alpha_smooth, kf$P_smooth, kf$P_cov_smooth, sens$dalpha_smooth, sens$dP_smooth, sens$dP_cov_smooth); drho <- dCs$dC1 - phi * dCs$dC2; dR_obs <- rowSums(sens$dP_smooth - 2 * resid * sens$dalpha_smooth)
    for (j in seq_len(p)) Ii[, bn[j], l] <- Ii[, bn[j], l] - (1 / q) * Xmat[, j] * (1 - phi^2) * drho
    Ii[, "a1", l] <- Ii[, "a1", l] - (1 / P1) * sens$dalpha_smooth[, 1]; Ii[, "P1", l] <- Ii[, "P1", l] - (1 / (2 * P1^2)) * (sens$dP_smooth[, 1] + 2 * (kf$alpha_smooth[, 1] - a1) * sens$dalpha_smooth[, 1])
    Ii[, "sigma_eta2", l] <- Ii[, "sigma_eta2", l] - (1 / (2 * q^2)) * (dCs$dC0 - 2 * phi * dCs$dC1 + phi^2 * dCs$dC2); Ii[, "sigma_eps2", l] <- Ii[, "sigma_eps2", l] - (1 / (2 * r^2)) * dR_obs }
  Ii }
scores_per_unit_free <- function(y, Xmat, mu, fit) {
  N <- nrow(y); Tt <- ncol(y); p <- ncol(Xmat); beta <- fit$beta; a1 <- fit$a1; P1 <- fit$P1; q <- fit$sigma_eta2; r <- fit$sigma_eps2; phi <- tanh(as.vector(Xmat %*% beta)); kf <- fit$kf
  Cs <- compute_C(kf$alpha_smooth, kf$P_smooth, kf$P_cov_smooth); rho <- Cs$C1 - phi * Cs$C2; resid <- y - mu - kf$alpha_smooth; R_obs <- rowSums(resid^2 + kf$P_smooth); R_tr <- Cs$C0 - 2 * phi * Cs$C1 + phi^2 * Cs$C2
  S <- cbind(Xmat * ((1 - phi^2) * rho / q), a1 = (kf$alpha_smooth[, 1] - a1) / P1, P1 = -1 / (2 * P1) + (kf$P_smooth[, 1] + (kf$alpha_smooth[, 1] - a1)^2) / (2 * P1^2), sigma_eta2 = -(Tt - 1) / (2 * q) + R_tr / (2 * q^2), sigma_eps2 = -Tt / (2 * r) + R_obs / (2 * r^2))
  colnames(S) <- c(colnames(Xmat), "a1", "P1", "sigma_eta2", "sigma_eps2"); S }
cluster_inference_free <- function(y, Xmat, mu, fit, cluster, G_df = NULL) {
  Ii <- info_per_unit_free(y, Xmat, mu, fit); S <- scores_per_unit_free(y, Xmat, mu, fit); p <- ncol(Xmat); pn <- dimnames(Ii)[[2]]; I <- apply(Ii, c(2, 3), sum); asym <- max(abs(I - t(I))) / max(abs(I)); I <- (I + t(I)) / 2
  ev <- eigen(I, symmetric = TRUE, only.values = TRUE)$values; Iinv <- solve(I); cl <- unique(cluster); G <- length(cl); Sg <- t(sapply(cl, function(g) colSums(S[cluster == g, , drop = FALSE]))); meat <- crossprod(Sg); V0 <- Iinv %*% meat %*% Iinv
  Vmd <- matrix(0, length(pn), length(pn)); for (g in cl) { idx <- which(cluster == g); Ig <- apply(Ii[idx, , , drop = FALSE], c(2, 3), sum); Ag <- solve(diag(length(pn)) - Ig %*% Iinv); v <- Ag %*% colSums(S[idx, , drop = FALSE]); Vmd <- Vmd + v %*% t(v) }; Vmd <- Iinv %*% Vmd %*% Iinv
  list(se_model = setNames(sqrt(diag(Iinv)), pn), se_cr0 = setNames(sqrt(diag(V0)), pn), se_md = setNames(sqrt(diag(Vmd)), pn), I_obs = I, score_sum = colSums(S), min_eig = min(ev), asymmetry = asym, G = G) }
