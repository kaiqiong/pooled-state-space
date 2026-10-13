# Covariance-only diagnostics at unchanged free-variance estimates.
# All coordinates are raw parameters: beta, a1, P1, q, r.
direct_information_free <- function(y, Xmat, mu, fit) {
  N <- nrow(y); Tt <- ncol(y); p <- ncol(Xmat)
  q <- fit$sigma_eta2; r <- fit$sigma_eps2; P1 <- fit$P1
  kf <- fit$kf; phi <- tanh(as.vector(Xmat %*% fit$beta))
  Cs <- compute_C(kf$alpha_smooth, kf$P_smooth, kf$P_cov_smooth)
  rho <- Cs$C1 - phi * Cs$C2
  w <- (1 - phi^2) * ((1 - phi^2) * Cs$C2 + 2 * phi * rho)
  pn <- c(colnames(Xmat), "a1", "P1", "sigma_eta2", "sigma_eps2")
  D <- matrix(0, p + 4L, p + 4L, dimnames = list(pn, pn))
  D[seq_len(p), seq_len(p)] <- crossprod(Xmat, Xmat * (w / q))
  D["a1", "a1"] <- N / P1
  S <- kf$P_smooth[, 1] + (kf$alpha_smooth[, 1] - fit$a1)^2
  A <- Cs$C0 - 2 * phi * Cs$C1 + phi^2 * Cs$C2
  R <- rowSums(kf$P_smooth + (y - mu - kf$alpha_smooth)^2)
  D["P1", "P1"] <- -N / (2 * P1^2) + sum(S) / P1^3
  D["sigma_eta2", "sigma_eta2"] <- -N * (Tt - 1) / (2 * q^2) + sum(A) / q^3
  D["sigma_eps2", "sigma_eps2"] <- -N * Tt / (2 * r^2) + sum(R) / r^3
  D["a1", "P1"] <- D["P1", "a1"] <- sum(kf$alpha_smooth[, 1] - fit$a1) / P1^2
  bq <- colSums(Xmat * ((1 - phi^2) * rho / q^2))
  D[seq_len(p), "sigma_eta2"] <- bq
  D["sigma_eta2", seq_len(p)] <- bq
  # Retain these cross terms even though their pooled sums can be tiny.
  (D + t(D)) / 2
}

covariance_from_information <- function(H, Sg, p) {
  H <- (H + t(H)) / 2
  ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  pd <- all(is.finite(ev)) && min(ev) > 0
  empty <- setNames(rep(NA_real_, ncol(H)), colnames(H))
  inv <- if (pd) tryCatch(solve(H), error = function(e) NULL) else NULL
  if (is.null(inv)) return(list(model = empty[seq_len(p)], cr0 = empty[seq_len(p)],
      model_all = empty, cr0_all = empty, pd = pd, finite = FALSE,
      min_eigenvalue = min(ev), condition_number = NA_real_))
  V0 <- inv %*% crossprod(Sg) %*% t(inv)
  good <- all(is.finite(inv)) && all(is.finite(V0)) &&
          all(diag(inv) > 0) && all(diag(V0) > 0)
  if (!good) return(list(model = empty[seq_len(p)], cr0 = empty[seq_len(p)],
      model_all = empty, cr0_all = empty, pd = pd, finite = FALSE,
      min_eigenvalue = min(ev), condition_number = max(ev) / min(ev)))
  sm <- sqrt(diag(inv)); sc <- sqrt(diag(V0))
  list(model = sm[seq_len(p)], cr0 = sc[seq_len(p)], model_all = sm,
       cr0_all = sc, pd = pd, finite = TRUE, min_eigenvalue = min(ev),
       condition_number = max(ev) / min(ev))
}

evaluate_variants <- function(H, D, Sg, p) {
  list(full = covariance_from_information(H, Sg, p),
       no_correction = covariance_from_information(D, Sg, p),
       nuisance_fixed = covariance_from_information(
         H[seq_len(p), seq_len(p), drop = FALSE],
         Sg[, seq_len(p), drop = FALSE], p))
}
