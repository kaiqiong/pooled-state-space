generate_covariates <- function(N, sex_prob = 0.8, age_biomarker_cor = 0.3) {
  Sex <- rbinom(N, 1, sex_prob)
  Sigma <- matrix(c(1, age_biomarker_cor, age_biomarker_cor, 1), 2, 2)
  Z <- matrix(rnorm(N * 2), N, 2) %*% chol(Sigma)
  cbind(Intercept = 1, Age = Z[, 1], Sex = Sex, Biomarker = Z[, 2])
}

generate_data_mv <- function(Xmat, beta, Tt, mu, a1, P1, sigma_eta2, sigma_eps2) {
  N <- nrow(Xmat)
  phi <- tanh(as.vector(Xmat %*% beta))
  alpha <- matrix(0, N, Tt); y <- matrix(0, N, Tt)
  alpha[, 1] <- rnorm(N, a1, sqrt(P1))
  y[, 1] <- mu + alpha[, 1] + rnorm(N, 0, sqrt(sigma_eps2))
  for (t in 2:Tt) {
    alpha[, t] <- phi * alpha[, t - 1] + rnorm(N, 0, sqrt(sigma_eta2))
    y[, t] <- mu + alpha[, t] + rnorm(N, 0, sqrt(sigma_eps2))
  }
  list(y = y, alpha = alpha, phi = phi)
}

