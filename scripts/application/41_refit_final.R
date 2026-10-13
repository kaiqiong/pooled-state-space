## Step 41: prepare the final-window design and the fixed-SD20 comparison fit.
## Primary boundary inference is produced separately by Step 51.
source(file.path(APP_ROOT, "R", "state_space_em.R"))
w <- readRDS("meal_windows_final_65min.rds")
sub <- read.csv(file.path(APP_DATA, "bio.csv"), check.names = FALSE)
sub$subject <- sprintf("CGMacros-%03d", sub$subject)
sub$group <- cut(sub[["A1c PDL (Lab)"]], breaks = c(-Inf, 5.7, 6.5, Inf), labels = c("healthy", "prediabetic", "diabetic"), right = FALSE)
rownames(sub) <- sub$subject
## Six aggregate historical initialization coefficients, preserved from Step 41.
## No historical participant cache is needed to run the analysis.
start_beta <- c(1.585644613623052, 0.053624330419135155, 0.079270609589127183,
                0.014654732116248322, 0.058572078290040568, 0.30150286379177538)
s <- sub[w$meta$subject, ]; stopifnot(!any(is.na(s$Age)))
Xmat <- cbind(Intercept = 1, Age = as.numeric(scale(s$Age)), Sex = as.numeric(s$Gender == "F"), BMI = as.numeric(scale(s$BMI)),
              Prediabetic = as.numeric(s$group == "prediabetic"), Diabetic = as.numeric(s$group == "diabetic"))
Y <- w$Y - w$meta$baseline_gl; subj <- w$meta$subject; G <- length(unique(subj)); q <- qt(0.975, G - 1); se2 <- 400
cat("final: N =", nrow(Y), " by group:", table(s$group), " windows/participant median", median(table(subj)), "\n")
fit <- em_regression_mv(Y, Xmat, 0, list(beta = start_beta, a1 = mean(Y[, 1]), P1 = max(var(Y[, 1]), 1), sigma_eta2 = 40, sigma_eps2 = se2), thres = 1e-7, max_iter = 4000, fix_sigma_eps2 = se2)
cr2 <- compute_SE_regression_mv_cluster_CR2(Y, Xmat, 0, fit, subj); cr0 <- compute_SE_regression_mv_cluster_robust(Y, Xmat, 0, fit, subj)
ph <- tanh(Xmat %*% fit$beta)
cmp <- data.frame(parameter = colnames(Xmat), new_est = fit$beta, se_model = cr2$se_model_based, se_cr0 = cr2$se_cr0, se_cr2 = cr2$se, lo_cr2 = fit$beta - q * cr2$se, hi_cr2 = fit$beta + q * cr2$se)
print(cmp, digits = 3, row.names = FALSE)
cat(sprintf("phi by group: %s | a1 %.1f P1 %.0f sigma_eta2 %.1f | loglik %.1f | iter %d | max leverage %.2f\n", paste(round(tapply(ph, s$group, mean), 4), collapse = " / "), fit$a1, fit$P1, fit$sigma_eta2, fit$loglik, fit$iter, max(cr2$leverage)))
## naive two-stage
nv <- naive_two_stage(Y, Xmat); cat("naive: ", round(nv$beta_hat, 3), "\n")
## collinear-as-missing sensitivity
Yna <- Y; Yna[w$flags != "kink"] <- NA; cat("collinear slots set NA:", sum(is.na(Yna)), sprintf("(%.1f%%)", 100 * mean(is.na(Yna))), "\n")
fna <- em_regression_mv_na(Yna, Xmat, 0, list(beta = fit$beta, a1 = fit$a1, P1 = fit$P1, sigma_eta2 = fit$sigma_eta2, sigma_eps2 = se2), thres = 1e-7, max_iter = 4000, fix_sigma_eps2 = se2)
cat("collinear-as-missing: beta", round(fna$beta, 3), " phi by group", round(tapply(tanh(Xmat %*% fna$beta), s$group, mean), 4), " sigma_eta2", round(fna$sigma_eta2, 1), " iter", fna$iter, "\n")
saveRDS(list(slot_time = w$slot_time, fit = fit, se_cr2 = cr2, se_cr0 = cr0, naive_fit = nv, fit_na = fna, Xmat = Xmat, Y = Y, Yna = Yna, meta = cbind(w$meta, s[, c("Age", "Gender", "BMI", "group")]), flags = w$flags, sd_eps_primary = 20, cmp = cmp), "FINAL_fitted_models_final.rds")
write.csv(cmp, "fixed_sd20_coefficients.csv", row.names = FALSE); cat("saved FINAL_fitted_models_final.rds\n")
