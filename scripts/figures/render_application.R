## Main Figure 3: CGMacros application. Uses locally generated episode-level data.
## Usage: Rscript scripts/figures/render_application.R [saved_input_directory] [output_directory]
##  (a) recovery trajectories by clinical group (thin lines: windows; bold: group mean)
##  (b) estimates over the measurement-error SD: 0 = boundary (primary, the maximum-likelihood value), 2-30 fixed (sensitivity)
##  (c) forest plot: joint model-based / joint CR0 / joint participant bootstrap (primary) / naive participant bootstrap
suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(dplyr); library(tidyr)})
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
stopifnot(length(script_arg) == 1L)
figure_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
here <- normalizePath(file.path(figure_dir, "../.."))
args <- commandArgs(trailingOnly=TRUE)
rd <- if (length(args)>=1) normalizePath(args[1]) else file.path(here, "results/application")
output_dir <- if (length(args)>=2) normalizePath(args[2]) else file.path(here,"results/figures")
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
f <- readRDS(file.path(rd, "PRIMARY_boundary.rds"))   # primary analysis: measurement-error variance at its boundary (0)

theme_paper <- theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey93"),
        legend.position = "bottom", plot.tag = element_text(face = "bold"))
grp_lab <- c(healthy = "Healthy", prediabetic = "Prediabetic", diabetic = "Type 2 diabetes")
grp_col <- c(Healthy = "#1b9e77", Prediabetic = "#d95f02", `Type 2 diabetes` = "#7570b3")

## ---- (a) trajectories ----
phi_hat <- as.vector(tanh(f$Xmat %*% f$fit$beta))
Y <- f$Y; N <- nrow(Y); Tt <- ncol(Y)
traj <- data.frame(window = rep(seq_len(N), Tt), t = rep((0:(Tt - 1)) * 5, each = N),
                   y = as.vector(Y), group = rep(grp_lab[as.character(f$meta$group)], Tt))
grp_mean <- traj %>% group_by(group, t) %>% summarise(y = mean(y), .groups = "drop")
phi_by_grp <- data.frame(group = grp_lab[names(tapply(phi_hat, f$meta$group, mean))],
                         phi = tapply(phi_hat, f$meta$group, mean))
phi_by_grp$lab <- sprintf("mean~hat(phi)==%.3f", phi_by_grp$phi)
set.seed(1); show_idx <- unlist(lapply(split(seq_len(N), traj$group[1:N]), function(ix) sample(ix, min(60, length(ix)))))
pA <- ggplot() +
  geom_line(data = traj[traj$window %in% show_idx, ], aes(t, y, group = window, colour = group), alpha = 0.18, linewidth = 0.3) +
  geom_line(data = grp_mean, aes(t, y, colour = group), linewidth = 1.3) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
  geom_text(data = phi_by_grp, aes(x = 64, y = 118, label = lab), parse = TRUE, hjust = 1, size = 3.2) +
  facet_wrap(~group) + scale_colour_manual(values = grp_col, guide = "none") +
  labs(x = "Minutes after the selected glucose peak", y = "Glucose above pre-meal baseline (mg/dL)") +
  coord_cartesian(ylim = c(-40, 125)) + theme_paper

## ---- (b) profile (safeguarded EM, v2): coefficients on top, phi below ----
prof <- read.csv(file.path(rd, "sigma_eps2_profile_final.csv"))
b0 <- setNames(f$fit$beta, colnames(f$Xmat))
prof <- bind_rows(data.frame(sd_eps = 0, phi_mean = mean(phi_hat), phi_min = min(phi_hat), phi_max = max(phi_hat), Age = b0["Age"], Sex = b0["Sex"], BMI = b0["BMI"], Prediabetic = b0["Prediabetic"], Diabetic = b0["Diabetic"], row.names = NULL),
                  prof[, c("sd_eps", "phi_mean", "phi_min", "phi_max", "Age", "Sex", "BMI", "Prediabetic", "Diabetic")])
cov_col <- c(Age = "#1f78b4", `Sex (female)` = "#e7298a", BMI = "#a6761d", Prediabetic = "#d95f02", `Type 2 diabetes` = "#7570b3")   # group contrasts keep the clinical-group colours of panel (a); other covariates avoid them
cov_lab <- c(Age = "Age", Sex = "Sex (female)", BMI = "BMI", Prediabetic = "Prediabetic", Diabetic = "Type 2 diabetes")
pl <- prof %>% select(sd_eps, Age, Sex, BMI, Prediabetic, Diabetic) %>%
  tidyr::pivot_longer(-sd_eps, names_to = "covariate", values_to = "coef") %>%
  mutate(covariate = factor(cov_lab[covariate], levels = names(cov_col)))
pB1 <- ggplot(pl, aes(sd_eps, coef, colour = covariate, shape = covariate)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey50") +
  geom_vline(xintercept = 0, linetype = 3) +
  geom_line() + geom_point(size = 1.8) +
  scale_colour_manual(values = cov_col) + scale_shape_manual(values = c(16, 17, 15, 3, 7)) +
  scale_x_continuous(breaks = c(0, 6, 10, 15, 20, 26, 30)) +
  labs(x = NULL, y = "Coefficient (tanh-link)", colour = NULL, shape = NULL) +
  theme_paper + theme(legend.position = "right", legend.text = element_text(size = 8), legend.key.height = grid::unit(10, "pt"))
pB2 <- ggplot(prof, aes(sd_eps)) +
  geom_ribbon(aes(ymin = phi_min, ymax = phi_max), fill = "grey50", alpha = 0.25) +
  geom_line(aes(y = phi_mean), colour = "grey30", linewidth = 0.9) + geom_point(aes(y = phi_mean), colour = "grey30", size = 1.8) +
  geom_vline(xintercept = 0, linetype = 3) +
  annotate("text", x = 0.6, y = 0.985, label = "primary (boundary)", size = 2.8, hjust = 0, colour = "grey20") +
  scale_x_continuous(breaks = c(0, 6, 10, 15, 20, 26, 30)) +
  labs(x = expression(paste(sigma[epsilon], " (mg/dL); 0 = boundary (primary)")),
       y = expression(paste(hat(phi)[i], " (mean, range)"))) +
  coord_cartesian(ylim = c(0.89, 1.0)) + theme_paper
pB <- (pB1 / pB2) + plot_layout(heights = c(1.4, 1))

## ---- (c) forest ----
covs <- c("Age", "Sex", "BMI", "Prediabetic", "Diabetic")
cov_lab <- c(Age = "Age", Sex = "Sex (female)", BMI = "BMI", Prediabetic = "Prediabetic", Diabetic = "Type 2 diabetes")
cj <- read.csv(file.path(rd, "primary_boundary_joint.csv")); cn <- read.csv(file.path(rd, "primary_boundary_naive.csv"))
lev <- c("Joint: model-based SE", "Joint: cluster CR0 sandwich", "Joint: participant bootstrap (primary)", "Naive two-stage: participant bootstrap")
fd <- bind_rows(
  data.frame(parameter = cj$parameter, est = cj$estimate, lo = cj$estimate - qnorm(.975) * cj$se_model, hi = cj$estimate + qnorm(.975) * cj$se_model, estimator = lev[1]),
  data.frame(parameter = cj$parameter, est = cj$estimate, lo = cj$cr0_lo_t, hi = cj$cr0_hi_t, estimator = lev[2]),
  data.frame(parameter = cj$parameter, est = cj$estimate, lo = cj$boot_pct_lo, hi = cj$boot_pct_hi, estimator = lev[3]),
  data.frame(parameter = cn$parameter, est = cn$estimate, lo = cn$boot_pct_lo, hi = cn$boot_pct_hi, estimator = lev[4])) %>%
  filter(parameter %in% covs) %>%
  mutate(parameter = factor(cov_lab[parameter], levels = rev(cov_lab)), estimator = factor(estimator, levels = lev))
lev_lab <- c("Joint: model-based SE", expression(paste("Joint: cluster CR0 sandwich, ", t[43])), "Joint: participant bootstrap, percentile (primary)", "Naive two-stage: participant bootstrap, percentile")
pC <- ggplot(fd, aes(x = est, y = parameter, colour = estimator, shape = estimator)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2, position = position_dodge(width = 0.75), linewidth = 0.55) +
  geom_point(position = position_dodge(width = 0.75), size = 2.2, fill = "white") +
  scale_colour_manual(values = c("grey55", "#6baed6", "#08306b", "#e31a1c"), labels = lev_lab) +
  scale_shape_manual(values = c(1, 21, 16, 17), labels = lev_lab) +
  labs(x = expression(paste("Effect on ", phi[i], " (tanh-link scale), 95% CI")), y = NULL, colour = NULL, shape = NULL) +
  theme_paper + theme(legend.direction = "vertical", legend.text = element_text(size = 9))


## Restore the original three-panel application figure, including variance sensitivity.
fig <- pA / (wrap_elements(full = pB) | pC) + plot_layout(heights = c(1, 1.25)) + plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")")
save_pair <- function(plot,stem,width,height) {
  ggsave(file.path(output_dir,paste0(stem,".pdf")),plot,width=width,height=height)
  ggsave(file.path(output_dir,paste0(stem,".png")),plot,width=width,height=height,dpi=300,bg="white")
}
save_pair(fig,"figure_3_application",10,8.2)

write.csv(prof,file.path(output_dir,"application_variance_profile.csv"),row.names=FALSE)
write.csv(fd,file.path(output_dir,"application_intervals.csv"),row.names=FALSE)
