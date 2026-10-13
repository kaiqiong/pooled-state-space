## Schematic for Section 2.1 (three panels: hierarchy, model terms, half-life), matched to the PRIMARY ANALYSIS of the application: episodes are 14 readings at
## 5-minute spacing (65 minutes) starting at the post-meal glucose peak, centered on the glucose at meal onset
## (the pre-meal reference, mu_i = 0 on the centered scale). Simulated values use the primary-fit parameters
## (mean initial excess a1 = 56 and phi = 0.93 / 0.96 as fitted; the spread of initial values and the innovation variance are
## reduced from the fitted values so that single episodes resemble the group-mean trajectories); a small illustrative observation error is drawn
## in panel B although the primary fit estimates that variance at zero. Every drawn episode starts at its observed maximum.
suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(dplyr)})
set.seed(7); Delta <- 5; Tt <- 14
## Illustrative draws are kept only if the first observed reading is the maximum of the episode, mirroring the
## peak anchoring of the application (the latent maximum need not coincide with the observed one).
sim_episode <- function(phi, a1 = 56, P1 = 100, q = 12, r = 30, mu = 0) { repeat { a <- numeric(Tt); a[1] <- abs(rnorm(1, a1, sqrt(P1))); for (t in 2:Tt) a[t] <- phi * a[t - 1] + rnorm(1, 0, sqrt(q)); Y <- mu + a + rnorm(Tt, 0, sqrt(r)); if (which.max(Y) == 1) break }
  data.frame(t = (seq_len(Tt) - 1) * Delta, X = mu + a, Y = Y) }
## ---- Panel A: two participants, three post-meal episodes each, on one recording day (centered scale) ----
eps <- list(); starts <- list(g1 = c(0.2, 1.75, 3.3), g2 = c(0.1, 1.65, 3.2)); k <- 0   # episode start (peak) times on a compressed axis; gaps not to scale
for (g in 1:2) for (j in 1:3) { k <- k + 1; e <- sim_episode(c(0.93, 0.96)[g]); e$g <- g; e$i <- k; e$x <- starts[[g]][j] + e$t / 60; eps[[k]] <- e }
ep <- bind_rows(eps); pc <- c(`1` = "#4d4d4d", `2` = "#2b7bba")   # participant colours deliberately differ from the clinical-group palette used in panel C and Figure 3
pA <- ggplot(ep, aes(x, Y, group = i, colour = factor(g))) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.4) + geom_line(linewidth = 0.5) + geom_point(size = 1.3) +
  facet_wrap(~g, ncol = 1, labeller = labeller(g = c(`1` = "list('participant '*g == 1, C[1] == '{1, 2, 3}')", `2` = "list('participant '*g == 2, C[2] == '{4, 5, 6}')"), .default = label_parsed)) +
  geom_text(data = ep %>% group_by(g, i) %>% summarise(x = min(x), y = 118, lab = sprintf("episode i = %d", first(i)), .groups = "drop"), aes(x, y, label = lab), hjust = 0, size = 3, colour = "grey20", inherit.aes = FALSE) +
  geom_text(data = ep %>% group_by(g, i) %>% summarise(x = min(x), y = Y[1], .groups = "drop"), aes(x, y, label = "peak"), hjust = 1.2, vjust = 0.3, size = 2.6, colour = "grey40", inherit.aes = FALSE) +
  geom_segment(data = data.frame(g = 1:2, x = c(3.3, 3.2)), aes(x = x, xend = x + Delta / 60, y = -22, yend = -22), inherit.aes = FALSE, colour = "grey30", arrow = arrow(ends = "both", length = unit(0.08, "cm"))) +
  geom_segment(data = data.frame(g = 1:2, x = c(3.3, 3.2)), aes(x = x, xend = x, y = -22, yend = 70), inherit.aes = FALSE, colour = "grey70", linetype = 3) +
  geom_segment(data = data.frame(g = 1:2, x = c(3.3, 3.2) + Delta / 60), aes(x = x, xend = x, y = -22, yend = 70), inherit.aes = FALSE, colour = "grey70", linetype = 3) +
  geom_text(data = data.frame(g = 1:2, x = c(3.3, 3.2) + Delta / 120), aes(x = x, y = -22, label = "Delta==5~min"), parse = TRUE, inherit.aes = FALSE, vjust = 1.6, size = 3.1, colour = "grey30") +
  annotate("text", x = 0.05, y = 0, label = "episode-specific reference: zero after centering", hjust = 0, vjust = 1.5, size = 2.9, colour = "grey30") +
  scale_colour_manual(values = pc, guide = "none") + scale_x_continuous(breaks = NULL, limits = c(0, 4.5)) + coord_cartesian(ylim = c(-40, 128)) +
  labs(x = NULL, y = expression(paste("Centered glucose, ", Y[it] - mu[i], " (mg/dL)")), title = "Repeated episodes within participants") +
  theme_bw(base_size = 11) + theme(panel.grid = element_blank(), strip.background = element_rect(fill = "grey93"), plot.title = element_text(size = 11), axis.ticks.x = element_blank())
## ---- Panel B: one episode on the original glucose scale ----
mu <- 96; set.seed(11); e <- sim_episode(0.94, mu = mu); ta <- 4; te <- which.max(abs(e$Y - e$X)[6:13]) + 5
pB <- ggplot(e, aes(t)) +
  geom_hline(yintercept = mu, linetype = 2, colour = "grey30") + annotate("text", x = 65, y = mu, label = "mu[i]*': pre-meal glucose (reference)'", parse = TRUE, hjust = 1, vjust = -0.5, size = 3.4) +
  geom_line(aes(y = X), colour = "#1f78b4", linewidth = 1) + geom_point(aes(y = Y), size = 2.3) +
  annotate("segment", x = e$t[ta], xend = e$t[ta], y = mu, yend = e$X[ta], colour = "#1f78b4", linewidth = 0.9, arrow = arrow(length = unit(0.15, "cm"))) +
  annotate("text", x = e$t[ta] - 1.5, y = mu + 0.45 * (e$X[ta] - mu), label = "alpha[it]==X[it]-mu[i]", parse = TRUE, hjust = 1, size = 3.6, colour = "#1f78b4") +
  annotate("segment", x = e$t[te], xend = e$t[te], y = e$X[te], yend = e$Y[te], colour = "#b2182b", linewidth = 0.9, arrow = arrow(length = unit(0.15, "cm"))) +
  annotate("text", x = e$t[te] + 1.5, y = e$Y[te], label = "epsilon[it]==Y[it]-X[it]", parse = TRUE, hjust = 0, vjust = if (e$Y[te] > e$X[te]) -0.4 else 1.3, size = 3.6, colour = "#b2182b") +
  annotate("text", x = e$t[1] + 3, y = e$Y[1], label = "Y[it]~'(readings)'", parse = TRUE, hjust = 0, vjust = 0.4, size = 3.4) +
  annotate("text", x = e$t[ta] + 1.5, y = e$X[ta], label = "X[it]~'(latent glucose)'", parse = TRUE, hjust = 0, vjust = -0.9, size = 3.4, colour = "#1f78b4") +
  annotate("label", x = 65, y = max(e$Y) + 4, label = "Y[it]==X[it]+epsilon[it]", parse = TRUE, hjust = 1, size = 3.8, label.size = 0, fill = "grey97") +
  annotate("label", x = 65, y = max(e$Y) - 6, label = "alpha[it]==phi[i]*alpha[i*','*t-1]+eta[it]", parse = TRUE, hjust = 1, size = 3.8, label.size = 0, fill = "grey97") +
  scale_x_continuous(breaks = seq(0, 65, 15)) + coord_cartesian(xlim = c(-1, 66), ylim = c(mu - 10, max(e$Y) + 10)) +
  labs(x = "Minutes after the selected peak", y = "Glucose (mg/dL)", title = "Observation and latent-state model") +
  theme_bw(base_size = 11) + theme(panel.grid = element_blank(), plot.title = element_text(size = 11))
## ---- Panel C: expected deviation a1 phi^(t-1) and the model-based half-life, at the three group-mean persistence
##      coefficients of the primary fit (0.930, 0.943, 0.961: half-lives 48, 59, 87 min); the curve is the model mean ----
phis <- c(Healthy = 0.930, Prediabetic = 0.943, `Type 2 diabetes` = 0.961); a1c <- 56
hl <- Delta * log(2) / (-log(phis)); gcol <- c(Healthy = "#1b9e77", Prediabetic = "#d95f02", `Type 2 diabetes` = "#7570b3")
cur <- do.call(rbind, lapply(names(phis), function(k) data.frame(group = k, s = seq(0, 95, length.out = 401), m = a1c * phis[k]^(seq(0, 95, length.out = 401) / Delta), obs = seq(0, 95, length.out = 401) <= 65)))
pts <- do.call(rbind, lapply(names(phis), function(k) data.frame(group = k, s = (seq_len(Tt) - 1) * Delta, m = a1c * phis[k]^(seq_len(Tt) - 1))))
hld <- data.frame(group = names(phis), h = hl, y = a1c / 2)
pC <- ggplot(cur, aes(s, m, colour = group)) + geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.4) +
  geom_hline(yintercept = a1c / 2, linetype = 3, colour = "grey50") +
  geom_line(aes(linetype = obs), linewidth = 0.9, show.legend = FALSE) + scale_linetype_manual(values = c(`TRUE` = "solid", `FALSE` = "dotted"), guide = "none") + geom_point(data = pts, size = 1.1) +
  annotate("rect", xmin = 65, xmax = 96, ymin = -Inf, ymax = Inf, fill = "grey92", alpha = 0.5) +
  geom_segment(data = hld, aes(x = h, xend = h, y = -4, yend = y, colour = group), linetype = 2, linewidth = 0.5) +
  geom_point(data = hld, aes(h, y, colour = group), size = 2.6, shape = 21, fill = "white", stroke = 1.1) +
  geom_text(data = hld, aes(h, -4, label = sprintf("%.0f", h), colour = group), vjust = 1.6, size = 3.1, show.legend = FALSE) +
  annotate("text", x = 1, y = a1c / 2, label = "a[1]/2", parse = TRUE, hjust = 0, vjust = -0.5, size = 3.4, colour = "grey30") +
  annotate("text", x = 2, y = 12, label = "h[i] == frac(Delta*log(2), -log(phi[i]))", parse = TRUE, hjust = 0, vjust = 0.5, size = 3.8) +
  geom_text(data = data.frame(group = names(phis), s = 95, m = c(55, 50, 45)), aes(s, m, label = sprintf("list('%s: '*phi[i] == '%.3f', h[i] == %.0f~min)", group, phis, hl), colour = group), parse = TRUE, hjust = 1, size = 2.9, show.legend = FALSE) +
  scale_colour_manual(values = gcol, guide = "none") +
  scale_x_continuous(breaks = c(0, 30, 65, 90)) + scale_y_continuous(breaks = c(0, a1c / 2, a1c), labels = c("0", "", expression(a[1]))) +
  coord_cartesian(xlim = c(0, 96), ylim = c(-9, 60)) +
  labs(x = "Minutes after the selected peak", y = expression(paste("Expected deviation ", a[1]*phi[i]^{t-1})), title = "Persistence and model-based half-life") +
  theme_bw(base_size = 11) + theme(panel.grid = element_blank(), plot.title = element_text(size = 11), legend.position = "none")
row2 <- (plot_spacer() | pC | plot_spacer()) + plot_layout(widths = c(0.22, 1, 0.22))
fig <- (pA | pB) / row2 + plot_layout(heights = c(1, 0.72)) + plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")")
script_arg <- grep("^--file=", commandArgs(), value=TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), "../.."))
args <- commandArgs(trailingOnly=TRUE)
out_dir <- if(length(args)) args[1] else file.path(root,"results/figures")
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
out <- file.path(out_dir,"figure_1_model_schematic")
ggsave(paste0(out,".pdf"),fig,width=10.5,height=7.8,device=cairo_pdf)
ggsave(paste0(out,".png"),fig,width=10.5,height=7.8,dpi=300,bg="white")
