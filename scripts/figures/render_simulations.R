# Figures 2 and S1. Saved-summary rendering; no fits are performed.
# Rscript scripts/figures/render_simulations.R [input_directory] [output_directory]
suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(dplyr); library(tidyr)})
script_arg <- grep("^--file=", commandArgs(), value=TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), "../.."))
args <- commandArgs(trailingOnly=TRUE)
data_dir <- if(length(args)>=1) normalizePath(args[1]) else file.path(root,"reference/figures/simulation")
out_dir <- if(length(args)>=2) args[2] else file.path(root,"results/figures")
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
reg4 <- c("Intercept", "Age", "Sex", "Biomarker")
reg6 <- c("Intercept", "Age", "Sex", "BMI", "Prediabetic", "Diabetic")
base_theme <- theme_bw(base_size = 10.5) + theme(
  panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
  panel.grid.major.x = element_line(colour = "grey92", linewidth = .25),
  plot.title = element_text(size = 11.5, face = "bold"),
  plot.subtitle = element_text(size = 9, margin = margin(b = 7)),
  axis.text = element_text(size = 9), axis.title = element_text(size = 10),
  strip.background = element_rect(fill = "grey95", colour = "grey75", linewidth = .3),
  strip.text = element_text(size = 9), legend.position = "bottom",
  legend.title = element_blank(), legend.text = element_text(size = 8.3),
  legend.key.width = grid::unit(11,"pt"), legend.key.height = grid::unit(9,"pt"),
  legend.margin = margin(0,0,0,0), legend.spacing.x = grid::unit(3,"pt"),
  plot.margin = margin(7,9,7,5))
method_order <- c("model", "cr0", "md", "bootstrap")
method_labels <- c(model = "Model-based", cr0 = "CR0", md = "MD-style", bootstrap = "Bootstrap")
method_colour <- c(model = "#666666", cr0 = "#D55E00", md = "#8D62A5", bootstrap = "#0072B2")
tau_labels <- c(`0` = "No heterogeneity", `0.4` = "Participant heterogeneity (tau = 0.4)")

# A: average RMSE across the same four regression coefficients, not root-average-MSE.
a0 <- read.csv(file.path(data_dir,"benchmark.csv"))
stopifnot(nrow(a0) == 36L, setequal(unique(a0$parameter),reg4))
a <- a0 %>% group_by(sigma_eps2,method) %>% summarise(mean_rmse = mean(rmse),.groups="drop")
map_a <- c(joint_EM="Joint",naive_two_stage="OLS two-stage",per_subject_kalman="State-space two-stage")
a$method <- factor(a$method,levels=names(map_a),labels=map_a)
pA <- ggplot(a,aes(factor(sigma_eps2),mean_rmse,group=method,colour=method,shape=method)) +
  geom_line(linewidth=.65) + geom_point(size=2.4) +
  scale_colour_manual(values=c("Joint"="#0072B2","OLS two-stage"="#D55E00","State-space two-stage"="#009E73")) +
  scale_shape_manual(values=c(16,17,15)) +
  scale_y_continuous(limits=c(0,.15),breaks=seq(0,.15,.05),expand=expansion(mult=c(0,.03))) +
  labs(title="(A) Pooling short episodes",subtitle="Independent episodes; M = 60 per variance",
       x=expression(paste("Observation variance ",sigma[epsilon]^2)),y="Mean coefficient RMSE") +
  guides(colour=guide_legend(nrow=1),shape=guide_legend(nrow=1)) + base_theme
write.csv(a,file.path(out_dir,"panel_A_derived.csv"),row.names=FALSE)

# B: full observed information versus the direct term only, at common fitted estimates
# with observation variance estimated rather than fixed during fitting.
b0 <- read.csv(file.path(data_dir,"information_ablation.csv"))
stopifnot(all(c("n_total","n_valid","method","variant") %in% names(b0)),
          all(b0$n_total == 500L),all(b0$n_valid == 500L))
b <- b0 %>% filter(tau==0,variant %in% c("full","no_correction"),method %in% c("model","cr0"))
stopifnot(nrow(b)==24L,setequal(unique(b$parameter),reg6),
          !anyDuplicated(b[c("variant","method","parameter")]))
b$group <- paste(b$method,ifelse(b$variant=="full","Full","Direct only"),sep="|")
b_order <- c("model|Full","model|Direct only","cr0|Full","cr0|Direct only")
b_labels <- c("Model: full", "Model: direct only", "CR0: full", "CR0: direct only")
b$position <- 5-match(b$group,b_order)+(match(b$parameter,reg6)-3.5)*.064
b$information <- factor(ifelse(b$variant=="full","Full","Direct only"),levels=c("Full","Direct only"))
pB <- ggplot(b,aes(coverage,position,colour=information,shape=information)) +
  geom_vline(xintercept=.95,linetype=2,colour="grey45",linewidth=.45) + geom_point(size=2.2,alpha=.9) +
  scale_colour_manual(values=c(Full="#0072B2",`Direct only`="#D55E00"),guide="none") +
  scale_shape_manual(values=c(16,17),guide="none") +
  scale_y_continuous(breaks=4:1,labels=b_labels,limits=c(.5,4.5)) +
  scale_x_continuous(limits=c(.63,.985),breaks=c(.65,.75,.85,.95)) +
  labs(title="(B) Contribution of the correction",subtitle="Variance estimated; tau = 0; M = 500",
       x="95% interval coverage",y=NULL) + base_theme
write.csv(b,file.path(out_dir,"panel_B_derived.csv"),row.names=FALSE)

# C/D use matched free-variance fits: the same 200 datasets per scenario
# for all four methods and all six regression coefficients.
pair <- read.csv(file.path(data_dir,"clustered_free_variance.csv"))
stopifnot(nrow(pair)==48L,all(pair$n==200L),setequal(unique(pair$parameter),reg6),
          setequal(unique(pair$method),method_order))
pair$scenario <- factor(pair$tau,levels=c(0,.4),labels=tau_labels)
pair$method <- factor(pair$method,levels=method_order)
pair$position <- 5-match(pair$method,method_order)+(match(pair$parameter,reg6)-3.5)*.075
pC <- ggplot(pair,aes(coverage,position,colour=method)) +
  geom_vline(xintercept=.95,linetype=2,colour="grey45",linewidth=.45) + geom_point(size=2.1,alpha=.9) +
  facet_wrap(~scenario,ncol=1) +
  scale_colour_manual(values=method_colour,guide="none") +
  scale_y_continuous(breaks=4:1,labels=method_labels,limits=c(.45,4.55)) +
  scale_x_continuous(limits=c(.24,1.005),breaks=c(.3,.5,.7,.95)) +
  labs(title="(C) Coverage under clustering",subtitle="Variance estimated; M = 200 per scenario",
       x="95% interval coverage",y=NULL) + base_theme
cr <- pair %>% filter(method=="cr0") %>% select(tau,parameter,cr0_length=mean_length)
d <- pair %>% left_join(cr,by=c("tau","parameter")) %>% mutate(length_ratio=mean_length/cr0_length)
stopifnot(all(abs(d$length_ratio[d$method=="cr0"]-1)<1e-12))
pD <- ggplot(d,aes(length_ratio,position,colour=method)) +
  geom_vline(xintercept=1,linetype=2,colour="grey45",linewidth=.45) + geom_point(size=2.1,alpha=.9) +
  facet_wrap(~scenario,ncol=1) +
  scale_colour_manual(values=method_colour,guide="none") +
  scale_y_continuous(breaks=4:1,labels=method_labels,limits=c(.45,4.55)) +
  scale_x_continuous(limits=c(.18,1.59),breaks=c(.25,.5,1,1.5)) +
  labs(title="(D) Corresponding interval lengths",subtitle="Same datasets and intervals as panel C",
       x="Mean interval length / CR0 mean length",y=NULL) + base_theme
write.csv(d,file.path(out_dir,"panels_CD_derived.csv"),row.names=FALSE)
fig <- (pA | pB) / (pC | pD) + plot_layout(heights=c(1,1.27))
ggsave(file.path(out_dir,"figure_2_simulations.pdf"),fig,width=8.1,height=7.15,device=cairo_pdf)
ggsave(file.path(out_dir,"figure_2_simulations.png"),fig,width=8.1,height=7.15,dpi=300,bg="white")

# Supplement: show all four regression coefficients, rather than averages hiding departures.
j0 <- read.csv(file.path(data_dir,"episode_size.csv"))
j <- j0 %>% filter(parameter %in% reg4)
stopifnot(nrow(j)==48L,all(j$n_se_ok==100L),setequal(unique(j$N),c(50,100,200,400)),
          setequal(unique(j$T),c(15,30,60)))
j$parameter <- factor(j$parameter,levels=reg4)
j$T <- factor(j$T,levels=c(15,30,60))
cols <- c(`15`="#0072B2",`30`="#D55E00",`60`="#009E73")
j_rmse <- ggplot(j,aes(N,rmse,colour=T,shape=T,group=T)) +
  geom_line(linewidth=.55) + geom_point(size=1.9) + facet_wrap(~parameter,nrow=1) +
  scale_x_log10(breaks=c(50,100,200,400)) +
  scale_colour_manual(values=cols,name="Readings per episode") +
  scale_shape_manual(values=c(16,17,15),name="Readings per episode") +
  labs(title="(A) Estimation error",x="Number of episodes",y="Coefficient RMSE") + base_theme +
  theme(legend.position="none",axis.text.x=element_text(size=8))
j_cov <- ggplot(j,aes(N,coverage,colour=T,shape=T,group=T)) +
  geom_hline(yintercept=.95,linetype=2,colour="grey45",linewidth=.4) +
  geom_line(linewidth=.55) + geom_point(size=1.9) + facet_wrap(~parameter,nrow=1) +
  scale_x_log10(breaks=c(50,100,200,400)) +
  scale_y_continuous(limits=c(.88,1.005),breaks=c(.9,.95,1)) +
  scale_colour_manual(values=cols,name="Readings per episode") +
  scale_shape_manual(values=c(16,17,15),name="Readings per episode") +
  labs(title="(B) Model-based interval coverage",x="Number of episodes",y="95% interval coverage") + base_theme +
  theme(axis.text.x=element_text(size=8),legend.title=element_text(size=9))
supp <- j_rmse / j_cov + plot_layout(heights=c(1,1.05))
ggsave(file.path(out_dir,"figure_S1_episode_size.pdf"),supp,width=8.2,height=5.0,device=cairo_pdf)
ggsave(file.path(out_dir,"figure_S1_episode_size.png"),supp,width=8.2,height=5.0,dpi=300,bg="white")
write.csv(j,file.path(out_dir,"supplement_size_derived.csv"),row.names=FALSE)
cat("Saved main and supplementary simulation figures in",out_dir,"\n")
