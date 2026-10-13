# Figures 4 and S3, from aggregate result summaries only.
# Rscript scripts/figures/render_references.R [aggregate_input_directory] [output_directory]
suppressPackageStartupMessages({library(ggplot2); library(patchwork); library(dplyr)})
script_arg <- grep("^--file=",commandArgs(),value=TRUE)
root <- normalizePath(file.path(dirname(sub("^--file=","",script_arg)),"../.."))
args <- commandArgs(trailingOnly=TRUE)
input_dir <- if(length(args)>=1) args[1] else file.path(root,"reference/figures/application")
output_dir <- if(length(args)>=2) args[2] else file.path(root,"results/figures")
dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
save_pair <- function(plot,stem,width,height) {
 ggsave(file.path(output_dir,paste0(stem,".pdf")),plot,width=width,height=height)
 ggsave(file.path(output_dir,paste0(stem,".png")),plot,width=width,height=height,dpi=300,bg="white")
}
grp_lab <- c(healthy="Healthy",prediabetic="Prediabetic",diabetic="Type 2 diabetes")
grp_col <- c(Healthy="#1b9e77",Prediabetic="#d95f02",`Type 2 diabetes`="#7570b3")
## Larger type for the standalone reference displays.
theme_paper <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey93"),
        legend.position = "bottom", plot.tag = element_text(face = "bold"))

ref_choices <- read.csv(file.path(input_dir,"fig_reference_choices.csv"))
stopifnot(nrow(ref_choices)==5L,all(ref_choices$B==2000L),all(ref_choices$interval=="Bootstrap"))
ref_choices$display <- factor(ref_choices$display,levels=rev(ref_choices$display))
ref_choices$interval <- factor(ref_choices$interval,levels="Bootstrap")
pReference <- ggplot(ref_choices, aes(est, display, colour=interval, shape=interval)) +
  geom_vline(xintercept=0, linetype=2, colour="grey55") +
  geom_errorbarh(aes(xmin=lo,xmax=hi),height=.16,linewidth=.65) +
  geom_point(size=2.5) +
  geom_text(aes(label=sprintf("%.3f",est)),vjust=-1.15,size=3.3,show.legend=FALSE) +
  scale_shape_manual(values=c(Bootstrap=16), guide="none") +
  scale_colour_manual(values=c(Bootstrap="#08306b"), guide="none") +
  scale_x_continuous(breaks=c(-.2,0,.2,.4),limits=c(-.29,.48)) +
  labs(x="Diabetes coefficient (participant-bootstrap 95% percentile interval, B = 2000)",y=NULL,colour=NULL,shape=NULL,tag="(b)") +
  theme_paper + theme(axis.text.y=element_text(size=10,lineheight=.95),
    legend.text=element_text(size=10),legend.key.width=grid::unit(12,"pt"),
    plot.margin=margin(5,10,5,5))

reference_choice_figure <- pReference + labs(tag=NULL) +
  theme(axis.text.y=element_text(size=11),axis.title=element_text(size=12),
        legend.text=element_text(size=11))
save_pair(reference_choice_figure,"figure_4_reference_choices",7.4,4.6)

## Noiseless two-reading illustration: the dotted lines are implied paths,
## not observed trajectories or fitted mean curves for the CGMacros data.
toy_panel <- function(reference,colour,tag) {
  phi <- (145-reference)/(150-reference)
  pred <- data.frame(t=seq(5,65,by=1))
  pred$y <- reference+(150-reference)*phi^(pred$t/5)
  readings <- data.frame(t=c(0,5),y=c(150,145))
  ggplot() +
    geom_hline(yintercept=reference,linetype=2,colour=colour,linewidth=.65) +
    geom_line(data=pred,aes(t,y),linetype=3,colour=colour,linewidth=.8) +
    geom_line(data=readings,aes(t,y),linewidth=.7) +
    geom_point(data=readings,aes(t,y),size=2.5) +
    annotate("text",x=c(0,6),y=c(153,145),label=c("150","145"),vjust=-.4,size=3.2) +
    annotate("segment",x=c(0,5),xend=c(0,5),y=reference,yend=c(150,145),
      colour="grey45",linewidth=.4) +
    annotate("label",x=c(-.5,6),y=c((reference+150)/2,(reference+145)/2),
      label=c(150-reference,145-reference),size=3,label.size=0,label.padding=grid::unit(.05,"lines")) +
    annotate("text",x=34,y=164,label=sprintf("atop(paste(phi, ' = ', frac(%d,%d), ' = %.2f'), t[1/2] == %.0f~min)",
      145-reference,150-reference,phi,5*log(2)/(-log(phi))),parse=TRUE,colour=colour,size=3.4) +
    annotate("text",x=64,y=reference,label=paste("Reference =",reference),
      hjust=1,vjust=1.4,colour=colour,size=3.1) +
    scale_x_continuous(breaks=c(0,20,40,60)) +
    coord_cartesian(xlim=c(-4,66),ylim=c(94,172)) +
    labs(x="Minutes",y="Glucose (mg/dL)",tag=tag) + theme_paper +
    theme(axis.text=element_text(size=9),axis.title=element_text(size=10))
}

reference_curves <- read.csv(file.path(input_dir,"fig_reference_descriptive_means.csv"))
reference_means <- read.csv(file.path(input_dir,"fig_reference_descriptive_levels.csv"))
stopifnot(sum(reference_means$n)==1171L)
reference_curves$group <- factor(reference_curves$group,levels=unname(grp_lab))
reference_means$group <- factor(reference_means$group,levels=unname(grp_lab))
pDescriptive <- ggplot(reference_curves,aes(t,y,colour=group)) +
  geom_hline(yintercept=0,colour="grey25",linewidth=.6) +
  geom_hline(data=reference_means,aes(yintercept=reference,colour=group),
    linetype=2,linewidth=.7) +
  geom_line(linewidth=1) + geom_point(size=1.5) +
  scale_colour_manual(values=grp_col) +
  scale_x_continuous(breaks=c(0,20,40,60)) +
  labs(x="Minutes after glucose peak",y="Mean glucose deviation\n(mg/dL)",
    colour=NULL,tag="(c)") +
  theme_paper + theme(legend.text=element_text(size=10),axis.title=element_text(size=10))
reference_figure <- ((toy_panel(100,"#1f78b4","(a)") | toy_panel(130,"#d95f02","(b)")) /
    pDescriptive) + plot_layout(heights=c(1,1))
save_pair(reference_figure,"figure_S3_reference_explanation",7.4,6.9)

