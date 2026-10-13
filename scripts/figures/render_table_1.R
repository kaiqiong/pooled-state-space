# Main Table 1 from checked numeric summaries. No simulation is run.
# Rscript scripts/figures/render_table_1.R [simulation_summary_directory] [boundary_summary_csv] [output_directory]
script_arg <- grep("^--file=",commandArgs(),value=TRUE)
script_dir <- dirname(normalizePath(sub("^--file=","",script_arg)))
root <- normalizePath(file.path(script_dir,"../.."))
args <- commandArgs(trailingOnly=TRUE)
d <- if(length(args)>=1) args[1] else file.path(root,"reference/figures/simulation")
boundary_file <- if(length(args)>=2) args[2] else file.path(root,"reference/tables/boundary_coverage.csv")
out <- if(length(args)>=3) args[3] else file.path(root,"results/tables")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
a <- read.csv(file.path(d,"benchmark.csv"))
b <- read.csv(file.path(d,"information_ablation.csv"))
c <- read.csv(file.path(d,"clustered_free_variance.csv"))
z <- read.csv(boundary_file)
stopifnot(nrow(a)==36L,all(b$n_total==500),all(b$n_valid==500),nrow(c)==48L,all(c$n==200),nrow(z)==12L,all(z$n_boot==1000L),all(z$B_boot==999L))
endrow <- " \\\\"
range3 <- function(x) {stopifnot(length(x)==6L,all(is.finite(x)));sprintf("%.3f--%.3f",min(x),max(x))}
range_bias <- function(x) {
 stopifnot(length(x)==6L,all(is.finite(x)))
 limits<-round(range(x),4);limits[limits==0]<-0
 values<-sprintf("%.4f",limits)
 values[limits<0]<-paste0("$",values[limits<0],"$")
 paste(values,collapse="--")
}
parameters<-c("Intercept","Age","Sex","BMI","Prediabetic","Diabetic")
methods<-c("model","cr0","md","bootstrap")
stopifnot(all(c("target","mean_est","bias_vs_target","emp_sd","n_fit")%in%names(z)),
          all(c("target","mean_estimate","bias","empirical_sd","empirical_sd_same200")%in%names(c)))
estimation_rows<-list()
estimation<-function(r,t) {
 if(r==0) {
  x<-z[z$tau==t,];stopifnot(nrow(x)==6L,setequal(x$parameter,parameters),all(x$n_fit==1000L))
  x<-x[match(parameters,x$parameter),]
  e<-data.frame(sigma_eps2=r,tau=t,M=x$n_fit,parameter=x$parameter,
                target=x$target,mean_estimate=x$mean_est,bias=x$bias_vs_target,
                empirical_sd=x$emp_sd)
 } else {
  all<-c[c$tau==t,];stopifnot(nrow(all)==24L,setequal(all$method,methods))
  x<-all[all$method=="model",];stopifnot(nrow(x)==6L,setequal(x$parameter,parameters))
  x<-x[match(parameters,x$parameter),]
  # The four uncertainty procedures share the same estimates and matched cohort.
  shared<-c("n","target","mean_estimate","bias","empirical_sd","empirical_sd_same200")
  for(m in methods) {
   other<-all[all$method==m,];stopifnot(nrow(other)==6L,setequal(other$parameter,parameters))
   other<-other[match(parameters,other$parameter),]
   for(v in shared)stopifnot(isTRUE(all.equal(x[[v]],other[[v]],tolerance=1e-12)))
  }
  stopifnot(all(x$n==200L),isTRUE(all.equal(x$empirical_sd,x$empirical_sd_same200,tolerance=1e-12)))
  e<-data.frame(sigma_eps2=r,tau=t,M=x$n,parameter=x$parameter,
                target=x$target,mean_estimate=x$mean_estimate,bias=x$bias,
                empirical_sd=x$empirical_sd)
 }
 stopifnot(all(is.finite(as.matrix(e[,c("M","target","mean_estimate","bias","empirical_sd")]))),
           all(e$empirical_sd>=0),max(abs(e$bias-(e$mean_estimate-e$target)))<1e-12)
 # empirical_sd uses the M-1 sample-variance denominator; empirical MSE uses M.
 e$rmse<-sqrt(e$bias^2+(e$M-1)/e$M*e$empirical_sd^2)
 e
}
pa <- character()
for(p in c("Intercept","Age","Sex","Biomarker")) {
 v<-unlist(lapply(c(.25,1,4),function(r) lapply(c("joint_EM","naive_two_stage","per_subject_kalman"),function(m){x<-a[a$parameter==p&a$sigma_eps2==r&a$method==m,];stopifnot(nrow(x)==1L);x})),recursive=FALSE)
 rm<-vapply(v,function(x)x$rmse,numeric(1));rb<-vapply(v,function(x)x$rel_bias_pct,numeric(1))
 pa<-c(pa,paste0(p," & ",paste(sprintf("%.3f",rm),collapse=" & "),endrow,"[-1.5pt]"),paste0(" & ",paste(sprintf("$(%.1f)$",rb),collapse=" & "),endrow))
 if(p!="Biomarker")pa<-c(pa,"\\addlinespace[1.5pt]")
}
pb<-character();labels<-c(full="Full information",no_correction="Correction $B_i$ omitted",nuisance_fixed="Nuisance uncertainty ignored")
for(t in c(0,.4))for(v in names(labels)){
 cells<-c(if(v=="full")format(t)else "",labels[[v]])
 for(m in c("model","cr0")){x<-b[b$tau==t&b$variant==v&b$method==m,];cells<-c(cells,range3(x$ratio_to_emp_sd),range3(x$coverage))}
 pb<-c(pb,paste0(paste(cells,collapse=" & "),endrow))
 if(t==0&&v=="nuisance_fixed")pb<-c(pb,"\\midrule")
}
pc<-character()
for(r in c(0,400))for(t in c(0,.4)){
 e<-estimation(r,t);estimation_rows[[length(estimation_rows)+1L]]<-e
 cells<-c(r,t,if(r==0)1000 else 200,range_bias(e$bias),range3(e$rmse))
 for(m in methods){
  vals<-if(r==0)z[z$tau==t,paste0("cov_",if(m=="bootstrap")"boot"else m)] else c$coverage[c$tau==t&c$method==m]
  cells<-c(cells,range3(vals))
 }
 pc<-c(pc,paste0(paste(cells,collapse=" & "),endrow))
}
layout<-readLines(file.path(script_dir,"table_1_layout.tex"),warn=FALSE)
inject<-function(lines,marker,rows){i<-match(marker,lines);stopifnot(!is.na(i));c(lines[seq_len(i-1)],rows,lines[seq.int(i+1,length(lines))])}
layout<-inject(layout,"@@PANEL_A@@",pa);layout<-inject(layout,"@@PANEL_B@@",pb);layout<-inject(layout,"@@PANEL_C@@",pc)
writeLines(layout,file.path(out,"table_1.tex"))
writeLines(pa,file.path(out,"table_1_A_rows.tex"));writeLines(pb,file.path(out,"table_1_B_rows.tex"));writeLines(pc,file.path(out,"table_1_C_rows.tex"))
write.csv(a,file.path(out,"table_1_A_data.csv"),row.names=FALSE)
write.csv(b,file.path(out,"table_1_B_data.csv"),row.names=FALSE)
write.csv(c,file.path(out,"table_1_C_positive_data.csv"),row.names=FALSE)
write.csv(z,file.path(out,"table_1_C_boundary_data.csv"),row.names=FALSE)
write.csv(do.call(rbind,estimation_rows),file.path(out,"table_1_C_estimation.csv"),row.names=FALSE)
cat("Generated Table 1 from numeric summaries in",out,"\n")
