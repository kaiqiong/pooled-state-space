# Regenerate numerical views of Supplementary Tables S1--S10 from summary CSVs.
# CSVs retain full precision; TeX tabular fragments are reviewable numerical views.
# Rscript scripts/figures/render_supplement_tables.R [summary_directory] [output_directory]
script_arg<-grep("^--file=",commandArgs(),value=TRUE)
root<-normalizePath(file.path(dirname(sub("^--file=","",script_arg)),"../.."))
args<-commandArgs(trailingOnly=TRUE)
d<-if(length(args)>=1)args[1]else file.path(root,"reference/simulation/tables")
out<-if(length(args)>=2)args[2]else file.path(root,"results/tables")
mode<-if(length(args)>=3)args[3]else "reference"
stopifnot(mode%in%c("reference","recomputed"))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
read<-function(n)read.csv(file.path(d,n),check.names=FALSE)
num<-function(x,k=3)formatC(round(x,k)+0,format="f",digits=k)
escape<-function(x){x<-gsub("_","\\_",as.character(x),fixed=TRUE);x<-gsub("%","\\%",x,fixed=TRUE);x}
save<-function(x,stem){
 write.csv(x,file.path(out,paste0(stem,".csv")),row.names=FALSE)
 header<-paste(escape(names(x)),collapse=" & ")
 rows<-apply(x,1,function(z)paste(escape(z),collapse=" & "))
 tex<-c(paste0("\\begin{tabular}{",paste(rep("l",ncol(x)),collapse=""),"}"),"\\toprule",paste0(header," \\\\"),"\\midrule",paste0(rows," \\\\"),"\\bottomrule","\\end{tabular}")
 writeLines(tex,file.path(out,paste0(stem,".tex")))
}
reg<-c("Intercept","Age","Sex","BMI","Prediabetic","Diabetic")
methods<-c("model","cr0","md","bootstrap")
# S1: use actual individual-parameter differences, rather than copied maxima.
a<-read("part_a_se_analytic_vs_numeric.csv")
v<-do.call(rbind,lapply(unique(a$model),function(m){x<-a[a$model==m,];data.frame(model=m,parameters=nrow(x),max_abs_se_difference=max(x$abs_diff))}))
for(pair in list(c("four_covariates","part_e1_mv_se_analytic_vs_numeric.csv"),c("interaction","part_h1_interaction_se_analytic_vs_numeric.csv"))){x<-read(pair[2]);v<-rbind(v,data.frame(model=pair[1],parameters=nrow(x),max_abs_se_difference=max(x$abs_diff)))}
if(file.exists(file.path(d,"fixed_variance_se_analytic_vs_numeric.csv"))){x<-read("fixed_variance_se_analytic_vs_numeric.csv");v<-rbind(v,data.frame(model="fixed_observation_variance",parameters=nrow(x),max_abs_se_difference=max(x$abs_diff)))}
save(v,"table_S1_validation_checks")
if(mode=="reference") {
 published<-read.csv(file.path(root,"reference/tables/table_S1_published.csv"),check.names=FALSE)
 save(published,"table_S1_validation")
} else save(v,"table_S1_validation")
# S2.
e<-read("part_e2_mv_coverage_summary.csv");save(e,"table_S2_base_coverage")
# S3: matched free-variance cohort, all methods on identical 200 datasets.
p<-read("paired_posvar200_free_variance_comparison.csv")
stopifnot(nrow(p)==48L,all(p$n==200L))
rowsA<-rowsB<-list();j<-0L
for(t in c(0,.4))for(k in reg){j<-j+1L;x<-p[p$tau==t&p$parameter==k,];x<-x[match(methods,x$method),];stopifnot(nrow(x)==4L,!anyNA(x$method))
 rowsA[[j]]<-data.frame(tau=t,coefficient=k,model=x$coverage[1],CR0=x$coverage[2],MD=x$coverage[3],bootstrap=x$coverage[4],bootstrap_minus_CR0=x$bootstrap_minus_coverage[2],paired_MCSE_CR0=x$paired_mcse[2],bootstrap_minus_MD=x$bootstrap_minus_coverage[3],paired_MCSE_MD=x$paired_mcse[3])
 rowsB[[j]]<-data.frame(tau=t,coefficient=k,model_length=x$mean_length[1],CR0_length=x$mean_length[2],MD_length=x$mean_length[3],bootstrap_length=x$mean_length[4],model_SE_ratio=x$se_empirical_sd_ratio[1],CR0_SE_ratio=x$se_empirical_sd_ratio[2],MD_SE_ratio=x$se_empirical_sd_ratio[3],bootstrap_SE_ratio=x$se_empirical_sd_ratio[4])
}
save(do.call(rbind,rowsA),"table_S3_A_matched_coverage");save(do.call(rbind,rowsB),"table_S3_B_matched_lengths")
# S4: each specification evaluated against its own working target.
f<-read("free_fixed_matched200_comparison.csv");rows<-list();j<-0L
for(t in c(0,.4))for(m in methods){j<-j+1L;x<-f[f$tau==t&f$method==m,];stopifnot(nrow(x)==6L);rows[[j]]<-data.frame(tau=t,method=m,free_min=min(x$coverage_free),free_max=max(x$coverage_free),fixed_min=min(x$coverage_fixed),fixed_max=max(x$coverage_fixed),max_absolute_coefficient_change=max(abs(x$coverage_free-x$coverage_fixed)))}
save(do.call(rbind,rows),"table_S4_variance_comparison")
# S5: common 1000-dataset boundary cohort. All original fitted estimates retained.
z<-read("part_m_boundary_v2_working_summary.csv")
stopifnot(nrow(z)==12L,all(z$n_fit==1000L),all(z$n_boot==1000L),all(z$B_boot==999L))
save(z[,c("tau","parameter","bias_vs_working","emp_sd","cov_model_working","cov_cr0_working","cov_md_working","cov_boot_working")],"table_S5_A_boundary_coverage")
cols<-c("tau","parameter","emp_sd",unlist(lapply(c("model","cr0","md","boot"),function(m)c(paste0(m,"_mean_length"),paste0("se_ratio_",m)))))
save(z[,cols],"table_S5_B_boundary_lengths")
# S6: comparison against CONDITIONAL generating coefficients for both methods.
x<-read("part_p_mixed_comparator_v2_summary.csv")
save(x[,c("tau","parameter","wi_bias","wi_sd","wi_cov_cr0","wi_cov_md","mx_bias","mx_sd","mx_se","mx_cov","rmse_ratio_mx_over_wi")],"table_S6_mixed_effects")
# S7: retain scenario labels N1--N6; definitions are in OUTPUTS.txt.
x<-read("part_n_final_summary.csv");x<-x[x$parameter%in%c("Prediabetic","Diabetic"),]
save(x[,c("scenario","parameter","truth","joint_mean","joint_sd","naive_mean","naive_sd","joint_reject","naive_reject","joint_cover","naive_cover")],"table_S7_differential_attenuation")
# S8.
x<-read("part_l_heavytail_summary.csv");e<-e[match(x$parameter,e$parameter),]
save(data.frame(parameter=x$parameter,truth=x$truth,Gaussian_bias=e$mean_est-e$truth,Gaussian_coverage=e$coverage_95,t3_bias=x$bias,t3_RMSE=x$rmse,t3_SE_ratio=x$mean_se/x$emp_sd,t3_coverage=x$coverage),"table_S8_heavy_tails")
# S9: elapsed times are hardware-dependent, unlike the statistical definitions.
x<-read("part_o_direct_optimization_v2.csv")
save(x[,c("design","M","em_median_iter","em_cap","em_median_secs","same_optimum","max_abs_dll","direct_cap","direct_median_secs","capped_same_optimum","capped_nonconv","capped_median_secs")],"table_S9_optimization")
# S10: validated best converged feasible cross-cap fits, all six coefficients.
x<-read("cap_sensitivity_simulation_ms.csv");ap<-read("cap_sensitivity_application_ms.csv")
rows<-list();j<-0L
for(s in unique(x$set))for(t in sort(unique(x$tau[x$set==s])))for(cap in c(.99,.995,.999)){
 q<-x[x$set==s&x$tau==t&x$cap==cap,];if(!nrow(q))next;j<-j+1L;stopifnot(nrow(q)==6L);di<-q[q$parameter=="Diabetic",]
 rows[[j]]<-data.frame(set=s,tau=t,cap=cap,n=di$n,original_fits_at_cap=di$cap_active_fits,bootstrap_refits_at_cap_percent=100*di$boot_constrained_frac,diabetes_median_bootstrap_length=di$medlen_boot,maximum_endpoint_change=max(q$max_abs_endpoint_shift_vs_995),CR0_min=min(q$cov_cr0),CR0_max=max(q$cov_cr0),bootstrap_min=min(q$cov_boot),bootstrap_max=max(q$cov_boot))
}
save(do.call(rbind,rows),"table_S10_cap_simulations")
save(ap[ap$parameter=="Diabetic",c("cap","estimate","loglik","status","max_phi","boot_lo","boot_hi","boot_constrained","boot_nonconv")],"table_S10_cap_application")
cat("Regenerated supplemental numerical table views in",out,"\n")
