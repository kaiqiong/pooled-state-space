source("R/simulation_design.R")
# Portable October 13 replication. Run from the repository root; see docs/SIMULATIONS.txt.
## Bounded derivative verification AT STORED FITS: no EM or optimization.
args <- commandArgs(trailingOnly=TRUE)
root <- getwd()
out <- file.path(root, 'results/simulation/free_variance_ablation/validation')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
setwd(root)
source('R/state_space_em.R')
source('R/free_variance_cluster_inference.R')
stopifnot(requireNamespace('numDeriv',quietly=TRUE))
f <- pss_simulation_design("positive")
Xmat<-f$Xmat;subj<-f$meta$subject;N<-nrow(Xmat);Tt<-ncol(f$Y);pn<-colnames(Xmat);p<-ncol(Xmat)
beta0<-setNames(f$truth$beta,pn);a1_0<-f$truth$a1;P1_0<-f$truth$P1;se2_0<-f$truth$sigma_eta2;r0<-400
## Literal data-generation operations from sim_clustered_freevar.R.
gen<-function(tau,X=Xmat,cl=subj){ug<-unique(cl);b<-setNames(rnorm(length(ug),0,tau),ug);phi<-tanh(as.vector(X%*%beta0)+b[cl]);n<-nrow(X);al<-rnorm(n,a1_0,sqrt(P1_0));y<-matrix(0,n,Tt);y[,1]<-al+rnorm(n,0,sqrt(r0));for(t in 2:Tt){al<-phi*al+rnorm(n,0,sqrt(se2_0));y[,t]<-al+rnorm(n,0,sqrt(r0))};y}
loglik_fun<-function(par,y){phi<-tanh(as.vector(Xmat%*%par[1:p]));kf<-kalman_run(y,phi,0,par[p+1],par[p+2],par[p+3],par[p+4]);Fv<-kf$P_pred+par[p+4];e<-y-kf$alpha_pred;-.5*sum(log(2*pi*Fv)+e^2/Fv)}
rows<-list();score_rows<-list();gradient_sensitivity<-list();t0<-Sys.time()
for(tau in c(0,.4)){
 fit<-readRDS(sprintf('results/simulation/freevar/tau%s/m001.rds',tau));stopifnot(fit$m==1L,fit$tau==tau,isTRUE(fit$conv))
 set.seed(410000+10000*round(10*tau)+1);y<-gen(tau)
 fit$kf<-kalman_run(y,tanh(as.vector(Xmat%*%fit$beta)),0,fit$a1,fit$P1,fit$sigma_eta2,fit$sigma_eps2)
 par<-setNames(c(fit$beta,fit$a1,fit$P1,fit$sigma_eta2,fit$sigma_eps2),c(pn,'a1','P1','sigma_eta2','sigma_eps2'))
 ci<-cluster_inference_free(y,Xmat,0,fit,subj)
 H<--numDeriv::hessian(loglik_fun,par,y=y)
 g<-numDeriv::grad(loglik_fun,par,y=y)
 for(step in c(.0001,.001,.01,.05)){gs<-numDeriv::grad(loglik_fun,par,y=y,method.args=list(d=step,r=6));gradient_sensitivity[[length(gradient_sensitivity)+1]]<-data.frame(tau=tau,relative_step=step,richardson_levels=6,gradient_max_abs_difference=max(abs(ci$score_sum-gs)),parameter=names(par)[which.max(abs(ci$score_sum-gs))])}
 gr<-numDeriv::grad(loglik_fun,par,y=y,method.args=list(d=.01,r=6))
 pooled<-compute_SE_regression_mv(y,Xmat,0,fit)$se
 sn<-sqrt(diag(solve(H)))
 vals<-data.frame(tau=tau,m=1,seed=410000+10000*round(10*tau)+1,original_iterations=fit$iter,reconstructed_loglik=loglik_fun(par,y),stored_loglik=fit$loglik,loglik_abs_difference=abs(loglik_fun(par,y)-fit$loglik),information_frobenius_relative=norm(ci$I_obs-H,'F')/norm(H,'F'),information_max_abs_difference=max(abs(ci$I_obs-H)),information_asymmetry=ci$asymmetry,analytic_min_eigenvalue=ci$min_eig,numerical_min_eigenvalue=min(eigen(H,symmetric=TRUE,only.values=TRUE)$values),gradient_max_abs_difference=max(abs(ci$score_sum-g)),gradient_refined_d01_r6_max_abs_difference=max(abs(ci$score_sum-gr)),analytic_score_max_abs=max(abs(ci$score_sum)),numerical_score_max_abs=max(abs(g)),model_se_vs_numerical_max_relative_all10=max(abs(ci$se_model-sn)/sn),model_se_vs_pooled_max_relative_beta=max(abs(ci$se_model[pn]-pooled[pn])/pooled[pn]),model_se_vs_saved_max_relative_all10=max(abs(ci$se_model-fit$se$model)/fit$se$model),cr0_se_vs_saved_max_relative_all10=max(abs(ci$se_cr0-fit$se$cr0)/fit$se$cr0),md_se_vs_saved_max_relative_all10=max(abs(ci$se_md-fit$se$md)/fit$se$md))
 stopifnot(all(is.finite(unlist(vals))),vals$loglik_abs_difference<1e-8)
 rows[[length(rows)+1]]<-vals
 score_rows[[length(score_rows)+1]]<-data.frame(tau=tau,parameter=names(par),analytic_score=ci$score_sum,numerical_gradient=g,absolute_difference=abs(ci$score_sum-g),numerical_gradient_refined=gr,absolute_difference_refined=abs(ci$score_sum-gr),analytic_model_se=ci$se_model,numerical_model_se=sn)
}
write.csv(do.call(rbind,gradient_sensitivity),file.path(out,'gradient_step_sensitivity.csv'),row.names=FALSE)
write.csv(do.call(rbind,rows),file.path(out,'saved_fit_derivative_checks.csv'),row.names=FALSE)
write.csv(do.call(rbind,score_rows),file.path(out,'saved_fit_derivative_components.csv'),row.names=FALSE)
capture.output({cat('Checks at original cached m001 parameters; no EM or optimization.\nElapsed seconds:',as.numeric(Sys.time()-t0,units='secs'),'\n');print(do.call(rbind,rows),digits=10);sessionInfo()},file=file.path(out,'saved_fit_derivative_checks.txt'))
cat('Saved-fit derivative check completed in',as.numeric(Sys.time()-t0,units='secs'),'seconds\n')
print(do.call(rbind,rows),digits=7)
