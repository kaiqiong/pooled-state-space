# Independent dense-Gaussian checks and a short estimation/resampling smoke test.
# Run from repository root. Outputs are synthetic and cannot be manuscript results.
source('R/state_space_em.R')
source('R/free_variance_cluster_inference.R')
stopifnot(requireNamespace('numDeriv', quietly=TRUE))
out <- 'results/validation/free_information'
dir.create(out, recursive=TRUE, showWarnings=FALSE)
rows <- list()
check <- function(name, error, tolerance) {
  rows[[length(rows)+1L]] <<- data.frame(check=name,error=error,tolerance=tolerance,pass=is.finite(error)&&error<tolerance)
}
dense <- function(th, x, Tt) {
  p <- length(x); phi <- tanh(sum(x*th[seq_len(p)]))
  a <- th[p+1]; P <- th[p+2]; q <- th[p+3]; r <- th[p+4]
  v <- numeric(Tt); v[1] <- P
  for (t in 2:Tt) v[t] <- phi^2*v[t-1]+q
  V <- outer(seq_len(Tt),seq_len(Tt),function(i,j) phi^abs(i-j)*v[pmin(i,j)])+diag(r,Tt)
  list(mean=a*phi^(0:(Tt-1)), V=V)
}
ll_unit <- function(th, y, x) {
  z <- dense(th,x,length(y)); L <- chol(z$V); e <- y-z$mean
  -0.5*(length(y)*log(2*pi)+2*sum(log(diag(L)))+sum(backsolve(L,e,transpose=TRUE)^2))
}
fit_at <- function(th,y,X) {
  p<-ncol(X)
  list(beta=th[seq_len(p)],a1=th[p+1],P1=th[p+2],sigma_eta2=th[p+3],sigma_eps2=th[p+4],
       kf=kalman_run(y,tanh(as.vector(X%*%th[seq_len(p)])),0,th[p+1],th[p+2],th[p+3],th[p+4]))
}
relerr <- function(a,b) max(abs(a-b)/(1+abs(b)))
set.seed(20261013)
X <- cbind(Intercept=1,x=rnorm(9)); th <- c(.4,-.2,.8,1.3,.7,.5)
y <- matrix(rnorm(9*7),9,7); fit <- fit_at(th,y,X)
Ii <- info_per_unit_free(y,X,0,fit); Si <- scores_per_unit_free(y,X,0,fit)
for (i in seq_len(nrow(y))) {
  ll <- function(z) ll_unit(z,y[i,],X[i,])
  check(paste0('episode_',i,'_Hessian'),relerr(Ii[i,,],-numDeriv::hessian(ll,th)),2e-5)
  check(paste0('episode_',i,'_score'),relerr(Si[i,],numDeriv::grad(ll,th)),1e-6)
}
check('information_symmetry',max(abs(Ii-aperm(Ii,c(1,3,2)))),1e-9)
N <- 120L; Tt <- 14L; X <- cbind(Intercept=1,x=rnorm(N)); cluster <- rep(seq_len(30),each=4)
y <- t(vapply(seq_len(N),function(i){z<-dense(th,X[i,],Tt);z$mean+as.vector(t(chol(z$V))%*%rnorm(Tt))},numeric(Tt)))
init <- list(beta=c(.3,-.1),a1=.5,P1=1,sigma_eta2=.6,sigma_eps2=.6)
fit <- em_regression_mv(y,X,0,init,thres=1e-7,max_iter=5000)
check('EM_converged',as.numeric(!fit$converged),.5)
thfit <- c(fit$beta,fit$a1,fit$P1,fit$sigma_eta2,fit$sigma_eps2)
ll <- function(z) sum(vapply(seq_len(N),function(i)ll_unit(z,y[i,],X[i,]),numeric(1)))
check('innovation_vs_dense_loglik',abs(fit$loglik-ll(thfit)),1e-8)
ci <- cluster_inference_free(y,X,0,fit,cluster)
check('full_information_vs_dense_Hessian',relerr(ci$I_obs,-numDeriv::hessian(ll,thfit)),3e-5)
check('finite_model_CR0_MD_SEs',as.numeric(any(!is.finite(c(ci$se_model,ci$se_cr0,ci$se_md)))),.5)
set.seed(20261014)
boots <- matrix(NA_real_,8,ncol(X))
for(b in seq_len(nrow(boots))) {
  picked <- sample(unique(cluster),length(unique(cluster)),replace=TRUE)
  ix <- unlist(lapply(picked,function(g)which(cluster==g)),use.names=FALSE)
  f <- em_regression_mv(y[ix,],X[ix,],0,fit,thres=1e-6,max_iter=5000)
  check(paste0('bootstrap_',b,'_converged'),as.numeric(!f$converged),.5)
  boots[b,] <- f$beta
}
write.csv(boots,file.path(out,'synthetic_bootstrap_estimates.csv'),row.names=FALSE)
checks <- do.call(rbind,rows)
write.csv(checks,file.path(out,'checks.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
print(checks,row.names=FALSE)
stopifnot(all(checks$pass))
cat('PASS:',nrow(checks),'independent free-variance and smoke checks.\n')
