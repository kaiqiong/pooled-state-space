# Portable October 13 replication. Run from the repository root; see docs/SIMULATIONS.txt.
# Validated summaries of saved simulations; no model fitting and no writes to raw data.
# The strict final path requires every planned ID, B=999, matching full-fit estimates,
# finite retained draws and converged diagnostics. Constrained converged refits are retained.
# Regenerated caches are checked against their cohort, IDs, B and fit diagnostics.
atomic_csv <- function(x,path) { z<-tempfile(tmpdir=dirname(path)); write.csv(x,z,row.names=FALSE); if(!file.rename(z,path)) stop('Cannot install ',path) }
atomic_rds <- function(x,path) { z<-tempfile(tmpdir=dirname(path)); saveRDS(x,z); if(!file.rename(z,path)) stop('Cannot install ',path) }
run_summary <- function(design, root, out, mode) {
  if(length(mode)!=1L || !mode %in% c('final','monitor')) stop('mode must be final or monitor')
  imp<-root; tb<-file.path(imp,'results/simulation/tables')
  src<-file.path(tb,if(design=='posvar') 'part_q_raw.rds' else 'part_m_boundary_v2_raw.rds')
  target_file<-file.path(tb,'expected_likelihood_targets.rds'); allraw<-readRDS(src); raw<-allraw$raw
  tg<-readRDS(target_file)[[design]]; pn<-names(tg$generating$beta); p<-length(pn)
  planned<-if(design=='posvar') 1:200 else 1:1000; expected_main<-if(design=='posvar') 500L else 1000L
  scans<-lapply(c('0','0.4'),function(key) {
    fs<-list.files(file.path(imp,'results/simulation/boot_b999',paste0(design,'_tau',key)), '^m[0-9]+\\.rds$',full.names=TRUE)
    id<-as.integer(sub('^m([0-9]+)\\.rds$','\\1',basename(fs)))
    if(anyDuplicated(id)) stop('Duplicate cached dataset IDs')
    list(key=key,files=fs,id=id,complete=setequal(id,planned))
  })
  if(mode=='final' && !all(vapply(scans,`[[`,logical(1),'complete'))) {
    stop('Final export refused: planned IDs 1:',max(planned),' required in BOTH scenarios; current counts ',
         paste(vapply(scans,function(z) length(intersect(z$id,planned)),integer(1)),collapse='/'),
         '. Use monitor for completion counts; no final outputs written.')
  }
  if(mode=='monitor') {
    dir.create(out,recursive=TRUE,showWarnings=FALSE)
    mon<-do.call(rbind,lapply(scans,function(z)data.frame(design=design,tau=as.numeric(z$key),planned=length(planned),completed=length(intersect(z$id,planned)),missing_ids=paste(setdiff(planned,z$id),collapse=';'),unexpected_ids=paste(setdiff(z$id,planned),collapse=';'))))
    atomic_csv(mon,file.path(out,paste0(design,'_bootstrap_monitor.csv')));print(mon);return(invisible(mon))
  }
  rows<-list(); long<-list(); paired<-list(); diagnostics<-list(); comparisons<-list(); input_files<-c(src,target_file)
  ci_stat<-function(lo,hi,target) {
    u<-is.finite(lo)&is.finite(hi); if(!any(u))return(c(n=0,coverage=NA,coverage_mcse=NA,mean_length=NA,length_mcse=NA,below_target=NA,above_target=NA))
    cv<-mean(lo[u]<=target&target<=hi[u]); len<-hi[u]-lo[u]
    c(n=sum(u),coverage=cv,coverage_mcse=sqrt(cv*(1-cv)/sum(u)),mean_length=mean(len),length_mcse=sd(len)/sqrt(length(len)),below_target=mean(hi[u]<target),above_target=mean(lo[u]>target))
  }
  for(sc in scans) {
    key<-sc$key; tau<-as.numeric(key); r<-raw[[key]]
    # Match through saved ID, never through an index after filtering failed fits.
    ids<-vapply(r,function(z)as.integer(z$m),integer(1)); if(anyDuplicated(ids)||!setequal(ids,seq_len(expected_main)))stop('Main dataset IDs invalid')
    names(r)<-as.character(ids)
    ok<-vapply(r,function(z)if(design=='posvar')is.list(z)&&identical(z$err,FALSE) else is.list(z)&&identical(z$fit_error,''),logical(1))
    rr<-r[ok]; idok<-as.integer(names(rr)); est<-t(vapply(rr,function(z)as.numeric(z$beta),numeric(p)))
    if(any(!is.finite(est)))stop('Nonfinite main estimate')
    target<-if(tau==0) tg$generating$beta else tg$beta_star; beta0<-tg$generating$beta
    lo_b<-hi_b<-matrix(numeric(),0,p); se_b<-matrix(numeric(),0,p); bid<-integer(); bdiag<-list()
      bid<-planned; lo_b<-hi_b<-se_b<-matrix(NA_real_,length(bid),p); cap_fraction<-numeric(length(bid))
      for(k in seq_along(bid)) {
        m<-bid[k]; fn<-sc$files[match(m,sc$id)]; z<-readRDS(fn)
        required<-c('design','tau','m','B','beta_full','draws','diag','n_unresolved','refit_dev_from_saved')
        if(!all(required%in%names(z))||!identical(z$design,design)||z$tau!=tau||z$m!=m||z$B!=999L)stop('Cache metadata mismatch: ',fn)
        if(!as.character(m)%in%names(rr))stop('Bootstrap main fit failed or missing: ',fn)
        if(max(abs(z$beta_full-r[[as.character(m)]]$beta))>1e-6||z$refit_dev_from_saved>1e-6)stop('Full-fit mismatch: ',fn)
        if(!identical(dim(z$draws),c(999L,as.integer(p)))||!identical(colnames(z$draws),pn)||nrow(z$diag)!=999L||!identical(as.integer(z$diag$b),1:999))stop('Draw/diagnostic layout invalid: ',fn)
        # A nonzero unresolved count requires explicit review before final export;
        # never silently compute final coverage from varying accepted B.
        if(any(!is.finite(z$diag$iter))||any(z$diag$iter<1)||any(!is.finite(z$diag$loglik))||any(!is.finite(z$beta_full)))stop('Nonfinite iteration/full-fit diagnostics: ',fn)
        if(z$n_unresolved!=0L||any(!is.finite(z$draws))||anyNA(z$diag$converged)||any(!z$diag$converged)||anyNA(z$diag$status)||any(!z$diag$status%in%c('interior','constrained'))||anyNA(z$diag$attempt)||any(!z$diag$attempt%in%1:2))stop('Unresolved or invalid bootstrap diagnostics: ',fn)
        lo_b[k,]<-apply(z$draws,2,quantile,.025,type=7);hi_b[k,]<-apply(z$draws,2,quantile,.975,type=7);se_b[k,]<-apply(z$draws,2,sd)
        cap_fraction[k]<-mean(z$diag$status=='constrained')
        bdiag[[k]]<-data.frame(tau=tau,m=m,B=999L,constrained=sum(z$diag$status=='constrained'),second_attempt=sum(z$diag$attempt==2),unresolved=z$n_unresolved,max_iter=max(z$diag$iter),refit_deviation=z$refit_dev_from_saved)
      }
      input_files<-c(input_files,sc$files[match(planned,sc$id)]);diagnostics[[key]]<-do.call(rbind,bdiag)
    
    oldids<-idok[vapply(rr,function(z)!is.null(z$boot_lo),logical(1))]
    # Historical diagnostic columns are retained for CSV compatibility only.
    # They never contribute to the final intervals or trigger pilot fitting.
    for(j in seq_len(p)) {
      row<-data.frame(tau=tau,parameter=pn[j],generating=beta0[j],target=target[j],target_working=target[j],n_fit=nrow(est),mean_est=mean(est[,j]),bias_vs_target=mean(est[,j])-target[j],bias_vs_working=mean(est[,j])-target[j],bias=mean(est[,j])-target[j],bias_mc_se=sd(est[,j])/sqrt(nrow(est)),bias_vs_generating=mean(est[,j])-beta0[j],emp_sd=sd(est[,j]),export_state='complete_cohort_B999',numerical_review=if(design=='boundary')'see_separate_multistart_and_cap_diagnostics'else'not_certified_by_exporter')
      hits<-list(); masks<-list(); mainlimits<-list()
      for(method in c('model','cr0','md')) {
        s<-vapply(rr,function(z)as.numeric(z$se[[method]][j]),numeric(1)); u<-is.finite(s)&s>0;crit<-if(method=='model')qnorm(.975)else qt(.975,43)
        lo<-est[,j]-crit*s;hi<-est[,j]+crit*s;lo[!u]<-hi[!u]<-NA; a<-ci_stat(lo,hi,target[j]); ag<-ci_stat(lo,hi,beta0[j])
        row[[paste0('n_',method)]]<-sum(u);row[[paste0('se_',method)]]<-mean(s[u]);row[[paste0('se_ratio_',method)]]<-mean(s[u])/sd(est[u,j]);row[[paste0('cov_',method)]]<-a['coverage'];row[[paste0('cov_',method,'_working')]]<-a['coverage'];row[[paste0('cov_',method,'_generating')]]<-ag['coverage']
        for(nm in c('coverage_mcse','mean_length','length_mcse','below_target','above_target'))row[[paste0(method,'_',nm)]]<-a[nm]
        long[[length(long)+1]]<-data.frame(design=design,tau=tau,parameter=pn[j],method=method,sample='all_main',target_type='working',as.list(a),mean_se=mean(s[u]),se_empirical_sd_ratio=mean(s[u])/sd(est[u,j]))
        mainlimits[[method]]<-list(lo=lo,hi=hi,s=s)
      }
      row$n_boot<-length(bid);row$B_boot<-if(length(bid))999L else NA_integer_
      row$se_boot<-row$cov_boot<-row$cov_boot_working<-row$cov_boot_generating<-NA_real_
      if(length(bid)) {
        a<-ci_stat(lo_b[,j],hi_b[,j],target[j]);ag<-ci_stat(lo_b[,j],hi_b[,j],beta0[j]);ix<-match(bid,idok)
        row$se_boot<-mean(se_b[,j]);row$cov_boot<-row$cov_boot_working<-a['coverage'];row$cov_boot_generating<-ag['coverage'];row$se_ratio_boot<-mean(se_b[,j])/sd(est[ix,j])
        for(nm in c('coverage_mcse','mean_length','length_mcse','below_target','above_target'))row[[paste0('boot_',nm)]]<-a[nm]
        long[[length(long)+1]]<-data.frame(design=design,tau=tau,parameter=pn[j],method='bootstrap',sample='planned_bootstrap',target_type='working',as.list(a),mean_se=mean(se_b[,j]),se_empirical_sd_ratio=mean(se_b[,j])/sd(est[ix,j]))
        limits<-lapply(mainlimits,function(z)list(lo=z$lo[ix],hi=z$hi[ix],s=z$s[ix]));limits$bootstrap<-list(lo=lo_b[,j],hi=hi_b[,j],s=se_b[,j])
        for(method in names(limits)) {
          v<-limits[[method]]; a<-ci_stat(v$lo,v$hi,target[j]);u<-is.finite(v$lo)&is.finite(v$hi)
          paired[[length(paired)+1]]<-data.frame(design=design,tau=tau,parameter=pn[j],method=method,sample=paste0('same_',length(bid),'_IDs'),target_type='working',as.list(a),mean_se=mean(v$s[u]),se_empirical_sd_ratio=mean(v$s[u])/sd(est[ix[u],j]))
        }
        for(method in c('model','cr0','md')) {
          v<-limits[[method]]; u<-is.finite(v$lo)&is.finite(v$hi);delta<-as.numeric(lo_b[u,j]<=target[j]&target[j]<=hi_b[u,j])-as.numeric(v$lo[u]<=target[j]&target[j]<=v$hi[u])
          comparisons[[length(comparisons)+1]]<-data.frame(tau=tau,parameter=pn[j],comparison=paste0('bootstrap_minus_',method),n=sum(u),coverage_difference=mean(delta),paired_mcse=sd(delta)/sqrt(length(delta)))
        }
      }
      row$n_boot_b99<-length(oldids);row$cov_boot_b99<-NA_real_
      if(length(oldids)) { old<-r[as.character(oldids)]; row$cov_boot_b99<-mean(vapply(old,function(z)z$boot_lo[j]<=target[j]&&target[j]<=z$boot_hi[j],logical(1))) }
      row$fit_errors<-sum(!ok);row$not_converged<-sum(!vapply(rr,function(z)if(design=='posvar')z$conv else z$converged,logical(1)))
      row$constrained<-sum(vapply(rr,function(z)z$status=='constrained',logical(1)));row$median_iter<-median(vapply(rr,function(z)as.numeric(z$iter),numeric(1)))
      rows[[length(rows)+1]]<-row
    }
  }
  df<-do.call(rbind,rows);rownames(df)<-NULL
  prefix<-if(design=='posvar')'part_q_summary_final'else'part_m_boundary_v2_working_summary'
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  atomic_csv(df,file.path(out,paste0(prefix,'.csv')));atomic_csv(do.call(rbind,long),file.path(out,paste0(prefix,'_intervals.csv')))
  if(length(paired))atomic_csv(do.call(rbind,paired),file.path(out,paste0(prefix,'_paired.csv')))
  if(length(comparisons))atomic_csv(do.call(rbind,comparisons),file.path(out,paste0(prefix,'_paired_differences.csv')))
  if(length(diagnostics))atomic_csv(do.call(rbind,diagnostics),file.path(out,paste0(prefix,'_diagnostics.csv')))
  manifest<-data.frame(path=input_files,md5=unname(tools::md5sum(input_files)))
  atomic_csv(manifest,file.path(out,paste0(prefix,'_inputs.csv')))
  meta<-list(cohort_complete=TRUE,numerical_review=if(design=='boundary')'Original fitted intervals retained; see the separately regenerated starting-value and cap diagnostics'else'export validates cache completeness and recorded diagnostics, not global optimality or interval validity',summary_file=paste0(prefix,'.csv'),summary_md5=unname(tools::md5sum(file.path(out,paste0(prefix,'.csv')))),source_checksums=tools::md5sum(file.path(imp,c('scripts/simulation/cluster_bootstrap.R','R/state_space_em.R','scripts/simulation/working_targets.R'))),design=design,state=mode,planned_ids=planned,bootstrap_B=999L,bootstrap_ids=planned,scenarios=c(0,.4),main_M=expected_main,zero_unresolved_required=TRUE,quantile_type=7L,generated_at=format(Sys.time(),tz='UTC',usetz=TRUE),input_manifest=manifest,configuration_hash_available_in_legacy_caches=FALSE)
  atomic_rds(meta,file.path(out,paste0(prefix,'_metadata.rds')))
  print(df[,c('tau','parameter','target','bias_vs_target','cov_model','cov_cr0','cov_md','cov_boot','n_boot','cov_boot_b99')],row.names=FALSE)
  invisible(df)
}
