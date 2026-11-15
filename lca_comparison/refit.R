invisible(Sys.setlocale('LC_ALL','English_United States.utf8'))
source('setup.R',encoding='UTF-8')
rows<-list();models<-list()
for(k in 3:4){
 cat('Independent starts for K=',k,'\n');flush.console()
 f<-fit(ld,k,nrep=200,seed=20261106+k)
 stopifnot(abs(sum(ll(f,ld))-f$llik)<1e-6,abs(f$aic-(-2*f$llik+2*f$npar))<1e-8)
 models[[as.character(k)]]<-f
 rows[[as.character(k)]]<-data.frame(K=k,LL=f$llik,q=f$npar,AIC=f$aic,BIC=f$bic,aBIC=-2*f$llik+f$npar*log((nrow(ld)+2)/24),best_replicates=sum(abs(f$attempts-f$llik)<1e-4),starts=length(f$attempts),iterations=f$numiter,min_class=min(f$P),entropy=1+sum(f$posterior*log(pmax(f$posterior,1e-300)))/(nrow(ld)*log(k)))
 write.csv(do.call(rbind,rows),file.path(out,'independent_refits.csv'),row.names=FALSE)
 saveRDS(models,file.path(out,'independent_refits.rds'))
 cat('K',k,'LL',format(f$llik,digits=15),'best repeats',rows[[as.character(k)]]$best_replicates,'\n');flush.console()
}
