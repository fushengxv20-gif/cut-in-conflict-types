invisible(Sys.setlocale('LC_ALL','English_United States.utf8'))
root<-dirname(normalizePath(sub('^--file=', '', grep('^--file=', commandArgs(), value=TRUE)[1]), winslash='/'))
library(jsonlite);library(poLCA)
js<-function(x,n)write_json(x,file.path(root,paste0(n,'.json')),auto_unbox=TRUE,pretty=TRUE,digits=12,na='null')
primary<-readRDS(file.path(root,'lca_primary.rds'));fit<-read.csv(file.path(root,'tables','lca_fit_primary.csv'))
# The one-class model has a unique solution and was fitted from one start.
fit$starts[1]<-1;fit$best_replicates[1]<-1;write.csv(fit,file.path(root,'tables','lca_fit_primary.csv'),row.names=FALSE,na='')
f<-primary[[which.min(fit$BIC)]]
nf<-read.csv(file.path(root,'tables','lca_fit_no_emotion.csv'));nf<-nf[nf$K<=3,];nf$starts[1]<-1;nf$best_replicates[1]<-1
write.csv(nf,file.path(root,'tables','lca_fit_no_emotion.csv'),row.names=FALSE,na='')
js(list(no_emotion=nf),'lca_sensitivity')
d<-read.csv(file.path(root,'final_analysis_data.csv'));ld<-data.frame(injury=d$injury,penalty=d$penalty,collision=d$collision+1,emotion=d$emotion_recorded+1)
pairs<-combn(names(ld),2,simplify=FALSE)
out<-lapply(pairs,function(v){ok<-complete.cases(ld[,v]);obs<-table(factor(ld[ok,v[1]],levels=1:ncol(f$probs[[v[1]]])),factor(ld[ok,v[2]],levels=1:ncol(f$probs[[v[2]]])));pred<-t(f$probs[[v[1]]])%*%diag(f$P)%*%f$probs[[v[2]]]*sum(ok);res<-(obs-pred)/sqrt(pred);list(pair=paste(v,collapse=' x '),n=sum(ok),discrepancy=sum(res^2),max_abs=max(abs(res)))})
js(list(pairs=out,selected_K=which.min(fit$BIC),no_emotion_selected_K=nf$K[which.min(nf$BIC)],no_emotion_candidate_note='K=1-3 only. K>=4 has more parameters than the 23 independent complete response-pattern probabilities and is excluded from class enumeration.'),'lca_diagnostics')
# Independent likelihood recomputation from the reported category probabilities.
contrib<-sapply(1:nrow(ld),function(i){pr<-f$P;for(v in names(ld)){if(!is.na(ld[i,v]))pr<-pr*f$probs[[v]][,ld[i,v]]};log(sum(pr))})
stopifnot(abs(sum(contrib)-f$llik)<1e-5,all(sapply(f$probs,function(x)max(abs(rowSums(x)-1)))<1e-8))
cat('Independent observed-data likelihood:',sum(contrib),'\n')
