invisible(Sys.setlocale('LC_ALL','English_United States.utf8'))
root <- dirname(normalizePath(sub('^--file=', '', grep('^--file=', commandArgs(), value=TRUE)[1]), winslash='/'))
library(brglm2); library(poLCA); library(jsonlite)
dir.create(file.path(root,'tables'),showWarnings=FALSE)
d <- read.csv(file.path(root,'case_data.csv'))
stopifnot(nrow(d)==1003,!anyDuplicated(d$case_id),all(d$injury %in% c(1:4,NA)))
# Collision: supplementary flag, or a principal process coded as vehicle or chain collision.
d$collision <- ifelse(d$collision_record %in% 1 | d$process %in% c(2,5),1L,
                      ifelse(is.na(d$process),NA_integer_,0L))
stopifnot(!anyNA(d$emotion), all(d$emotion %in% c(0, 1)))
d$emotion_recorded <- as.integer(d$emotion)
d$sex <- factor(d$sex1,1:2,c('Male','Female'))
d$age_group <- factor(d$age1,1:4,c('G1','G2','G3','G4'))
d$time_group <- relevel(factor(d$time,1:3,c('Night','Day','Evening')),'Day')
d$road_group <- factor(d$road,1:5,c('Urban','Expressway','National','Residential','Rural'))
d$vehicle_group <- factor(d$vehicle1,1:4,c('Car','TwoWheeler','BusTruck','Special'))
d$type <- factor(d$conflict,1:5,paste0('T',1:5))
d$base_sample <- complete.cases(d[,c('sex1','age1')])
d$adjusted_sample <- complete.cases(d[,c('sex1','age1','time','road','vehicle1')])
write.csv(d,file.path(root,'final_analysis_data.csv'),row.names=FALSE,na='')
js <- function(x,name) write_json(x,file.path(root,paste0(name,'.json')),pretty=TRUE,auto_unbox=TRUE,digits=12,na='null')
tabout <- function(x,name) write.csv(x,file.path(root,'tables',paste0(name,'.csv')),row.names=FALSE,na='')
wilson <- function(x,n) {z<-qnorm(.975);p<-x/n;den<-1+z*z/n;mid<-(p+z*z/(2*n))/den;h<-z*sqrt(p*(1-p)/n+z*z/(4*n*n))/den;c(p=p,lo=mid-h,hi=mid+h)}
inj <- table(d$type,d$injury)
injtab <- do.call(rbind,lapply(1:5,function(t){n<-sum(inj[t,]);data.frame(type=paste0('T',t),n=n,no=inj[t,1],minor=inj[t,2],serious=inj[t,3],death=inj[t,4],mean=sum(inj[t,]*(1:4))/n,t(wilson(inj[t,3],n)),t(wilson(inj[t,4],n)))}))
names(injtab)[8:13]<-c('serious_p','serious_lo','serious_hi','death_p','death_lo','death_hi')
tabout(injtab,'injury')
chi <- function(tab,seed=20261019) {
 a<-suppressWarnings(chisq.test(tab));set.seed(seed);mc<-chisq.test(tab,simulate.p.value=TRUE,B=100000)
 list(n=sum(tab),chi2=unname(a$statistic),df=unname(a$parameter),p=a$p.value,mc_p=mc$p.value,min_expected=min(a$expected),small_expected=sum(a$expected<5),cells=length(tab),counts=unclass(tab),pearson=unclass(a$residuals),adjusted=unclass(a$stdres),holm=matrix(p.adjust(2*pnorm(-abs(a$stdres)),method='holm'),nrow=nrow(tab)),cramers_v=sqrt(unname(a$statistic)/(sum(tab)*min(nrow(tab)-1,ncol(tab)-1))))
}
tests<-list(injury=chi(inj),time=chi(table(d$time,d$type)),road=chi(table(d$road,d$type)),sex_cause=chi(table(d$sex1,d$cause)),sex_process=chi(table(d$sex1,d$process)))
js(tests,'association_results')
desc<-lapply(c('sex1','sex2','age1','age2','vehicle1','vehicle2','time','road','cause','process','injury','penalty','emotion_recorded','collision','conflict'),function(nm){tb<-table(d[[nm]]);data.frame(variable=nm,code=names(tb),count=as.integer(tb),valid_n=sum(tb),percent=100*as.integer(tb)/sum(tb),missing=sum(is.na(d[[nm]])))})
tabout(do.call(rbind,desc),'descriptives')
fitreg<-function(form,data,name){
 # Identical predictor patterns are collapsed into multinomial count rows; the likelihood is unchanged.
 predictors<-all.vars(form)[-1]
 groups<-interaction(data[,predictors,drop=FALSE],drop=TRUE,lex.order=TRUE)
 idx<-match(levels(groups),as.character(groups));g<-data[idx,predictors,drop=FALSE]
 Y<-unclass(table(groups,factor(data$type,levels=paste0('T',1:5))))
 g$Y<-I(Y);ff<-form;ff[[2]]<-quote(Y)
 f<-brmultinom(ff,data=g,ref=5,type='AS_mean',maxit=300,epsilon=1e-8)
 s<-summary(f);b<-coef(f);se<-s$standard.errors
 out<-do.call(rbind,lapply(1:nrow(b),function(i) data.frame(model=name,type=rownames(b)[i],term=colnames(b),B=b[i,],SE=se[i,],OR=exp(b[i,]),lower=exp(b[i,]-qnorm(.975)*se[i,]),upper=exp(b[i,]+qnorm(.975)*se[i,]),p=2*pnorm(-abs(b[i,]/se[i,])),n=nrow(data))))
 tabout(out,paste0('regression_',name))
 cat('Regression',name,'N=',nrow(data),'converged=',f$converged,'\n');out
}
reg<-fitreg(type~sex+age_group+time_group+road_group+vehicle_group,d[d$adjusted_sample,],'categorical_adjusted')
ld<-data.frame(injury=d$injury,penalty=d$penalty,collision=d$collision+1,emotion=d$emotion_recorded+1)
fl<-cbind(injury,penalty,collision,emotion)~1
fitset<-function(ld,formula,prefix,nrep=100){
 models<-list();rows<-list()
 usedvars<-all.vars(formula)
 levels_n<-sapply(ld[,usedvars,drop=FALSE],function(x)length(unique(na.omit(x))))
 maxK<-min(5,floor(prod(levels_n)/(sum(levels_n-1)+1)))
 for(k in seq_len(maxK)){
 set.seed(20261019+k)
 f<-poLCA(formula,ld,nclass=k,maxiter=10000,tol=1e-8,na.rm=FALSE,nrep=if(k==1)1 else nrep,verbose=FALSE,calc.se=TRUE)
 models[[k]]<-f
 rows[[k]]<-data.frame(K=k,LL=f$llik,q=f$npar,AIC=f$aic,BIC=f$bic,aBIC=-2*f$llik+f$npar*log((nrow(ld)+2)/24),best_replicates=sum(abs(f$attempts-f$llik)<1e-4),starts=length(f$attempts),iterations=f$numiter,min_class=min(f$P),entropy=if(k==1)NA else 1+sum(f$posterior*log(pmax(f$posterior,1e-300)))/(nrow(ld)*log(k)))
 stopifnot(abs(f$aic-(-2*f$llik+2*f$npar))<1e-6)
 cat('LCA',prefix,'K=',k,'LL=',f$llik,'BIC=',f$bic,'\n')
 }
 saveRDS(models,file.path(root,paste0('lca_',prefix,'.rds')))
 tabout(do.call(rbind,rows),paste0('lca_fit_',prefix))
 list(models=models,fit=do.call(rbind,rows))
}
lc<-fitset(ld,fl,'primary')
k<-which.min(lc$fit$BIC);f<-lc$models[[k]]
# Classes are labelled by increasing probability of criminal penalty, then collision.
ord<-order(f$probs$penalty[,2],f$probs$collision[,2])
prof<-do.call(rbind,lapply(names(f$probs),function(v){m<-f$probs[[v]][ord,,drop=FALSE];data.frame(indicator=v,category=rep(1:ncol(m),each=k),class=rep(paste0('C',1:k),ncol(m)),probability=as.vector(m))}))
post<-f$posterior[,ord,drop=FALSE];colnames(post)<-paste0('C',1:k)
tabout(prof,'lca_profiles');tabout(data.frame(class=paste0('C',1:k),proportion=f$P[ord],modal_n=tabulate(max.col(post),nbins=k),mean_posterior=sapply(1:k,function(c)mean(post[max.col(post)==c,c]))),'lca_classes')
tabout(data.frame(case_id=d$case_id,type=d$type,post,modal_class=paste0('C',max.col(post))),'lca_posterior')
cross<-aggregate(post,list(type=d$type),mean);tabout(cross,'lca_type_weighted')
js(list(selected_K=k,fit=lc$fit,profiles=prof,proportion=f$P[ord],entropy=lc$fit$entropy[k]),'lca_results')
# Refit without emotional expression to assess overlap with interaction coding.
lc3<-fitset(ld,cbind(injury,penalty,collision)~1,'no_emotion',50)
writeLines(capture.output(sessionInfo()),file.path(root,'R_session_info.txt'))
cat('ALL ANALYSES COMPLETED\n')
