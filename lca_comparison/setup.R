invisible(Sys.setlocale('LC_ALL','English_United States.utf8'))
root <- normalizePath(getwd(), winslash='/')
library(poLCA)
library(jsonlite)
src <- file.path(root,'../reproducibility')
out <- root
dir.create(out,showWarnings=FALSE,recursive=TRUE)
d <- read.csv(file.path(src,'final_analysis_data.csv'))
ld <- data.frame(injury=d$injury,penalty=d$penalty,collision=d$collision+1L,emotion=d$emotion_recorded+1L)
fl <- cbind(injury,penalty,collision,emotion)~1
old <- readRDS(file.path(src,'lca_primary.rds'))
js <- function(x,n) write_json(x,file.path(out,n),pretty=TRUE,auto_unbox=TRUE,digits=15,na='null')
fit <- function(dat,k,nrep=20,seed=1,initial=NULL) {
 set.seed(seed)
 poLCA(fl,dat,nclass=k,nrep=nrep,maxiter=20000,tol=1e-9,verbose=FALSE,calc.se=FALSE,na.rm=FALSE,probs.start=initial)
}
ll <- function(f,dat) {
 z<-matrix(rep(f$P,each=nrow(dat)),nrow(dat))
 for(v in names(dat)) {ok<-!is.na(dat[[v]]); z[ok,]<-z[ok,,drop=FALSE]*t(f$probs[[v]][,dat[ok,v],drop=FALSE])}
 log(rowSums(z))
}
align <- function(f,ref) {
 perm <- function(x) if(length(x)==1) list(x) else unlist(lapply(x,function(a) lapply(perm(setdiff(x,a)),function(b)c(a,b))),recursive=FALSE)
 pp<-perm(seq_along(f$P))
 dist<-sapply(pp,function(p) sum(sapply(names(f$probs),function(v) sum((f$probs[[v]][p,,drop=FALSE]-ref$probs[[v]])^2))))
 pp[[which.min(dist)]]
}
