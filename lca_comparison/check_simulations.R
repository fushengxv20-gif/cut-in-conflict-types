invisible(Sys.setlocale('LC_ALL','English_United States.utf8'))
source('setup.R',encoding='UTF-8')
rows<-list()
for(i in c(0,17,44)) {
 z<-read.csv(file.path(out,paste0('simulation_',i,'.csv')))
 for(k in 3:4){
  f<-fit(z,k,20,202611090+i+k)
  rows[[length(rows)+1]]<-data.frame(replicate=i,K=k,LL=f$llik)
  cat('Crosscheck',i,k,format(f$llik,digits=15),'\n');flush.console()
  write.csv(do.call(rbind,rows),file.path(out,'simulation_poLCA_crosschecks.csv'),row.names=FALSE)
 }
}
