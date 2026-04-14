##############################################################################
##############################################################################
######                                                                  ######
######            Deposit Compeition and Financial Fragility            ###### 
######            Calculate the Default Probability from CDS            ######
######                          03/03/2014                              ######
######           By Mark Egan, Ali Hortacsu and Gregor Matvos           ######
######                                                                  ######
##############################################################################
##############################################################################

setwd("")
### Parameters ###
libor <<- 0.03
recovery <<- 0.4
haz<- 0.02
cds<-0.012424885

### Functions ###
fn.cds.pv <- function(haz,cds){
  time<-seq(1,5,1)
  d.time<-seq(0.5,4.5,1)
  s.prob<-(1-haz)^time
  sell.pv.1<-sum(cds*s.prob*exp(-libor*time))
  sell.pv.2<-sum(cds*0.5*((1-haz)^(time-1)*haz)*exp(-libor*d.time))
  buy.pv<- sum((1-recovery)*((1-haz)^(time-1)*haz)*exp(-libor*d.time))
  return(((sell.pv.1+sell.pv.2-buy.pv)*100)^2) 
}
fn.cds.to.haz <- function(cds){
  res <-optimize(f=fn.cds.pv, interval=c(0,1), cds=cds,tol = 10e-8)
  return(res$minimum)
}

### Results ###
cds<- t(seq(0,0.80,0.0001))
haz<- apply(cds,2, fn.cds.to.haz)
results<-data.frame(cbind(t(cds),haz))
colnames(results)<-c("cds","hazard")
write.table(results,"cds to hazard.csv",sep=",", col.names=TRUE,row.names=FALSE, append=FALSE)
