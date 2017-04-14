rm(list=ls(all=TRUE))
gc()

library(lattice)
library(mcmcplots)
library(rjags)
library(coda)
library(RColorBrewer)
#load.module("dic")?
#library(dclone)?


load("Y_1HM.Rdata")
load("Y_PHM.Rdata")
load("Dens.Rdata")

annee <- 20
surfaceERR <- 1907.83

data_scorff <- list(
year=annee,
SRR=surfaceERR,
sw=c(509,784,694,467,555,252,305,310,582,229,1119,452,890,483,284,260,764,386,378,599),
NBEGGsw=3485,
PROPFEMsw=0.45,
Ysw=round(Y_1HM$SCORFF[paste(seq(1994,2013,1))]),
SW=c(80,96,84,70,27,94,46,41,22,59,66,129,104,94,88,116,65,217,145,103),
NBEGGSW=5569,
PROPFEMSW=0.8,
YSW=round(Y_PHM$SCORFF[paste(seq(1994,2013,1))]),
PARR=c(9926,13410,8866,6445,10850,20210,3169,17710,13120,24470,13890,15280,25830,21570,17230,21600,33670,32270,25210,25970)/surfaceERR
)
EGG <- (data_scorff$SW*data_scorff$PROPFEMSW*data_scorff$NBEGGSW + data_scorff$sw*data_scorff$PROPFEMsw*data_scorff$NBEGGsw)/data_scorff$SRR
stock.pred <- c(seq(0,1000,10),2000)

init_scorff <- list(
a=0.15,
Rmax=13,
Fsw=rep(0.2,annee),
FSW=rep(0.2,annee),
n=2,
theta=data_scorff$PARR/EGG
)


monitor <- c("lg.e","Rmax","a","Fsw","FSW","n","theta","alphapred","thetapred","thetaquant","PARRpred","muPARRpred","tauPARRpred")


nchains <- 1
nadapt <- 10000
nburnin <- 10000
niter <- 100000
nthin <- 10


jm <- jags.model(file="Survie_Scorff_normal_JAGS.txt",data=data_scorff,inits=init_scorff, n.chains=nchains, n.adapt=nadapt)
update(jm, n.iter=nburnin)
chains <- coda.samples(model=jm,variable.names=monitor,n.iter=niter,thin=nthin)

Mat <- as.matrix(chains)

x = "Rmax"
Rmax.mcmc <- Mat[,which(substr(colnames(Mat),1,nchar(x)+1)==x)]
Rmax.med <- median(Rmax.mcmc)

x = "a"
a.mcmc <- Mat[,which(substr(colnames(Mat),1,nchar(x)+1)==x)]
a.med <- median(a.mcmc)

x = "n"
n.mcmc <- Mat[,which(substr(colnames(Mat),1,nchar(x)+1)==x)]
n.med <- median(n.mcmc)


S50max.mcmc <- Rmax.mcmc/a.mcmc
S50max.med <- median(S50max.mcmc)

x <- "lg.e"
lg.e.mcmc <- Mat[,which(substr(colnames(Mat),1,nchar(x)+1)==paste(x,"[",sep=""))]
lg.e.med <- apply(lg.e.mcmc,2,median)

x <- "Fsw"
Fsw.mcmc <- Mat[,which(substr(colnames(Mat),1,nchar(x)+1)==paste(x,"[",sep=""))]
Fsw.med <- apply(Fsw.mcmc,2,median)

x <- "FSW"
FSW.mcmc <- Mat[,which(substr(colnames(Mat),1,nchar(x)+1)==paste(x,"[",sep=""))]
FSW.med <- apply(FSW.mcmc,2,median)

x <- "theta"
theta.mcmc <- Mat[,which(substr(colnames(Mat),1,nchar(x)+1)==paste(x,"[",sep=""))]
theta.med <- apply(theta.mcmc,2,median)

x <- "alphapred"
alphapred.mcmc <- Mat[,which(substr(colnames(Mat),1,nchar(x)+1)==paste(x,"[",sep=""))]
alphapred.med <- apply(alphapred.mcmc,2,median)

betapred.mcmc <- n.mcmc-alphapred.mcmc
betapred.med <- apply(betapred.mcmc,2,median)

x <- "thetapred"
thetapred.mcmc <- Mat[,which(substr(colnames(Mat),1,nchar(x)+1)==paste(x,"[",sep=""))]
thetapred.med <- apply(thetapred.mcmc,2,median)

x <- "thetaquant"
thetaquant.mcmc <- Mat[,which(substr(colnames(Mat),1,nchar(x)+1)==paste(x,"[",sep=""))]
thetaquant.med <- apply(thetaquant.mcmc,2,median)

x <- "PARRpred"
PARRpred.mcmc <- Mat[,which(substr(colnames(Mat),1,nchar(x)+1)==paste(x,"[",sep=""))]
PARRpred.med <- apply(PARRpred.mcmc,2,median)

x <- "muPARRpred"
muPARRpred.mcmc <- Mat[,which(substr(colnames(Mat),1,nchar(x)+1)==paste(x,"[",sep=""))]
muPARRpred.med <- apply(muPARRpred.mcmc,2,median)

x <- "tauPARRpred"
tauPARRpred.mcmc <- Mat[,which(substr(colnames(Mat),1,nchar(x)+1)==paste(x,"[",sep=""))]
tauPARRpred.med <- apply(tauPARRpred.mcmc,2,median)

PARRpred.0.975 <- apply(PARRpred.mcmc,2,quantile,probs=0.975)
PARRpred.0.025 <- apply(PARRpred.mcmc,2,quantile,probs=0.025)

PARRquant.mcmc <- matrix(rep(0,nrow(thetaquant.mcmc),ncol(thetaquant.mcmc)),nrow=nrow(thetaquant.mcmc),ncol=ncol(thetaquant.mcmc))
for(t in 1:annee)
{
PARRquant.mcmc[,t] <- thetaquant.mcmc[,t]*EGG[t]
}

q <-c()
PARR.ranked <- c()
EGG.rank <- order(EGG)
for (t in 1:annee)
{
PARR.ranked[EGG.rank[t]] <- data_scorff$PARR[t]
}

for(t in 1:annee)
{
	for (i in 1:100)
	{
	if(round(quantile(PARRquant.mcmc[,t],probs=i*0.01),0)==round(PARR.ranked[t],0))
	q[t] <- i*0.01
	}
}


plot(sort(q),1:(annee)/(annee),xlim=c(0,1),xlab="quantile", ylab="frequence")
boxplot(PARRquant.mcmc)
stock.pred <- c(seq(0,1000,10),2000)


muthetapred.mcmc <- matrix(rep(0,nrow(alphapred.mcmc),ncol(alphapred.mcmc)),nrow=nrow(alphapred.mcmc),ncol=ncol(alphapred.mcmc))
for (i in 1 : ncol(muthetapred.mcmc))
{
muthetapred.mcmc[,i] <- alphapred.mcmc[,i] / n.mcmc
}
muthetapred.med <- apply(muthetapred.mcmc,2,mean) 
thetapred.0.975 <- apply(thetapred.mcmc,2,quantile,probs=0.975)
thetapred.0.025 <- apply(thetapred.mcmc,2,quantile,probs=0.025)

P50.mcmc <- matrix(rep(0,length(a.mcmc)*ncol(PARRpred.mcmc)),nrow=length(a.mcmc),ncol=ncol(PARRpred.mcmc))
P30.mcmc <- matrix(rep(0,length(a.mcmc)*ncol(PARRpred.mcmc)),nrow=length(a.mcmc),ncol=ncol(PARRpred.mcmc))
P70.mcmc <- matrix(rep(0,length(a.mcmc)*ncol(PARRpred.mcmc)),nrow=length(a.mcmc),ncol=ncol(PARRpred.mcmc))

for (i in 1 :ncol(PARRpred.mcmc))
{
  for (j in 1 : length(a.mcmc))
  {
    if((PARRpred.mcmc[j,i]-0.3*Rmax.mcmc[j]) < 0) P30.mcmc[j,i] <- 1 
    if((PARRpred.mcmc[j,i]-0.3*Rmax.mcmc[j]) >= 0) P30.mcmc[j,i] <- 0   
    if((PARRpred.mcmc[j,i]-0.5*Rmax.mcmc[j]) < 0) P50.mcmc[j,i] <- 1 
    if((PARRpred.mcmc[j,i]-0.5*Rmax.mcmc[j]) >= 0) P50.mcmc[j,i] <- 0 
    if((PARRpred.mcmc[j,i]-0.7*Rmax.mcmc[j]) < 0) P70.mcmc[j,i] <- 1 
    if((PARRpred.mcmc[j,i]-0.7*Rmax.mcmc[j]) >= 0) P70.mcmc[j,i] <- 0 
  }
}
P30.mean <- apply(P30.mcmc,2,mean)
P50.mean <- apply(P50.mcmc,2,mean)
P70.mean <- apply(P70.mcmc,2,mean)


alpha.mcmc <- matrix(rep(0,length(a.mcmc)*annee),nrow=length(a.mcmc),ncol=annee)
for (t in 1 : annee)
{
alpha.mcmc[,t] <- n.mcmc/((1/a.mcmc)+(EGG[t]/Rmax.mcmc))
}

beta.mcmc <- n.mcmc - alpha.mcmc
boxplot(alpha.mcmc,outliers=F)
boxplot(beta.mcmc)


#----------------------------------------
# Plot mcmc chains for the main variables
# ---------------------------------------
year<- 1:annee
pred <- 1:ncol(PARRpred.mcmc)
x <- 1:1000

pdf("Normal_conv.pdf")
  
  par(mfrow=c(4,4))
  var="lg.e"
  for (y in 1:annee)
  {
  traceplot(chains[,paste("lg.e[",year[y],"]",sep="")],ylab=paste(var), main = paste("lg.e",1994+year[y],sep=" "))
  densplot(chains[,paste("lg.e[",year[y],"]",sep="")],ylab=paste(var), main = paste("lg.e",1994+year[y],sep=" "))
  }
  
  par(mfrow=c(4,4))
  var="Fsw"
  for (y in 1:annee)
  {
  traceplot(chains[,paste("Fsw[",year[y],"]",sep="")],ylab=paste(var), main = paste("Fsw",1994+year[y],sep=" "))
  densplot(chains[,paste("Fsw[",year[y],"]",sep="")],ylab=paste(var), main = paste("Fsw",1994+year[y],sep=" "))
  }
  
  par(mfrow=c(4,4))
  var="FSW"
  for (y in 1:annee)
  {
  traceplot(chains[,paste("FSW[",year[y],"]",sep="")],ylab=paste(var), main = paste("FSW",1994+year[y],sep=" "))
  densplot(chains[,paste("FSW[",year[y],"]",sep="")],ylab=paste(var), main = paste("FSW",1994+year[y],sep=" "))
  }
  
  par(mfrow=c(4,4))
  var="theta"
  for (y in 1:annee)
  {
  traceplot(chains[,paste("theta[",year[y],"]",sep="")],ylab=paste(var), main = paste("theta",1994+year[y],sep=" "))
  densplot(chains[,paste("theta[",year[y],"]",sep="")],ylab=paste(var), main = paste("theta",1994+year[y],sep=" "))
  }

  par(mfrow=c(3,2))

  var="n"
  
	traceplot(chains[,"n"],ylab=paste(var))   
	densplot(chains[,"n"],ylab=paste(var),col="red",xlim=c(0,1000))
	curve(dunif(x,0,1000),from=0,to=1000,n=500, add=T)


  var="a"
  
	traceplot(chains[,"a"],ylab=paste(var))    
	densplot(chains[,"a"],ylab=paste(var),col="red",xlim=c(0,1))
	curve(dbeta(x,1,1),from=0,to=1,n=500, add=T)
	

  var="Rmax"
  
    traceplot(chains[,"Rmax"],ylab=paste(var))
    #PRIORS VS POSTERIOR 
    densplot(chains[,"Rmax"],ylab=paste(var), xlim=c(10,100),col="red")
    curve(dunif(x,0,100),from=0,to=100,n=500, add=T)

dev.off()

pdf("Normal_pred.pdf")

par(mfrow=c(4,4))

var="PARRpred"
for (y in pred)
{
  traceplot(chains[,paste("PARRpred[",y,"]",sep="")],ylab=paste(var), main = paste("PARRpred S=",stock.pred[y],sep=" "))
  densplot(chains[,paste("PARRpred[",y,"]",sep="")],ylab=paste(var), main = paste("PARRpred S=",stock.pred[y],sep=" "))
}

par(mfrow=c(4,4))

var="thetapred"
for (y in pred)
{
  traceplot(chains[,paste("thetapred[",y,"]",sep="")],ylab=paste(var), main = paste("thetapred S=",stock.pred[y],sep=" "))
  densplot(chains[,paste("thetapred[",y,"]",sep="")],ylab=paste(var), main = paste("thetapred S=",stock.pred[y],sep=" "))
}

dev.off()

### ------------------- ###
###  Analyse graphique  ###
### ------------------- ###
pdf("Normal_graph_interet.pdf")

par(mfrow=c(1,1))

## --post. jointe a_Rmax-- ##

plot(a.mcmc,Rmax.mcmc,main=paste("Joint post."),xlab="a",ylab="Rmax",axes=F)
axis(1,pos=0)
axis(2,las=2)

## ------R?sidus------ ##

boxplot(lg.e.mcmc,axes=F,main="Res. log")
axis(2,las=2)
axis(1,at=c(seq(1,20,5),20),labels=paste(c(seq(1995,2014,5),2014)))


## ---Courbes B-H--- ##

plot(0,0,type="n",xlim=c(0,2000), ylim=c(0,60), xlab="", ylab="",axes=F, main="SR relationship")
polygon(c(stock.pred,stock.pred[length(stock.pred):1]),
c(PARRpred.0.975,PARRpred.0.025[length(PARRpred.0.025):1]), col="grey", border = NA)
points(stock.pred,muthetapred.med*stock.pred,type="l")
points(EGG, data_scorff$PARR, pch=19)
axis(1)
mtext("Egg /m^2",1,2)
axis(2, las=2)
mtext("Parr /m^2",2,2)

plot(0,0,type="n",xlim=c(0,2000), ylim=c(0,1), xlab="", ylab="",axes=F, main="BH survival")
polygon(c(stock.pred,stock.pred[length(stock.pred):1]),
c(thetapred.0.975,thetapred.0.025[length(thetapred.0.025):1]), col="grey", border = NA)
points(stock.pred,muthetapred.med,type="l")
points(EGG, theta.med, pch=19)
axis(1,pos=0)
mtext("Egg /m?",1,2)
axis(2, las=2,pos=0)
mtext("Survie egg-to-parr",2,2)

## --Analyse risque-- ##

colors <- c(brewer.pal(3,"Greens")[3],brewer.pal(3,"Greys")[3],brewer.pal(3,"Reds")[3])

plot(stock.pred,P30.mean,ylim=c(0,1),type="l",axes=F,xlab="",ylab="",col=colors[1], main="Risk analysis")
points(stock.pred,P50.mean,col=colors[2],type="l")
points(stock.pred,P70.mean,col=colors[3],type="l")

axis(1,pos=0)
mtext("Stock",1,2)
axis(2,pos=0,las=2)
mtext("P(R < x*Rmax)",2,2)
legend("topright",fill=colors,legend=c("x = 0.3","x = 0.5","x = 0.7"),bty="n")

dev.off()

