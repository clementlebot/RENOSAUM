rm(list=ls(all=TRUE))
gc()

setwd("M:/Stage M2/Mod?le SCORFF/SR Scorff")

library(rjags)
library(coda)
#load.module("dic")?
#library(dclone)?


load("Y_1HM.Rdata")
load("Y_PHM.Rdata")
load("Dens.Rdata")

annee <- 21

data_scorff <- list(
year=annee,
SRR=1907.83,
sw=c(509,784,694,467,555,252,305,310,582,229,1119,452,890,483,284,260,764,386,378,599,729),
NBEGGsw=3485,
PROPFEMsw=0.45,
Ysw=round(Y_1HM$SCORFF[paste(seq(1994,2014,1))]),
SW=c(80,96,84,70,27,94,46,41,22,59,66,129,104,94,88,116,65,217,145,103,139),
NBEGGSW=5569,
PROPFEMSW=0.8,
YSW=round(Y_PHM$SCORFF[paste(seq(1994,2014,1))]),
PARR=Dens$SCORFF[paste(seq(1995,2015,1))]
)

EGG <- (data_scorff$SW*data_scorff$PROPFEMSW*data_scorff$NBEGGSW + data_scorff$sw*data_scorff$PROPFEMsw*data_scorff$NBEGGsw)/data_scorff$SRR

init_scorff <- list(
tau=1,
a=0.15,
Rmax=40,
PARRpred=c(seq(0.5,50,0.5),100),
Fsw=rep(0.2,annee),
FSW=rep(0.2,annee)
)



monitor <- c("tau","e.s","a","Rmax","PARRpred","Fsw","FSW")


nchains <- 1
nadapt <- 10000
nburnin <- 10000
niter <- 100000
nthin <- 10


jm <- jags.model(file="SR_Scorff_test_JAGS.txt",data=data_scorff,inits=init_scorff, n.chains=nchains, n.adapt=nadapt)
update(jm, n.iter=nburnin)
chains <- coda.samples(model=jm,variable.names=monitor,n.iter=niter,thin=nthin)

gelman.diag(chains)


Mat <- as.matrix(chains)

head(Mat)
colnames(Mat)

x <- "tau"
tau.mcmc <- Mat[,which(substr(colnames(Mat),1,nchar(x)+1)==x)]
tau.med <- median(tau.mcmc)

x = "Rmax"
Rmax.mcmc <- Mat[,which(substr(colnames(Mat),1,nchar(x)+1)==x)]
Rmax.med <- median(Rmax.mcmc)

x = "a"
a.mcmc <- Mat[,which(substr(colnames(Mat),1,nchar(x)+1)==x)]
a.med <- median(a.mcmc)


S50max.mcmc <- Rmax.mcmc/a.mcmc
S50max.med <- median(S50max.mcmc)

x <- "e.s"
e.s.mcmc <- Mat[,which(substr(colnames(Mat),1,nchar(x)+1)==paste(x,"[",sep=""))]
e.s.med <- apply(e.s.mcmc,2,median)

x <- "Fsw"
Fsw.mcmc <- Mat[,which(substr(colnames(Mat),1,nchar(x)+1)==paste(x,"[",sep=""))]
Fsw.med <- apply(Fsw.mcmc,2,median)

x <- "FSW"
FSW.mcmc <- Mat[,which(substr(colnames(Mat),1,nchar(x)+1)==paste(x,"[",sep=""))]
FSW.med <- apply(FSW.mcmc,2,median)

x <- "PARRpred"
PARRpred.mcmc <- Mat[,which(substr(colnames(Mat),1,nchar(x)+1)==paste(x,"[",sep=""))]
PARRpred.med <- apply(PARRpred.mcmc,2,median)

muPARR.mcmc <- matrix(rep(0,length(a.mcmc)*annee),nrow=length(a.mcmc),ncol=annee)
for (i in 1 : annee)
{
muPARR.mcmc[,i] <- log(EGG[i]/(1/a.mcmc+EGG[i]/Rmax.mcmc))
}
muPARR.med <- apply(muPARR.mcmc,2,median)

exp()

muPARRpred.mcmc <- matrix(rep(0,length(a.mcmc)*101),nrow=length(a.mcmc),ncol=101)
for (i in 1 :100)
{
muPARRpred.mcmc[,i] <- log((i*10)/(1/a.mcmc+(i*10)/Rmax.mcmc))
}
muPARRpred.mcmc[,101] <- log(2000/(1/a.mcmc+2000/Rmax.mcmc))

muPARRpred.med <- apply(muPARRpred.mcmc,2,median)
muPARRpred.0.975 <- apply(muPARRpred.mcmc,2,quantile,probs=0.975)
muPARRpred.0.025 <- apply(muPARRpred.mcmc,2,quantile,probs=0.025)

P50.mcmc <- matrix(rep(0,length(a.mcmc)*101),nrow=length(a.mcmc),ncol=101)

for (i in 1 :101)
{
	for (j in 1 : length(a.mcmc))
	{
	if((PARRpred.mcmc[j,i]-0.5*Rmax.mcmc[j]) > 0) P50.mcmc[j,i] <- 1 
	if((PARRpred.mcmc[j,i]-0.5*Rmax.mcmc[j]) <= 0) P50.mcmc[j,i] <- 0 
	}
}

P50.moy <- apply(P50.mcmc,2,mean)


#----------------------------------------
# Plot mcmc chains for the main variables
# ---------------------------------------
year<- seq(1,annee,5)
n.plot <- 4
x <- 1:1000

# Marine survival and maturing probability
# ----------------------------------------
  
  var="e.s"
  par(mfrow=c(n.plot/2,n.plot))
  for (y in 1:n.plot)
  {
  traceplot(chains[,paste("e.s[",year[y],"]",sep="")],ylab=paste(var), main = paste("e.s",", year= ",year[y],sep=""))
  densplot(chains[,paste("e.s[",year[y],"]",sep="")],ylab=paste(var), main = paste("e.s",", year= ",year[y],sep=""))
  }

  var="PARRpred"
  par(mfrow=c(n.plot/2,n.plot))
  for (y in 1:n.plot)
  {
  traceplot(chains[,paste("PARRpred[",year[y],"]",sep="")],ylab=paste(var), main = paste("PARRpred",", year= ",year[y],sep=""))
  densplot(chains[,paste("PARRpred[",year[y],"]",sep="")],ylab=paste(var), main = paste("PARRpred",", year= ",year[y],sep=""))
  }

  var="Fsw"
  par(mfrow=c(n.plot/2,n.plot))
  for (y in 1:n.plot)
  {
  traceplot(chains[,paste("Fsw[",year[y],"]",sep="")],ylab=paste(var), main = paste("Fsw",", year= ",year[y],sep=""))
  densplot(chains[,paste("Fsw[",year[y],"]",sep="")],ylab=paste(var), main = paste("Fsw",", year= ",year[y],sep=""))
  }

  var="FSW"
  par(mfrow=c(n.plot/2,n.plot))
  for (y in 1:n.plot)
  {
  traceplot(chains[,paste("FSW[",year[y],"]",sep="")],ylab=paste(var), main = paste("FSW",", year= ",year[y],sep=""))
  densplot(chains[,paste("FSW[",year[y],"]",sep="")],ylab=paste(var), main = paste("FSW",", year= ",year[y],sep=""))
  }

  var="a"
  
	traceplot(chains[,"a"],ylab=paste(var), main = "a")    
	densplot(chains[,"a"],ylab=paste(var), main = "a",col="red",xlim=c(0,1))
	curve(dbeta(x,1,1),from=0,to=1,n=500, add=T)


  var="Rmax"
  
    traceplot(chains[,"Rmax"],ylab=paste(var), main = "Rmax")
    #PRIORS VS POSTERIOR 
    densplot(chains[,"Rmax"],ylab=paste(var), main = "Rmax", xlim=c(10,100),col="red")

plot(a.mcmc,Rmax.mcmc)

quantile(a.mcmc,probs = c(0.05,0.5,0.95))

  var="tau"
  
    traceplot(chains[,"tau"],ylab=paste(var), main = "tau")
    #PRIORS VS POSTERIOR 
    densplot(chains[,"tau"],ylab=paste(var), main = "tau",col="red")
    curve(dgamma(x,0.01,0.01),from=0,to=1,n=500,add=T)


par(mfrow=c(1,1),mar=c(5,5,0,0))
plot(0,0,type="n",xlim=c(0,2000), ylim=c(0,60), xlab="", ylab="",axes=F)
polygon(c(seq(0,1000,10),2000,2000,sort(seq(0,1000,10),decreasing=T)),
c(0,exp(muPARRpred.0.975),sort(exp(muPARRpred.0.025),decreasing=T),0), col="grey", border = NA)
points(c(seq(0,1000,10),2000),c(0,exp(muPARRpred.med)),type="l")
points(EGG, data_scorff$PARR, pch=19)
axis(1)
mtext("Egg /m?",1,3)
axis(2, las=2)
mtext("Parr /m?",2,4)

windows()
plot()




