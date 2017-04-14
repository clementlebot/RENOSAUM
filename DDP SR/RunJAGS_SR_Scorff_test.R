rm(list=ls(all=TRUE))
gc()

setwd("M:/Stage M2/Modele SCORFF/SR")

library(rjags)
library(coda)
#load.module("dic")?
#library(dclone)?


load("Y_1HM.Rdata")
load("Y_PHM.Rdata")
load("Dens.Rdata")

annee <- 20
surfaceERR <- 1907.83

data_scorff <- list(
year=annee,
SRR=1907.83,
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
tau=1,
a=0.15,
Rmax=10,
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
PARRpred.0.975 <- apply(PARRpred.mcmc,2,quantile,probs=0.975)
PARRpred.0.025 <- apply(PARRpred.mcmc,2,quantile,probs=0.025)

muPARR.mcmc <- matrix(rep(0,length(a.mcmc)*annee),nrow=length(a.mcmc),ncol=annee)
for (i in 1 : annee)
{
muPARR.mcmc[,i] <- log(EGG[i]/(1/a.mcmc+EGG[i]/Rmax.mcmc))
}
muPARR.med <- apply(muPARR.mcmc,2,median)

muPARRpred.mcmc <- matrix(rep(0,length(a.mcmc)*101),nrow=length(a.mcmc),ncol=101)
for (i in 1 :100)
{
muPARRpred.mcmc[,i] <- log((i*10)/(1/a.mcmc+(i*10)/Rmax.mcmc))
}
muPARRpred.mcmc[,101] <- log(2000/(1/a.mcmc+2000/Rmax.mcmc))
muPARRpred.med <- apply(muPARRpred.mcmc,2,median)


P50.mcmc <- matrix(rep(0,length(a.mcmc)*101),nrow=length(a.mcmc),ncol=101)

for (i in 1 :101)
{
	for (j in 1 : length(a.mcmc))
	{
	if((PARRpred.mcmc[j,i]-0.5*Rmax.mcmc[j]) > 0) P50.mcmc[j,i] <- 1 
	if((PARRpred.mcmc[j,i]-0.5*Rmax.mcmc[j]) <= 0) P50.mcmc[j,i] <- 0 
	}
}
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


#----------------------------------------
# Plot mcmc chains for the main variables
# ---------------------------------------
year<- 1:annee
pred <- 1:ncol(PARRpred.mcmc)
x <- 1:1000

pdf("Conv.pdf")
  var="e.s"
  par(mfrow=c(4,4)) 
  for (y in year)
  {
  traceplot(chains[,paste("e.s[",year[y],"]",sep="")],ylab=paste(var), main = paste("e.s",1994+year[y],sep=" "))
  densplot(chains[,paste("e.s[",year[y],"]",sep="")],ylab=paste(var), main = paste("e.s",1994+year[y],sep=" "))
  }

  var="Fsw"
  par(mfrow=c(4,4)) 
  for (y in year)
  {
  traceplot(chains[,paste("Fsw[",year[y],"]",sep="")],ylab=paste(var), main = paste("Fsw",1994+year[y],sep=" "))
  densplot(chains[,paste("Fsw[",year[y],"]",sep="")],ylab=paste(var), main = paste("Fsw",1994+year[y],sep=" "))
  }

  var="FSW"
  par(mfrow=c(4,4)) 
  for (y in year)
  {
  traceplot(chains[,paste("FSW[",year[y],"]",sep="")],ylab=paste(var), main = paste("FSW",1994+year[y],sep=" "))
  densplot(chains[,paste("FSW[",year[y],"]",sep="")],ylab=paste(var), main = paste("FSW",1994+year[y],sep=" "))
  }

  par(mfrow=c(3,2))
  
  var="tau"
  
  traceplot(chains[,"tau"],ylab=paste(var), main = "tau")
  #PRIORS VS POSTERIOR 
  densplot(chains[,"tau"],ylab=paste(var), main = "tau",col="red")
  curve(dgamma(x,0.01,0.01),from=0,to=7,n=500,add=T)
  
  
  var="a"
  
	traceplot(chains[,"a"],ylab=paste(var), main = "a")    
	densplot(chains[,"a"],ylab=paste(var), main = "a",col="red",xlim=c(0,1))
	curve(dbeta(x,1,1),from=0,to=1,n=500, add=T)


  var="Rmax"
  
    traceplot(chains[,"Rmax"],ylab=paste(var), main = "Rmax")
    #PRIORS VS POSTERIOR 
    densplot(chains[,"Rmax"],ylab=paste(var), main = "Rmax", xlim=c(10,100),col="red")
    curve(dunif(x,0,100),from=0,to=100,n=500, add=T)

dev.off()

pdf("Pred.pdf")

    var="PARRpred"
    par(mfrow=c(4,4)) 
    for (y in pred)
    {
      traceplot(chains[,paste("PARRpred[",y,"]",sep="")],ylab=paste(var), main = paste("PARRpred S=",stock.pred[y],sep=""))
      densplot(chains[,paste("PARRpred[",y,"]",sep="")],ylab=paste(var), main = paste("PARRpred",", year= ",stock.pred[y],sep=""))
    }
dev.off()

pdf("Graphe_interet.pdf")
    
par(mfrow=c(2,2),mar=c(5,5,2,2))
plot(a.mcmc,Rmax.mcmc,main=paste("Joint post."),xlab="a",ylab="Rmax",axes=F)
axis(1,pos=0)
axis(2,las=2)

boxplot(e.s.mcmc,axes=F,main=paste("Res. log"))
axis(2,at=-2:1,labels=-2:1,las=2,pos=0)
axis(1,at=c(seq(1,20,5),20),labels=paste(c(seq(1995,2014,5),2014)))

plot(0,0,type="n",xlim=c(0,2000), ylim=c(0,40), xlab="", ylab="",axes=F,main="SR curve")
polygon(c(stock.pred,stock.pred[length(stock.pred):1]),
c(0,PARRpred.0.975,PARRpred.0.025[length(PARRpred.0.025):1],0), col="grey", border = NA)
points(stock.pred,c(0,exp(muPARRpred.med)),type="l")
points(EGG, data_scorff$PARR, pch=19)
axis(1)
mtext("Egg /m^2",1,3)
axis(2, las=2)
mtext("Parr /m^2",2,4)


### ------------------- ###
### Diagramme de risque ###
### ------------------- ###

colors <- c(brewer.pal(3,"Greens")[3],brewer.pal(3,"Greys")[3],brewer.pal(3,"Reds")[3])

plot(c(seq(10,1000,10),2000),P30.mean,ylim=c(0,1),type="l",axes=F,xlab="",ylab="",col=colors[1], main="Risk analysis")
points(c(seq(10,1000,10),2000),P50.mean,col=colors[2],type="l")
points(c(seq(10,1000,10),2000),P70.mean,col=colors[3],type="l")

axis(1,pos=0)
mtext("Stock",1,3)
axis(2,pos=0,las=2)
mtext("P(R < x*Rmax)",2,3)
legend("topright",fill=colors,legend=c("x = 0.3","x = 0.5","x = 0.7"),bty="n")

dev.off()


