
rm(list=ls())

library(mnormt)
library(lme4)
library(MASS)
library(mvtnorm)
library(REEMtree)
library(heavy)
library(nlme)
library(rpart)
library(rpart.plot)
library(ggplot2)
library(dplyr)
library(plyr)
library(lattice)
library(Metrics)
library(data.table)
source("Heavy_Algorithm.R")

#### Data Simulation Process-7  ####

#### From Multivariate Normal Distribution with Outlier #####

p <- 4  # number of fixed effect 
q <- 2 # number of random effect 


n <- 100 # number of person / group / cluster
m <- 50  # number of repeated measurements
subject <- rep(c(1:n),each=m) #subject id's
N<- n*m # number of observations

D=matrix(c(0.25,0,0.125,0.25),q,q)  # variance of random effect
I = diag(m)  # variance of errors


## simulate data from Multivariate Normal Distribution
set.seed(22)
b <- rmnorm(n, rep(0, q), D)
e <- rmnorm(n, rep(0, m), I)

head(b) 
head(e)


set.seed(12)
X1<-matrix(runif(N,min=0,max=10), nrow=n, ncol=m)
X2<-matrix(runif(N,min=0,max=10), nrow=n, ncol=m)
X3<-matrix(runif(N,min=0,max=10), nrow=n, ncol=m)

mu_1=-20 
mu_2=-10 
mu_3= 10 
mu_4= 20




y <- matrix(0,n,m)

for( i in 1:n) {
  
  for( j in 1:m) {
  
  if(X1[i,j] <=5 & X2[i,j]<= 5 ) {
    
    y[i,j]=mu_1 +  cbind(1,X1[i,j] )%*%b[i,] +  e[i,j]
    
  } else if(X1[i,j] <=5 & X2[i,j]> 5 ) {
    
    y[i,j]=mu_2 + cbind(1,X1[i,j] )%*%b[i,] +   e[i,j]  
    
  } else if (X1[i,j] > 5 & X2[i,j]<= 5 ) {
    
    y[i,j]=mu_3 +  cbind(1,X1[i,j] )%*%b[i,] +  e[i,j] 
    
  } else if(X1[i,j] > 5 & X2[i,j] > 5 ) {
    
    y[i,j]=mu_4 +  cbind(1,X1[i,j] )%*%b[i,] +  e[i,j] 
    
  } else {
    
    
  }
  
  
  }
}


X1<-cbind(as.factor(1:n),X1)
colnames(X1)<-c("id",1:m)
X1<-as.data.table(X1)
head(X1)
X1<-melt(X1,id.vars="id", variable.name = "variable",value.name = "X1")
head(X1)
X1<-X1[,-2]
head(X1)
X1<-ddply(X1, "id")
head(X1,10)


X2<-cbind(as.factor(1:n),X2)
colnames(X2)<-c("id",1:m)
X2<-as.data.table(X2)
head(X2)
X2<-melt(X2,id.vars="id", variable.name = "variable",value.name = "X2")
head(X2)
X2<-X2[,-2]
head(X2)
X2<-ddply(X2, "id")
head(X2,10)

X3<-cbind(as.factor(1:n),X3)
colnames(X3)<-c("id",1:m)
X3<-as.data.table(X3)
head(X3)
X3<-melt(X3,id.vars="id", variable.name = "variable",value.name = "X3")
head(X3)
X3<-X3[,-2]
head(X3)
X3<-ddply(X3, "id")
head(X3,10)

y<-cbind(as.factor(1:n),y)
colnames(y)<-c("id",1:m)
y<-as.data.table(y)
head(y)
y<-melt(y,id.vars="id", variable.name = "variable",value.name = "y")
head(y)
y<-y[,-2]
head(y)
y<-ddply(y, "id")
head(y,10)


simdat<-cbind(y,X1,X2,X3)
head(simdat)
simdat<-simdat[,-c(3,5,7)]
head(simdat)

#Adding Outliers
simdat$y[496]<-simdat$y[496]+100
simdat$y[497]<-simdat$y[497]+100
simdat$y[498]<-simdat$y[498]+100
simdat$y[499]<-simdat$y[499]+100
simdat$y[500]<-simdat$y[500]+100

####evolution of simulated data from Multivariate Normal Distribution with Outlier ####

reem_model<-REEMtree(y ~X1+X2+X3 , data=simdat ,random=~1|id,cv=TRUE)
reem_model
rpart.plot(reem_model$Tree,digits = 5,roundint=FALSE, main="Decision Tree of \n REEMtree Model")
reem_model$BetweenMatrix
reem_model$EffectModel
reem_predict<-fitted(reem_model)
rmse(simdat$y,reem_predict)
mean((simdat$y - reem_predict)^2)   #mean square error
AIC(reem_model)
BIC(reem_model)



reem_model_heavy<-REEMtree_heavy_mode(y ~X1+X2+X3 , data=simdat ,random=~1, groups=~id,cv=TRUE,family_h=Student(df = 5))
reem_model_heavy
rpart.plot(reem_model_heavy$Tree,digits = 5,roundint=FALSE, main="Decision Tree of \n Heavy_REEMtree Model")
reem_model_heavy$BetweenMatrix
reem_model_heavy$EffectModel
heavy_predict<-reem_model_heavy$fitted
rmse(simdat$y,heavy_predict)
mean((simdat$y - reem_model_heavy$fitted)^2)  #mean square error
reem_model_heavy$AIC
reem_model_heavy$BIC



lme_model<-lme(y ~X1+X2+X3 , data = simdat,random=~1|id)
lme_model
lme_model$logLik
VarCorr(lme_model)  ### estimated variance of errors
lme_predict<-fitted(lme_model)
rmse(simdat$y,fitted(lme_model))
mean((simdat$y - fitted(lme_model))^2)    #mean square error
AIC(lme_model)
BIC(lme_model)



