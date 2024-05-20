
rm(list = ls())

library(REEMtree)
library(heavy)
library(nlme)
library(rpart)
library(rpart.plot)
library(ggplot2)
library(dplyr)
library(Metrics)
library(qpcR)
library(lattice)
library(ALA)
source("Heavy_Algorithm.R")


#### Descriptive Statistics and Data Preparation for Lung_Function_Growth ####
data(fev1)


fev1$fev<-exp(fev1$logFEV1) #back transformation of logarithmic of fev to normal

fev1$fev[1]<-fev1$fev[1]+10 #add noise 

fev1$fev[1994]<-fev1$fev[1994]+10  #add noise 

fev1<-as.data.frame(fev1)

dim(fev1)
head(fev1)
summary(fev1)

### Lattice Graphs 20 by 20 ###

#1:20

index<-1:20
sample1 <- fev1[fev1$id %in% index, ] 
sample1<-as.data.frame(sample1)

sample1$id<-as.factor(sample1$id)

xyplot(fev~age | id, 
       data=sample1, 
       panel=function(x,y){
         panel.xyplot(x, y)
         panel.lmline(x,y)
       },  as.table=T)


#21:40

index<-21:40
sample2 <- fev1[fev1$id %in% index, ] 
sample2<-as.data.frame(sample2)

sample2$id<-as.factor(sample2$id)

xyplot(fev~age | id, 
       data=sample2, 
       panel=function(x,y){
         panel.xyplot(x, y)
         panel.lmline(x,y)
       },  as.table=T)

#41:60

index<-41:60
sample3 <- fev1[fev1$id %in% index, ] 
sample3<-as.data.frame(sample3)

sample3$id<-as.factor(sample3$id)

xyplot(fev~age | id, 
       data=sample3, 
       panel=function(x,y){
         panel.xyplot(x, y)
         panel.lmline(x,y)
       },  as.table=T)

#61:80

index<-61:80
sample4 <- fev1[fev1$id %in% index, ] 
sample4<-as.data.frame(sample4)

sample4$id<-as.factor(sample4$id)

xyplot(fev~age | id, 
       data=sample4, 
       panel=function(x,y){
         panel.xyplot(x, y)
         panel.lmline(x,y)
       },  as.table=T)

#81:100

index<-81:100
sample5 <- fev1[fev1$id %in% index, ] 
sample5<-as.data.frame(sample5)

sample5$id<-as.factor(sample5$id)

xyplot(fev~age | id, 
       data=sample5, 
       panel=function(x,y){
         panel.xyplot(x, y)
         panel.lmline(x,y)
       },  as.table=T)

##### Model construction for Lung_Function_Growth  #### 


#### REEMTree Model ####
reem_model<-REEMtree(fev ~height+age , data=fev1 ,random=~1|id,cv=TRUE)
reem_model
rpart.plot(reem_model$Tree,digits = 5,roundint=FALSE, main="Decision Tree of \n REEMtree Model")

reem_model$BetweenMatrix
reem_model$RandomEffects
reem_model$EffectModel
reem_predict<-fitted(reem_model) #predicted values
rmse(fev1$fev,reem_predict) # root mean square error
mean((fev1$fev - reem_predict)^2)   #mean square error
AIC(reem_model)
BIC(reem_model)

# predicted vs actual values plot of REEMTree Model
reem_actualvspredict<-ggplot(fev1, aes(x=reem_predict, y=fev)) + 
  geom_point() +
  geom_abline(intercept=0, slope=1) +
  labs(x='Fitted Values', y='Actual Values', title='Fitted vs. Actual Values for REEMTree')

reem_actualvspredict

#### Heavy_REEMTree Model ####

set.seed(8)
reem_model_heavy<-REEMtree_heavy_mode(fev ~height+age  , data=fev1 ,random=~1, groups=~id,cv=TRUE)
reem_model_heavy
rpart.plot(reem_model_heavy$Tree,digits = 5,roundint=FALSE, main="Decision Tree of \n Heavy_REEMtree Model")

reem_model_heavy$BetweenMatrix
reem_model_heavy$RandomEffects
reem_model_heavy$EffectModel
heavy_predict<-reem_model_heavy$fitted #predicted values
rmse(fev1$fev,heavy_predict) # root mean square error
mean((fev1$fev - reem_model_heavy$fitted)^2)  #mean square error
reem_model_heavy$AIC
reem_model_heavy$BIC

# predicted vs actual values plot of Heavy_REEMTree Model
heavy_actualvsfitted<-ggplot(fev1, aes(x=heavy_predict, y=fev)) + 
  geom_point() +
  geom_abline(intercept=0, slope=1) +
  labs(x='Fitted Values', y='Actual Values', title='Fitted vs. Actual Values for Heavy_REEMTree')

heavy_actualvsfitted

#### LME Model (LMM) ####

lme_model<-lme(fev ~age+height , data = fev1,random=~1|id)
lme_model
lme_model$logLik
VarCorr(lme_model)  ### estimated variance of errors
lme_predict<-fitted(lme_model) #predicted values
rmse(fev1$fev,fitted(lme_model)) # root mean square error
mean((fev1$fev - fitted(lme_model))^2)    #mean square error
AIC(lme_model)
BIC(lme_model)

# predicted vs actual values plot of LME Model (LMM)
lme_actualvspredict<-ggplot(fev1, aes(x=fitted(lme_model), y=fev)) + 
  geom_point() +
  geom_abline(intercept=0, slope=1) +
  labs(x='Fitted Values', y='Actual Values', title='Fitted vs. Actual Values for LME')

lme_actualvspredict





