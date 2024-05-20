

rm(list = ls())

library(heavy)
library(rpart.plot)
library(REEMtree)
library(nlme)
library(rpart)
library(ggplot2)
library(dplyr)
library(Metrics)
library(qpcR)
library(lattice)
source("Heavy_Algorithm.R")



#### Descriptive Statistics and Data Preparation for Alcohol_Usage ####

alcohol1 <- read.table("https://stats.idre.ucla.edu/stat/r/examples/alda/data/alcohol1_pp.txt", header=T, sep=",")
alcohol1$alcuse[1]<-alcohol1$alcuse[1]+10 #add noise 
alcohol1$id<-as.factor(alcohol1$id)
dim(alcohol1)
head(alcohol1)
summary(alcohol1)
table(alcohol1$id)
table(alcohol1$coa)
table(alcohol1$male)

### Lattice Graphs 20 by 20 ###

#1:20

index<-1:20
sample1 <- alcohol1[alcohol1$id %in% index, ] 
sample1<-as.data.frame(sample1)

sample1$id<-as.factor(sample1$id)

xyplot(alcuse~age | id, 
       data=sample1, 
       panel=function(x,y){
         panel.xyplot(x, y)
         panel.lmline(x,y)
       },  as.table=T)

#21:40

index<-21:40
sample2 <- alcohol1[alcohol1$id %in% index, ] 
sample2<-as.data.frame(sample2)

sample2$id<-as.factor(sample2$id)

xyplot(alcuse~age | id, 
       data=sample2, 
       panel=function(x,y){
         panel.xyplot(x, y)
         panel.lmline(x,y)
       },  as.table=T)


#41:60
index<-41:60
sample3 <- alcohol1[alcohol1$id %in% index, ] 
sample3<-as.data.frame(sample3)

sample3$id<-as.factor(sample3$id)

xyplot(alcuse~age | id, 
       data=sample3, 
       panel=function(x,y){
         panel.xyplot(x, y)
         panel.lmline(x,y)
       },  as.table=T)


#61:80
index<-61:80
sample4 <- alcohol1[alcohol1$id %in% index, ] 
sample4<-as.data.frame(sample4)

sample4$id<-as.factor(sample4$id)

xyplot(alcuse~age | id, 
       data=sample4, 
       panel=function(x,y){
         panel.xyplot(x, y)
         panel.lmline(x,y)
       },  as.table=T)


#81:100
index<-81:100
sample5 <- alcohol1[alcohol1$id %in% index, ] 
sample5<-as.data.frame(sample5)

sample5$id<-as.factor(sample5$id)

xyplot(alcuse~age | id, 
       data=sample5, 
       panel=function(x,y){
         panel.xyplot(x, y)
         panel.lmline(x,y)
       },  as.table=T)


##### Model construction for Alcohol_Usage  #### 

#### REEMTree Model ####

REEMresult_alch<-REEMtree(alcuse ~ age + coa + male +peer , data=alcohol1, random=~1|id,cv=TRUE)
REEMresult_alch
rpart.plot(REEMresult_alch$Tree,digits = 5,roundint=FALSE, main="Decision Tree of \n REEMtree Model")


REEMresult_alch$BetweenMatrix
REEMresult_alch$RandomEffects
REEMresult_alch$EffectModel
reem_predict_alch<-fitted(REEMresult_alch) #predicted values
rmse(alcohol1$alcuse,reem_predict_alch) # root mean square error
mean((alcohol1$alcuse - reem_predict_alch)^2)   #mean square error
AIC(REEMresult_alch)
BIC(REEMresult_alch)

# predicted vs actual values plot of REEMTree Model 

reem_actualvspredict_alch<-ggplot(alcohol1, aes(x=reem_predict_alch, y=alcuse)) + 
  geom_point() +
  geom_abline(intercept=0, slope=1) +
  labs(x='Fitted Values', y='Actual Values', title='Fitted vs. Actual Values for REEMTree')

reem_actualvspredict_alch


#### Heavy_REEMTree Model ####
set.seed(22)
REEMresult_heavy_alch<-REEMtree_heavy_mode(alcuse~age + coa + male +peer ,data=alcohol1,random= ~1, groups= ~id,cv=TRUE)
REEMresult_heavy_alch
rpart.plot(REEMresult_heavy_alch$Tree,digits = 5,roundint=FALSE, main="Decision Tree of \n Heavy_REEMtree Model")


REEMresult_heavy_alch$BetweenMatrix
REEMresult_heavy_alch$RandomEffects
REEMresult_heavy_alch$EffectModel
heavy_predict_alch<-REEMresult_heavy_alch$fitted #predicted values
rmse(alcohol1$alcuse,heavy_predict_alch) # root mean square error
mean((alcohol1$alcuse - REEMresult_heavy_alch$fitted)^2)  #mean square error
REEMresult_heavy_alch$AIC
REEMresult_heavy_alch$BIC

# predicted vs actual values plot of Heavy_REEMTree Model
heavy_actualvsfitted_alch<-ggplot(alcohol1, aes(x=heavy_predict_alch, y=alcuse)) + 
  geom_point() +
  geom_abline(intercept=0, slope=1) +
  labs(x='Fitted Values', y='Actual Values', title='Fitted vs. Actual Values for Heavy_REEMTree')

heavy_actualvsfitted_alch

#### LME Model (LMM) ####

lme_model_alch<-lme(alcuse ~age + coa + male +peer , data = alcohol1,random=~1|id)
lme_model_alch
lme_model_alch$logLik


VarCorr(lme_model_alch)  ### estimated variance of errors

lme_predict_alch<-fitted(lme_model_alch) #predicted values
rmse(alcohol1$alcuse,fitted(lme_model_alch)) # root mean square error
mean((alcohol1$alcuse - fitted(lme_model_alch))^2)    #mean square error
AIC(lme_model_alch)
BIC(lme_model_alch)

# predicted vs actual values plot of LME Model (LMM)

lme_actualvspredict_alch<-ggplot(alcohol1, aes(x=lme_predict_alch, y=alcuse)) + 
  geom_point() +
  geom_abline(intercept=0, slope=1) +
  labs(x='Predicted Values', y='Actual Values', title='Predicted vs. Actual Values')

lme_actualvspredict_alch



