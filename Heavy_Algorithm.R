#This script includes function named REEMtree_heavy_mode codes and its helper functions.

#REEMtree_heavy_mode is adapted by Ibrahim Hakki Erduran from REEMtree function by Simonoff & Sela.


### Function: tree
# This is the generic (S3-type) of function that returns the tree associated
# with an object (such as a RE-EM tree or a FE-EM tree).
tree <- function(object,...){
  UseMethod("tree")
}



### Function: REEMtree_heavy_mode



# Inputs:
#	formula - a formula (as in the usual LM/rpart calls)
#	data - data (as in the usual LM/rpart calls)
#	random - the formula describing the random effects (as in heavyLme)
#	groups - a vector containing the group identifier for each observation
#	subset - subset of the data to use for fitting
#	initialRandomEffects [0] - a vector of initial values for random effects
#	ErrorTolerance [0.001]
#	MaxIterations [1000]
# 	likelihoodCheck [TRUE] - should the likelihood be used to check for
#		convergence?  (if not, the random effects are checked instead)
#	verbose [FALSE] - print intermediate trees

#	### These options pertain to the RPART part of estimation
#	tree.control [rpart.control] - controls to be passed through to rpart
#	cv [TRUE] - Should cross-validation be used?
#	cpmin [0.0001] - complexity parameter used in building a tree before cross-validation
#	cpcv [0.01] - complexity used for pruning in a cross-validated tree
#	no.SE [0] - number of standard errors used in pruning (0 if unused)


#	### These options pertain to the heavyLME part of estimation
#	heavy_control [heavy.control(fix.shape = TRUE)] - controls to be passed through to heavyLme
#	method ["REML"] - "ML" or "REML", depending on whether the effects should be estimated with maximum likelihood or restricted maximum likelihood
#	correlation [NULL] - an option CorStruct object describing the within-group correlation structure
# family_h -the error distribution to be used in the model. By default the Student-t distribution with 4 degrees of freedom is considered.

REEMtree_heavy_mode <- function(formula, data, random, subset=NULL, initialRandomEffects=rep(0,TotalObs), 
                       ErrorTolerance=0.001, MaxIterations=1000, verbose=FALSE, tree.control=rpart.control(), 
                       cv=TRUE, cpmin = 0.001, no.SE =1,
                       heavy_control=heavy.control(fix.shape = TRUE),groups,method="REML",family_h=Student(df = 4)){
  TotalObs <- dim(data)[1]
  
  originaldata <- data
  
  # Subset the data if necessary
  if(identical(subset, NULL)){ 
    subs <- rep(TRUE, dim(data)[1])
  } else {
    subs <- subset
  }
  
  # Parse formula
  Predictors <- paste(attr(terms(formula),"term.labels"),collapse="+")
  TargetName <- formula[[2]]
  # Remove the name of the data frame if necessary
  if(length(TargetName)>1) 
    TargetName <-TargetName[3]
  if(verbose) print(paste("Target variable: ", TargetName))
  
  Target <- data[,toString(TargetName)]
  
  
  # Condition that indicates the loop has not converged or 
  # run out of iterations
  ContinueCondition <- TRUE
  
  iterations <- 0
  
  # Get initial values
  AdjustedTarget <- Target - initialRandomEffects
  oldlik <- -Inf
  
  # Make a new data frame to include all the new variables
  newdata <- data
  newdata[, "SubsetVector"] <- subs
  
 
  
  while(ContinueCondition){
    
    # Current values of variables
    newdata[,"AdjustedTarget"] <- AdjustedTarget
    iterations <- iterations+1
    
    
    
    # Compute current tree
    if (cv) {
      tree1 <- rpart(formula(paste(c("AdjustedTarget", Predictors),
                                   collapse = "~")), data = newdata, subset = subs,
                     method = "anova", control = rpart.control(cp=cpmin))
      if (nrow(tree1$cptable)==1){
        tree <- tree1}
      else {
        cventry <- which.min(tree1$cptable[, "xerror"])
        if (no.SE == 0){
          cpcv <- tree1$cptable[cventry, "CP"]
          tree <- prune(tree1, cp=cpcv)}
        else {
          xerrorcv <- tree1$cptable[cventry, "xerror"]
          sexerrorcv <- xerrorcv + tree1$cptable[cventry, "xstd"] * no.SE
          cpcvse <- tree1$cptable[which.max(tree1$cptable[, "xerror"] <= sexerrorcv), "CP"]
          tree <- prune(tree1, cp=cpcvse)}
      }
    }
    else {
      tree <- rpart(formula(paste(c("AdjustedTarget", Predictors),
                                  collapse = "~")), data = newdata, subset = subs,
                    method = "anova", control = tree.control)
    }
    if(verbose) print(tree)
    
    ## Estimate New Random Effects and Errors using heavyLme
    # Get variables that identify the node for each observation
    newdata[,"nodeInd"] <- 0
    newdata[subs,"nodeInd"] <- tree$where
   
    # Fit linear model with nodes as predictors (we use the original target so likelihoods are comparable)
    # Check that the fitted tree has at least two nodes.
    if(min(tree$where)==max(tree$where)){
      heavyfit <- heavyLme(formula(paste(c(toString(TargetName),1), collapse="~")), data=newdata, random=random, 
                      control=heavy_control,groups=groups,family=family_h)
      
      
      
    } else {
      heavyfit <- heavyLme(formula(paste(c(toString(TargetName),"as.factor(nodeInd)"), collapse="~")), data=newdata, random=random, 
                           control=heavy_control, groups=groups,family=family_h)
      
    }
    adjtarg <-  unique(cbind(tree$where, heavyfit$Fitted[,2]))
    tree$frame[adjtarg[,1],]$yval <- adjtarg[,2]
    
    
    if(verbose){
      print(heavyfit)
      print(paste("Estimated Error Variance = ", heavyfit$scale))
      print("Estimated Random Effects Variance = ")
    }
    
    # Get the likelihood to check on convergence
    newlik <- heavyfit$logLik
    if(verbose) print(paste("Log likelihood: ", newlik))
    
    ContinueCondition <- (newlik-oldlik>ErrorTolerance & iterations < MaxIterations)
    oldlik <- newlik
  

    
    # Extract random effects to make the new adjusted target
   AllEffects <- heavyfit$Resid[,1]-heavyfit$Resid[,dim(heavyfit$Resid)[2]]
   AdjustedTarget[subs] <- Target[subs] - AllEffects
    

    
  }
  
  residuals <- rep(NA, length=length(Target))
  residuals[subs] <- Target[subs]-heavyfit$Fitted[,2]
  
  attr(residuals, "label") <- NULL
  
  
  adjtarg <- unique(cbind(tree$where, heavyfit$Fitted[,2]))
  tree$frame[adjtarg[,1],]$yval <- adjtarg[,2]
  
Aic<- -2*newlik+(2*length(formula) ) + (2*length(formula) *(length(formula) +1)/(dim(data)[1]-length(formula) -1)) # Akaike information criterion calculation
Bic<- -2 * newlik + length(formula)  * log(dim(data)[1])   # Bayesian information criterion  calculation
 
  result <- list(Tree=tree, EffectModel=heavyfit, RandomEffects=heavyfit$ranef,
                 BetweenMatrix=heavyfit$theta,
                 ErrorVariance=heavyfit$scale^2, data=data, logLik=newlik,
                 IterationsUsed=iterations, Formula=formula, Random=random, Subset=subs,
                 ErrorTolerance=ErrorTolerance, 
                 residuals=residuals, method=method, cv=cv, heavy_control=heavy_control, tree.control=tree.control,
                fitted=heavyfit$Fitted[,2],adjusted_tar=AdjustedTarget,AIC=Aic,BIC=Bic)
  class(result) <- "REEMtree_heavy"
  
  return(result)
}



















































































































































