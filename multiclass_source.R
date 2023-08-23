rm(list = ls())
gc(reset = TRUE)
library("reticulate")
use_python("/usr/bin/python3.10",required = TRUE)
reticulate::py_config()
#repl_python()
##Now we will start our project
#1. Loading packages
#R packages
library(cellWise)
library(GSE)
library(data.table)
library(mvtnorm)
library(parallel)
library(foreach)
library(doParallel)
library(ggplot2)
library(tictoc)
library(beepr)
library(xtable)
#Python Packages
np = import("numpy")
pd = import("pandas")
mat = import("matplotlib")
mat$use('WebAgg')
plt = import("matplotlib.pyplot")
sk = import("sklearn")
ski = import("sklearnex")
ski$patch_sklearn()
import_main(convert = TRUE)
import_builtins(convert = TRUE)
##########
#Defining the Global Variables


depth = function(x,mu,S){
  z = 1/(1+sqrt(mahalanobis(x,mu,S)))
  return(z)
}

ncores = detectCores(logical = FALSE)-2

rowindex = function(d){
  rownames(d) = 1:nrow(d)
  return(d)
}
ColNames = function(d){
  names(d) = c("Mean","Sd")
  return(d)
}
Depth = function(est){
  z = apply(Train[,-ncol(Train)],FUN = function(x)depth(x,mu = est$Mean,S=est$Sd),MARGIN = 1)
  return(z)
}
TDepth = function(est){
  z = apply(Test[,-ncol(Test)],FUN = depth,mu = est$Mean,S=est$Sd,MARGIN = 1)
  return(z)
}




Data = function(Train = data.table(),Test = data.table()){
  
  D = split(Train,factor(Train$label))
  D = lapply(D,FUN = function(data)subset(data,select = -label))
  D = lapply(D,FUN = rowindex)
  D = lapply(D,FUN = function(d)apply(d,FUN = as.numeric,MARGIN = c(1,2)))
  
  
  
  nr = lapply(D,FUN = function(data)list(apply(data,MARGIN = 2,FUN = mean),cov(data)))
  nr = lapply(nr,FUN = ColNames)
  
  mcd = lapply(D, FUN = function(d)cellMCD(d,alpha = 0.4))
  mcd = lapply(mcd,function(est)list(est$mu,est$S))
  mcd = lapply(mcd,FUN = ColNames)
  
  tsgs = lapply(D, FUN = TSGS)
  tsgs = lapply(tsgs,function(est)list(getLocation(est),getScatter(est)))
  tsgs = lapply(tsgs,FUN = ColNames)
  
  di = lapply(D, FUN = DI)
  di = lapply(di,function(est)list(est$center,est$cov))
  di = lapply(di,FUN = ColNames)
  
  X_Original_Depth = lapply(nr,FUN = Depth)
  X_Original_Depth = data.table(data.frame(X_Original_Depth))
  
  X_MCD_Depth = lapply(mcd,FUN = Depth)
  X_MCD_Depth = data.table(data.frame(X_MCD_Depth))
  
  X_TSGS_Depth = lapply(tsgs,FUN = Depth)
  X_TSGS_Depth = data.table(data.frame(X_TSGS_Depth))
  
  X_DI_Depth = lapply(di,FUN = Depth)
  X_DI_Depth = data.table(data.frame(X_DI_Depth))
  
  
  
  
  
  TD = split(Test,factor(Test$Tlabel))
  TD = lapply(TD,FUN = function(data)subset(data,select = -Tlabel))
  TD = lapply(TD,FUN = function(d)apply(d,FUN = as.numeric,MARGIN = c(1,2)))
  
  TX_Original_Depth = lapply(nr,FUN = TDepth)
  TX_Original_Depth = data.table(data.frame(TX_Original_Depth))
  
  TX_MCD_Depth = lapply(mcd,FUN = TDepth)
  TX_MCD_Depth = data.table(data.frame(TX_MCD_Depth))
  
  TX_TSGS_Depth = lapply(tsgs,FUN = TDepth)
  TX_TSGS_Depth = data.table(data.frame(TX_TSGS_Depth))
  
  TX_DI_Depth = lapply(di,FUN = TDepth)
  TX_DI_Depth = data.table(data.frame(TX_DI_Depth))
  
  L = list(X_MCD_Depth,TX_MCD_Depth,X_TSGS_Depth,TX_TSGS_Depth,X_DI_Depth,TX_DI_Depth,
           X_Original_Depth,TX_Original_Depth,Train[,-ncol(Train)],Test[,-ncol(Test)])
  return(L)
}
# da = Data(Train,Test)
# np_train = da[[1]]
# np_test = da[[2]]



##****Fitting Machine Learning Models*****
##*
##'[Defining Function]
##*

Model = function(np_train,np_test,priors)
{
  np_train <- np$asarray(np_train, dtype = "float64")
  np_test <- np$asarray(np_test, dtype = "float64")
  
  
  Value = numeric(6)
  ##LDA(reduced Data)
  LDA = sk$discriminant_analysis$LinearDiscriminantAnalysis(priors = np$array(priors))
  LDA_fit = LDA$fit(np$array(np_train),np$array(label))
  LDA_predict = LDA$predict(np$array(np_test))
  Value[1] = 1-LDA$score(np$array(np_test),np$array(Tlabel))
  ##############################################################################
  
  ##############################################################################
  ##QDA(reduced data)
  QDA = sk$discriminant_analysis$QuadraticDiscriminantAnalysis(priors = np$array(priors))
  QDA_fit = QDA$fit(np$array(np_train),np$array(label))
  QDA_predict = QDA$predict(np$array(np_test))
  Value[2] = 1-QDA$score(np$array(np_test),np$array(Tlabel))
  #############################################################################
  ##KNN(reduced data)
  CV = sk$model_selection$GridSearchCV
  score = sk$model_selection$cross_val_score
  KF = sk$model_selection$KFold
  
  tic()
  KNN = sk$neighbors$KNeighborsClassifier
  K = list("n_neighbors" = as.integer(seq(1,25,1)))
  KNN_CV = CV(KNN(n_jobs = as.integer(-1)),K,cv = KF(n_splits = as.integer(10),shuffle = TRUE),n_jobs = as.integer(-1))
  KNN_Fit = KNN_CV$fit(np$array(np_train),np$array(label))
  KNN_Fit$best_params_
  KNN_predict = KNN_Fit$predict(np$array(np_test))
  Value[3] = 1-KNN_CV$score(np$array(np_test),np$array(Tlabel))
  toc()
  ############################################################################
  tic()
  SVML = sk$svm$SVC
  param_SVCL = list("C" = as.numeric(seq(0.1,5,0.5)),"degree" = list(as.integer(1)))
  SVML = CV(SVML(),param_SVCL,cv = KF(n_splits = as.integer(10),shuffle = TRUE),n_jobs = as.integer(-1))
  SVML_Fit = SVML$fit(np$array(np_train),np$array(label))
  SVML_predict = SVML$predict(np$array(np_test))
  Value[4] = 1-SVML$score(np$array(np_test),np$array(Tlabel))
  toc()
  ############################################################################
  ##SVM radial(reduced data)
  tic()
  SVMR = sk$svm$SVC
  param_SVCR = list("C" = as.numeric(seq(0.1,5,0.5)),"kernel" = list("rbf"),gamma = as.numeric(seq(0.1,5,0.5)))
  SVMR = CV(SVMR(),param_SVCR,cv = KF(n_splits = as.integer(10),shuffle = TRUE),n_jobs = as.integer(-1))
  SVMR_Fit = SVMR$fit(np$array(np_train),np$array(label))
  SVMR_predict = SVMR$predict(np$array(np_test))
  Value[5] = 1-SVMR$score(np$array(np_test),np$array(Tlabel))
  toc()
  ############################################################################
  ##SVM polynomial(reduced data)
  tic()
  SVMN = sk$svm$SVC(max_iter = as.integer(100))
  param_SVCN = list("kernel" = list("poly"),"C" = as.numeric(seq(0.1,2,0.5)),"degree" = as.integer(seq(2,5,1)))
  SVMN = CV(SVMN,param_SVCN,cv = KF(n_splits = as.integer(10),shuffle = TRUE),n_jobs = as.integer(-1))
  SVMN_Fit = SVMN$fit(np$array(np_train),np$array(label))
  SVMN_predict = SVMN$predict(np$array(np_test))
  Value[6] = 1-SVMN$score(np$array(np_test),np$array(Tlabel))
  toc()
  return(Value)
}
#Model(np_train,np_test,priors)
