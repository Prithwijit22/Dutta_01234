rm(list = ls())
gc(reset = TRUE)
source("generatet.R",local = TRUE)
library("reticulate")
use_python("/usr/bin/python3.10",required = TRUE)
reticulate::py_config()
#repl_python()
##Now we will start our project
#1. Loading packages
#R packages
library(cellWise)
library(ggplot2)
library(robustbase)
library(lattice)
#library(caret)
library(GSE)
library(mvtnorm)
library(parallel)
library(foreach)
library(doParallel)
library(tictoc)
library(beepr)
library(RPushbullet)
library(xtable)
#Python Packages
np = import("numpy")
pd = import("pandas")
sci = import("scipy")
mat = import("matplotlib")
mat$use('WebAgg')
plt = import("matplotlib.pyplot")
sk = import("sklearn")
hm = import("hummingbird.ml")
import_main(convert = TRUE)
import_builtins(convert = TRUE)
##########




label = c(replicate(n_train,"1"),replicate(n_train,"2"))
Tlabel = c(replicate(n_test,"1"),replicate(n_test,"2"))

distance = function(x,mu,S){sqrt(mahalanobis(x,mu,S))}  ##Calculating Mahalanobis Distance
depth = function(x,mu,S){1/(1+sqrt(mahalanobis(x,mu,S)))}  ##Calculating Mahalanobis Depth
ncores = detectCores(logical = FALSE)-2
tic()
##****Fitting Machine Learning Models*****
##*
##'[Preparing the Dataset]
##*
Data = function(DataMatrix1,DataMatrix2,TDataMatrix1,TDataMatrix2)
{
  
  #1.Train Data(For Class - 1)
  registerDoParallel(4)
  tic()
  Mu1_Mean = apply(DataMatrix1,FUN = mean,MARGIN = 2)
  Sigma1_S = cov(DataMatrix1)
  #cellMap(DataMatrix1,rowlabels = 1:n_train,nrowsinblock = 4,drawCircles = FALSE)
  DataMatrix1 = apply(DataMatrix1,FUN = as.numeric,MARGIN = 2)
  estim1_MCD = cellMCD(DataMatrix1,alpha = 0.6)
  Mu1_MCD = estim1_MCD$mu              ##Estimate of Mu
  Sigma1_MCD = estim1_MCD$S           ##Estimate of Sigma
  typeof(Sigma1_MCD[1,1])
  #2SGS
  estim1_2SGS = TSGS(DataMatrix1)
  Mu1_2SGS = getLocation(estim1_2SGS)              ##Estimate of Mu
  Sigma1_2SGS = getScatter(estim1_2SGS)           ##Estimate of Sigma
  #DI
  estim1_DI = DI(DataMatrix1)
  Mu1_DI = estim1_DI$center              ##Estimate of Mu
  Sigma1_DI = estim1_DI$cov           ##Estimate of Sigma
  
  
  
  
  
  ##2.Train Data(For class-2)
  Mu2_Mean = apply(DataMatrix2,FUN = mean,MARGIN = 2)
  Sigma2_S = cov(DataMatrix2)
  #cellMap(DataMatrix2,rowlabels = 1:n_train,nrowsinblock = 4,drawCircles = FALSE)
  DataMatrix2 = apply(DataMatrix2,FUN = as.numeric,MARGIN = 2)
  estim2_MCD = cellMCD(DataMatrix2,alpha = 0.6)
  Mu2_MCD = estim2_MCD$mu              ##Estimate of Mu
  Sigma2_MCD = estim2_MCD$S           ##Estimate of Sigma
  #2SGS
  estim2 = TSGS(DataMatrix2)
  Mu2_2SGS = getLocation(estim2)              ##Estimate of Mu
  Sigma2_2SGS = getScatter(estim2)          ##Estimate of Sigma
  #DI
  estim2_DI = DI(DataMatrix2)
  Mu2_DI = estim2_DI$center              ##Estimate of Mu
  Sigma2_DI = estim2_DI$cov           ##Estimate of Sigma
  
  
  
  
  
  DataMatrix = rbind(DataMatrix1,DataMatrix2)
  DataMatrix = apply(DataMatrix,FUN = as.numeric,MARGIN = 2)
  
  #Original
  X1_Original = apply(DataMatrix,FUN = distance,mu = Mu1_Mean,S = Sigma1_S,MARGIN = 1)
  X2_Original = apply(DataMatrix,FUN = distance,mu = Mu2_Mean,S = Sigma2_S,MARGIN = 1)
  X_Original = cbind(X1_Original,X2_Original)
  X_Original = apply(X_Original,FUN = as.numeric,MARGIN = 2)
  plot(X_Original,col = label,pch = 16)    ##Final plot based on the training data
  
  #MCD
  X1_MCD = apply(DataMatrix,FUN = distance,mu = Mu1_MCD,S = Sigma1_MCD,MARGIN = 1)
  X2_MCD = apply(DataMatrix,FUN = distance,mu = Mu2_MCD,S = Sigma2_MCD,MARGIN = 1)
  X_MCD = cbind(X1_MCD,X2_MCD)
  X_MCD = apply(X_MCD,FUN = as.numeric,MARGIN = 2)
  plot(X_MCD,col = label,pch = 16)    ##Final plot based on the training data
  #2SGS
  X1_2SGS = apply(DataMatrix,FUN = distance,mu = Mu1_2SGS,S = Sigma1_2SGS,MARGIN = 1)
  X2_2SGS = apply(DataMatrix,FUN = distance,mu = Mu2_2SGS,S = Sigma2_2SGS,MARGIN = 1)
  X_2SGS = cbind(X1_2SGS,X2_2SGS)
  X_2SGS = apply(X_2SGS,FUN = as.numeric,MARGIN = 2)
  plot(X_2SGS,col = label,pch = 16)    ##Final plot based on the training data
  #DI
  
  X1_DI = apply(DataMatrix,FUN = distance,mu = Mu1_DI,S = Sigma1_DI,MARGIN = 1)
  X2_DI = apply(DataMatrix,FUN = distance,mu = Mu2_DI,S = Sigma2_DI,MARGIN = 1)
  X_DI = cbind(X1_DI,X2_DI)
  X_DI = apply(X_DI,FUN = as.numeric,MARGIN = 2)
  #plot(X_DI,col = label,pch = 16)    ##Final plot based on the training data
  
  
  
  
  
  
  
  
  
  
  #Original_Depth
  X1_Original_Depth = apply(DataMatrix,FUN = depth,mu = Mu1_Mean,S = Sigma1_S,MARGIN = 1)
  X2_Original_Depth = apply(DataMatrix,FUN = depth,mu = Mu2_Mean,S = Sigma2_S,MARGIN = 1)
  X_Original_Depth = cbind(X1_Original_Depth,X2_Original_Depth)
  X_Original_Depth = apply(X_Original_Depth,FUN = as.numeric,MARGIN = 2)
  plot(X_Original_Depth,col = label,pch = 16)    ##Final plot based on the training data_Depth
  
  
  
  #MCD_Depth
  X1_MCD_Depth = apply(DataMatrix,FUN = depth,mu = Mu1_MCD,S = Sigma1_MCD,MARGIN = 1)
  X2_MCD_Depth = apply(DataMatrix,FUN = depth,mu = Mu2_MCD,S = Sigma2_MCD,MARGIN = 1)
  X_MCD_Depth = cbind(X1_MCD_Depth,X2_MCD_Depth)
  X_MCD_Depth = apply(X_MCD_Depth,FUN = as.numeric,MARGIN = 2)
  plot(X_MCD_Depth,col = label,pch = 16)    ##Final plot based on the training data
  
  
  #2SGS_Depth
  X1_2SGS_Depth = apply(DataMatrix,FUN = depth,mu = Mu1_2SGS,S = Sigma1_2SGS,MARGIN = 1)
  X2_2SGS_Depth = apply(DataMatrix,FUN = depth,mu = Mu2_2SGS,S = Sigma2_2SGS,MARGIN = 1)
  X_2SGS_Depth = cbind(X1_2SGS_Depth,X2_2SGS_Depth)
  X_2SGS_Depth = apply(X_2SGS_Depth,FUN = as.numeric,MARGIN = 2)
  #plot(X_2SGS_Depth,col = label,pch = 16)    ##Final plot based on the training data
  
  #DI_Depth
  X1_DI_Depth = apply(DataMatrix,FUN = depth,mu = Mu1_DI,S = Sigma1_DI,MARGIN = 1)
  X2_DI_Depth = apply(DataMatrix,FUN = depth,mu = Mu2_DI,S = Sigma2_DI,MARGIN = 1)
  X_DI_Depth = cbind(X1_DI_Depth,X2_DI_Depth)
  X_DI_Depth = apply(X_DI_Depth,FUN = as.numeric,MARGIN = 2)
  #plot(X_DI_Depth,col = label,pch = 16)    ##Final plot based on the training data
  
  
  
  
  
  ##Test Data(For Class-1)
  #cellMap(TDataMatrix1,rowlabels = 1:n_test,nrowsinblock = 10,drawCircles = FALSE)
  TDataMatrix1 = apply(TDataMatrix1,FUN = as.numeric,MARGIN = 2)
  
  
  
  
  ##2.Test Data(For class-2)
  #cellMap(TDataMatrix2,rowlabels = 1:n_test,nrowsinblock = 10,drawCircles = FALSE)
  TDataMatrix2 = apply(TDataMatrix2,FUN = as.numeric,MARGIN = 2)
  
  
  
  
  TDataMatrix = rbind(TDataMatrix1,TDataMatrix2)
  TDataMatrix = apply(TDataMatrix,FUN = as.numeric,MARGIN = 2)
  
  #Original
  TX1_Original = apply(TDataMatrix,FUN = distance,mu = Mu1_Mean,S = Sigma1_S,MARGIN = 1)
  TX2_Original = apply(TDataMatrix,FUN = distance,mu = Mu2_Mean,S = Sigma2_S,MARGIN = 1)
  TX_Original = cbind(TX1_Original,TX2_Original)
  TX_Original = apply(TX_Original,FUN = as.numeric,MARGIN = 2)
  plot(TX_Original,col = Tlabel,pch = 16)    ##Final plot based on the training data
  
  
  #MCD
  TX1_MCD = apply(TDataMatrix,FUN = distance,mu = Mu1_MCD,S = Sigma1_MCD,MARGIN = 1)
  TX2_MCD = apply(TDataMatrix,FUN = distance,mu = Mu2_MCD,S = Sigma2_MCD,MARGIN = 1)
  TX_MCD = cbind(TX1_MCD,TX2_MCD)
  TX_MCD = apply(TX_MCD,FUN = as.numeric,MARGIN = 2)
  plot(TX_MCD,col = Tlabel,pch = 16)
  #2SGS
  TX1_2SGS = apply(TDataMatrix,FUN = distance,mu = Mu1_2SGS,S = Sigma1_2SGS,MARGIN = 1)
  TX2_2SGS = apply(TDataMatrix,FUN = distance,mu = Mu2_2SGS,S = Sigma2_2SGS,MARGIN = 1)
  TX_2SGS = cbind(TX1_2SGS,TX2_2SGS)
  TX_2SGS = apply(TX_2SGS,FUN = as.numeric,MARGIN = 2)
  plot(TX_2SGS,col = Tlabel,pch = 16)
  #DI
  TX1_DI = apply(TDataMatrix,FUN = distance,mu = Mu1_DI,S = Sigma1_DI,MARGIN = 1)
  TX2_DI = apply(TDataMatrix,FUN = distance,mu = Mu2_DI,S = Sigma2_DI,MARGIN = 1)
  TX_DI = cbind(TX1_DI,TX2_DI)
  TX_DI = apply(TX_DI,FUN = as.numeric,MARGIN = 2)
  plot(TX_DI,col = Tlabel,pch = 16)
  
  #Original_Depth
  TX1_Original_Depth = apply(TDataMatrix,FUN = depth,mu = Mu1_Mean,S = Sigma1_S,MARGIN = 1)
  TX2_Original_Depth = apply(TDataMatrix,FUN = depth,mu = Mu2_Mean,S = Sigma2_S,MARGIN = 1)
  TX_Original_Depth = cbind(TX1_Original_Depth,TX2_Original_Depth)
  TX_Original_Depth = apply(TX_Original_Depth,FUN = as.numeric,MARGIN = 2)
  plot(TX_Original_Depth,col = Tlabel,pch = 16) 
  
  #MCD_Depth
  TX1_MCD_Depth = apply(TDataMatrix,FUN = depth,mu = Mu1_MCD,S = Sigma1_MCD,MARGIN = 1)
  TX2_MCD_Depth = apply(TDataMatrix,FUN = depth,mu = Mu2_MCD,S = Sigma2_MCD,MARGIN = 1)
  TX_MCD_Depth = cbind(TX1_MCD_Depth,TX2_MCD_Depth)
  TX_MCD_Depth = apply(TX_MCD_Depth,FUN = as.numeric,MARGIN = 2)
  plot(TX_MCD_Depth,col = Tlabel,pch = 16)
  
  #2SGS_Depth
  TX1_2SGS_Depth = apply(TDataMatrix,FUN = depth,mu = Mu1_2SGS,S = Sigma1_2SGS,MARGIN = 1)
  TX2_2SGS_Depth = apply(TDataMatrix,FUN = depth,mu = Mu2_2SGS,S = Sigma2_2SGS,MARGIN = 1)
  TX_2SGS_Depth = cbind(TX1_2SGS_Depth,TX2_2SGS_Depth)
  TX_2SGS_Depth = apply(TX_2SGS_Depth,FUN = as.numeric,MARGIN = 2)
  plot(TX_2SGS_Depth,col = Tlabel,pch = 16)
  
  #DI_Depth
  TX1_DI_Depth = apply(TDataMatrix,FUN = depth,mu = Mu1_DI,S = Sigma1_DI,MARGIN = 1)
  TX2_DI_Depth = apply(TDataMatrix,FUN = depth,mu = Mu2_DI,S = Sigma2_DI,MARGIN = 1)
  TX_DI_Depth = cbind(TX1_DI_Depth,TX2_DI_Depth)
  TX_DI_Depth = apply(TX_DI_Depth,FUN = as.numeric,MARGIN = 2)
  plot(TX_DI_Depth,col = Tlabel,pch = 16)
  L = list(X_MCD,TX_MCD,X_2SGS,TX_2SGS,X_DI,TX_DI,X_Original,TX_Original,X_MCD_Depth,TX_MCD_Depth,
           X_2SGS_Depth,TX_2SGS_Depth,X_DI_Depth,TX_DI_Depth,X_Original_Depth,TX_Original_Depth,
           DataMatrix,TDataMatrix)
  toc()
  return(L)
}


mcd = function(data)
{
  x = cellMCD(data,alpha = 0.6)
  mu = x$mu
  S = x$S
  return(list(mu,S))
}


tsgs = function(data)
{
  x = TSGS(data)
  mu = getLocation(x)
  S = getScatter(x)
  return(list(mu,S))
}


di= function(data)
{
  x = DI(data)
  mu = x$center
  S = x$cov
  return(list(mu,S))
}








##****Fitting Machine Learning Models*****
##*
##'[Defining Function]
##*

Model = function(Train,Test,method = "usual")
{
  registerDoParallel(cores = ncores)
  Value = numeric(7)
  
  m_mu1 = mcd(Train[1:n_train,])[[1]]
  m_mu2 = mcd(Train[(n_train+1):(2*n_train),])[[1]]
  m_S = mcd(Train)[[2]]
  m_S1 = mcd(Train[1:n_train,])[[2]]
  m_S2 = mcd(Train[(n_train+1):(2*n_train),])[[2]]
  
  t_mu1 = tsgs(Train[1:n_train,])[[1]]
  t_mu2 = tsgs(Train[(n_train+1):(2*n_train),])[[1]]
  t_S = tsgs(Train)[[2]]
  t_S1 = tsgs(Train[1:n_train,])[[2]]
  t_S2 = tsgs(Train[(n_train+1):(2*n_train),])[[2]]
  
  
  
  d_mu1 = di(Train[1:n_train,])[[1]]
  d_mu2 = di(Train[(n_train+1):(2*n_train),])[[1]]
  d_S = di(Train)[[2]]
  d_S1 = di(Train[1:n_train,])[[2]]
  d_S2 = di(Train[(n_train+1):(2*n_train),])[[2]]  
  
  
  
  
  ##LDA(reduced Data)
  if (method == "usual")
  {
    LDA = sk$discriminant_analysis$LinearDiscriminantAnalysis()
    LDA_fit = LDA$fit_transform(np$array(Train),np$array(label))
    LDA_predict = LDA$predict(np$array(Test))
    Value[1] = 1-LDA$score(Test,Tlabel)
  }
  else if (method == "mcd")
  {
    LDA_predict = ifelse(mahalanobis(Test,m_mu1,m_S)>mahalanobis(Test,m_mu2,m_S),"2","1")
    Value[1] = 1-sum(Tlabel==LDA_predict)/(2*n_test)
  }
  else if (method == "tsgs")
  {
    LDA_predict = ifelse(mahalanobis(Test,t_mu1,t_S)>mahalanobis(Test,t_mu2,t_S),"2","1")
    Value[1] = 1-sum(Tlabel==LDA_predict)/(2*n_test)
  }
  else if (method == "di")
  {
    LDA_predict = ifelse(mahalanobis(Test,d_mu1,d_S)>mahalanobis(Test,d_mu2,d_S),"2","1")
    Value[1] = 1-sum(Tlabel==LDA_predict)/(2*n_test)
  }
  ##############################################################################
  
  ##############################################################################
  ##QDA(reduced data)
  if(method == "usual")
  {
    QDA = sk$discriminant_analysis$QuadraticDiscriminantAnalysis()
    QDA_fit = QDA$fit(np$array(Train),np$array(label))
    QDA_predict = QDA$predict(np$array(Test))
    Value[2] = 1-QDA$score(Test,Tlabel)
  }
  else if(method == "mcd")
  {
    QDA_predict = ifelse(mahalanobis(Test,m_mu1,m_S1)>mahalanobis(Test,m_mu2,m_S2),"2","1")
    Value[2] = 1-sum(Tlabel==QDA_predict)/(2*n_test)
  }
  else if(method == "tsgs")
  {
    QDA_predict = ifelse(mahalanobis(Test,t_mu1,t_S1)>mahalanobis(Test,t_mu2,t_S2),"2","1")
    Value[2] = 1-sum(Tlabel==QDA_predict)/(2*n_test)
  }
  else if(method == "di")
  {
    QDA_predict = ifelse(mahalanobis(Test,d_mu1,d_S1)>mahalanobis(Test,d_mu2,d_S2),"2","1")
    Value[2] = 1-sum(Tlabel==QDA_predict)/(2*n_test)
  }
  #############################################################################
  ##KNN(reduced data)
  tic()
  CV = sk$model_selection$GridSearchCV
  score = sk$model_selection$cross_val_score
  KF = sk$model_selection$KFold
  
  KNN = sk$neighbors$KNeighborsClassifier
  K = list("n_neighbors" = as.integer(seq(1,25,1)))
  KNN_CV = CV(KNN(),K,cv = KF(n_splits = as.integer(10),shuffle = TRUE))
  KNN_Fit = KNN_CV$fit(np$array(Train),np$array(label))
  KNN_Fit$best_params_
  KNN_predict = KNN_Fit$predict(np$array(Test))
  Value[3] = 1-KNN_CV$score(Test,Tlabel)
  toc()
  #############################################################################
  ##Random Forest(reduced data)
  tic()
  RF = sk$ensemble$RandomForestClassifier
  B = list("n_estimators"=as.integer(seq(20,50,5)))
  RF_CV = CV(RF(),B,cv = KF(n_splits = as.integer(10),shuffle = TRUE))
  RF_Fit = RF_CV$fit(np$array(Train),np$array(label))
  RF_Fit$best_params_
  RF_predict = RF_CV$predict(np$array(Test))
  Value[4] = 1-RF_CV$score(Test,Tlabel)
  toc()
  ############################################################################
  tic()
  SVML = sk$svm$SVC
  param_SVCL = list("C" = as.numeric(seq(0.1,10,0.5)),"degree" = list(as.integer(1)))
  SVML = CV(SVML(),param_SVCL,cv = KF(n_splits = as.integer(10),shuffle = TRUE))
  SVML_Fit = SVML$fit(np$array(Train),np$array(label))
  SVML_predict = SVML$predict(np$array(Test))
  Value[5] = 1-SVML$score(Test,Tlabel)
  toc()
  ############################################################################
  ##SVM radial(reduced data)
  tic()
  SVMR = sk$svm$SVC
  param_SVCR = list("C" = as.numeric(seq(0.1,10,0.5)),"kernel" = list("rbf"))
  SVMR = CV(SVMR(),param_SVCR,cv = KF(n_splits = as.integer(10),shuffle = TRUE))
  SVMR_Fit = SVMR$fit(np$array(Train),np$array(label))
  SVMR_predict = SVMR$predict(np$array(Test))
  Value[6] = 1-SVMR$score(Test,Tlabel)
  toc()
  ############################################################################
  ##SVM polynomial(reduced data)
  tic()
  SVMN = sk$svm$SVC
  param_SVCN = list("C" = as.numeric(seq(0.1,10,0.5)),"degree" = as.integer(seq(2,5,1)))
  SVMN = CV(SVMN(),param_SVCN,cv = KF(n_splits = as.integer(10),shuffle = TRUE))
  SVMN_Fit = SVMN$fit(np$array(Train),np$array(label))
  SVMN_predict = SVMN$predict(np$array(Test))
  Value[7] = 1-SVMN$score(Test,Tlabel)
  toc()
  return(Value)
}
tic()

Summary = function(DataMatrix1,DataMatrix2,TDataMatrix1,TDataMatrix2)
{
Final = np$zeros(list(as.integer(7),as.integer(9),as.integer(r)))
cl = makePSOCKcluster(2)
registerDoParallel(cl,cores = 2)
foreach(i=1:r) %do%
  {
    
    
    D = Data(DataMatrix1,DataMatrix2,TDataMatrix1,TDataMatrix2)
    print(paste("Number of iterations:",i))
    Final[,1,i] = Model(D[[1]],D[[2]],method = "mcd")
    Final[,2,i] = Model(D[[3]],D[[4]],method = "tsgs")
    Final[,3,i] = Model(D[[5]],D[[6]],method = "di")
    Final[,4,i] = Model(D[[7]],D[[8]])
    Final[,5,i] = Model(D[[9]],D[[10]],method = "mcd")
    Final[,6,i] = Model(D[[11]],D[[12]],method = "tsgs")
    Final[,7,i] = Model(D[[13]],D[[14]],method = "di")
    Final[,8,i] = Model(D[[15]],D[[16]])
    Final[,9,i] = Model(D[[17]],D[[18]])
    print(paste("Number of iterations:",i))
    
    
  }
stopCluster(cl)
toc()
beep(3)

ME = matrix(0,14,9)
Mean = np$mean(Final,axis = as.integer(2))
Sd = np$round(np$std(Final,axis = as.integer(2)),as.integer(4))
Mean[1,]

for(i in 1:7)
{
  ME[2*i-1,] = Mean[i,]
  ME[2*i,] = Sd[i,]
}
rownames(ME) = c("LDA","","QDA","","KNN","","RF","","SVML","","SVMR","","SVMN","")
colnames(ME) = c("MCD","2SGS","DI","Original","MCD","2SGS","DI","Original","Original")
print(ME)
print(xtable(ME,digits = 4))
return(ME)
}
pbPost("note","Complilation Done")

