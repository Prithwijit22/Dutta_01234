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
library(GSE)
library(mvtnorm)
library(parallel)
library(foreach)
library(doParallel)
library(LaplacesDemon)
library(tictoc)
library(beepr)
library(xtable)
#Python Packages
np = import("numpy")
pd = import("pandas")
sk = import("sklearn")
import_main(convert = TRUE)
import_builtins(convert = TRUE)
##########
#Defining the Global Variables

loc = 2
sca = 5
n_train = 400
d_train = 25
percent = 0.20
Mu_train1 = numeric(d_train)
Sigma_train1 = diag(d_train)
Mu_train2 = numeric(d_train)+loc
Sigma_train2 = sca*diag(d_train)
n_test = 1000
d_test = d_train
Mu_test1 = numeric(d_test)
Sigma_test1 = diag(d_test)
Mu_test2 = numeric(d_test)+loc
Sigma_test2 = sca*diag(d_test)
##Creating the
r = 1
gamma = 10
df = 5


label = c(replicate(n_train,"1"),replicate(n_train,"2"))
Tlabel = c(replicate(n_test,"1"),replicate(n_test,"2"))

#distance = function(x,mu,S){sqrt(mahalanobis(x,mu,S))}  ##Calculating Mahalanobis Distance
depth = function(x,mu,S){1/(1+sqrt(mahalanobis(x,mu,S)))}
est = function(x,estimate){x[which(is.na(x))] = estimate[which(is.na(x))]
return(x)
}
ncores = detectCores(logical = FALSE)-2
tic()
##****Fitting Machine Learning Models*****
##*
##'[Preparing the Dataset]
##*
Data = function()
{
  #1.Train Data(For Class - 1)
  registerDoParallel(4)
  tic()
  DataMatrix1 = generatet(n = n_train,df = df,mu = Mu_train1,S = Sigma_train1,
                          perout = percent,gamma = gamma,outlierType = "cellwisePlain")$X
  Mu1_Mean = apply(DataMatrix1,FUN = mean,MARGIN = 2)
  Mu1_Mean
  Sigma1_S = cov(DataMatrix1)
  #cellMap(DataMatrix1,rowlabels = 1:n_train,nrowsinblock = 4,drawCircles = FALSE)
  DataMatrix1 = apply(DataMatrix1,c(1,2),as.numeric)
  estim1_MCD = cellMCD(DataMatrix1,alpha = 0.6)
  D1_MCD = estim1_MCD$preds
  Mu1_MCD = estim1_MCD$mu              ##Estimate of Mu
  Sigma1_MCD = estim1_MCD$S           ##Estimate of Sigma
  
  #2SGS
  estim1_2SGS = TSGS(DataMatrix1)
  Mu1_2SGS = getLocation(estim1_2SGS)              ##Estimate of Mu
  Sigma1_2SGS = getScatter(estim1_2SGS)           ##Estimate of Sigma
  D1_2sgs = getFiltDat(estim1_2SGS)
  #na = matrix(NA,nrow = nrow(D),ncol = ncol(D))
  D1_2sgs = t(apply(D1_2sgs,FUN = est,estimate = Mu1_2SGS,MARGIN = 1))
  
  #DI
  estim1_DI = DI(DataMatrix1)
  D1_DI = estim1_DI$Ximp
  Mu1_DI = estim1_DI$center              ##Estimate of Mu
  Sigma1_DI = estim1_DI$cov           ##Estimate of Sigma
  
  
  
  
  
  ##2.Train Data(For class-2)
  DataMatrix2 = rmvt(n = n_train,df = df,mu = Mu_train2,S = Sigma_train2)
  Mu2_Mean = apply(DataMatrix2,FUN = mean,MARGIN = 2)
  Mu2_Mean
  Sigma2_S = cov(DataMatrix2)
  #cellMap(DataMatrix2,rowlabels = 1:n_train,nrowsinblock = 4,drawCircles = FALSE)
  DataMatrix2 = apply(DataMatrix2,c(1,2),as.numeric)
  estim2_MCD = cellMCD(DataMatrix2,alpha =0.6)
  D2_MCD = estim2_MCD$preds
  Mu2_MCD = estim2_MCD$mu              ##Estimate of Mu
  Sigma2_MCD = estim2_MCD$S           ##Estimate of Sigma
  #2SGS
  estim2 = TSGS(DataMatrix2)
  Mu2_2SGS = getLocation(estim2)              ##Estimate of Mu
  Sigma2_2SGS = getScatter(estim2)          ##Estimate of Sigma
  D2_2sgs = getFiltDat(estim2)
  D2_2sgs = t(apply(D2_2sgs,FUN = est,estimate = Mu2_2SGS,MARGIN = 1))
  
  
  #DI
  estim2_DI = DI(DataMatrix2)
  D2_DI = estim2_DI$Ximp
  Mu2_DI = estim2_DI$center              ##Estimate of Mu
  Sigma2_DI = estim2_DI$cov           ##Estimate of Sigma
  
  D_MCD = rbind(D1_MCD,D2_MCD)
  D_MCD = as.data.frame(scale(D_MCD,center = TRUE,scale = TRUE))
  D_MCD = apply(D_MCD,c(1,2),FUN = as.numeric)
  D_2sgs = rbind(D1_2sgs,D2_2sgs)
  D_2sgs = as.data.frame(scale(D_2sgs,center = TRUE,scale = TRUE))
  D_2sgs = apply(D_2sgs,c(1,2),FUN = as.numeric)
  D_DI = rbind(D1_DI,D2_DI)
  D_DI = as.data.frame(scale(D_DI,center = TRUE,scale = TRUE))
  D_DI = apply(D_DI,c(1,2),FUN = as.numeric)
  
  
  
  DataMatrix = rbind(DataMatrix1,DataMatrix2)
  DataMatrix = apply(DataMatrix,c(1,2),FUN = as.numeric)
  D = as.data.frame(scale(DataMatrix,center = TRUE,scale = TRUE)) 
  
  
  
  
  #Original
  X1_Original_Depth = apply(D,FUN = depth,mu = Mu1_Mean,S = Sigma1_S,MARGIN = 1)
  X2_Original_Depth = apply(D,FUN = depth,mu = Mu2_Mean,S = Sigma2_S,MARGIN = 1)
  X_Original_Depth = cbind(X1_Original_Depth,X2_Original_Depth)
  X_Original_Depth = apply(X_Original_Depth,FUN = as.numeric,MARGIN = c(1,2))
  plot(X_Original_Depth,col = label,pch = 16)    ##Final plot based on the training data_Depth
  
  
  
  #MCD_Depth
  X1_MCD_Depth = apply(D_MCD,FUN = depth,mu = Mu1_MCD,S = Sigma1_MCD,MARGIN = 1)
  X2_MCD_Depth = apply(D_MCD,FUN = depth,mu = Mu2_MCD,S = Sigma2_MCD,MARGIN = 1)
  X_MCD_Depth = cbind(X1_MCD_Depth,X2_MCD_Depth)
  X_MCD_Depth = apply(X_MCD_Depth,FUN = as.numeric,MARGIN = c(1,2))
  plot(X_MCD_Depth,col = label,pch = 16)    ##Final plot based on the training data
  
  
  #2SGS_Depth
  X1_2SGS_Depth = apply(D_2sgs,FUN = depth,mu = Mu1_2SGS,S = Sigma1_2SGS,MARGIN = 1)
  X2_2SGS_Depth = apply(D_2sgs,FUN = depth,mu = Mu2_2SGS,S = Sigma2_2SGS,MARGIN = 1)
  X_2SGS_Depth = cbind(X1_2SGS_Depth,X2_2SGS_Depth)
  X_2SGS_Depth = apply(X_2SGS_Depth,FUN = as.numeric,MARGIN = c(1,2))
  plot(X_2SGS_Depth,col = label,pch = 16)    ##Final plot based on the training data
  
  #DI_Depth
  X1_DI_Depth = apply(D_DI,FUN = depth,mu = Mu1_DI,S = Sigma1_DI,MARGIN = 1)
  X2_DI_Depth = apply(D_DI,FUN = depth,mu = Mu2_DI,S = Sigma2_DI,MARGIN = 1)
  X_DI_Depth = cbind(X1_DI_Depth,X2_DI_Depth)
  X_DI_Depth = apply(X_DI_Depth,FUN = as.numeric,MARGIN = c(1,2))
  plot(X_DI_Depth,col = label,pch = 16)    ##Final plot based on the training data
  
  Mu_Mean = (Mu1_Mean+Mu2_Mean)/2
  
  
  
  ##Test Data(For Class-1)
  TDataMatrix1 = generatet(n = n_test,df = df,mu = Mu_test1,S = Sigma_test1,
                           perout = percent,gamma = gamma,outlierType = "cellwisePlain")$X
  #cellMap(TDataMatrix1,rowlabels = 1:n_test,nrowsinblock = 10,drawCircles = FALSE)
  TDataMatrix1 = apply(TDataMatrix1,c(1,2),as.numeric)
  
  ##2.Test Data(For class-2)
  TDataMatrix2 = rmvt(n = n_test,df =df,mu = Mu_test2,S = Sigma_test2)
  #cellMap(TDataMatrix2,rowlabels = 1:n_test,nrowsinblock = 10,drawCircles = FALSE)
  TDataMatrix2 = apply(TDataMatrix2,c(1,2), as.numeric)
  
  TDataMatrix = rbind(TDataMatrix1,TDataMatrix2)
  TDataMatrix = apply(TDataMatrix,c(1,2),as.numeric)
  TD = as.data.frame(scale(TDataMatrix,center = TRUE,scale = TRUE)) 
  
  Testim_MCD = cellMCD(TDataMatrix,alpha = 0.6)
  TD_MCD = Testim_MCD$preds
  TD_MCD = as.data.frame(scale(TD_MCD,center = TRUE,scale = TRUE))
  TD_MCD = apply(TD_MCD,c(1,2),FUN = as.numeric)  
  #2SGS
  Testim_2SGS = TSGS(TDataMatrix)
  TD_2sgs = getFiltDat(Testim_2SGS)
  TD_2sgs = t(apply(TD_2sgs,FUN = est,estimate = Mu_Mean,MARGIN = 1))
  TD_2sgs = as.data.frame(scale(TD_2sgs,center = TRUE,scale = TRUE))
  TD_2sgs = apply(TD_2sgs,c(1,2),FUN = as.numeric)  
  #DI
  Testim_DI = DI(TDataMatrix)
  TD_DI = Testim_DI$Ximp
  TD_DI = as.data.frame(scale(TD_DI,center = TRUE,scale = TRUE))
  TD_DI = apply(TD_DI,c(1,2),FUN = as.numeric)
  
  
  
  #Original_Depth
  TX1_Original_Depth = apply(TD,FUN = depth,mu = Mu1_Mean,S = Sigma1_S,MARGIN = 1)
  TX2_Original_Depth = apply(TD,FUN = depth,mu = Mu2_Mean,S = Sigma2_S,MARGIN = 1)
  TX_Original_Depth = cbind(TX1_Original_Depth,TX2_Original_Depth)
  TX_Original_Depth = apply(TX_Original_Depth,FUN = as.numeric,MARGIN = c(1,2))
  plot(TX_Original_Depth,col = Tlabel,pch = 16) 
  
  #MCD_Depth
  TX1_MCD_Depth = apply(TD_MCD,FUN = depth,mu = Mu1_MCD,S = Sigma1_MCD,MARGIN = 1)
  TX2_MCD_Depth = apply(TD_MCD,FUN = depth,mu = Mu2_MCD,S = Sigma2_MCD,MARGIN = 1)
  TX_MCD_Depth = cbind(TX1_MCD_Depth,TX2_MCD_Depth)
  TX_MCD_Depth = apply(TX_MCD_Depth,FUN = as.numeric,MARGIN = c(1,2))
  plot(TX_MCD_Depth,col = Tlabel,pch = 16)
  
  #2SGS_Depth
  TX1_2SGS_Depth = apply(TD_2sgs,FUN = depth,mu = Mu1_2SGS,S = Sigma1_2SGS,MARGIN = 1)
  TX2_2SGS_Depth = apply(TD_2sgs,FUN = depth,mu = Mu2_2SGS,S = Sigma2_2SGS,MARGIN = 1)
  TX_2SGS_Depth = cbind(TX1_2SGS_Depth,TX2_2SGS_Depth)
  TX_2SGS_Depth = apply(TX_2SGS_Depth,FUN = as.numeric,MARGIN = c(1,2))
  plot(TX_2SGS_Depth,col = Tlabel,pch = 16)
  
  #DI_Depth
  TX1_DI_Depth = apply(TD_DI,FUN = depth,mu = Mu1_DI,S = Sigma1_DI,MARGIN = 1)
  TX2_DI_Depth = apply(TD_DI,FUN = depth,mu = Mu2_DI,S = Sigma2_DI,MARGIN = 1)
  TX_DI_Depth = cbind(TX1_DI_Depth,TX2_DI_Depth)
  TX_DI_Depth = apply(TX_DI_Depth,FUN = as.numeric,MARGIN = c(1,2))
  plot(TX_DI_Depth,col = Tlabel,pch = 16)
  L = list(X_MCD_Depth,TX_MCD_Depth,X_2SGS_Depth,TX_2SGS_Depth,X_DI_Depth,TX_DI_Depth,
           X_Original_Depth,TX_Original_Depth,DataMatrix,TDataMatrix)
  toc()
  return(L)
}




##****Fitting Machine Learning Models*****
##*
##'[Defining Function]
##*

Model = function(Train,Test)
{
  
  registerDoParallel(cores = ncores)
  Value = numeric(7)
  ##LDA(reduced Data)
  LDA = sk$discriminant_analysis$LinearDiscriminantAnalysis()
  LDA_fit = LDA$fit_transform(np$array(Train),np$array(label))
  LDA_predict = LDA$predict(np$array(Test))
  Value[1] = 1-LDA$score(Test,Tlabel)
  ##############################################################################
  
  ##############################################################################
  ##QDA(reduced data)
  QDA = sk$discriminant_analysis$QuadraticDiscriminantAnalysis()
  QDA_fit = QDA$fit(np$array(Train),np$array(label))
  QDA_predict = QDA$predict(np$array(Test))
  Value[2] = 1-QDA$score(Test,Tlabel)
  #############################################################################
  ##KNN(reduced data)
  CV = sk$model_selection$GridSearchCV
  score = sk$model_selection$cross_val_score
  KF = sk$model_selection$KFold
  
  tic()
  KNN = sk$neighbors$KNeighborsClassifier
  K = list("n_neighbors" = as.integer(seq(1,25,1)))
  KNN_CV = CV(KNN(n_jobs = as.integer(-1)),K,cv = KF(n_splits = as.integer(10),shuffle = TRUE),n_jobs = as.integer(-1))
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
  RF_CV = CV(RF(n_jobs = as.integer(-1)),B,cv = KF(n_splits = as.integer(10),shuffle = TRUE),n_jobs = as.integer(-1))
  RF_Fit = RF_CV$fit(np$array(Train),np$array(label))
  RF_Fit$best_params_
  RF_predict = RF_CV$predict(np$array(Test))
  Value[4] = 1-RF_CV$score(Test,Tlabel)
  toc()
  ############################################################################
  tic()
  SVML = sk$svm$SVC
  param_SVCL = list("C" = as.numeric(seq(0.1,10,0.5)),"degree" = list(as.integer(1)))
  SVML = CV(SVML(),param_SVCL,cv = KF(n_splits = as.integer(10),shuffle = TRUE),n_jobs = as.integer(-1))
  SVML_Fit = SVML$fit(np$array(Train),np$array(label))
  SVML_predict = SVML$predict(np$array(Test))
  Value[5] = 1-SVML$score(Test,Tlabel)
  toc()
  ############################################################################
  ##SVM radial(reduced data)
  tic()
  SVMR = sk$svm$SVC
  param_SVCR = list("C" = as.numeric(seq(0.1,10,0.5)),"kernel" = list("rbf"))
  SVMR = CV(SVMR(),param_SVCR,cv = KF(n_splits = as.integer(10),shuffle = TRUE),n_jobs = as.integer(-1))
  SVMR_Fit = SVMR$fit(np$array(Train),np$array(label))
  SVMR_predict = SVMR$predict(np$array(Test))
  Value[6] = 1-SVMR$score(Test,Tlabel)
  toc()
  ############################################################################
  ##SVM polynomial(reduced data)
  tic()
  SVMN = sk$svm$SVC
  param_SVCN = list("C" = as.numeric(seq(0.1,10,0.5)),"degree" = as.integer(seq(2,5,1)))
  SVMN = CV(SVMN(),param_SVCN,cv = KF(n_splits = as.integer(10),shuffle = TRUE),n_jobs = as.integer(-1))
  SVMN_Fit = SVMN$fit(np$array(Train),np$array(label))
  SVMN_predict = SVMN$predict(np$array(Test))
  Value[7] = 1-SVMN$score(Test,Tlabel)
  toc()
  return(Value)
}
tic()
Final = np$zeros(list(as.integer(7),as.integer(5),as.integer(r)))
cl = makePSOCKcluster(2)
registerDoParallel(cl,cores = 2)
foreach(i=1:r) %do%
  {
    
    
    D = Data()
    print(paste("Number of iterations:",i))
    Final[,1,i] = Model(D[[1]],D[[2]])
    Final[,2,i] = Model(D[[3]],D[[4]])
    Final[,3,i] = Model(D[[5]],D[[6]])
    Final[,4,i] = Model(D[[7]],D[[8]])
    Final[,5,i] = Model(D[[9]],D[[10]])
    print(paste("Number of iterations:",i))
    
    
  }
stopCluster(cl)
toc()
beep(4)

ME = matrix(0,14,5)
Mean = np$mean(Final,axis = as.integer(2))
Sd = np$round(np$std(Final,axis = as.integer(2)),as.integer(4))
Mean[1,]

for(i in 1:7)
{
  ME[2*i-1,] = Mean[i,]*100
  ME[2*i,] = Sd[i,]
}
rownames(ME) = c("LDA","","QDA","","KNN","","RF","","SVML","","SVMR","","SVMN","")
colnames(ME) = c("MCD","2SGS","DI","Non-Robust","NoTransformation")
ME = apply(ME,MARGIN = c(1,2),as.numeric)
print(ME)
xtable(ME,digits = 5)
pbPost("note","Complilation Done")


