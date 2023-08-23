rm(list = ls())
gc(reset = TRUE)

##Now we will start our project
#1. Loading packages
#R packages
library(cellWise)
library(ggplot2)
library(caret)
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
library(xtable)


##########
#Defining the Global Variables
loc = 1
sca = 1
n_train = 400
d_train = 25
percent = 0.10
Mu_train1 = numeric(d_train)-loc
Sigma_train1 = diag(d_train)
Mu_train2 = numeric(d_train)+loc
Sigma_train2 = sca*diag(d_train)
n_test = 1000
d_test = d_train
Mu_test1 = numeric(d_test)-loc
Sigma_test1 = diag(d_test)
Mu_test2 = numeric(d_test)+loc
Sigma_test2 = sca*diag(d_test)
##Creating the
r = 5
gamma1 = (1:10)*10
gamma2 = (1:10)*10


gamma = expand.grid(gamma1,gamma2)

a1 = seq(0.1,1,length.out = 20)
a2 = seq(0.1,1,length.out = 20)
a3 = seq(0.1,1,length.out = 20)
a = expand.grid(a1,a2)
b = expand.grid(a1,a2,a3)

label = c(replicate(n_train,"1"),replicate(n_train,"2"))
Tlabel = c(replicate(n_test,"1"),replicate(n_test,"2"))

### Calculation for the Bayes Classifier
B = numeric(100)

for(i in 1:100)
{
  class_1 = rmvnorm(n_test,mean = Mu_test1,sigma = Sigma_test1)
  class_2 = rmvnorm(n_test,mean = Mu_test2,sigma = Sigma_test2)
  class = rbind(class_1,class_2)
  class = apply(X = class, MARGIN = c(1,2), FUN = as.numeric)
  

  
  Bayes = ifelse(dmvnorm(class,mean = Mu_test1,sigma = Sigma_test1)>
                   dmvnorm(class,mean = Mu_test2,sigma = Sigma_test2),"1","2")
  
  B[i] = mean(Bayes!=Tlabel)
}

mean(B)
sd(B)



DD = function(Train, label,a1, a2 = NULL, a3 = NULL)
{
  z0 = Train[,1]
  ##Degree = 1
  d1 = poly(Train[,2],degree = 1,raw = TRUE)
  z1 = d1%*%a1
  ypred1 = ifelse(z0>z1,"1","2")
  mis1 = data.frame(ypred1!=label)
  mis1 = apply(mis1,FUN = sum,MARGIN = 2)/length(label)
  
  ##Degree = 2
  a = expand.grid(a1,a2)
  d2 = poly(Train[,2],degree = 2,raw = TRUE)
  z2 = d2%*%t(a)
  ypred2 = ifelse(z0>z2,"1","2")
  mis2 = data.frame(ypred2!=label)
  mis2 = apply(mis2,FUN = sum,MARGIN = 2)/length(label)
  min(mis2)
  ##Degree = 3
  b = expand.grid(a1,a2,a3)
  d3 = poly(Train[,2],degree = 3,raw = TRUE)
  z3 = d3%*%t(b)
  ypred3 = ifelse(z0>z3,"1","2")
  mis3 = data.frame(ypred3!=label)
  mis3 = apply(mis3,FUN = sum,MARGIN = 2)/length(label)
  min_mis = c(min(mis1),min(mis2),min(mis3))
  pos = c(which.min(mis1),which.min(mis2),which.min(mis3))
  best_degree = which.min(min_mis)
  G = list(mis1 = mis1,mis2 = mis2,mis3 = mis3)
  return(G)
}



Predict = function(Test,Tlabel,degree,a)
{
  z0 = Test[,1]
  z = poly(Test[,2],degree,raw = TRUE)
  z1 = z%*%(as.numeric(a))
  predictions = ifelse(z0>z1,"1","2")
  mis = sum(predictions!=Tlabel)/length(Tlabel)
  return(mis)
}



num_folds = 10

KFoldCrossValidation = function(Train,label,num_folds,a1,a2,a3){
  # Create the folds using createFolds function
  folds = createFolds(label, k = num_folds)
  
  Mis1 = matrix(0,num_folds,length(a1))
  Mis2 = matrix(0,num_folds,length(a1)*length(a2))
  Mis3 = matrix(0,num_folds,length(a1)*length(a2)*length(a3))
  # Iterate over the folds
  for (fold in 1:num_folds) {
    # Get the training and testing indices for the current fold
    train_indices = folds[[fold]]
    
    # Create the training and testing datasets
    train_data = Train[train_indices, ]
    train_label = label[train_indices]
    
    # Train the model using the training set
    model = DD(Train = train_data,label = train_label,a1 = a1,a2  =a1,a3 = a1)
    
    # Calculate the performance metric for the current fold
    Mis1[fold,] = model$mis1
    Mis2[fold,] = model$mis2
    Mis3[fold,] = model$mis3
  }
  
  # Calculate the average performance metric across all folds
  M1 = apply(Mis1,FUN = mean,MARGIN = 2)
  M2 = apply(Mis2,FUN = mean,MARGIN = 2)
  M3 = apply(Mis3,FUN = mean,MARGIN = 2)
  
  # Calculate the misclassification rate and find the best degree
  min_mis = c(min(M1),min(M2),min(M3))
  pos = c(which.min(M1),which.min(M2),which.min(M3))
  best_degree = which.min(min_mis)
  
  return(list(pos = pos[best_degree], best_degree = best_degree))
}






Model = function(Train,Test,label = label,Tlabel = Tlabel){
  G = KFoldCrossValidation(Train = Train,label = label,num_folds = 5,a1 = a1,a2 = a2,a3 = a3)
  best_pos = G[[1]]
  best_degree =G[[2]]
  if(best_degree==1)
    best_params = a1[best_pos]
  if(best_degree==2)
    best_params = a[best_pos,]
  if(best_degree==3)
    best_params = b[best_pos,]
  best_params
  pre = Predict(Test,Tlabel,degree = best_degree,a = best_params)
  return(pre)
}



distance = function(x,mu,S){sqrt(mahalanobis(x,mu,S))}  ##Calculating Mahalanobis Distance
depth = function(x,mu,S){1/(1+sqrt(mahalanobis(x,mu,S)))}
est = function(x,estimate){x[which(is.na(x))] = estimate[which(is.na(x))]
return(x)
}
ncores = detectCores(logical = FALSE)-2

MC = matrix(0,nrow = length(gamma1),length(gamma2))
TS = matrix(0,nrow = length(gamma1),length(gamma2))
DI = matrix(0,nrow = length(gamma1),length(gamma2))
NR = matrix(0,nrow = length(gamma1),length(gamma2))
Ori = matrix(0,nrow = length(gamma1),length(gamma2))

myCluster <- makeCluster(2, type = "FORK") # why "FORK"?
registerDoParallel(myCluster)

tic()
##****Fitting Machine Learning Models*****
##*
##'[Preparing the Dataset]
##*

foreach(i=1:length(gamma1))%:%
  foreach(j=1:length(gamma2))%do%{
Data = function()
{
  DataMatrix1 = generateData(n = n_train,d  =d_train,mu = Mu_train1,S = Sigma_train1,
                             perout = percent,gamma = gamma1[i],outlierType = "cellwisePlain")$X
  Mu1_Mean = apply(DataMatrix1,FUN = mean,MARGIN = 2)
  Mu1_Mean
  Sigma1_S = cov(DataMatrix1)
  #cellMap(DataMatrix1,rowlabels = 1:n_train,nrowsinblock = 4,drawCircles = FALSE)
  DataMatrix1 = apply(DataMatrix1,c(1,2),as.numeric)
  estim1_MCD = cellMCD(DataMatrix1)
  D1_MCD = estim1_MCD$preds
  Mu1_MCD = estim1_MCD$mu              ##Estimate of Mu
  Sigma1_MCD = estim1_MCD$S           ##Estimate of Sigma
  
  #2SGS
  estim1_2SGS = TSGS(DataMatrix1)
  Mu1_2SGS = getLocation(estim1_2SGS)              ##Estimate of Mu
  Sigma1_2SGS = getScatter(estim1_2SGS)           ##Estimate of Sigma
  #DI
  estim1_DI = DI(DataMatrix1)
  Mu1_DI = estim1_DI$center              ##Estimate of Mu
  Sigma1_DI = estim1_DI$cov           ##Estimate of Sigma
  
  
  
  ##2.Train Data(For class-2)
  DataMatrix2 = generateData(n = n_train,d  =d_train,mu = Mu_train2,S = Sigma_train2,
                             perout = percent,gamma = gamma2[j],outlierType = "cellwisePlain")$X
  Mu2_Mean = apply(DataMatrix2,FUN = mean,MARGIN = 2)
  Sigma2_S = cov(DataMatrix2)
  #cellMap(DataMatrix2,rowlabels = 1:n_train,nrowsinblock = 4,drawCircles = FALSE)
  DataMatrix2 = apply(DataMatrix2,c(1,2),as.numeric)
  estim2_MCD = cellMCD(DataMatrix2)
  D2_MCD = estim2_MCD$preds
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
  DataMatrix = apply(DataMatrix,c(1,2),FUN = as.numeric)
  
  
  
  #Original
  X1_Original_Depth = apply(DataMatrix,FUN = depth,mu = Mu1_Mean,S = Sigma1_S,MARGIN = 1)
  X2_Original_Depth = apply(DataMatrix,FUN = depth,mu = Mu2_Mean,S = Sigma2_S,MARGIN = 1)
  X_Original_Depth = cbind(X1_Original_Depth,X2_Original_Depth)
  X_Original_Depth = apply(X_Original_Depth,FUN = as.numeric,MARGIN = c(1,2))   ##Final plot based on the training data_Depth
  
  
  #MCD_Depth
  X1_MCD_Depth = apply(DataMatrix,FUN = depth,mu = Mu1_MCD,S = Sigma1_MCD,MARGIN = 1)
  X2_MCD_Depth = apply(DataMatrix,FUN = depth,mu = Mu2_MCD,S = Sigma2_MCD,MARGIN = 1)
  X_MCD_Depth = cbind(X1_MCD_Depth,X2_MCD_Depth)
  X_MCD_Depth = apply(X_MCD_Depth,FUN = as.numeric,MARGIN = c(1,2))    ##Final plot based on the training data
  
  
  #2SGS_Depth
  X1_2SGS_Depth = apply(DataMatrix,FUN = depth,mu = Mu1_2SGS,S = Sigma1_2SGS,MARGIN = 1)
  X2_2SGS_Depth = apply(DataMatrix,FUN = depth,mu = Mu2_2SGS,S = Sigma2_2SGS,MARGIN = 1)
  X_2SGS_Depth = cbind(X1_2SGS_Depth,X2_2SGS_Depth)
  X_2SGS_Depth = apply(X_2SGS_Depth,FUN = as.numeric,MARGIN = c(1,2))##Final plot based on the training data
  
  #DI_Depth
  X1_DI_Depth = apply(DataMatrix,FUN = depth,mu = Mu1_DI,S = Sigma1_DI,MARGIN = 1)
  X2_DI_Depth = apply(DataMatrix,FUN = depth,mu = Mu2_DI,S = Sigma2_DI,MARGIN = 1)
  X_DI_Depth = cbind(X1_DI_Depth,X2_DI_Depth)
  X_DI_Depth = apply(X_DI_Depth,FUN = as.numeric,MARGIN = c(1,2))  ##Final plot based on the training data
  
  
  
  
  
  ##Test Data(For Class-1)
  #cellMap(TDataMatrix1,rowlabels = 1:n_test,nrowsinblock = 10,drawCircles = FALSE)
  TDataMatrix1 = generateData(n = n_test,d  =d_test,mu = Mu_test1,S = Sigma_test1,
                              perout = percent,gamma = gamma1[i],outlierType = "cellwisePlain")$X
  TDataMatrix1 = apply(TDataMatrix1,c(1,2),as.numeric)
  
  ##2.Test Data(For class-2)
  #cellMap(TDataMatrix2,rowlabels = 1:n_test,nrowsinblock = 10,drawCircles = FALSE)
  TDataMatrix2 = generateData(n = n_test,d  =d_test,mu = Mu_test2,S = Sigma_test2,
                              perout = percent,gamma = gamma2[j],outlierType = "cellwisePlain")$X
  TDataMatrix2 = apply(TDataMatrix2,c(1,2), as.numeric)
  
  
  
  
  
  TDataMatrix = rbind(TDataMatrix1,TDataMatrix2)
  TDataMatrix = apply(TDataMatrix,c(1,2),as.numeric)
  
  
  #Original_Depth
  TX1_Original_Depth = apply(TDataMatrix,FUN = depth,mu = Mu1_Mean,S = Sigma1_S,MARGIN = 1)
  TX2_Original_Depth = apply(TDataMatrix,FUN = depth,mu = Mu2_Mean,S = Sigma2_S,MARGIN = 1)
  TX_Original_Depth = cbind(TX1_Original_Depth,TX2_Original_Depth)
  TX_Original_Depth = apply(TX_Original_Depth,FUN = as.numeric,MARGIN = c(1,2))
  
  #MCD_Depth
  TX1_MCD_Depth = apply(TDataMatrix,FUN = depth,mu = Mu1_MCD,S = Sigma1_MCD,MARGIN = 1)
  TX2_MCD_Depth = apply(TDataMatrix,FUN = depth,mu = Mu2_MCD,S = Sigma2_MCD,MARGIN = 1)
  TX_MCD_Depth = cbind(TX1_MCD_Depth,TX2_MCD_Depth)
  TX_MCD_Depth = apply(TX_MCD_Depth,FUN = as.numeric,MARGIN = c(1,2))
  
  #2SGS_Depth
  TX1_2SGS_Depth = apply(TDataMatrix,FUN = depth,mu = Mu1_2SGS,S = Sigma1_2SGS,MARGIN = 1)
  TX2_2SGS_Depth = apply(TDataMatrix,FUN = depth,mu = Mu2_2SGS,S = Sigma2_2SGS,MARGIN = 1)
  TX_2SGS_Depth = cbind(TX1_2SGS_Depth,TX2_2SGS_Depth)
  TX_2SGS_Depth = apply(TX_2SGS_Depth,FUN = as.numeric,MARGIN = c(1,2))
  
  #DI_Depth
  TX1_DI_Depth = apply(TDataMatrix,FUN = depth,mu = Mu1_DI,S = Sigma1_DI,MARGIN = 1)
  TX2_DI_Depth = apply(TDataMatrix,FUN = depth,mu = Mu2_DI,S = Sigma2_DI,MARGIN = 1)
  TX_DI_Depth = cbind(TX1_DI_Depth,TX2_DI_Depth)
  TX_DI_Depth = apply(TX_DI_Depth,FUN = as.numeric,MARGIN = c(1,2))
  L = list(X_MCD_Depth,TX_MCD_Depth,X_2SGS_Depth,TX_2SGS_Depth,X_DI_Depth,TX_DI_Depth,
           X_Original_Depth,TX_Original_Depth,DataMatrix,TDataMatrix)
  toc()
  return(L)
}


M = matrix(0,r,5)
V = Data()
print(paste("The value of gamma1 is:", gamma1[i]))
print(paste("The value of gamma2 is:", gamma2[j]))
for(k in 1:r){
  M[k,1] = Model(Train = V[[1]],Test = V[[2]],label = label,Tlabel = Tlabel)
  M[k,2] = Model(V[[3]],V[[4]],label = label,Tlabel = Tlabel)
  M[k,3] = Model(V[[5]],V[[6]],label = label,Tlabel = Tlabel)
  M[k,4] = Model(V[[7]],V[[8]],label = label,Tlabel = Tlabel)
  M[k,5] = Model(V[[9]],V[[10]],label = label,Tlabel = Tlabel)
  print(paste("The Number of iteration:", k))
}
m = apply(M,FUN = mean,MARGIN = 2)

MC[i,j] = m[1]
TS[i,j] = m[2]
DI[i,j] = m[3]
NR[i,j] = m[4]
Ori[i,j] = m[5]

}

stopCluster(myCluster)
toc()







library(rgl)


# Create a new rgl scene
open3d()

# Plot the first surface
surface3d(gamma1, gamma2, 100*MC, color = "darkblue", alpha = 0.9)

# Plot the second surface
surface3d(gamma1, gamma2, 100*TS, color = "darkgreen", alpha = 0.9)
surface3d(gamma1, gamma2, 100*DI, color = "darkorange", alpha = 0.9)
surface3d(gamma1, gamma2, 100*NR, color = "darkviolet", alpha = 0.9)
surface3d(gamma1, gamma2, 100*Ori, color = "darkgrey", alpha = 0.9)

light3d(specular = "white", ambient = "#00000000")

# Add labels and title to the plot
axes3d()
title3d("3D Surfaces for Different Classifiers",xlab = "gamma1", ylab = "gamma2", zlab = "ME")

# Rotate the plot
view3d(theta = -0.5, phi = -80,fov = 60, zoom = 0.7)

rglwidget()
close3d()




MC = readRDS("MC.rds")
TS = readRDS("TS.rds")
DI = readRDS("DI.rds")
NR= readRDS("NR.rds")
Ori = readRDS("Ori.rds")