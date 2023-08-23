generatet =function(n, df,mu, S, perout, gamma,
             outlierType = "casewise", seed = NULL)
{
  if (!outlierType %in% c("casewise", "cellwisePlain", "cellwiseStructured", 
                          "both")) {
    stop("outlierType should be one of \"casewise\", \"cellwisePlain\",\n    \"cellwiseStructured\" or \"both\"")
  }
  if (is.numeric(seed)) {
    set.seed(seed)
  }
  X = LaplacesDemon::rmvt(n,mu,S,df)
  Sigma = S
  d = length(mu)
  indcells <- c()
  indrows <- c()
  if (perout > 0) {
    if (outlierType == "casewise") {
      replacement <- eigen(Sigma)$vectors[, d]/sqrt(mahalanobis(eigen(Sigma)$vectors[, 
                                                                                     d], mu, Sigma))
      replacement <- gamma * replacement * sqrt(d) * d
      ind <- seq_len(floor(perout * n))
      X[ind, ] <- matrix(replacement, nrow = length(ind), 
                         ncol = d, byrow = TRUE)
      indrows <- ind
    }
    else if (outlierType == "cellwisePlain") {
      ind <- replicate(d, sample(seq_len(n), perout * n, 
                                 replace = FALSE))
      ind <- as.vector(t(t(ind) + n * (0:(d - 1))))
      X[ind] <- gamma
      indcells <- ind
    }
    else if (outlierType == "cellwiseStructured") {
      ind <- replicate(d, sample(seq_len(n), perout * n, 
                                 replace = FALSE))
      ind <- as.vector(t(t(ind) + n * (0:(d - 1))))
      W <- array(0, dim(X))
      W[ind] <- 1
      for (i in seq_len(n)) {
        continds <- which(W[i, ] == 1)
        if (length(continds) > 0) {
          eigen_out <- eigen(Sigma[continds, continds])$vectors
          replacement <- eigen_out[, length(continds)]/sqrt(mahalanobis(eigen_out[, 
                                                                                  length(continds)], mu[continds], Sigma[continds, 
                                                                                                                         continds]))
          X[i, continds] <- replacement * gamma * sqrt(length(continds))
          indcells <- c(indcells, i + (continds - 1) * 
                          n)
        }
      }
    }
    else if (outlierType == "both") {
      replacement <- eigen(Sigma)$vectors[, d]/sqrt(mahalanobis(eigen(Sigma)$vectors[, 
                                                                                     d], mu, Sigma))
      replacement <- gamma * replacement * d * sqrt(d)
      ind <- seq_len(floor(perout/2 * n))
      X[ind, ] <- matrix(replacement, nrow = length(ind), 
                         ncol = d, byrow = TRUE)
      indrows <- ind
      startind <- (floor(perout/2 * n) + 1)
      ind <- replicate(d, sample(startind:n, ceiling(perout/2 * 
                                                       n), replace = FALSE))
      ind <- as.vector(t(t(ind) + n * (0:(d - 1))))
      W <- array(0, dim(X))
      W[ind] <- 1
      for (i in startind:n) {
        continds <- which(W[i, ] == 1)
        if (length(continds) > 0) {
          eigen_out <- eigen(Sigma[continds, continds])$vectors
          replacement <- eigen_out[, length(continds)]/sqrt(mahalanobis(eigen_out[, 
                                                                                  length(continds)], mu[continds], Sigma[continds, 
                                                                                                                         continds]))
          X[i, continds] <- replacement * gamma * sqrt(length(continds))
          indcells <- c(indcells, i + (continds - 1) * 
                          n)
        }
      }
    }
  }
  return(list(X = X, indcells = indcells, indrows = indrows))
}




###################################################################################



generateC =function(n,mu, S, perout, gamma,
                    outlierType = "casewise", seed = NULL)
{
  if (!outlierType %in% c("casewise", "cellwisePlain", "cellwiseStructured", 
                          "both")) {
    stop("outlierType should be one of \"casewise\", \"cellwisePlain\",\n    \"cellwiseStructured\" or \"both\"")
  }
  if (is.numeric(seed)) {
    set.seed(seed)
  }
  X = LaplacesDemon::rmvc(n,mu,S)
  Sigma = S
  d = length(mu)
  indcells <- c()
  indrows <- c()
  if (perout > 0) {
    if (outlierType == "casewise") {
      replacement <- eigen(Sigma)$vectors[, d]/sqrt(mahalanobis(eigen(Sigma)$vectors[, 
                                                                                     d], mu, Sigma))
      replacement <- gamma * replacement * sqrt(d) * d
      ind <- seq_len(floor(perout * n))
      X[ind, ] <- matrix(replacement, nrow = length(ind), 
                         ncol = d, byrow = TRUE)
      indrows <- ind
    }
    else if (outlierType == "cellwisePlain") {
      ind <- replicate(d, sample(seq_len(n), perout * n, 
                                 replace = FALSE))
      ind <- as.vector(t(t(ind) + n * (0:(d - 1))))
      X[ind] <- gamma
      indcells <- ind
    }
    else if (outlierType == "cellwiseStructured") {
      ind <- replicate(d, sample(seq_len(n), perout * n, 
                                 replace = FALSE))
      ind <- as.vector(t(t(ind) + n * (0:(d - 1))))
      W <- array(0, dim(X))
      W[ind] <- 1
      for (i in seq_len(n)) {
        continds <- which(W[i, ] == 1)
        if (length(continds) > 0) {
          eigen_out <- eigen(Sigma[continds, continds])$vectors
          replacement <- eigen_out[, length(continds)]/sqrt(mahalanobis(eigen_out[, 
                                                                                  length(continds)], mu[continds], Sigma[continds, 
                                                                                                                         continds]))
          X[i, continds] <- replacement * gamma * sqrt(length(continds))
          indcells <- c(indcells, i + (continds - 1) * 
                          n)
        }
      }
    }
    else if (outlierType == "both") {
      replacement <- eigen(Sigma)$vectors[, d]/sqrt(mahalanobis(eigen(Sigma)$vectors[, 
                                                                                     d], mu, Sigma))
      replacement <- gamma * replacement * d * sqrt(d)
      ind <- seq_len(floor(perout/2 * n))
      X[ind, ] <- matrix(replacement, nrow = length(ind), 
                         ncol = d, byrow = TRUE)
      indrows <- ind
      startind <- (floor(perout/2 * n) + 1)
      ind <- replicate(d, sample(startind:n, ceiling(perout/2 * 
                                                       n), replace = FALSE))
      ind <- as.vector(t(t(ind) + n * (0:(d - 1))))
      W <- array(0, dim(X))
      W[ind] <- 1
      for (i in startind:n) {
        continds <- which(W[i, ] == 1)
        if (length(continds) > 0) {
          eigen_out <- eigen(Sigma[continds, continds])$vectors
          replacement <- eigen_out[, length(continds)]/sqrt(mahalanobis(eigen_out[, 
                                                                                  length(continds)], mu[continds], Sigma[continds, 
                                                                                                                         continds]))
          X[i, continds] <- replacement * gamma * sqrt(length(continds))
          indcells <- c(indcells, i + (continds - 1) * 
                          n)
        }
      }
    }
  }
  return(list(X = X, indcells = indcells, indrows = indrows))
}


##########################################################################


generateL =function(n,mu, S, perout, gamma,
                    outlierType = "casewise", seed = NULL)
{
  if (!outlierType %in% c("casewise", "cellwisePlain", "cellwiseStructured", 
                          "both")) {
    stop("outlierType should be one of \"casewise\", \"cellwisePlain\",\n    \"cellwiseStructured\" or \"both\"")
  }
  if (is.numeric(seed)) {
    set.seed(seed)
  }
  X = LaplacesDemon::rmvl(n,mu,Sigma = S)
  Sigma = S
  d = length(mu)
  indcells <- c()
  indrows <- c()
  if (perout > 0) {
    if (outlierType == "casewise") {
      replacement <- eigen(Sigma)$vectors[, d]/sqrt(mahalanobis(eigen(Sigma)$vectors[, 
                                                                                     d], mu, Sigma))
      replacement <- gamma * replacement * sqrt(d) * d
      ind <- seq_len(floor(perout * n))
      X[ind, ] <- matrix(replacement, nrow = length(ind), 
                         ncol = d, byrow = TRUE)
      indrows <- ind
    }
    else if (outlierType == "cellwisePlain") {
      ind <- replicate(d, sample(seq_len(n), perout * n, 
                                 replace = FALSE))
      ind <- as.vector(t(t(ind) + n * (0:(d - 1))))
      X[ind] <- gamma
      indcells <- ind
    }
    else if (outlierType == "cellwiseStructured") {
      ind <- replicate(d, sample(seq_len(n), perout * n, 
                                 replace = FALSE))
      ind <- as.vector(t(t(ind) + n * (0:(d - 1))))
      W <- array(0, dim(X))
      W[ind] <- 1
      for (i in seq_len(n)) {
        continds <- which(W[i, ] == 1)
        if (length(continds) > 0) {
          eigen_out <- eigen(Sigma[continds, continds])$vectors
          replacement <- eigen_out[, length(continds)]/sqrt(mahalanobis(eigen_out[, 
                                                                                  length(continds)], mu[continds], Sigma[continds, 
                                                                                                                         continds]))
          X[i, continds] <- replacement * gamma * sqrt(length(continds))
          indcells <- c(indcells, i + (continds - 1) * 
                          n)
        }
      }
    }
    else if (outlierType == "both") {
      replacement <- eigen(Sigma)$vectors[, d]/sqrt(mahalanobis(eigen(Sigma)$vectors[, 
                                                                                     d], mu, Sigma))
      replacement <- gamma * replacement * d * sqrt(d)
      ind <- seq_len(floor(perout/2 * n))
      X[ind, ] <- matrix(replacement, nrow = length(ind), 
                         ncol = d, byrow = TRUE)
      indrows <- ind
      startind <- (floor(perout/2 * n) + 1)
      ind <- replicate(d, sample(startind:n, ceiling(perout/2 * 
                                                       n), replace = FALSE))
      ind <- as.vector(t(t(ind) + n * (0:(d - 1))))
      W <- array(0, dim(X))
      W[ind] <- 1
      for (i in startind:n) {
        continds <- which(W[i, ] == 1)
        if (length(continds) > 0) {
          eigen_out <- eigen(Sigma[continds, continds])$vectors
          replacement <- eigen_out[, length(continds)]/sqrt(mahalanobis(eigen_out[, 
                                                                                  length(continds)], mu[continds], Sigma[continds, 
                                                                                                                         continds]))
          X[i, continds] <- replacement * gamma * sqrt(length(continds))
          indcells <- c(indcells, i + (continds - 1) * 
                          n)
        }
      }
    }
  }
  return(list(X = X, indcells = indcells, indrows = indrows))
}






###########################################################################
generateSN = function(n,d,location,scale,alpha,perout, gamma,
                      outlierType = "casewise", seed = NULL)
  {
    if (!outlierType %in% c("casewise", "cellwisePlain", "cellwiseStructured", 
                            "both")) {
      stop("outlierType should be one of \"casewise\", \"cellwisePlain\",\n    \"cellwiseStructured\" or \"both\"")
    }
    if (is.numeric(seed)) {
      set.seed(seed)
    }
    x = sn::rsn(n*d,xi = location,omega = scale,alpha = alpha)
    X = matrix(x,nrow = n,ncol = d)
    mu = numeric(d)+location
    Sigma = scale*(diag(d))
    
    indcells <- c()
    indrows <- c()
    if (perout > 0) {
      if (outlierType == "casewise") {
        replacement <- eigen(Sigma)$vectors[, d]/sqrt(mahalanobis(eigen(Sigma)$vectors[, 
                                                                                       d], mu, Sigma))
        replacement <- gamma * replacement * sqrt(d) * d
        ind <- seq_len(floor(perout * n))
        X[ind, ] <- matrix(replacement, nrow = length(ind), 
                           ncol = d, byrow = TRUE)
        indrows <- ind
      }
      else if (outlierType == "cellwisePlain") {
        ind <- replicate(d, sample(seq_len(n), perout * n, 
                                   replace = FALSE))
        ind <- as.vector(t(t(ind) + n * (0:(d - 1))))
        X[ind] <- gamma
        indcells <- ind
      }
      else if (outlierType == "cellwiseStructured") {
        ind <- replicate(d, sample(seq_len(n), perout * n, 
                                   replace = FALSE))
        ind <- as.vector(t(t(ind) + n * (0:(d - 1))))
        W <- array(0, dim(X))
        W[ind] <- 1
        for (i in seq_len(n)) {
          continds <- which(W[i, ] == 1)
          if (length(continds) > 0) {
            eigen_out <- eigen(Sigma[continds, continds])$vectors
            replacement <- eigen_out[, length(continds)]/sqrt(mahalanobis(eigen_out[, 
                                                                                    length(continds)], mu[continds], Sigma[continds, 
                                                                                                                           continds]))
            X[i, continds] <- replacement * gamma * sqrt(length(continds))
            indcells <- c(indcells, i + (continds - 1) * 
                            n)
          }
        }
      }
      else if (outlierType == "both") {
        replacement <- eigen(Sigma)$vectors[, d]/sqrt(mahalanobis(eigen(Sigma)$vectors[, 
                                                                                       d], mu, Sigma))
        replacement <- gamma * replacement * d * sqrt(d)
        ind <- seq_len(floor(perout/2 * n))
        X[ind, ] <- matrix(replacement, nrow = length(ind), 
                           ncol = d, byrow = TRUE)
        indrows <- ind
        startind <- (floor(perout/2 * n) + 1)
        ind <- replicate(d, sample(startind:n, ceiling(perout/2 * 
                                                         n), replace = FALSE))
        ind <- as.vector(t(t(ind) + n * (0:(d - 1))))
        W <- array(0, dim(X))
        W[ind] <- 1
        for (i in startind:n) {
          continds <- which(W[i, ] == 1)
          if (length(continds) > 0) {
            eigen_out <- eigen(Sigma[continds, continds])$vectors
            replacement <- eigen_out[, length(continds)]/sqrt(mahalanobis(eigen_out[, 
                                                                                    length(continds)], mu[continds], Sigma[continds, 
                                                                                                                           continds]))
            X[i, continds] <- replacement * gamma * sqrt(length(continds))
            indcells <- c(indcells, i + (continds - 1) * 
                            n)
          }
        }
      }
    }
    return(list(X = X, indcells = indcells, indrows = indrows))
  }



