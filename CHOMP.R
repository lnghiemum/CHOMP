library(Matrix)
library(glmnet)

computeProjection <- function(one_beta) {
  # compute the projection matrix corresponding to one vector or one matrix
  one_beta = as.matrix(one_beta)
  a = min(sqrt(apply(one_beta,2,function(j) sum(j^2))))
  if (rcond(one_beta)<1e-10) 
    return(NA)
  else {
    # None of the columns equal to 0
    Pr = one_beta %*% solve(crossprod(one_beta)) %*% t(one_beta)
    return(Pr)
  }
}

getKernelMatrix <- function(X, y, method = c("sir", "save", "phdy"), nslices = 10){
  # get Kernel Matrix that corresponds to inverse regression methods
  n <- length(y)
  p <- ncol(X)
  order_y <- order(y) 
  slice_size <- floor(n/nslices);
  
  # Divide the feature matrix into slices
  slicedData <- array(0, dim = c(nslices, slice_size, p))
  for (h in 1:nslices){
    begin_index = slice_size*(h-1)+1
    end_index = slice_size*h
    slicedData[h, , ] = X[order_y[begin_index:end_index], ] 
  }
  
  M <- list()
  
  if ("sir" %in% method){
    tmp <- apply(slicedData, c(1, 3), mean, na.rm = TRUE)
    Msir <- cov(tmp)
    M <- c(M, Msir = list(Msir))
    rm(tmp)
  } 
  
  if ("save" %in% method){
    tmp <- lapply(1:nslices, function(oneslice) {
      Q <- cov(X) - cov(slicedData[oneslice, ,])
      Q %*% Q
    })
    Msave <- Reduce("+", tmp)/nslices
    M <- c(M, Msave = list(Msave))
    rm(tmp)
  } 
  
  if ('phdy' %in% method){
    y = y - mean(y)
    Mphdy <- t(X*y) %*% X/n
    M <- c(M, Mphdy = list(Mphdy))
  } 
  return(M)
}

ICProjection = function(fit, beta_ref_PIC, tau){
  # compute the PIC. 
  # fit: CHOMP/Adaptive CHOMP fit object, on which we use PIC to select the tuning parameter and final estimator
  # beta_ref_PIC: an estimator, typically the unpenalized estimator, that is consistent
  # tau: model complexity term
  
  p <- nrow(fit$beta)
  G = beta_ref_PIC %*% solve(crossprod(beta_ref_PIC)) %*% t(beta_ref_PIC)
  
  frobnorm <- apply(fit$beta, 2, function(x){
    tmp <- computeProjection(x)
    if (is.na(sum(tmp))) return(NA) else return(norm(tmp - G, "F")^2)
  })
  
  maxbetafit <- apply(abs(fit$beta),2,max)
  ICseq <- p*frobnorm + tau*fit$df
  ICseq[fit$df == 0 | maxbetafit < 1e-10] <- max(ICseq) 
  betaIC1 = coef(fit, s = fit$lambda[which.min(ICseq)])[-1]
  
  return(list(betaIC1 = betaIC1,
              IC = min(ICseq), 
              frobnorm = frobnorm,
              ICseq = ICseq))
}

AdaptCHOMP.fit = function(R, kappa, beta_ref_weight = NULL, beta_ref_PIC, tau){
  # R = t(L) transpose of the Cholesky factor; 
  # kappa: pseudo-response computed corresponding to each inverse regression method
  # beta_ref_weight: the initial estimator that is used for defining adaptive weight; if not supplied, then (unadaptive) CHOMP is fit
  # beta_ref_PIC: the initial estimator that is used to compute PIC 
  # tau: penalty for PIC
  
  p <- dim(R)[1]
  R <- as.matrix(R)
  kappa <- as.numeric(kappa)
  
  if (is.null(beta_ref_weight)){
    beta_ref_weight <- rep(1, p)
  }
    
  beta_ref_weight[beta_ref_weight == 0] <- 1e-8
  m_adapt_old = glmnet(x = R, y = kappa, penalty.factor = 1/abs(beta_ref_weight),
                       intercept = FALSE)
  max_lambda = m_adapt_old$lambda[1]
  # Grid of tuning parameters 
  lambda_seq = exp(seq(max(log(max_lambda),1e-6), log(max(1e-5*max_lambda, 1e-5)),
                     len=100))
  # Fit the adaptive CHOMP estimates
  m_adapt <- glmnet(x = R, y = kappa, penalty.factor = 1/abs(beta_ref_weight),
                   intercept = FALSE,
                   lambda = lambda_seq)
  
  ICnew <- ICProjection(fit = m_adapt, beta_ref_PIC = beta_ref_PIC, tau = tau)
  betaICnew <- ICnew$betaIC1
  betas <- betaICnew
  return(list(betas = betas, IC = ICnew$IC))
}

AdaptCHOMPwithPIC <- function(X, y, method = "sir", nslices = 10, d, gamma.pow = 2, 
                              adaptive = TRUE){
  # X : a n \times p design matrix; y : n \times 1 outcome
  # method: base method, either "sir", "save", or "phd"
  # nslices: require for sir and save
  # d: number of dimensions
  # gamma.pow: the power for defining the Adaptive CHOMP weight estimator, default is 2
  # adaptive: whether adaptive CHOMP should be fitted, if FALSE, unweighted CHOMP is fit
  
  SIR.kernel <- getKernelMatrix(X, y, method = method, nslices = nslices) 
  eta <- eigen(SIR.kernel$Msir)$vectors[, 1:d, drop = FALSE]
  unpen.fit <- solve(cov(X), eta)
  
  # ChOMP and Adaptive CHOMP fit
  R <- chol(cov(X)) # R is an upper triangular matrix
  L <- t(R)
  kappa <- forwardsolve(L, eta)
  # We fit CHOMP for each dimension separately
  if (adaptive){
    fit <- sapply(1:d, function(one_dim){
      AdaptCHOMP.fit(R, kappa[, one_dim, drop = FALSE], beta_ref_weight = NULL, beta_ref_PIC = unpen.fit[, one_dim, drop = FALSE], tau = log(p))
    })
  }
  else{
    fit <- sapply(1:d, function(one_dim){
      AdaptCHOMP.fit(R, kappa[, one_dim, drop = FALSE], beta_ref_weight = unpen.fit[, one_dim, drop = FALSE]^gamma.pow, beta_ref_PIC = unpen.fit[, one_dim, drop = FALSE], tau = log(p))  
    })
  }
  # Adaptive CHOMP with weight defined based on SIR.fit  
  
  estBeta <- as(Reduce(cbind, fit['betas', ]), 'sparseMatrix')
  colnames(estBeta) <- paste("Dim", 1:d)
  return(estBeta)
}
# compute the unpenalized SIR, assuming d = 2 being known

