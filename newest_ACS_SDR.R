library(Matrix)
library(MASS)
require(snowfall)
sfLibrary(glmnet)
library(e1071)
library(lattice)
library(parallel)
library(doParallel)
library(snowfall)
library(stringr)
library(lattice)
library(ggplot2)
library(glmnet)
library(doParallel)
library(energy)
library(msda)

AR <- function(rho, p){
  m <- matrix(0, p, p)
  for (i in 1:p){
    for (j in 1:p){
      m[i,j] <- rho**(abs(i-j))
    }
  }
  return(m)
}

# Cut small values in a matrix to zero. 
cut_mat <- function(Beta, thrd, rank){
  l <- length(Beta)
  for (i in 1:l){
    if(is.null(Beta[[i]])) next
    mat <- as.matrix(Beta[[i]])
    nobs <- nrow(mat)
    nvars <- ncol(mat)
    r <- rank[i]
    if(r == 0){
      Beta[[i]] <- matrix(0, nobs, nvars)
    }else{
      vec <- as.vector(mat)
      vec[abs(vec) < thrd] <- 0
      Beta[[i]] <- matrix(vec, nobs, nvars)
    }
  }
  return(Beta)
}

# Evaluation based on distance correlation. (Szekely et al., 2007)
eval_dc <- function(Beta, x, y){
  if(!is.list(Beta)){Beta <- list(Beta)}
  l <- length(Beta)
  result <- sapply(seq_len(l), function(i){
    if(is.null(Beta[[i]])){
      NA
    }else{
      mat <- as.matrix(Beta[[i]])
      dcor(x %*% mat, y)
    }
  })
  return(result)
}

# Compute M and U matrices from observation data.
#########
# Input:
# x: n x p observation matrix for predictor.
# y: n-dimensional observation vector for response.
# yclass: Discretized response taking values in 1,...,H.
# type: Specifying the specific SEAS method. "sir" means SEAS-SIR, "intra" means SEAS-Intra and "pfc" means SEAS-PFC.
# FUN: the user-specified function f in SEAS-PFC. The default is f(y) = (y, y^2, y^3).
# categorical: A logical value indicating whether y is categorical.
# H: The number of slices. The default value is 5.
MU <- function(x, y, yclass=NULL, type='sir', FUN = NULL, categorical = FALSE, H = 5){
  if(is.null(yclass)){ # Construct the discretized response
    if(categorical == FALSE){
      ybreaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H), na.rm=TRUE))
      yclass <- cut(y, breaks = ybreaks, include.lowest = TRUE, labels = FALSE)
      nclass <- as.integer(length(unique(yclass)))
    }
    else if(categorical == TRUE){
      yclass <- y
    }
  }
  cls <- sort(unique(yclass))
  nclass <- length(cls)
  nobs <- as.integer(dim(x)[1])
  nvars <- as.integer(dim(x)[2])
  prior <- sapply(cls, function(i){mean(yclass == i)})
  mu <- colMeans(x)
  x_c <- x - matrix(mu, nobs, nvars, byrow = TRUE) # centered predictor
  M <- crossprod(x_c/sqrt(nobs)) # sample covariance of X
  
  if(type == 'sir'){
    U <- matrix(0, nvars, nclass)
    for (i in 1:nclass){
      U[, i] <- colMeans(x_c[yclass == cls[i],, drop=FALSE])
    }
  }else if(type == 'intra'){
    y_c <- y - mean(y)
    U <- matrix(0, nvars, nclass)
    lb <- quantile(y_c, 0.1)[[1]]
    ub <- quantile(y_c, 0.9)[[1]]
    y_c <- sapply(y_c, cut_func, lb = lb, ub = ub)
    for (i in 1:nclass){
      y_copy <- y_c
      y_copy[yclass!=cls[i]] <- 0
      U[, i] <- (1/nobs) * t(x_c) %*% (y_copy - mean(y_copy))
    }
  }else if(type == 'pfc'){
    if(is.null(FUN)) Fmat <- cbind(y, y^2, y^3) # the default function
    else Fmat <- t(sapply(y, FUN))
    Fmat_mean <- colMeans(Fmat)
    Fmat_c <- Fmat - matrix(Fmat_mean, NROW(Fmat), NCOL(Fmat), byrow = TRUE) # centered function f
    lb <- apply(Fmat_c, 2, quantile, 0.1)
    ub <- apply(Fmat_c, 2, quantile, 0.9)
    for(i in 1:NCOL(Fmat_c)){
      Fmat_c[,i] <- sapply(Fmat_c[,i], cut_func, lb[i], ub[i])
    }
    U <- (1/nobs)*(t(x_c) %*% Fmat_c)
  }
  list(M = M, U = U, nclass = nclass, prior=prior)
}

# Cut extreme values in the samples. This function is used in MU function.
cut_func <- function(x, lb, ub){
  if(x < lb){
    return(lb)
  } else if(x > ub){
    return(ub)
  } else{
    return(x)
  }
}

# Estimate the rank of a matrix.
rank_func <- function(B, thrd){
  d <- svd(B)$d
  r <- sum(d >= thrd)
  return(r)
}

# Subspace distance, defined in (19)
subspace <- function(A,B){
  if(is.vector(A)) A <- as.matrix(A)
  if(is.vector(B)) A <- as.matrix(B)
  Pa <- qr.Q(qr(A))
  Pa <- Pa %*% t(Pa)
  Pb <- qr.Q(qr(B))
  Pb <- Pb %*% t(Pb)
  d <- dim(A)[2]
  return(norm(Pa-Pb, type="F")/sqrt(2*d))
}

## ------------------------------------------------ ##
## The utility functions imported from R package 'msda'.
## These functions are used in 'msda' and 'cv.msda' functions.
formatoutput <- function(fit, maxit, pmax, p, H) {
  nalam <- fit$nalam
  ntheta <- fit$ntheta[seq(nalam)]
  nthetamax <- max(ntheta)
  lam <- fit$alam[seq(nalam)]
  theta_vec <- fit$theta
  errmsg <- err(fit$jerr, maxit, pmax)  ### error messages from fortran
  switch(paste(errmsg$n), `1` = stop(errmsg$msg, call. = FALSE), `-1` = cat(errmsg$msg))
  if(nthetamax > 0){
    ja <- fit$itheta[seq(nthetamax)]
    theta <- lapply(seq_len(nalam), function(i){
      tmp <- theta_vec[(pmax * H * (i-1) + 1):(pmax * H * i)]
      a <- matrix(tmp, pmax, H, byrow = TRUE)[seq(nthetamax), , drop = FALSE]
      theta_i <- matrix(0, p, H)
      theta_i[ja,] <- a
      theta_i
    })
  }
  else{
    theta <- lapply(seq(nalam), function(x){matrix(0, p, H)})
  }
  list(theta = theta, lambda = lam)
}

err <- function(n, maxit, pmax) {
  if (n == 0) 
    msg <- ""
  if (n > 0) {
    # fatal error
    if (n < 7777) 
      msg <- "Memory allocation error; contact package maintainer"
    if (n == 10000) 
      msg <- "All penalty factors are <= 0"
    n <- 1
    msg <- paste("in the fortran code -", msg)
  }
  if (n < 0) {
    # non fatal error
    if (n > -10000) 
      msg <- paste("Convergence for ", -n, "th lambda value not reached after maxit=", maxit, " iterations; solutions for larger lambdas returned.\n", sep = "")
    if (n < -10000) 
      msg <- paste("Number of nonzero coefficients along the path exceeds pmax=", pmax, " at ", -n - 10000, "th lambda value; solutions for larger lambdas returned.\n", sep = "")
    if (n < -20000) 
      msg <- paste("Number of nonzero coefficients along the path exceeds dfmax=", pmax, " at ", -n - 20000, "th lambda value; solutions for larger lambdas returned.\n", sep = "")
    n <- -1
  }
  list(n = n, msg = msg)
}

lamfix <- function(lam){
  llam <- log(lam)
  if(length(llam) >= 3){lam[1] <- exp(2 * llam[2] - llam[3])}
  lam
}

seas <- function(x = NULL, y = NULL, yclass = NULL, d = NULL, categorical=FALSE, H=5, type = 'sir', M = NULL, U = NULL, nobs = NULL, lam1 = NULL, lam2 = NULL, gamma = NULL, lam1_fac=seq(1.0,0.01, length.out = 10), lam2_fac=seq(0.01,0.5, length.out = 10), FUN = NULL, eps = 1e-3, maxit = 1e+3, ...){
  
  if(is.null(M) || is.null(U)){ # Generate M and U matrices
    if(missing(x) || missing(y)) stop("Missing x or y.")
    if(is.data.frame(x)) x <- as.matrix(x)
    if(is.null(yclass)){
      if(categorical == FALSE){
        ybreaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H), na.rm=TRUE))
        yclass <- cut(y, breaks = ybreaks, include.lowest = TRUE, labels = FALSE)
        nclass <- as.integer(length(unique(yclass)))
      }
      else if(categorical == TRUE){
        yclass <- y
      }
    }
    if(any(table(yclass) < 5)) warning(sprintf("The sample size of class %d is less than 5\n", which(table(yclass) < 5)))
    if(is.null(gamma)){
      gamma <- c(10,30,50)
    }
    if(is.null(lam1) || is.null(lam2)){ # Automatically generate the tuning parameter sequence using the function cv.msda.
      fit_1 <- cv.msda(x, y, yclass = yclass, type=type, nlambda=10, lambda.factor=0.5, nfolds = 5, FUN = FUN, maxit=1e3)
      M <- fit_1$M  # The M matrix based on the full data.
      U <- fit_1$U  # The U matrix based on the full data.
      id_max_msda <- fit_1$id
      lam1_max_msda <- fit_1$lam_max  # The optimal lambda from msda
      beta_msda <- as.matrix(fit_1$beta)  # The optimal matrix from msda
      if(is.null(lam1)) lam1 <- (lam1_max_msda)*lam1_fac
      if(is.null(lam2)) lam2 <- svd(beta_msda)$d[1] * matrix(gamma, ncol = 1) %*% matrix(lam2_fac, nrow = 1)
      if (all(lam2 == 0)){
        lam2 <- 0
        warning("The automatically generated lambda 2 is zero, no nuclear norm penalty is imposed.")
      }
    }else{
      MU_out <- MU(x, y, yclass, type, FUN)
      M <- MU_out$M
      U <- MU_out$U
    }
    nobs <- as.integer(dim(x)[1])
    nvars <- as.integer(dim(x)[2])
  }
  else{
    if(is.null(lam1) || is.null(lam2) || is.null(gamma)) stop("Sequences lam1, lam2 or gamma is missing.")
    if(is.null(nobs)) stop("Missing nobs.")
    nvars <- NCOL(M)
  }
  
  ## Error code
  code <- 0
  
  if(is.vector(lam1) && (length(lam1) == 1) && (lam1 == 0) && is.vector(lam2) && (length(lam2) == 1) && (lam2 == 0)){ # For degenerate case where lambda1 = lambda2 = 0, return B = M^{-1} U directly.
    B <- solve(M) %*% U
    if(is.null(d)) beta <- svd(B)$u
    else if(d == 0) beta <- matrix(0, nrow(Bnew), ncol(Bnew))
    else beta <- svd(B)$u[,1:d,drop=FALSE]
    vec <- as.vector(beta)
    vec[abs(vec) < 1e-3] <- 0
    beta <- matrix(vec, nrow(beta), ncol(beta))
    rank <- NCOL(beta)
    output <- list(beta = beta, B = B, rank = rank, lam1 = lam1, lam2 = lam2, code = code)
  }
  else{
    # Fit with admm function
    fit <- admm(M, U, nobs, nvars, lam1, lam2, gamma, eps, maxit, d, ...)
    B_l <- fit$B
    beta_l <- fit$beta
    if (all(sapply(beta_l, is.null))){
      code <- 1
      warning("No converged results returned.")
      return(list(beta = beta_l, code = code))
    }
    rank_l <- fit$rank
    s_l <- fit$s
    step_l <- fit$step
    time_l <- fit$time
    if(length(B_l) == 1){
      B_l = B_l[[1]]; beta_l = beta_l[[1]]; rank_l = rank_l[[1]]; s_l = s_l[[1]]; step_l = step_l[[1]]; time_l = time_l[[1]]
    }
    output <- list(beta = beta_l, B = B_l, rank = rank_l, s = s_l, lam1 = lam1, lam2 = lam2, gamma = gamma, step = step_l, time = time_l, code = code)
  }
  output
}

cv.seas <- function(x, y, yclass = NULL, d = NULL, categorical=FALSE, H=5, type = 'sir', lambda.factor=0.5, nlambda=10, nfolds = 5, foldid = NULL, lam1 = NULL, lam2 = NULL, gamma = NULL, lam1_fac=seq(1.0,0.01, length.out = 10), lam2_fac=seq(0.01,0.5, length.out = 10), plot = FALSE, FUN = NULL, eps = 1e-3, maxit = 1e+3, trace.it = FALSE, ...){
  # The inputs and outputs are similar to the ones in 'seas' functions. Only the different ones are listed below.
  # Inputs:
  # =======
  # nfolds: The number of folds in the cross-validation.
  # plot: If TRUE, (1) plot the evaluation for each tuning parameter in 'msda' function; (2) in each cross-validation data fold, plot the evaluation for each tuning parameter in 'seas' function.
  # trace.it: If TRUE, print the process of cross-validation.
  # 
  # Outputs:
  # ========
  # Refer to the outputs in 'seas' function.
  
  start_time <- Sys.time()
  if(is.data.frame(x)) x <- as.matrix(x)
  if(is.null(yclass)){
    if(categorical == FALSE){
      ybreaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H), na.rm=TRUE))
      yclass <- cut(y, breaks = ybreaks, include.lowest = TRUE, labels = FALSE)
    }
    else if(categorical == TRUE){
      yclass <- y
    }
  }
  if(any(table(yclass) < 5)) warning(sprintf("The sample size of class %d is less than 5\n", which(table(yclass) < 5)))
  
  M <- U <- M_fold <- U_fold <- NULL
  nobs <- dim(x)[1]
  
  if(is.null(foldid)){
    ord <- order(y)
    y <- y[ord]
    yclass <- yclass[ord]
    x <- x[ord,]
    if (nfolds < 3) stop("nfolds must be larger than 3")
    if (nfolds > nobs) stop("nfolds is larger than the sample size")
    count <- as.numeric(table(yclass))
    foldid <- c()
    for(cnt in count){
      foldid <- c(foldid, sample(rep(seq(nfolds), length = cnt)))
    }
  }else{nfolds <- length(unique(foldid))}
  
  if(is.null(gamma)){
    gamma <- c(10,30,50)
  }
  if(is.null(lam1) || is.null(lam2)){ # Automatically generate the tuning parameter sequence using the function cv.msda.
    fit_1 <- cv.msda(x, y, yclass = yclass, type=type, nlambda=nlambda, lambda.factor=lambda.factor, foldid = foldid, FUN = FUN, maxit=1e3, plot = plot)
    M <- fit_1$M  # The M matrix based on the full data.
    U <- fit_1$U  # The U matrix based on the full data.
    M_fold <- fit_1$M_fold  # The M matrix based on each four out of five folds.
    U_fold <- fit_1$U_fold  # The U matrix based on each four out of five folds.
    id_max_msda <- fit_1$id
    lam1_max_msda <- fit_1$lam_max
    beta_msda <- as.matrix(fit_1$beta)
    if(is.null(lam1)) lam1 <- (lam1_max_msda)*lam1_fac
    if(is.null(lam2)) lam2 <- svd(beta_msda)$d[1] * matrix(gamma, ncol = 1) %*% matrix(lam2_fac, nrow = 1)
    if (all(lam2 == 0)){
      lam2 <- 0
      warning("The automatically generated lambda 2 is zero, no nuclear norm penalty is imposed.")
    }
  }
  n1 <- length(lam1)
  n2 <- ifelse(is.null(dim(lam2)), length(lam2), dim(lam2)[2])
  n3 <- length(gamma)
  
  # The number of errors
  nerr <- 0
  code <- 0
  
  end_time <- Sys.time()
  time1 <- difftime(end_time, start_time, units = "secs")
  ## Record time1: estimate M and U matrices
  
  out_all <- lapply(1:nfolds, function(k){ # Cross-validation
    if(trace.it) cat(sprintf("Fold: %d/%d\n", k, nfolds))
    x_val <- x[foldid==k,,drop=FALSE]
    y_val <- y[foldid==k]
    
    if(is.null(M_fold) || is.null(U_fold)){
      x_train <- x[foldid!=k,,drop=FALSE]
      y_train <- y[foldid!=k]
      yclass_train <- yclass[foldid!=k]
      # Fit with seas function
      fit_fold <- seas(x_train, y_train, yclass = yclass_train, type = type, FUN = FUN, lam1 = lam1, lam2 = lam2, gamma = gamma, eps = eps, maxit = maxit, d = d)
    }
    else{
      fit_fold <- seas(M = M_fold[[k]], U = U_fold[[k]], nobs = sum(foldid!=k), lam1 = lam1, lam2 = lam2, gamma = gamma, eps = eps, maxit = maxit, d = d)
    }
    
    err <- 0
    beta_l <- fit_fold$beta
    rank_l <- fit_fold$rank
    step_l <- fit_fold$step
    time_l <- fit_fold$time
    
    eval_fold <- eval_dc(beta_l, x_val, y_val)  # The evaluation: distance correlation.
    ind <- which(sapply(beta_l, is.null))
    rank_l[ind] <- -1
    eval_fold[ind] <- min(eval_fold, na.rm = TRUE)
    
    if(plot){ # If true, plot the evaluation for each tuning parameter
      dat <- data.frame(x = 1:length(eval_fold), y = eval_fold, rank = as.factor(rank_l))
      g <- ggplot(dat, aes(x = x, y = y, col = rank))+
        geom_point(size = 1)+
        labs(title=paste0("Fold ", k), x="", y="Distance correlation")+
        theme_bw()
      #print(g)
    }
    out <- list(eval_fold, err)
    out
  })
  
  # Combine the evaluations from each fold.
  eval_all <- do.call(rbind, lapply(out_all, "[[", 1))
  errs <- do.call(c, lapply(out_all, "[[", 2))
  nerr <- sum(errs)
  
  if((nerr != 0) && (nerr != nfolds)){
    code <- 3
    warning(paste0("No converged results returned in", nerr, "folds."))
  }else if(nerr == nfolds){
    code <- 4
    warning("No converged results returned in any fold.")
    return(list(beta = NULL, code = code))
  }
  
  if(is.vector(eval_all)){
    eval_all <- as.matrix(eval_all)
  }
  
  # Compute the cross-validation mean and standard error.
  cvm <- colMeans(eval_all, na.rm=TRUE)
  cvsd <- sqrt(colMeans(scale(eval_all, cvm, FALSE)^2, na.rm = TRUE)/(nfolds-1))
  
  # Select the best tuning parameter.
  id_max <- which.max(cvm)
  id_lam1 <- ceiling(id_max/(n2*n3))
  id_gamma <- ceiling((id_max-(id_lam1-1)*(n2*n3))/n2)
  id_lam2 <- id_max-(id_lam1-1)*(n2*n3)-(id_gamma-1)*n2
  lam1_max <- lam1[id_lam1]
  gamma_max <- gamma[id_gamma]
  lam2_max <- ifelse(is.null(dim(lam2)), lam2[id_lam2], lam2[id_gamma,id_lam2])
  
  start_time <- Sys.time()
  # Refit with the selected tuning parameters.
  if(is.null(M) || is.null(U)){
    fit <- seas(x, y, yclass = yclass, type = type, FUN = FUN, lam1 = lam1_max, lam2 = lam2_max, gamma = gamma_max, eps = eps, maxit = maxit, d = d, ...)
  }else{
    fit <- seas(M = M, U = U, nobs = NROW(x), lam1 = lam1_max, lam2 = lam2_max, gamma = gamma_max, eps = eps, maxit = maxit, d = d, ...)
  }
  
  if(fit$code != 0){
    code <- 5
    warning("The estimated beta is null.")
    return(list(beta = NULL, code = code))
  }
  
  B <- fit$B
  beta <- fit$beta
  rank <- fit$rank
  
  end_time <- Sys.time()
  time2 <- difftime(end_time, start_time, units = "secs") # We do not include the time for tuning parameter selection.
  # Record time: one run with the selected tuning parameter
  
  time <- time1 + time2
  
  output <- list(beta = beta, B = B, rank = rank, eval = eval_all, id_lam1=id_lam1, id_lam2 = id_lam2, id_gamma = id_gamma, lam1 = lam1, lam2 = lam2, gamma = gamma, lam1_max = lam1_max, lam2_max = lam2_max, gamma_max = gamma_max, code = code, time = time)
  output
}

# admm algorithm function
admm <- function(M, U, nobs, nvars, lam1, lam2, gam, eps=1e-3, maxit=1e+3, d = NULL, ...){
  # Inputs:
  # =======
  # M: The M matrix in optimization problem. It is the sample covariance of predictor for SEAS-SIR, SEAS-Intra, and SEAS-PFC.
  # U: The U matrix in optimization problem.
  # nobs: The number of observations.
  # nvars: The number of predictors.
  # lam1, lam2, gamma: The user-specified sequences of tuning parameter lambda_1, lambda_2 and gamma.
  # eps: The tolerance of convergence in ADMM algorithm. The value is passed to 'admm' function.
  # maxit: The maximal iterations in ADMM algorithm. The value is passed to 'admm' function.
  # d: The true structural dimension. The default is NULL.
  # 
  # Outputs:
  # ========
  # beta: A list containing the estimated basis matrices of central subspace.
  # B: A list containing estimated B matrices.
  # rank: A vector containing the estimated ranks.
  # s: A vector containing the estimated sparsity levels.
  # step: A vector containing the number of iterations to converge for each tuning parameter.
  # nlam: The number of converged matrices.
  
  # since the user is required to provide lam1, then set flmin=1
  if(is.null(dim(U)))
  {
    U<-cbind(U,rep(0,length(U)))
  }
  opts <- list(...)
  if(is.null(opts$nlam)) opts$nlam <- as.integer(1)
  if(is.null(opts$H)) opts$H <- as.integer(dim(U)[2])
  if(is.null(opts$nvars)) opts$nvars <- as.integer(nvars)
  if(is.null(opts$pf)) opts$pf <- as.double(rep(1, nvars))
  if(is.null(opts$dfmax)) opts$dfmax <- as.integer(nobs)
  if(is.null(opts$pmax)) opts$pmax <- as.integer(min(nobs*2+20, nvars))
  if(is.null(opts$flmin)) opts$flmin <- as.double(1)
  if(is.null(opts$eps_inner)) opts$eps_inner <- as.double(1e-04)
  if(is.null(opts$maxit_inner)) opts$maxit_inner <- as.integer(1e+6)
  if(is.null(opts$sml)) opts$sml <- as.double(1e-6)
  if(is.null(opts$verbose)) opts$verbose <- as.integer(FALSE)
  if(is.null(opts$nalam)) opts$nalam <- integer(1)
  if(is.null(opts$theta)) opts$theta <- double(opts$pmax * opts$H * opts$nlam)
  if(is.null(opts$itheta)) opts$itheta <- integer(opts$pmax)
  if(is.null(opts$ntheta)) opts$ntheta <- integer(opts$nlam)
  if(is.null(opts$alam)) opts$alam <- double(opts$nlam)
  if(is.null(opts$npass)) opts$npass <- integer(1)
  if(is.null(opts$jerr)) opts$jerr <- integer(1)
  
  M0 <- M
  U0 <- U
  n1 <- length(lam1)
  n2 <- ifelse(is.null(dim(lam2)), length(lam2), ncol(lam2))
  n3 <- length(gam)
  nparams <- n1*n2*n3
  
  # The following lists save the corresponding objects for each tuning parameter, 
  B_l <- vector("list", nparams) # estimated B matrix
  beta_l <- vector("list", nparams) # estimated basis matrix
  step_l <- rep(NA_integer_, nparams) # iterations
  time_l <- rep(NA_real_, nparams) # running time
  rank_l <- rep(NA_integer_, nparams) # estimated rank
  s_l <- rep(NA_integer_, nparams) # estimated sparsity level
  
  # Count the number of converged matrices.
  nlam_cvg <- 0
  
  for(i in 1:n1){
    lambda1 <- as.double(lam1[i])
    
    for(j in 1:n3){
      gamma <- gam[j]
      
      for(k in 1:n2){
        lambda2 <- ifelse(is.null(dim(lam2)), lam2[k], lam2[j,k])
        
        M <- M0 + gamma*diag(rep(1,ncol(M0)), ncol(M0),ncol(M0))
        
        # Initialize three matrices
        Bold <- matrix(0,dim(U0)[1], dim(U0)[2])
        Cold <- matrix(0,dim(U0)[1], dim(U0)[2])
        etaold <- matrix(0,dim(U0)[1], dim(U0)[2])
        
        # The MAIN loop of admm method
        step <- 0    
        start_time <- Sys.time()
        
        repeat{
          step <- step + 1
          
          # Update B
          U <- U0 - etaold + gamma * Cold
          out_B <- updateB(M, U, lambda1, opts)
          Bnew <- out_B$Bnew
          jerr <- out_B$jerr
          if(jerr != 0) break
          
          # Update C
          Cnew <- updateC(Bnew, lambda2, gamma, etaold)
          
          # Update eta (omega in SEAS algorithm)
          etanew <- etaold + gamma * (Bnew - Cnew)
          
          # Code 1: success
          if(max(abs(Bnew - Cnew)) < eps){
            jerr <- 1
            break
          }
          # Code 404: then maximal iteration is reached
          if(step > maxit){
            jerr <- 404
            warning('Maximal iteration is reached.')
            break
          }
          Bold <- Bnew
          Cold <- Cnew
          etaold <- etanew
        }# End of repeat 
        end_time <- Sys.time()  # The time for each repeat
        time <- difftime(end_time, start_time, units = "secs")
        # Code < -10000: non-sparse matrix
        if(jerr < -10000){
          break
        }
        # Code 1: success, save the matrix and the related information.
        if(jerr==1){
          index <- (i-1)*n2*n3 + (j-1)*n2 + k
          nlam_cvg <- nlam_cvg + 1
          B_l[[index]] <- Bnew
          step_l[index] <- step
          time_l[index] <- time
          if(is.null(d)) rank <- rank_func(Cnew, thrd = eps)
          else rank <- d
          rank_l[index] <- rank
          # Cut and select the left singular vector of Bnew
          if(rank == 0){
            beta <- matrix(0, nrow(Bnew), ncol(Bnew))
          }else{
            tmp <- svd(Bnew)$u[,1:rank, drop = FALSE]
            vec <- as.vector(tmp)
            vec[abs(vec) < eps] <- 0
            beta <- matrix(vec, nrow(tmp), ncol(tmp))
          }
          beta_l[[index]] <- beta
          var_ind <- apply(beta, 1, function(x){any(x!=0)})
          s_l[index] <- sum(var_ind)
        }
      }# End of lambda2
      if(jerr < -10000) break
    }# End of gam
  }# End of lambda1
  return(list(beta = beta_l, B = B_l, rank = rank_l, s = s_l, step = step_l, time = time_l, nlam = nlam_cvg))
}

# Update B matrix in ADMM algorithm. Use the group-wise coordinate descent algorithm from 'msda' R package.
updateB <- function(M, U, lambda1, opts){
  U <- t(U)
  fit <- .Fortran("msda", obj = opts$nlam, opts$H, opts$nvars, as.double(M), as.double(U), opts$pf, opts$dfmax, opts$pmax, opts$nlam, opts$flmin, lambda1, opts$eps_inner, opts$maxit_inner, opts$sml, opts$verbose, nalam = opts$nalam, theta = opts$theta, itheta = opts$itheta, ntheta = opts$ntheta, alam = opts$alam, npass = opts$npass, jerr = opts$jerr)
  if(fit$jerr != 0){return(list(Bnew = NULL, jerr = fit$jerr))} # Code: non-zero, abnormal result.
  outlist <- formatoutput(fit, opts$maxit_inner, opts$pmax, opts$nvars, opts$H)
  Bnew <- as.matrix(outlist$theta[[1]])
  list(Bnew = Bnew, jerr = fit$jerr)
}

# Update C matrix in ADMM algorithm.
updateC <- function(Bnew, lambda2, gamma, etaold){
  Btemp <- Bnew + 1/gamma * etaold
  svd_B <- svd(Btemp)
  lamtemp <- pmax(0, svd_B$d-lambda2/gamma)
  Cnew <- svd_B$u %*% diag(lamtemp, nrow = length(lamtemp), ncol = length(lamtemp)) %*% t(svd_B$v)
  Cnew
}

# ------------------ revised functions from 'msda' package ---------------------- #
# Revise 'msda' function to accommodate other forms of M and U matrices.
# Some inputs are similar to the ones in 'seas' function. Please refer to 'msda' package documentation for more details of the arguments in 'msda' function.
# 
# Outputs:
# ========
# lambda: The tuning parameter sequence.
# theta: The list of estimated matrix.
# M: The M matrix from samples. It is the sample covariance of predictor for SEAS-SIR, SEAS-Intra, and SEAS-PFC.
# U: The U matrix from samples, which depends on the argument type.
# rank: The list of estimated rank for each matrix.
msda <- function(x, y, yclass=NULL, categorical=FALSE, H=5, type='sir', FUN = NULL, lambda.factor=NULL, nlambda=100, lambda=NULL, dfmax=NULL, pmax=NULL, pf=NULL, M = NULL, U = NULL, nobs=NULL, nclass=NULL, eps=1e-04, maxit=1e+06, sml=1e-06, verbose=FALSE, perturb=NULL){
  if(is.null(M) || is.null(U)){ # Generate M and U matrices
    if(missing(x) || missing(y)) stop("Missing x or y.")
    if(is.data.frame(x)) x <- as.matrix(x)
    if(is.null(yclass)){
      if(categorical == FALSE){
        ybreaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H), na.rm=TRUE))
        yclass <- cut(y, breaks = ybreaks, include.lowest = TRUE, labels = FALSE)
      }
      else if(categorical == TRUE){
        yclass <- y
      }
    }
    if(any(table(yclass) < 5)) warning(sprintf("The sample size of class %d is less than 5\n", which(table(yclass) < 5)))
    nclass <- as.integer(length(unique(yclass)))
    MU_out <- MU(x, y, yclass, type, FUN)
    M <- MU_out$M
    U <- MU_out$U
    nobs <- as.integer(dim(x)[1])
    nvars <- as.integer(dim(x)[2])
  }
  else{
    if(is.null(nobs)) stop("Missing nobs.")
    if(is.null(nclass)) stop("Missing nclass.")
    nvars <- NCOL(M)
  }
  
  if(is.null(lambda.factor)) lambda.factor <- ifelse((nobs - nclass)<=nvars, 0.2, 1e-03)
  if(is.null(dfmax)) dfmax <- nobs
  if(is.null(pmax)) pmax <- min(dfmax*2 + 20, nvars)
  if(is.null(pf)) pf <- rep(1, nvars)
  if (!is.null(perturb)) 
    diag(M) <- diag(M) + perturb
  H <- as.integer(dim(U)[2])
  ## parameter setup
  if (length(pf) != nvars) 
    stop("The size of penalty factor must be same as the number of input variables")
  maxit <- as.integer(maxit)
  verbose <- as.integer(verbose)
  sml <- as.double(sml)
  pf <- as.double(pf)
  eps <- as.double(eps)
  dfmax <- as.integer(dfmax)
  pmax <- as.integer(pmax)
  ## lambda setup
  nlam <- as.integer(nlambda)
  if (is.null(lambda)) {
    if (lambda.factor >= 1)
      stop("lambda.factor should be less than 1")
    flmin <- as.double(lambda.factor)
    ulam <- double(1)  #ulam=0 if lambda is missing
  } else {
    # flmin=1 if user define lambda
    flmin <- as.double(1)
    if (any(lambda < 0))
      stop("lambdas should be non-negative")
    ulam <- as.double(rev(sort(lambda)))
    nlam <- as.integer(length(lambda))
  }
  ## call Fortran core
  fit <- .Fortran("msda", obj = double(nlam), H, nvars, as.double(M), as.double(t(U)), pf, dfmax, pmax, nlam, flmin, ulam, eps, maxit, sml, verbose, nalam = integer(1), theta = double(pmax * H * nlam), itheta = integer(pmax), ntheta = integer(nlam),alam = double(nlam), npass = integer(1), jerr = integer(1))
  ## output
  outlist <- formatoutput(fit, maxit, pmax, nvars, H)
  rank <- rep(NA_integer_, length(outlist$theta))
  for (i in 1:length(outlist$theta)){
    if(!is.null(outlist$theta[[i]])){
      rank[i] <- rank_func(outlist$theta[[i]], thrd = 1e-3)
    }
  }
  if(is.null(lambda))
    outlist$lambda <- lamfix(outlist$lambda)
  outlist <- list(lambda = outlist$lambda, theta = outlist$theta, M = M, U = U, rank = rank)
  class(outlist) <- c("msda")
  outlist
}

# Revise 'cv.msda' function to accommodate other forms of M and U matrices. We also add the optional argument 'fold' to pass the user-specified folds index.
# Some inputs are similar to the ones in 'cv.seas' function. Please refer to 'msda' package documentation for more details of arguments used in 'cv.msda' function.
# 
# Outputs:
# ========
# beta: The optimal estimated matrix.
# id: The index of the optimal tuning parameter.
# lambda: The lambda sequence.
# lam_max: The optimal tuning parameter.
# rank: The rank of the optimal estimated matrix.
# M: The M matrix based on the full data. It is the sample covariance of predictor for SEAS-SIR, SEAS-Intra, and SEAS-PFC.
# U: The U matrix based on the full data, which depends on the argument type.
# M_fold: The M matrix list based on each cross-validation data fold.
# U_fold: The U matrix list based on each cross-validation data fold.
cv.msda <- function(x, y, yclass=NULL, categorical=FALSE, H=5, type='sir', lambda.factor=NULL, nlambda=100, nfolds=5, foldid = NULL, lambda = NULL, FUN = NULL, maxit = 1e3, plot = FALSE){
  if(is.data.frame(x)) x <- as.matrix(x)
  if(is.null(yclass)){
    if(categorical == FALSE){
      ybreaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H), na.rm=TRUE))
      yclass <- cut(y, breaks = ybreaks, include.lowest = TRUE, labels = FALSE)
    }
    else if(categorical == TRUE){
      yclass <- y
    }
  }
  if(any(table(yclass) < 5)) warning(sprintf("The sample size of class %d is less than 5\n", which(table(yclass) < 5)))
  nobs <- nrow(x)
  nclass <- length(unique(yclass))
  if(is.null(lambda.factor)) lambda.factor <- ifelse((nobs - nclass)<=nvars, 0.2, 1e-03)
  # Fit the model on the full data, obtain the lambda sequence.
  fit <- msda(x, y, yclass = yclass, type = type, lambda.factor = lambda.factor, nlambda = nlambda, lambda = lambda, FUN = FUN, maxit=maxit)
  lambda <- fit$lambda
  beta_l <- fit$theta
  M <- fit$M
  U <- fit$U
  rank_l <- fit$rank
  beta_l <- cut_mat(beta_l, 1e-3, rank_l)
  
  # Cross-validation
  if(is.null(foldid)){
    ord <- order(y)
    y <- y[ord]
    yclass <- yclass[ord]
    x <- x[ord,]
    count <- as.numeric(table(yclass))
    foldid <- c()
    for(cnt in count){
      foldid <- c(foldid, sample(rep(seq(nfolds), length = cnt)))
    }
  }
  else{
    nfolds <- length(unique(foldid))
  }
  
  cv_out <- lapply(1:nfolds, function(k){
    x_train <- x[foldid!=k,,drop=FALSE]
    x_val <- x[foldid==k,,drop=FALSE]
    y_train <- y[foldid!=k]
    y_val <- y[foldid==k]
    yclass_train <- yclass[foldid!=k]
    
    fit_fold <- msda(x_train, y_train, yclass_train, type = type, lambda.factor=lambda.factor, nlambda=nlambda, lambda = lambda, FUN = FUN, maxit=maxit)
    M_fold <- fit_fold$M
    U_fold <- fit_fold$U
    beta_fold <- fit_fold$theta
    rank_fold <- fit_fold$rank
    beta_fold <- cut_mat(beta_fold, 1e-3, rank_fold)
    
    # return evaluation of each fold
    eval_fold <- eval_dc(beta_fold, x_val, y_val)
    if(length(eval_fold) != length(lambda)){
      eval_fold <- c(eval_fold, rep(NA, length(lambda) - length(eval_fold)))
    }
    list(eval = eval_fold, M = M_fold, U = U_fold)
  })
  
  eval_all <- do.call(rbind, lapply(cv_out, "[[", 1))
  M_fold <- lapply(cv_out, "[[", 2)
  U_fold <- lapply(cv_out, "[[", 3)
  if(is.vector(eval_all)){
    eval_all <- t(as.matrix(eval_all))
  }
  
  ## No matrix is converged in any fold
  if(all(is.na(eval_all))) return(NULL)
  
  print(eval_all)
  print(typeof(eval_all))
  cvm <- colMeans(eval_all, na.rm = TRUE)
  # The optimal lambda1
  id_max <- which.max(cvm)
  lam_max <- lambda[id_max]
  beta <- as.matrix(beta_l[[id_max]])
  
  # Recalculate the rank
  rank <- rank_func(beta, thrd = 1e-3)
  
  if(plot){ # If TRUE, plot the cv evaluation for each tuning parameter.
    dat <- data.frame(x = 1:length(cvm), y = cvm)
    g <- ggplot(dat, aes(x = x, y = y))+
      geom_point(size = 1)+
      xlab("")+
      ylab("Distance correlation")+
      theme_bw()
    #print(g)
  }
  
  list(beta = beta, id = id_max, lambda = lambda, lam_max = lam_max, rank = rank, M = M, U = U, M_fold = M_fold, U_fold = U_fold)
}

matpower <- function(a,alpha){
  small <- .00000001
  if (length(c(a))==1) {return(a^(alpha))} else {
    p1<-nrow(a)
    eva<-eigen(a)$values
    eve<-eigen(a)$vectors
    eve<-eve/t(matrix((diag(t(eve)%*%eve)^0.5),p1,p1))
    index<-(1:p1)[abs(eva)>small]
    evai<-eva
    evai[index]<-(eva[index])^(alpha)
    ai<-eve%*%diag(evai,length(evai))%*%t(eve)
    return(ai)}
}

proj<-function(v)
{return(v%*%matpower(t(v)%*%v,-1)%*%t(v))}

##standarize x
stand <- function(x){
  n<-nrow(x)
  p<-ncol(x)
  xb <- apply(x, 2, mean)
  xb <- t(matrix(xb, p, n))
  x1 <- x - xb
  sigma <- t(x1) %*% (x1)/(n-1)
  eva <- eigen(sigma)$values
  eve <- eigen(sigma)$vectors
  sigmamrt <- eve %*% diag(1/sqrt(eva)) %*% t(eve)
  z <- sigmamrt %*% t(x1)
  return(t(z))
}

##deal with list
deal_list<-function(x)
{
  return(as.vector(x)[-1])
}

allzero<-function(X)  
{
  return(all(X==0))
}

sumzero<-function(X)
{
  return(length(which(X!=0)))
}

thresholding<-function(X,threshold)
{
  return(all(sum(X^2)<=threshold))
}

##slicing Y
slicing<-function(y,H) 
{
  n<-length(y)
  ord=order(y)
  y=y[ord]
  c=as.integer(n/H)
  if (length(levels(as.factor(y)))>H)
  {
    ytilde<-rep(0,H+1)
    ytilde[1]<-min(y)
    for (h in 1:(H-1))
    {
      ytilde[h+1]<-y[h*c+1]
    }  
  }
  if (length(levels(as.factor(y)))<=H)
  {
    H <- length(levels(as.factor(y)))
    ytilde<-rep(0,H+1)
    ytilde[1]=min(y)
    for (h in 1:(H-1))
    {
      ytilde[h+1]<-min(y[y>ytilde[h]])
    }
  } 
  ytilde[H+1]=max(y)+1
  prop<-rep(1,H)
  for (i in 1:H)
  {
    prop[i] = sum((y >= ytilde[i])&(y < ytilde[i+1]))/n
  }
  res<-list()
  res$H<-H
  res$ytilde<-ytilde
  res$prop<-prop
  return(res)
}

mhat_dr<-function(x,y,H){
  n<-nrow(x)
  p<-ncol(x)
  z<-stand(x)
  dy <- slicing(y,H)
  H <- dy$H
  ytilde <- dy$ytilde
  prop <- dy$prop
  ind<-matrix(0,n,H)
  zbar<-matrix(0,p,H)
  for (j in 1:(H-1))
  {
    ind[,j]<-((y >= ytilde[j])&(y < ytilde[j+1]))
    zbar[,j]<- (t(z)%*%(ind[,j]))/sum(ind[,j])
  }
  ind[,H]<-(y >= ytilde[H])
  zbar[,H]<- (t(z)%*%(ind[,H]))/sum(ind[,H])
  A<-matrix(0,p,p)
  B<-matrix(0,p,p)
  C<-0
  for (q in 1:H)
  {
    Z<-(t(z))[,ind[,q]==1]-zbar[,q]  
    A<-A + prop[q]*((Z%*%t(Z)/(sum(ind[,q])-1)+zbar[,q]%*%t(zbar[,q]))%*%  
                      (Z%*%t(Z)/(sum(ind[,q])-1)+zbar[,q]%*%t(zbar[,q])) - diag(1,p))
    B<-B + sqrt(prop[j])*(zbar[,q]%*%t(zbar[,q]))
    C<-C + sqrt(prop[j])*(t(zbar[,q])%*%zbar[,q])
  }
  C<-as.vector(C)
  M<-2*A + 2*(B%*%B) + 2*B*C
  return(M)
}

dr_kernel<-function(x,y,H)
{
  n <- dim(x)[1]
  p <-dim(x)[2]
  dy <- slicing(y,H)
  H <- dy$H
  ytilde <- dy$ytilde
  prop <- dy$prop
  xxbar<-numeric(0)
  freq<-rep(0,H)
  if(H>=3)
  {
    for (j in 2:(H-1))
    {
      for(i in 1:(j-1))
      {
        indj<-((y >= ytilde[j])&(y < ytilde[j+1]))
        indi<-((y >= ytilde[i])&(y < ytilde[i+1]))
        freq[j]=length(which(indj==TRUE))
        freq[i]=length(which(indi==TRUE))
        Sj<-matrix(0,ncol=n,nrow=freq[j])
        Sj[,indj]<-diag(freq[j])
        Si<-matrix(0,ncol=n,nrow=freq[i])
        Si[,indi]<-diag(freq[i])
        xxbar<-cbind(xxbar,sqrt(prop[j]*prop[i])*(2*x-1/prop[i]*t(Si)%*%Si%*%x-1/prop[j]*t(Sj)%*%Sj%*%x
                                                  +n*t(Si)%*%rep(1/freq[i],freq[i])%*%t(rep(1/freq[j],freq[j]))%*%Sj%*%x
                                                  +n*t(Sj)%*%rep(1/freq[j],freq[j])%*%t(rep(1/freq[i],freq[i]))%*%Si%*%x))
      }
    }
  }
  for(i in 1:(H-1))
  {
    indj<-(y >= ytilde[H])
    indi<-((y >= ytilde[i])&(y < ytilde[i+1]))
    freq[H]=length(which(indj==TRUE))
    freq[i]=length(which(indi==TRUE))
    Sj<-matrix(0,ncol=n,nrow=freq[H])
    Sj[,indj]<-diag(freq[H])
    Si<-matrix(0,ncol=n,nrow=freq[i])
    Si[,indi]<-diag(freq[i])
    xxbar<-cbind(xxbar,sqrt(prop[H]*prop[i])*(2*x-1/prop[i]*t(Si)%*%Si%*%x-1/prop[H]*t(Sj)%*%Sj%*%x
                                              +t(Si)%*%rep(1/freq[i],freq[i])%*%t(rep(1/freq[H],freq[H]))%*%Sj%*%x
                                              +t(Sj)%*%rep(1/freq[H],freq[H])%*%t(rep(1/freq[i],freq[i]))%*%Si%*%x))
  }
  return(xxbar)
}

save_kernel<-function(x,y,H)
{
  n <- dim(x)[1]
  p <-dim(x)[2]
  dy <- slicing(y,H)
  H <- dy$H
  ytilde <- dy$ytilde
  prop <- dy$prop
  freq<-prop*n
  xxbar<-numeric(0)
  for (j in 1:(H-1))
  {
    ind<-((y >= ytilde[j])&(y < ytilde[j+1]))
    freq[j]=length(which(ind==TRUE))
    S<-matrix(0,ncol=n,nrow=freq[j])
    S[,ind]<-diag(freq[j])-matrix(1/freq[j],freq[j],freq[j])
    xxbar<-cbind(xxbar,sqrt(prop[j])*(x-((n-1)/(freq[j]-1))*t(S)%*%S%*%x))
  } 
  ind<-(y >= ytilde[H])
  freq[H]=length(which(ind==TRUE))
  S<-matrix(0,ncol=n,nrow=freq[H])
  S[,ind]<-diag(freq[H])-matrix(1/freq[H],freq[H],freq[H])
  xxbar<-cbind(xxbar,sqrt(prop[H])*(x-((n-1)/(freq[H]-1))*t(S)%*%S%*%x))
  return(xxbar)
}

sir_kernel<-function(x,y,H)
{
  n <- dim(x)[1]
  p <-dim(x)[2]
  dy <- slicing(y,H)
  H <- dy$H
  ytilde <- dy$ytilde
  prop <- dy$prop
  freq<-prop*n
  xxbar<-numeric(0)
  for (j in 1:(H-1))
  {
    ind<-((y >= ytilde[j])&(y < ytilde[j+1]))
    freq[j]<-length(which(ind==TRUE))
    S<-matrix(0,ncol=n,nrow=freq[j])
    S[,ind]<-diag(freq[j])
    xxbar<-cbind(xxbar,sqrt(prop[j])*n*(rep(1/n,n)-t(S)%*%rep(1/freq[j],freq[j])))
  } 
  ind<-(y >= ytilde[H])
  freq[H]<-length(which(ind==TRUE))
  S<-matrix(0,ncol=n,nrow=freq[H])
  S[,ind]<-diag(freq[H])
  xxbar<-cbind(xxbar,sqrt(prop[H])*n*(rep(1/n,n)-t(S)%*%rep(1/freq[H],freq[H])))
  return(xxbar)
}

Sparse_vector<-function(Y,design_matrix,nfolds=5,lambda=NULL)
{
  if(is.null(lambda))
  {
    fit=cv.glmnet(x=design_matrix,y=Y,nfolds=nfolds,standardize=FALSE, type.measure="mae",nlambda=100)
  }else
  {
    fit=cv.glmnet(x=design_matrix,y=Y,nfolds=nfolds,lambda=lambda,standardize=FALSE,type.measure="mae",nlambda=100)
  }
  return(as.numeric(coef(fit))[-1])
}

choose_weight<-function(x,rho=0.25)                
{
  x<-as.numeric(x)
  return(sum(x^2)^(-rho))
}

forward_column_selection <- function(X, threshold=1e-12, D=15) {
  n <- ncol(X)
  selected_cols <- c()
  col_norms <- apply(X, 2, function(col) sqrt(sum(col^2)))
  
  for (i in 1:D) {
    max_norm_index <- which.max(col_norms)
    max_norm_value <- col_norms[max_norm_index]
    if (max_norm_value < threshold) {
      break
    }
    selected_cols <- c(selected_cols, max_norm_index)
    selected_col <- X[, max_norm_index]
    for (j in 1:n) {
      if (!j %in% selected_cols) {
        projection <- sum(X[, j] * selected_col) / sum(selected_col^2) * selected_col
        X[, j] <- X[, j] - projection
        col_norms[j] <- sqrt(sum(X[, j]^2))
      }
    }
    col_norms[max_norm_index] <- 0
  }
  
  return(selected_cols)
}

save_lambda<-function(x,y,H)
{
  n<-nrow(x)
  p<-ncol(x)
  sig<-cov(x)
  xc<- t(t(x)-apply(x,2,mean))
  dy <- slicing(y,H)
  H <- dy$H
  ytilde <- dy$ytilde
  prop <- dy$prop
  xxbar<-numeric(0)
  for (j in 1:(H-1))
  {
    ind<-((y >= ytilde[j])&(y < ytilde[j+1]))
    xxbar<-cbind(xxbar,(sig-cov(xc[ind,])))
  } 
  ind<-(y >= ytilde[H])
  xxbar<-cbind(xxbar,(sig-cov(xc[ind,])))
  return(xxbar)
}

sir_lambda<-function(x,y,H)
{
  n<-nrow(x)
  p<-ncol(x)
  sig<-cov(x)
  xc<- t(t(x)-apply(x,2,mean))
  dy <- slicing(y,H)
  H <- dy$H
  ytilde <- dy$ytilde
  prop <- dy$prop
  xxbar<-numeric(0)
  for (j in 1:(H-1))
  {
    ind<-((y >= ytilde[j])&(y < ytilde[j+1]))
    xxbar<-cbind(xxbar,(apply(xc[ind,],2,mean)))
  } 
  ind<-(y >= ytilde[H])
  xxbar<-cbind(xxbar,(apply(xc[ind,],2,mean)))
  return(xxbar)
}


ADMM_estimation<-function(x, y, H1, H2, index, method, nfolds=5)
{
  p=ncol(x)
  n=nrow(x)
  ord <- order(y)
  y <- y[ord]
  x <- x[ord,]
  lam2=0
  lam1=0.96^(-50:50)
  gam=c(10,30,50)
  if(length(table(y))<nrow(x)/2)
  {yclass=y
  yclass_sir=y
  }else
  {
    ybreaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H1), na.rm=TRUE))
    yclass <- cut(y, breaks = ybreaks, include.lowest = TRUE, labels = FALSE)
    ybreaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H2), na.rm=TRUE))
    yclass_sir <- cut(y, breaks = ybreaks, include.lowest = TRUE, labels = FALSE)
  }
  #yclass_sir<-yclass_sir[ord]
  #yclass <- yclass[ord]
  count <- as.numeric(table(yclass))
  foldid <- c()
  for(cnt in count){
    foldid <- c(foldid, sample(rep(seq(nfolds), length = cnt)))}
  num=p*H1
  index_trans<-as.vector(t(matrix(c(1:num),p,H1)))
  out_all <- lapply(1:nfolds, function(k)
  { # Cross-validation
    x_val <- x[foldid==k,,drop=FALSE]
    y_val <- y[foldid==k]
    x_train <- x[foldid!=k,,drop=FALSE]
    y_train <- y[foldid!=k]
    yclass_train <- yclass[foldid!=k]
    yclass_sir_train<-yclass_sir[foldid!=k]
    # Fit with seas function
    M=cov(x_train)
    if(method=="save")
    {
      U=(save_lambda(x_train,yclass_train,H1)[,index_trans])[,index]
    }
    if(method=="ens")
    {
      U=(cbind(save_lambda(x_train,yclass_train,H1)[,index_trans],sir_lambda(x_train,yclass_sir_train,H2)))[,index]
    }
    if(method=="sir")
    {
      U=sir_lambda(x_train,yclass_sir_train,H2)
    }
    fit_fold <- admm(M, U, p,p, lam1=lam1, lam2=lam2, gam=gam, eps=1e-3, maxit=1e+3)
    beta_l <- fit_fold$beta
    B_l<-fit_fold$B
    rank_l <- fit_fold$rank
    step_l <- fit_fold$step
    time_l <- fit_fold$time
    eval_fold <- eval_dc(beta_l, x_val, y_val)  # The evaluation: distance correlation.
    ind <- which(sapply(beta_l, is.null))
    rank_l[ind] <- -1
    eval_fold[ind] <- min(eval_fold, na.rm = TRUE)
    return(eval_fold)
  })
  eval_all <- do.call(rbind, lapply(out_all, deal_list))
  cvm <- colMeans(eval_all, na.rm=TRUE)
  gam_max=c(10,30,50)[((which.max(cvm)+2)%%3+1)]
  lam1_max=lam1[(which.max(cvm)-1)%/%3+1]
  lam2_max <- 0
  M=cov(x)
  if(method=="save")
  {
    U=(save_lambda(x,y,H1)[,index_trans])[,index]
  }
  if(method=="ens")
  {
    U=cbind(save_lambda(x,y,H1)[,index_trans],sir_lambda(x,y,H2))[,index]
  }
  if(method=="sir")
  {
    U=sir_lambda(x,y,H2)
  }
  fit<-admm(M, U, p,p, lam1_max, lam2_max, gam_max, eps=1e-3, maxit=1e+3)
  if(length(index)==1)
  {
    return((fit$B[[1]])[,1])
  }else
  {
    return(fit$B[[1]])
  }
}

ADMM_estimation_rough<-function(x, y, H1, H2, index, method, nfolds=5)
{
  p=ncol(x)
  lam2=0
  lam1=seq(0.2,2,0.2)
  gam=c(10,30,50)
  if(length(table(y))<nrow(x)/2)
  {yclass=y
  yclass_sir=y
  }else
  {
    ybreaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H1), na.rm=TRUE))
    yclass <- cut(y, breaks = ybreaks, include.lowest = TRUE, labels = FALSE)
    ybreaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H2), na.rm=TRUE))
    yclass_sir <- cut(y, breaks = ybreaks, include.lowest = TRUE, labels = FALSE)
  }
  ord <- order(y)
  y <- y[ord]
  yclas_sir<-yclass_sir[ord]
  yclass <- yclass[ord]
  x <- x[ord,]
  count <- as.numeric(table(yclass))
  foldid <- c()
  for(cnt in count){
    foldid <- c(foldid, sample(rep(seq(nfolds), length = cnt)))}
  num=p*H1
  index_trans<-as.vector(t(matrix(c(1:num),p,H1)))
  out_all <- lapply(1:nfolds, function(k)
  { # Cross-validation
    x_val <- x[foldid==k,,drop=FALSE]
    y_val <- y[foldid==k]
    x_train <- x[foldid!=k,,drop=FALSE]
    y_train <- y[foldid!=k]
    yclass_train <- yclass[foldid!=k]
    yclass_sir_train<-yclass_sir[foldid!=k]
    # Fit with seas function
    M=cov(x_train)
    if(method=="save")
    {
      U=(save_lambda(x_train,yclass_train,H1)[,index_trans])[,index]
    }
    if(method=="ens")
    {
      U=(cbind(save_lambda(x_train,yclass_train,H1)[,index_trans],sir_lambda(x_train,yclass_sir_train,H2)))[,index]
    }
    time=Sys.time()
    if(method=="sir")
    {
      U=sir_lambda(x_train,yclass_sir_train,H2)
    }
    Sys.time()-time
    fit_fold <- admm(M, U, p,p, lam1=lam1, lam2=lam2, gam=gam, eps=1e-3, maxit=1e+3,d=1)
    beta_l <- fit_fold$beta
    rank_l <- fit_fold$rank
    step_l <- fit_fold$step
    time_l <- fit_fold$time
    eval_fold <- eval_dc(beta_l, x_val, y_val)  # The evaluation: distance correlation.
    ind <- which(sapply(beta_l, is.null))
    rank_l[ind] <- -1
    eval_fold[ind] <- min(eval_fold, na.rm = TRUE)
    return(eval_fold)
  })
  eval_all <- do.call(rbind, lapply(out_all, deal_list))
  cvm <- colMeans(eval_all, na.rm=TRUE)
  max_value<-max(cvm)
  index_cvm<-tail(which(cvm == max_value), n = 1)
  gam_max=c(10,30,50)[((index_cvm+2)%%3+1)]
  lam1_max=lam1[(index_cvm-1)%/%3+1]
  lam2_max <- 0
  M=cov(x)
  if(method=="save")
  {
    U=(save_lambda(x,y,H1)[,index_trans])[,index]
  }
  if(method=="ens")
  {
    U=cbind(save_lambda(x,y,H1)[,index_trans],sir_lambda(x,y,H2))[,index]
  }
  if(method=="sir")
  {
    U=sir_lambda(x,y,H2)
  }
  fit<-admm(M, U, p,p, lam1_max, lam2_max, gam_max, eps=1e-3, maxit=1e+3)
  if(length(index)==1)
  {
    return((fit$B[[1]])[,1])
  }else
  {
    return(fit$B[[1]])
  }
}

ADMM_estimation_weight<-function(x, y, H1, H2, weight, index, method, nfolds=5)
{
  p=ncol(x)
  n=nrow(x)
  ord <- order(y)
  y <- y[ord]
  x <- x[ord,]
  lam2=0
  lam1=0.96^(-30:70)*max(weight)
  gam=c(10,30,50)
  if(length(table(y))<nrow(x)/2)
  {yclass=y
  yclass_sir=y
  }else
  {
    ybreaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H1), na.rm=TRUE))
    yclass <- cut(y, breaks = ybreaks, include.lowest = TRUE, labels = FALSE)
    ybreaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H2), na.rm=TRUE))
    yclass_sir <- cut(y, breaks = ybreaks, include.lowest = TRUE, labels = FALSE)
  }
  count <- as.numeric(table(yclass))
  foldid <- c()
  for(cnt in count){
    foldid <- c(foldid, sample(rep(seq(nfolds), length = cnt)))}
  num=p*H1
  index_trans<-as.vector(t(matrix(c(1:num),p,H1)))
  out_all <- lapply(1:nfolds, function(k)
  { # Cross-validation
    x_val <- x[foldid==k,,drop=FALSE]
    y_val <- y[foldid==k]
    x_train <- x[foldid!=k,,drop=FALSE]
    y_train <- y[foldid!=k]
    yclass_train <- yclass[foldid!=k]
    yclass_sir_train<-yclass_sir[foldid!=k]
    # Fit with seas function
    M=cov(x_train)*outer(weight,weight)
    if(method=="save")
    {
      U=weight*(save_lambda(x_train,yclass_train,H1)[,index_trans])[,index]
    }
    if(method=="ens")
    {
      U=weight*(cbind(save_lambda(x_train,yclass_train,H1)[,index_trans],sir_lambda(x_train,yclass_sir_train,H2)))[,index]
    }
    if(method=="sir")
    {
      U=weight*sir_lambda(x_train,yclass_sir_train,H2)
    }
    fit_fold <- admm(M, U, p,p, lam1=lam1, lam2=lam2, gam=gam, eps=1e-3, maxit=1e+3)
    beta_l <- fit_fold$beta
    B_l<-fit_fold$B
    rank_l <- fit_fold$rank
    step_l <- fit_fold$step
    time_l <- fit_fold$time
    eval_fold <- eval_dc(beta_l, x_val, y_val)  # The evaluation: distance correlation
    return(eval_fold)
  })
  eval_all <- do.call(rbind, lapply(out_all, deal_list))
  cvm <- colMeans(eval_all, na.rm=TRUE)
  gam_max=c(10,30,50)[((which.max(cvm)+2)%%3+1)]
  lam1_max=lam1[(which.max(cvm)-1)%/%3+1]
  lam2_max <- 0
  M=cov(x)*outer(weight,weight)
  if(method=="save")
  {
    U=weight*(save_lambda(x,y,H1)[,index_trans])[,index]
  }
  if(method=="ens")
  {
    U=weight*cbind(save_lambda(x,y,H1)[,index_trans],sir_lambda(x,y,H2))[,index]
  }
  if(method=="sir")
  {
    U=weight*sir_lambda(x,y,H2)
  }
  fit<-admm(M, U, p,p, lam1_max, lam2_max, gam_max, eps=1e-3, maxit=1e+3)
  weight[weight==0]=Inf
  if(length(index)==1)
  {
    return((fit$B[[1]])[,1]/weight)
  }else
  {
    return(fit$B[[1]]/weight)
  }
}

Four_estimate_c<-function(X,Y,H1, H2,nfolds=5, D=10, lambda=NULL)
{
  p<-ncol(as.matrix(X))
  n<-nrow(as.matrix(X))
  kernel_matrix<-save_kernel(X, Y, H1)
  num=p*H1
  index_trans<-as.vector(t(matrix(c(1:num),p,H1)))
  kernel_matrix<-kernel_matrix[,index_trans]
  sir_kernel_matrix<-sir_kernel(X,Y,H2)
  time1=Sys.time()
  ##column estimation #######
  sfInit(parallel = TRUE,cpus=15)
  sfLibrary(glmnet)
  sparse_estimate<-sfApply(kernel_matrix,margin=2,fun=Sparse_vector,design_matrix=X,nfolds=nfolds,lambda=lambda)
  sfStop()
  index<-c()
  for(i in 1:dim(sparse_estimate)[2])
  {
    if(all(sparse_estimate[,i]==0))
    {
      index=append(index,i)
    }
  }
  index<-c(1:num)[-index]
  if(length(index)==0)
  {
    index=sample(1:num,1)
  }
  print(index)
  for(t in index)
  {
    sparse_estimate[,t]=ADMM_estimation_rough(X,Y,H1,H2,t,method="save")
  }
  index<-c()
  for(i in 1:dim(sparse_estimate)[2])
  {
    if(all(sparse_estimate[,i]==0))
    {
      index=append(index,i)
    }
  }
  index<-c(1:num)[-index]
  if(length(index)==0)
  {
    index=sample(1:num,1)
  }
  print(index)
  sparse_save_column=sparse_estimate[,index]
  sir_estimate<-matrix(0,p,ncol(sir_kernel_matrix))
  sfInit(parallel = TRUE,cpus=15)
  sfLibrary(glmnet)
  sir_estimate<- sfApply(sir_kernel_matrix,margin=2,fun=Sparse_vector,design_matrix=X,nfolds=nfolds,lambda=seq(0.01,0.2,0.01))
  sfStop()
  sparse_ens_column<-cbind(sparse_save_column,sir_estimate)
  if(length(index)==1)
  {
    sparse_save_column_index<-1
  }else
  {
    sparse_save_column_index<-forward_column_selection(sparse_save_column)
  }
  sparse_ens_column_index<-forward_column_selection(sparse_ens_column)
  sparse_save_column_index<-index[sparse_save_column_index]
  sparse_ens_column_index<-c(index,c((p*H1+1):((p*H1+H2))))[sparse_ens_column_index]
  if(length(sparse_ens_column_index)==0)
  {
    sparse_ens_column_index=sample(1:num,1)
  }
  print(sparse_save_column_index)
  print(sparse_ens_column_index)
  save_fit<-ADMM_estimation(X,Y,H1,H2,sparse_save_column_index,method="save")
  ens_fit<-ADMM_estimation(X,Y,H1,H2,sparse_ens_column_index,method="ens")
  save_column<-save_fit
  ens_column<-ens_fit
  weight_save<-1/apply(as.matrix(save_column),MARGIN=1, choose_weight)
  weight_ens<-1/apply(as.matrix(ens_column),MARGIN=1, choose_weight)
  #lam1_save<-save_fit$lam1
  #lam1_ens<-ens_fit$lam1
  if(allzero(weight_save))
  {
    save_column<-matrix(0,nrow(as.matrix(save_column)),ncol(as.matrix(save_column)))
  }else
  {
    save_column<-ADMM_estimation_weight(X,Y,H1,H2,weight_save,sparse_save_column_index,method="save")
  }
  if(allzero(weight_ens))
  {
    save_column<-matrix(0,nrow(as.matrix(ens_column)),ncol(as.matrix(ens_column)))
  }else
  {
    ens_column<-ADMM_estimation_weight(X,Y,H1,H2,weight_ens,sparse_ens_column_index,method="ens")
  }
  #solution_equal<-Matrix_estimation(data=t(t(X)/weight),ytilde=kernel_matrix,lambda=lambda)
  #result<-do.call(cbind,(lapply(solution_equal,deal_list)))/weight
  return(list(save_column=save_column,ens_column=ens_column))
}

#norm(proj(beta)-proj(svd(save_column)$u[,1:2]),type="2")
#norm(proj(beta)-proj(svd(ens_column)$u[,1:2]),type="2")

order_determination<-function(X,Y,H=2,r=20)
{
  p<-ncol(as.matrix(X))
  n<-nrow(as.matrix(X))
  s<-5
  if(p==1)
  {
    return(1)
  }
  hnzero<-rep(0,p+1)
  lambda<-eigen(mhat_dr(x=X,y=Y,H=H))$values
  lambda[p+1]=0
  for(i in 1:r)
  {
    augX<-cbind(X,matrix(rnorm(s*n,sd=0.5),n,s))
    V=eigen(mhat_dr(x=augX,y=Y,H=H))
    for(j in 2:(p+1))
    {
      hnzero[j]=hnzero[j]+sum(V$vectors[(p+1):(p+s),j-1]^2)/r
    }
  }
  phi=rep(0,p)
  phi[1]=lambda[1]/(1+lambda[1])
  phi[2:(p+1)]=lambda[2:(p+1)]/sum(lambda)+hnzero[2:(p+1)] 
  order_stat=which.min(phi)-1
  return(order_stat)
}


Model5<-function(X)
{
  p=ncol(X)
  n=nrow(X)
  beta<-rep(0,p)
  beta[1:5]<-rep(1,5)
  prob<-1/(1+exp((X%*%beta)^2-rep(t(beta)%*%cov(X)%*%beta*pchisq(0.5,1),n)))
  Y=rbinom(n,1,prob)
  return(Y)
}

Model6<-function(X)
{
  p=ncol(X)
  beta1<-rep(0,p)
  beta2<-rep(0,p)
  beta1[1:3]<-1
  beta2[(p-2):p]<-1
  Y=0.4*(X%*%beta1)^2+3*abs((X%*%beta2))^0.5+0.2*rnorm(n)
  return(Y)
}

Model7<-function(X)
{
  p=ncol(X)
  beta=rep(0,p)
  beta[1:2]=1
  Y=exp(X%*%beta)+rnorm(n)
  return(Y)
}

Model8<-function(X)
{
  p=ncol(X)
  beta1<-rep(0,p)
  beta2<-rep(0,p)
  beta1[1:2]<-1
  beta2[(p-1):p]<-1
  Y=(2*as.integer(sign(X%*%beta2)>0)-1)*(sign(X%*%beta1)*abs(X%*%beta1)^(1/3)+0.5)+0.2*rnorm(n)
  return(Y)
}

Exp_model<-function(X,idx)
{
  if(idx==5)
  {
    return(Model5(X))
  }
  if(idx==6)
  {
    return(Model6(X))
  }
  if(idx==7)
  {
    return(Model7(X))
  }
  if(idx==8)
  {
    return(Model8(X))
  }
}

beta_index<-function(x,p)
{
  if(x==5)
  {
    beta<-rep(0,p)
    beta[1:5]<-rep(1,5)
  }
  if(x==6)
  {
    beta1<-rep(0,p)
    beta2<-rep(0,p)
    beta1[1:3]<-1
    beta2[(p-2):p]<-1
    beta=cbind(beta1,beta2)
  }
  if(x==7)
  {
    beta=rep(0,p)
    beta[1:2]=1
  }
  if(x==8)
  {
    beta1<-rep(0,p)
    beta2<-rep(0,p)
    beta1[1:2]<-1
    beta2[(p-1):p]<-1
    beta=cbind(beta1,beta2)
  }
  return(beta)
}

cov_decay<-function(p,rho)
{
  cov1<-diag(p)
  for(i in 1:p)
  {
    for(j in 1:p)
    {
      cov1[i,j]=rho^(abs(i-j))
    }
  }
  return(cov1)
}

cov_dense_biway<-function(p,rho,a)
{
  cov<-diag(p)
  for(i in 1:p)
  {
    for(j in 1:p)
    {
      if((i<=a)||(j<=a)||(i>=(p-a+1))||(j>=(p-a+1)))
      {
        cov[i,j]=rho^abs(i-j)
      }else
        if(i==j)
        {
          cov[i,j]=1
        }
      else
      {
        cov[i,j]=rho
      }
    }
  }
  return(cov)
}

cov_dense_oneway<-function(p,rho,a)
{
  cov<-diag(p)
  for(i in 1:p)
  {
    for(j in 1:p)
    {
      if((i<=a)||(j<=a))
      {
        cov[i,j]=rho^abs(i-j)
      }else
        if(i==j)
        {
          cov[i,j]=1
        }
      else
      {
        cov[i,j]=rho
      }
    }
  }
  return(cov)
}

p=200
n=200
cov=cov_decay(p,rho=0.5)
#cov=cov_dense_oneway(p,rho=0.5,a=2)
time=0
count=0
supp_true<-c(1,2,199,200)
model_index=8
cov_index=1
beta<-beta_index(model_index,p)
#beta<-cbind(beta1,beta2)
H1=5
H2=5
Cv_save_column=0
Cv_ens_column=0
ICv_save_column=0
ICv_ens_column=0
col_save_column=0
col_ens_column=0
hatd_save_column=0
hatd_ens_column=0
loss2_save_column=rep(0,50)
loss2_ens_column=rep(0,50)
for(k in 1:100)
{
  X<-mvrnorm(n=n,rep(0,p),Sigma=cov)
  Y<-Exp_model(X,model_index)
  X=scale(X,scale=FALSE)
  result<-Four_estimate_c(X,Y,H1,H2)
  save_column<-cbind(result$save_column,rep(0,p))
  ens_column<-cbind(result$ens_column,rep(0,p))
  supp_save_column=which(abs(svd(save_column)$u[,1])>1e-6)
  supp_ens_column=which(abs(svd(ens_column)$u[,1])>1e-6)
  col_save_column= col_save_column+ncol(save_column)-1
  col_ens_column=col_ens_column+ncol(ens_column)-1
  col_save_column= col_save_column+ncol(save_column)-1
  col_ens_column=col_ens_column+ncol(ens_column)-1
  if(length(supp_save_column)==0)
  {
    supp_save_column=c(100)
  }
  else
  {
    Cv_save_column=Cv_save_column+length(intersect(supp_save_column,supp_true))
    ICv_save_column=ICv_save_column+length(setdiff(supp_save_column,supp_true))
  }
  if(length(supp_ens_column)==0)
  {
    supp_ens_column=c(100)
  }
  else
  {
    Cv_ens_column = Cv_ens_column + length(intersect(supp_ens_column, supp_true))
    ICv_ens_column = ICv_ens_column + length(setdiff(supp_ens_column, supp_true))
  }
  d_save_column=order_determination(X=X[,supp_save_column],Y=Y,H=5)
  if(d_save_column>ncol(as.matrix(save_column)))
  {
    d_save_column=ncol(as.matrix(save_column))
  }
  outcome_save_column=svd(save_column)$u[,1:d_save_column]
  loss2_save_column[k]=norm(proj(beta)-proj(outcome_save_column),type="2")
  hatd_save_column=hatd_save_column+d_save_column
  print(mean(loss2_save_column[1:k]))
  print(hatd_save_column/k)
  print(Cv_save_column/k)
  print(ICv_save_column/k)
  print(col_save_column/k)
  d_ens_column = order_determination(X = X[, supp_ens_column], Y = Y, H = 5)
  if(d_ens_column>ncol(as.matrix(ens_column)))
  {
    d_ens_column=ncol(as.matrix(ens_column))
  }
  outcome_ens_column = svd(ens_column)$u[, 1:d_ens_column]
  loss2_ens_column[k] = norm(proj(beta) - proj(outcome_ens_column), type = "2")
  hatd_ens_column=hatd_ens_column+d_ens_column
  print(mean(loss2_ens_column[1:k]))
  print(hatd_ens_column/k)
  print(Cv_ens_column/k)
  print(ICv_ens_column/k)
  print(col_ens_column/k)
}

print(sd(loss2_save_column)/sqrt(100))
print(sd(loss2_ens_column)/sqrt(100))

p=200
n=200
cov=cov_decay(p,rho=0.5)
#cov=cov_dense_oneway(p,rho=0.5,a=2)
time=0
count=0
supp_true<-c(1,2)
model_index=7
cov_index=1
beta<-beta_index(model_index,p)
#beta<-cbind(beta1,beta2)
H1=10
H2=10
Cv_a8=0
Cv_a10=0
ICv_a8=0
ICv_a10=0
loss2_a8=rep(0,200)
loss2_a10=rep(0,200)
for(k in 1:200)
{
  X<-mvrnorm(n=n,rep(0,p),Sigma=cov)
  Y<-Exp_model(X,model_index)
  X=scale(X,scale=FALSE)
  result=cv.seas(X,Y,H=5,d=1)$beta
  supp=which(result!=0)
  loss2_a8[k]=norm(proj(beta) - proj(result), type = "2")
  Cv_a8 = Cv_a8+ length(intersect(supp, supp_true))
  ICv_a8 = ICv_a8 + length(setdiff(supp, supp_true))
}

cov=cov_dense_oneway(p,rho=0.5,a=2)


args <- commandArgs(trailingOnly = TRUE)
output_dir <- args[1]

# 确保输出目录存在
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

data=c(mean(loss2_save_column[1:k]),hatd_save_column/k,Cv_save_column/k,ICv_save_column/k,col_save_column/k,sd(loss2_save_column)/sqrt(100),mean(loss2_ens_column[1:k]),hatd_ens_column/k,Cv_ens_column/k,ICv_ens_column/k,col_ens_column/k,sd(loss2_ens_column)/sqrt(100))
# 生成CSV文件
output_file <- file.path(output_dir, paste0("Model(", model_index, ")Sig(", cov_index, ")p(", p, ")H10.csv"))
write.csv(data, output_file)
