### Functions for Generalized Linear Models with Peeling and Instruments (GAMPI)
require('glmnet'); require('mltools'); require('igraph')

### GLMTLP
glmtlp <- function(X,y,lambda,tau = .5,beta_initial=NULL,intercept=TRUE,family="gaussian"){
  
  if (!intercept){ ### Fit with No intercept
    ### Initial estimate 
    if (is.numeric(beta_initial)){
      beta = beta_initial
    } else if (is.character(beta_initial)){
      if (beta_initial == "GLM"){
        beta = glm(y~X-1,family = family)$coefficients
      } else if (beta_initial == "Lasso"){
        beta = glmnet(X,y,intercept = F,standardize = F,lambda = lambda,family = family)$beta
      }
    } else if (is.null(beta_initial)){
      ### Default is the Lasso estimate
      beta = glmnet(X,y,intercept = F,standardize = F,lambda = lambda,family = family)$beta
    }
    
    beta_m = as.matrix(as.vector(beta),nrow = q)
    maxIter = 100; iter = 1; Tol = 1e-5
    obj1 = 1; obj0 = 0; obj = c(); obj.optimal = Inf; beta.optimal = NULL
    while (iter <= maxIter && abs(obj1-obj0) > Tol) {
      iter = iter + 1
      obj0 = obj1
      
      ###### Solve the subproblem
      if (sum(abs(beta_m) <= tau) == 0){
        beta = glmnet(X,y,intercept = F,standardize = F,lambda = 0,family = family)$beta
      } else {
        beta = glmnet(X,y,intercept = F,standardize = F,lambda = lambda/tau,penalty.factor = as.numeric(abs(beta_m) <= tau),family = family)$beta
      }
      obj1 = objective(X,y,beta,lambda/tau * as.numeric(abs(beta_m) <= tau),family = family)
      beta_m = beta
      obj = c(obj,obj1)
      
      if (obj1 == Inf){
        if (iter == 2){
          beta.optimal = beta
        }
        break
      }
      # store optimal solution
      if (obj1 < obj.optimal){
        obj.optimal = obj1; beta.optimal = beta;
      }
    }
    
    beta = as.vector(beta.optimal)
    return(list(beta=beta,obj=obj))
    
  } else { ### Fit With intercept
    ### Initial estimate 
    if (is.numeric(beta_initial)){
      beta = beta_initial
    } else if (is.character(beta_initial)){
      if (beta_initial == "GLM"){
        beta = glm(y~X,family = family)$coefficients[-1]
      } else if (beta_initial == "Lasso"){
        beta = glmnet(X,y,intercept = T,standardize = F,lambda = lambda,family = family)$beta
      }
    } else if (is.null(beta_initial)){
      ### Default is the Lasso estimate
      beta = glmnet(X,y,intercept = T,standardize = F,lambda = lambda,family = family)$beta
    }
    
    beta_m = as.matrix(as.vector(beta),nrow = q)
    maxIter = 100; iter = 1; Tol = 1e-5
    obj1 = 1; obj0 = 0; obj = c(); obj.optimal = Inf; beta.optimal = NULL; beta0.optimal = NULL;
    while (iter <= maxIter && abs(obj1-obj0) > Tol) {
      iter = iter + 1
      obj0 = obj1
      
      ###### Solve the subproblem
      if (sum(abs(beta_m) <= tau) == 0){
        fit = glmnet(X,y,intercept = T,standardize = F,lambda = 0,family = family)
        beta = fit$beta
        beta0 = fit$a0
      } else {
        fit = glmnet(X,y,intercept = T,standardize = F,lambda = lambda/tau,penalty.factor = as.numeric(abs(beta_m) <= tau),family = family)
        beta = fit$beta
        beta0 = fit$a0
      }
      obj1 = objective(cbind(1,X),y,c(beta0,as.vector(beta)),c(0,lambda/tau * as.numeric(abs(beta_m) <= tau)),family=family)
      beta_m = beta
      obj = c(obj,obj1)
      
      if (obj1 == Inf){
        if (iter == 2){
          beta.optimal = beta; beta0.optimal = beta0;
        }
        break
      }
      # store optimal solution
      if (obj1 < obj.optimal){
        obj.optimal = obj1; beta.optimal = beta; beta0.optimal = beta0;
      }
    }
    
    beta = as.vector(beta.optimal); beta0 = beta0.optimal
    return(list(beta=beta,beta0=beta0,obj=obj))
  }
}

### Compute objective function
objective <- function(X,y,beta,lambda,family="gaussian"){
  n = dim(X)[1]
  if (family == "gaussian"){
    sum((y - X%*% beta)^2)/(2*n) + sum(lambda * abs(beta)) #lambda is a vector
  } else if (family == "binomial"){
    sum(log(1 + exp( X%*% beta)) - y *X%*% beta)/n + sum(lambda * abs(beta))
  } else if (family == "poisson"){
    sum( exp( X%*% beta) - y *X%*% beta)/n + sum(lambda * abs(beta))
  }
}

### refit truncated (top K) GLM
trun_K_glm <- function(X,y,beta_tlp,K_list = 1:ncol(X),intercept=TRUE,family="gaussian"){
  
  n = nrow(X)
  q = ncol(X)
  
  if (!intercept){
    beta_matrix = matrix(rep(0,q*length(K_list)), ncol = length(K_list))
    fitted_values = matrix(rep(0,n*length(K_list)), ncol = length(K_list))
    
    if (sum(abs(beta_tlp)) != 0){
      for (s in 1:length(K_list)){
        K = K_list[s]
        if (K<q){
          indice = order(abs(beta_tlp), decreasing = T)[K+1]
          ### find nonzero index
          act_set = which(abs(beta_tlp) > abs(beta_tlp[indice]))
          if (length(act_set)==0){
            act_set = order(abs(beta_tlp), decreasing = T)[1:K]
          }
        } else if (K==q){
          act_set = which(beta_tlp != 0)
        }
        ### refit
        fit = glm(y~X[,act_set]-1,family = family)
        beta_refit = fit$coefficients
        beta = rep(0,q)
        beta[act_set] = beta_refit
        beta_matrix[,s] = beta
        fitted_values[,s] = fit$fitted.values
      } 
    }
    return(list(beta_matrix=beta_matrix,fitted_values=fitted_values))
  } else { # Fit with intercept
    
    beta_matrix = matrix(rep(0,q*length(K_list)), ncol = length(K_list))
    beta0_vec = rep(0,length(K_list))
    fitted_values = matrix(rep(0,n*length(K_list)), ncol = length(K_list))
    
    if (sum(abs(beta_tlp)) != 0){
      for (s in 1:length(K_list)){
        K = K_list[s]
        if (K<q){
          indice = order(abs(beta_tlp), decreasing = T)[K+1]
          ### find nonzero index
          act_set = which(abs(beta_tlp) > abs(beta_tlp[indice]))
          if (length(act_set)==0){
            act_set = order(abs(beta_tlp), decreasing = T)[1:K]
          }
        } else if (K==q){
          act_set = which(beta_tlp != 0)
        }
        ### refit
        fit = glm(y~X[,act_set],family = family)
        beta_refit = fit$coefficients[-1]
        beta = rep(0,q)
        beta[act_set] = beta_refit
        beta_matrix[,s] = beta
        beta0_vec[s] = fit$coefficients[1]
        fitted_values[,s] = fit$fitted.values
      } 
    }
    return(list(beta_matrix=beta_matrix,beta0_vec=beta0_vec,fitted_values=fitted_values))
  }
}

### GLMTLP with weights
glmtlp_weighted <- function(X,y,lambda,tau = .5,beta_initial=NULL,intercept=TRUE,is_penalize = rep(1,ncol(X)),family="gaussian"){
  
  if (!intercept){ ### Fit with No intercept
    ### Initial estimate 
    if (is.numeric(beta_initial)){
      beta = beta_initial
    } else if (is.character(beta_initial)){
      if (beta_initial == "GLM"){
        beta = glm(y~X-1,family = family)$coefficients
      } else if (beta_initial == "Lasso"){
        beta = glmnet(X,y,intercept = F,standardize = F,lambda = lambda,family = family,penalty.factor = is_penalize)$beta
      }
    } else if (is.null(beta_initial)){
      ### Default is the Lasso estimate
      beta = glmnet(X,y,intercept = F,standardize = F,lambda = lambda,family = family,penalty.factor = is_penalize)$beta
    }
    
    beta_m = as.matrix(as.vector(beta),nrow = q)
    maxIter = 100; iter = 1; Tol = 1e-5
    obj1 = 1; obj0 = 0; obj = c(); obj.optimal = Inf; beta.optimal = NULL
    while (iter <= maxIter && abs(obj1-obj0) > Tol) {
      iter = iter + 1
      obj0 = obj1
      
      ###### Solve the subproblem
      if (sum(as.numeric(abs(beta_m) <= tau) * is_penalize) == 0){
        beta = glmnet(X,y,intercept = F,standardize = F,lambda = 0,family = family)$beta
      } else {
        beta = glmnet(X,y,intercept = F,standardize = F,lambda = lambda/tau,penalty.factor = as.numeric(abs(beta_m) <= tau) * is_penalize,family = family)$beta
      }
      obj1 = objective(X,y,beta,lambda/tau * as.numeric(abs(beta_m) <= tau) * is_penalize,family=family)
      beta_m = beta
      obj = c(obj,obj1)
      
      if (obj1 == Inf){
        if (iter == 2){
          beta.optimal = beta
        }
        break
      }
      # store optimal solution
      if (obj1 < obj.optimal){
        obj.optimal = obj1; beta.optimal = beta;
      }
    }
    
    beta = as.vector(beta.optimal)
    return(list(beta=beta,obj=obj))
    
  } else { ### Fit With intercept
    ### Initial estimate 
    if (is.numeric(beta_initial)){
      beta = beta_initial
    } else if (is.character(beta_initial)){
      if (beta_initial == "GLM"){
        beta = glm(y~X,family = family)$coefficients[-1]
      } else if (beta_initial == "Lasso"){
        beta = glmnet(X,y,intercept = T,standardize = F,lambda = lambda,family = family,penalty.factor = is_penalize)$beta
      }
    } else if (is.null(beta_initial)){
      ### Default is the Lasso estimate
      beta = glmnet(X,y,intercept = T,standardize = F,lambda = lambda,family = family,penalty.factor = is_penalize)$beta
    }
    
    beta_m = as.matrix(as.vector(beta),nrow = q)
    maxIter = 100; iter = 1; Tol = 1e-5
    obj1 = 1; obj0 = 0; obj = c(); obj.optimal = Inf; beta.optimal = NULL; beta0.optimal = NULL;
    while (iter <= maxIter && abs(obj1-obj0) > Tol) {
      iter = iter + 1
      obj0 = obj1
      
      ###### Solve the subproblem
      if (sum(as.numeric(abs(beta_m) <= tau) * is_penalize) == 0){
        fit = glmnet(X,y,intercept = T,standardize = F,lambda = 0,family = family)
        beta = fit$beta
        beta0 = fit$a0
      } else {
        fit = glmnet(X,y,intercept = T,standardize = F,lambda = lambda/tau,penalty.factor = as.numeric(abs(beta_m) <= tau) * is_penalize,family = family)
        beta = fit$beta
        beta0 = fit$a0
      }
      obj1 = objective(cbind(1,X),y,c(beta0,as.vector(beta)),c(0,lambda/tau * as.numeric(abs(beta_m) <= tau)) * c(0,is_penalize),family=family )
      beta_m = beta
      obj = c(obj,obj1)
      
      if (obj1 == Inf){
        if (iter == 2){
          beta.optimal = beta; beta0.optimal = beta0;
        }
        break
      }
      # store optimal solution
      if (obj1 < obj.optimal){
        obj.optimal = obj1; beta.optimal = beta; beta0.optimal = beta0;
      }
    }
    
    beta = as.vector(beta.optimal); beta0 = beta0.optimal
    return(list(beta=beta,beta0=beta0,obj=obj))
  }
}

trun_K_glm_weighted <- function(X,y,beta_tlp,K_list = 1:ncol(X),intercept=TRUE,is_penalize = rep(1,ncol(X)),family="gaussian"){
  
  n = nrow(X)
  q = ncol(X)
  beta_tlp[is_penalize==0] = max(abs(beta_tlp)) + 1; # make unpenalized coefs the largest ones so that they are always included in the model

  if (!intercept){
    beta_matrix = matrix(rep(0,q*length(K_list)), ncol = length(K_list))
    fitted_values = matrix(rep(0,n*length(K_list)), ncol = length(K_list))
    
    if (sum(abs(beta_tlp)) != 0){
      for (s in 1:length(K_list)){
        K = K_list[s]
        if (K<q){
          indice = order(abs(beta_tlp), decreasing = T)[K+1]
          ### find nonzero index
          act_set = which(abs(beta_tlp) > abs(beta_tlp[indice]))
          if (length(act_set)==0){
            act_set = order(abs(beta_tlp), decreasing = T)[1:K]
          }
        } else if (K==q){
          act_set = which(beta_tlp != 0)
        }
        ### refit
        fit = glm(y~X[,act_set]-1,family = family)
        beta_refit = fit$coefficients
        beta = rep(0,q)
        beta[act_set] = beta_refit
        beta_matrix[,s] = beta
        fitted_values[,s] = fit$fitted.values
      } 
    }
    return(list(beta_matrix=beta_matrix,fitted_values=fitted_values))
  } else { # Fit with intercept
    
    beta_matrix = matrix(rep(0,q*length(K_list)), ncol = length(K_list))
    beta0_vec = rep(0,length(K_list))
    fitted_values = matrix(rep(0,n*length(K_list)), ncol = length(K_list))
    
    if (sum(abs(beta_tlp)) != 0){
      for (s in 1:length(K_list)){
        K = K_list[s]
        if (K<q){
          indice = order(abs(beta_tlp), decreasing = T)[K+1]
          ### find nonzero index
          act_set = which(abs(beta_tlp) > abs(beta_tlp[indice]))
          if (length(act_set)==0){
            act_set = order(abs(beta_tlp), decreasing = T)[1:K]
          }
        } else if (K==q){
          act_set = which(beta_tlp != 0)
        }
        ### refit
        fit = glm(y~X[,act_set],family = family)
        beta_refit = fit$coefficients[-1]
        beta = rep(0,q)
        beta[act_set] = beta_refit
        beta_matrix[,s] = beta
        beta0_vec[s] = fit$coefficients[1]
        fitted_values[,s] = fit$fitted.values
      } 
    }
    return(list(beta_matrix=beta_matrix,beta0_vec=beta0_vec,fitted_values=fitted_values))
  }
}


############################ Tuning parameter selection
### Parameter selection for GLMTLP using EBIC
BIC.unconstrained_glmtlp <- function(X,y,lambdas=NULL,K_list = 1:(ncol(X)-1),taus = c(0.5,1,2),intercept=TRUE,nlambda = 10,lambda.min=0.01,gamma = .5,family="gaussian"){
  n = dim(X)[1]
  p = dim(X)[2]
  
  if (is.null(lambdas)){
    lambda_max = max(abs(t(X) %*% y / n))
    lambdas = seq(lambda.min,lambda_max,length = nlambda)
  }
  
  nlambda = length(lambdas)
  nK = length(K_list)
  ntau = length(taus)
  
  dims <- c(nlambda,nK,ntau)
  BIC = array(rep(Inf,nlambda*nK*ntau),dims) # nlambda * K * ntau
  
  if (!intercept){
    for (j in 1:nlambda){
      for (m in 1:ntau){
        
        result = glmtlp(X,y,lambda = lambdas[j],tau = taus[m],intercept = intercept,family=family)
        
        beta = result$beta
        K_list_feasible = K_list[K_list <= sum(beta!=0)]
        
        if (length(K_list_feasible) > 0){
          beta_matrix = trun_K_glm(X,y,beta,K_list_feasible,intercept = intercept,family=family)$beta_matrix   ### q * K_list
          
          for (l in K_list_feasible){
            idx = which(K_list_feasible == l)
            idx_Klist = which(K_list == l)
            
            ### Calculate EBIC
            BIC[j,idx_Klist,m] = BIC_criterion(X,y,beta = beta_matrix[,idx],gamma = gamma,family=family)
          }
        }
      }
    }
  } else { # Fit with intercept
    for (j in 1:nlambda){
      for (m in 1:ntau){
        
        result = glmtlp(X,y,lambda = lambdas[j],tau = taus[m],intercept = intercept,family=family)
        
        beta = result$beta
        K_list_feasible = K_list[K_list <= sum(beta!=0)]
        
        if (length(K_list_feasible) > 0){
          trun_glm_result = trun_K_glm(X,y,beta,K_list_feasible,intercept = intercept,family=family)   ### q * K_list
          beta_matrix = trun_glm_result$beta_matrix
          beta0_vec = trun_glm_result$beta0_vec
          
          for (l in K_list_feasible){
            idx = which(K_list_feasible == l)
            idx_Klist = which(K_list == l)
            
            ### Calculate EBIC
            BIC[j,idx_Klist,m] = BIC_criterion(X,y,beta0=beta0_vec[idx],beta = beta_matrix[,idx],gamma = gamma,family=family)
          }
        }
      }
    }
  }
  
  optpar.min.indice = which.min(BIC)
  optpar.comb = arrayInd(optpar.min.indice,dim(BIC))
  
  return(list(optpar.lambda = lambdas[optpar.comb[1]], optpar.K = K_list[optpar.comb[2]], optpar.tau = taus[optpar.comb[3]],BIC=BIC))
  
}

### Extended BIC: http://www3.stat.sinica.edu.tw/statistica/oldpdf/A22n26.pdf
BIC_criterion <- function(X,y,beta0 = 0,beta,gamma = .5,family="gaussian"){
  n = dim(X)[1]; s = sum(beta!=0); p = length(beta); 
  
  if (family == "gaussian"){
    n * log(sum((y - (X%*% beta + beta0))^2)) - n * log(n) +  s * log(n) + 2 * gamma *  s * log(p)
  } else if (family == "binomial"){
    2 * sum(log(1 + exp( X%*% beta + beta0 )) - y * (X%*% beta + beta0)) +  s * log(n) + 2 * gamma *  s * log(p)
  } else if (family == "poisson"){
    2 * sum( exp( X%*% beta + beta0 ) - y * (X%*% beta + beta0)) +  s * log(n) + 2 * gamma *  s * log(p)
  } 
}

### Parameter selection for weighted GLMTLP using EBIC (used in the second step of GAMPI)
BIC.unconstrained_glmtlp_weighted <- function(X,y,lambdas=NULL,K_list = 1:(ncol(X)-1),taus = c(0.5,1,2),intercept=TRUE,is_penalize = rep(1,ncol(X)),nlambda = 10,lambda.min=0.01,gamma = .5,family="gaussian"){
  n = dim(X)[1]
  p = dim(X)[2]
  
  if (is.null(lambdas)){
    lambda_max = max(abs(t(X) %*% y / n))
    lambdas = seq(lambda.min,lambda_max,length = nlambda)
  }
  
  nlambda = length(lambdas)
  nK = length(K_list)
  ntau = length(taus)
  
  dims <- c(nlambda,nK,ntau)
  BIC = array(rep(Inf,nlambda*nK*ntau),dims) # nlambda * K * ntau
  
  if (!intercept){
    for (j in 1:nlambda){
      for (m in 1:ntau){
        
        result = glmtlp_weighted(X,y,lambda = lambdas[j],tau = taus[m],intercept = intercept,is_penalize=is_penalize,family=family)
        
        beta = result$beta
        K_list_feasible = K_list[K_list <= sum(beta!=0)]
        
        if (length(K_list_feasible) > 0){
          beta_matrix = trun_K_glm_weighted(X,y,beta,K_list_feasible,intercept = intercept,is_penalize=is_penalize,family=family)$beta_matrix   ### q * K_list
          
          for (l in K_list_feasible){
            idx = which(K_list_feasible == l)
            idx_Klist = which(K_list == l)
            
            ### Calculate EBIC
            BIC[j,idx_Klist,m] = BIC_criterion(X,y,beta = beta_matrix[,idx],gamma = gamma,family=family)
          }
        }
      }
    }
  } else { # Fit with intercept
    for (j in 1:nlambda){
      for (m in 1:ntau){
        
        result = glmtlp_weighted(X,y,lambda = lambdas[j],tau = taus[m],intercept = intercept,is_penalize=is_penalize,family=family)
        
        beta = result$beta
        K_list_feasible = K_list[K_list <= sum(beta!=0)]
        
        if (length(K_list_feasible) > 0){
          trun_glm_result = trun_K_glm_weighted(X,y,beta,K_list_feasible,intercept = intercept,is_penalize=is_penalize,family=family)   ### q * K_list
          beta_matrix = trun_glm_result$beta_matrix
          beta0_vec = trun_glm_result$beta0_vec
          
          for (l in K_list_feasible){
            idx = which(K_list_feasible == l)
            idx_Klist = which(K_list == l)
            
            ### Calculate EBIC
            BIC[j,idx_Klist,m] = BIC_criterion(X,y,beta0=beta0_vec[idx],beta = beta_matrix[,idx],gamma = gamma,family=family)
          }
        }
      }
    }
  }
  
  optpar.min.indice = which.min(BIC)
  optpar.comb = arrayInd(optpar.min.indice,dim(BIC))
  
  return(list(optpar.lambda = lambdas[optpar.comb[1]], optpar.K = K_list[optpar.comb[2]], optpar.tau = taus[optpar.comb[3]],BIC=BIC))
}


### Metrics (GLM deviance) for cross validation
GLM_deviance <- function(yhat,y,family="gaussian"){
  
  if (family == "gaussian"){
    D = sum((y-yhat)^2)
  } else if (family == "binomial"){
    ### Binomial deviance function
    eps = 1e-5
    yhat = exp(yhat) / ( 1 + exp(yhat) )
    y[y==0] = eps; y[y==1] = 1-eps; yhat[yhat == 0] = eps; yhat[yhat == 1] = 1-eps;
    D = 2 * sum(y * log(y / yhat) + (1-y) * log((1-y) / (1-yhat)))
    # or without having eps
    # D = -2 * sum(y * log(yhat) + (1-y) * log(1-yhat))
  } else if (family == "poisson"){
    ### Poisson deviance function
    eps = 1e-5
    yhat = exp(yhat)
    # https://peijin.medium.com/the-poisson-deviance-for-regression-d469b56959ce
    y[y==0] = eps; yhat[yhat == 0] = eps;
    D = 2 * sum(y * log(y / yhat)  - (y - yhat) )
    # or
    # D = -2 * sum(  y * log(yhat)  - yhat )
  }
  return(D)
}

### Parameter selection for GLMTLP using cross validation
cv.unconstrained_glmtlp_K_tau <- function(X,y,lambdas=NULL,fold = 5,type.measure = "deviance",K_list = 1:(ncol(X)-1),taus = c(0.5,1,2),intercept=TRUE,nlambda = 10,lambda.min=0.01,family="gaussian"){
  n = dim(X)[1]
  p = dim(X)[2]
  
  if (is.null(lambdas)){
    lambda_max = max(abs(t(X) %*% y / n))
    lambdas = seq(lambda.min,lambda_max,length = nlambda)
  }
  
  shuffle_index = sample(1:n,n)
  nlambda = length(lambdas)
  nK = length(K_list)
  ntau = length(taus)
  
  dims <- c(nlambda,nK,ntau,fold)
  CV_errs = array(rep(Inf,nlambda*fold*nK*ntau),dims) # nlambda * nfold * K * ntau
  
  if (!intercept){
    for (i in 1:fold){
      
      indice = shuffle_index[((i-1)*n/fold + 1):(i*n/fold)]
      Xtrain = X[-indice,]; ytrain = y[-indice]
      Xtest = X[indice,]; ytest = y[indice]
      
      for (j in 1:nlambda){
        for (m in 1:ntau){
          
          result = glmtlp(Xtrain,ytrain,lambda = lambdas[j],tau = taus[m],intercept = intercept,family=family)
          
          beta = result$beta
          K_list_feasible = K_list[K_list <= sum(beta!=0)]
          
          if (length(K_list_feasible) > 0){
            beta_matrix = trun_K_glm(Xtrain,ytrain,beta,K_list_feasible,intercept = intercept,family=family)$beta_matrix   ### q * K_list
            
            for (l in K_list_feasible){
              idx = which(K_list_feasible == l)
              idx_Klist = which(K_list == l)
              
              yhat = Xtest %*% beta_matrix[,idx]
              
              ### Calculate error
              if (type.measure == "class"){
                if (family == "binomial"){
                  yhat = exp(yhat) / ( 1 + exp(yhat) )
                  yhat = (yhat > 0.5)
                  CV_errs[j,idx_Klist,m,i] = sum(ytest != yhat)
                }
              } else if (type.measure == "deviance"){
                # GLM deviance
                CV_errs[j,idx_Klist,m,i] = GLM_deviance(yhat,ytest,family=family)
              }
            }
          }
        }
      }
    }
  } else { # Fit with intercept
    
    for (i in 1:fold){
      
      indice = shuffle_index[((i-1)*n/fold + 1):(i*n/fold)]
      Xtrain = X[-indice,]; ytrain = y[-indice]
      Xtest = X[indice,]; ytest = y[indice]
      
      for (j in 1:nlambda){
        for (m in 1:ntau){
          
          result = glmtlp(Xtrain,ytrain,lambda = lambdas[j],tau = taus[m],intercept = intercept,family=family)
          
          beta = result$beta
          K_list_feasible = K_list[K_list <= sum(beta!=0)]
          
          if (length(K_list_feasible) > 0){
            
            trun_glm_result = trun_K_glm(Xtrain,ytrain,beta,K_list_feasible,intercept = intercept,family=family)   ### q * K_list
            beta_matrix = trun_glm_result$beta_matrix
            beta0_vec = trun_glm_result$beta0_vec
            
            for (l in K_list_feasible){
              idx = which(K_list_feasible == l)
              idx_Klist = which(K_list == l)
              
              yhat = Xtest %*% beta_matrix[,idx] + beta0_vec[idx]
              
              ### Calculate error
              if (type.measure == "class"){
                if (family == "binomial"){
                  yhat = exp(yhat) / ( 1 + exp(yhat) )
                  yhat = (yhat > 0.5)
                  CV_errs[j,idx_Klist,m,i] = sum(ytest != yhat)
                }
              } else if (type.measure == "deviance"){
                # GLM deviance
                CV_errs[j,idx_Klist,m,i] = GLM_deviance(yhat,ytest,family=family)
              }
            }
          }
        }
      }
    }
  }
  
  
  CV_err = apply(CV_errs,c(1,2,3),mean)
  
  optpar.min.indice = which.min(CV_err)
  optpar.comb = arrayInd(optpar.min.indice, dim(CV_err))
  
  # one SE rule
  SE = sqrt(apply(CV_errs,c(1,2,3),var)/fold)
  minSE = CV_err[which.min(CV_err)]+SE[which.min(CV_err)]
  CV_err_K = CV_err[optpar.comb[1], ,optpar.comb[3]]
  minSEx = which(CV_err_K<=minSE)[1]
  optlam.SE = K_list[minSEx]

  return(list(optpar.lambda = lambdas[optpar.comb[1]], optpar.K.min = K_list[optpar.comb[2]], optpar.K.1se = K_list[minSEx], optpar.tau = taus[optpar.comb[3]],CV_errs=CV_errs))
  
}

### Parameter selection for weighted GLMTLP using cross validation (used in the second step of GAMPI)
cv.unconstrained_glmtlp_K_tau_weighted <- function(X,y,lambdas=NULL,fold = 5,type.measure = "deviance",K_list = 1:(ncol(X)-1),taus = c(0.5,1,2),intercept=TRUE,is_penalize = rep(1,ncol(X)),nlambda = 10,lambda.min=0.01,family="gaussian"){
  n = dim(X)[1]
  p = dim(X)[2]
  
  if (is.null(lambdas)){
    lambda_max = max(abs(t(X) %*% y / n))
    lambdas = seq(lambda.min,lambda_max,length = nlambda)
  }
  
  shuffle_index = sample(1:n,n)
  nlambda = length(lambdas)
  nK = length(K_list)
  ntau = length(taus)
  
  dims <- c(nlambda,nK,ntau,fold)
  CV_errs = array(rep(Inf,nlambda*fold*nK*ntau),dims) # nlambda * nfold * K * ntau
  
  if (!intercept){
    for (i in 1:fold){
      
      indice = shuffle_index[((i-1)*n/fold + 1):(i*n/fold)]
      Xtrain = X[-indice,]; ytrain = y[-indice]
      Xtest = X[indice,]; ytest = y[indice]
      
      for (j in 1:nlambda){
        for (m in 1:ntau){
          
          result = glmtlp_weighted(Xtrain,ytrain,lambda = lambdas[j],tau = taus[m],intercept = intercept,is_penalize=is_penalize,family=family)
          
          beta = result$beta
          K_list_feasible = K_list[K_list <= sum(beta!=0)]
          
          if (length(K_list_feasible) > 0){
            beta_matrix = trun_K_glm_weighted(Xtrain,ytrain,beta,K_list_feasible,intercept = intercept,is_penalize = is_penalize,family=family)$beta_matrix   ### q * K_list
            
            for (l in K_list_feasible){
              idx = which(K_list_feasible == l)
              idx_Klist = which(K_list == l)
              
              yhat = Xtest %*% beta_matrix[,idx]
              
              ### Calculate error
              if (type.measure == "class"){
                if (family == "binomial"){
                  yhat = exp(yhat) / ( 1 + exp(yhat) )
                  yhat = (yhat > 0.5)
                  CV_errs[j,idx_Klist,m,i] = sum(ytest != yhat)
                }
              } else if (type.measure == "deviance"){
                # GLM deviance
                CV_errs[j,idx_Klist,m,i] = GLM_deviance(yhat,ytest,family=family)
              }
            }
          }
        }
      }
    }
  } else { # Fit with intercept
    
    for (i in 1:fold){
      
      indice = shuffle_index[((i-1)*n/fold + 1):(i*n/fold)]
      Xtrain = X[-indice,]; ytrain = y[-indice]
      Xtest = X[indice,]; ytest = y[indice]
      
      for (j in 1:nlambda){
        for (m in 1:ntau){
          
          result = glmtlp_weighted(Xtrain,ytrain,lambda = lambdas[j],tau = taus[m],intercept = intercept,is_penalize=is_penalize,family=family)
          
          beta = result$beta
          K_list_feasible = K_list[K_list <= sum(beta!=0)]
          
          if (length(K_list_feasible) > 0){
            
            trun_glm_result = trun_K_glm_weighted(Xtrain,ytrain,beta,K_list_feasible,intercept = intercept,is_penalize = is_penalize,family=family)   ### q * K_list
            beta_matrix = trun_glm_result$beta_matrix
            beta0_vec = trun_glm_result$beta0_vec
            
            for (l in K_list_feasible){
              idx = which(K_list_feasible == l)
              idx_Klist = which(K_list == l)
              
              yhat = Xtest %*% beta_matrix[,idx] + beta0_vec[idx]
              
              ### Calculate error
              if (type.measure == "class"){
                if (family == "binomial"){
                  yhat = exp(yhat) / ( 1 + exp(yhat) )
                  yhat = (yhat > 0.5)
                  CV_errs[j,idx_Klist,m,i] = sum(ytest != yhat)
                }
              } else if (type.measure == "deviance"){
                # GLM deviance
                CV_errs[j,idx_Klist,m,i] = GLM_deviance(yhat,ytest,family=family)
              }
            }
          }
        }
      }
    }
  }
  
  
  CV_err = apply(CV_errs,c(1,2,3),mean)
  
  optpar.min.indice = which.min(CV_err)
  optpar.comb = arrayInd(optpar.min.indice, dim(CV_err))
  
  # one SE rule
  SE = sqrt(apply(CV_errs,c(1,2,3),var)/fold)
  minSE = CV_err[which.min(CV_err)]+SE[which.min(CV_err)]
  CV_err_K = CV_err[optpar.comb[1], ,optpar.comb[3]]
  minSEx = which(CV_err_K<=minSE)[1]
  optlam.SE = K_list[minSEx]

  return(list(optpar.lambda = lambdas[optpar.comb[1]], optpar.K.min = K_list[optpar.comb[2]], optpar.K.1se = K_list[minSEx], optpar.tau = taus[optpar.comb[3]],CV_errs=CV_errs))
  
}


### Peeling algorithm
peeling <- function(V){
  ### V: q by p, intervention by outcomes
  V0 = V; q0 = dim(V0)[1];
  V = V[apply(V,1,sum)!=0,]
  
  q = dim(V)[1];   p = dim(V)[2]
  V_t = V
  current_row_index = 1:q
  current_column_index = 1:p
  
  interv_relation = matrix(rep(0,q*p),ncol=p) # q by p
  ances_relation = matrix(rep(0,p*p),ncol=p) # p by p, parent by child relationship
  other_interv_relation = matrix(rep(0,q*p),ncol=p) # q by p
  
  prev_remove_index = c()
  
  while ( length(V_t) != 0){
    if (length(V_t) == 1 || is.vector(V_t)){
      V_t = matrix(V_t,nrow = length(current_row_index), ncol = length(current_column_index))
    }
    non_zero = apply(V_t!=0,1,sum)   ### row sum of non-zero entries, for each instrumental variable 
    min_l0_norm_row_index = which(non_zero == min(non_zero))  ### find instrumental variable with min l0 norm
    max_node_store = c()
    
    for (l in min_l0_norm_row_index){
      max_node = which(abs(V_t[l, ]) == max(abs(V_t[l, ]))); ### find the largest absolute value element index of the lth row
      
      ### X_l intervenes on Y_j
      interv_relation[current_row_index[l],current_column_index[max_node]] = 1;
      
      for (s in 1:p){
        if ( (! s %in%  current_column_index)  & V[current_row_index[l],s] != 0){
          ances_relation[current_column_index[max_node],s] = 1; # previously removed node that is a child of current node
          
          ### Find all the children of the previously removed node
          child = which(ances_relation[s,] == 1)
          ### Identify all ancestral relations
          ances_relation[current_column_index[max_node],child] = 1;
          
        }   
      }
      
      max_node_store = c(max_node_store, max_node);
      
      ### Find all intervention relations
      child = c()
      for (r in max_node){
        child_tmp = which(ances_relation[current_column_index[r],] == 1)
        child = c(child,child_tmp)
      }
      child = unique(child)
      
      other_interv_relation[current_row_index[l],child] = 1
    }
    
    C_t = c()
    A_t_complement = setdiff(1:length(current_row_index), min_l0_norm_row_index);
    B_t_complement = setdiff(1:length(current_column_index),max_node_store);
    
    for (m in A_t_complement){
      if (sum( abs(V_t[m, B_t_complement])) == 0){  # rows with all zeros
        C_t = c(C_t, m)
      }
    }
    
    for (r in C_t){
      non_zero_index = which(V_t[r,]!=0)
      interv_relation[current_row_index[r],current_column_index[non_zero_index]] = 1;
    }
    
    ### Update V_t matrix
    prev_remove_index = current_column_index[max_node_store];
    current_row_index = current_row_index[setdiff(1:length(current_row_index),c(min_l0_norm_row_index,C_t))];
    current_column_index = current_column_index[setdiff(1:length(current_column_index),max_node_store)];        
    
    V_t = V[current_row_index,current_column_index];
  }
  
  
  ### transform intervention relation back to the original intervention matrix
  interv_relation0 = matrix(rep(0,q0*p),ncol=p) # q0 by p
  interv_relation0[apply(V0,1,sum)!=0,] = interv_relation
  interv_relation = interv_relation0
  
  other_interv_relation0 = matrix(rep(0,q0*p),ncol=p) # q0 by p
  other_interv_relation0[apply(V0,1,sum)!=0,] = other_interv_relation
  other_interv_relation = other_interv_relation0
  
  return(list(interv_relation = interv_relation,ances_relation = ances_relation,other_interv_relation = other_interv_relation))
  
}



############################################# Evaluation metrics
####################### Evaluate the performance of methods
library(mltools)
library(igraph)

########## SHD
# codes from https://github.com/cran/pcalg/blob/master/R/pcalg.R
shd_input_adj_matrix <- function(m1,m2){
  ## Purpose: Compute Structural Hamming Distance between graphs g1 and g2
  ## ----------------------------------------------------------------------
  ## Arguments:
  ## - g1, g2: Input graphs
  ## (graph objects;connectivity matrix where m[x,y]=1 iff x->1
  ## and m[x,y]=m[y,x]=1 iff x-y; pcAlgo-objects)
  ## ----------------------------------------------------------------------
  ## Author: Markus Kalisch, Date:  1 Dec 2006, 17:21
  
  ## Idea: Transform g1 into g2
  ## Transform g1 and g2 into adjacency matrices
  ## m1 and m2 are adjacency matrices
  
  shd <- 0
  ## Remove superfluous edges from g1
  s1 <- m1 + t(m1)
  s2 <- m2 + t(m2)
  s1[s1 == 2] <- 1
  s2[s2 == 2] <- 1
  ds <- s1 - s2
  ind <- which(ds > 0)
  m1[ind] <- 0
  shd <- shd + length(ind)/2
  ## Add missing edges to g1
  ind <- which(ds < 0)
  m1[ind] <- m2[ind]
  shd <- shd + length(ind)/2
  ## Compare Orientation
  d <- abs(m1-m2)
  ## return
  shd + sum((d + t(d)) > 0)/2
}

# report evaluation metrics 
output_metrics <- function(true_adjacency_matrix,est_adjacency_matrix){
  
  p = dim(true_adjacency_matrix)[1]
  
  true_adjacency_matrix_no_diag = true_adjacency_matrix[-seq(1,p^2,p+1)]
  true_edge = which(true_adjacency_matrix_no_diag != 0)
  
  est_adjacency_matrix_no_diag = est_adjacency_matrix[-seq(1,p^2,p+1)]
  pred_edge = which(est_adjacency_matrix_no_diag != 0)
  
  false_edge = setdiff(1:(p^2-p),true_edge)
  pred_null_edge = setdiff(1:(p^2-p),pred_edge)
  
  TPR = length(intersect(pred_edge,true_edge))/length(true_edge)
  
  FPR = (length(intersect(pred_edge,false_edge)))/(length(false_edge))
  
  retrieved <- length(pred_edge)
  
  precision <- length(intersect(pred_edge,true_edge)) / retrieved
  
  recall <- length(intersect(pred_edge,true_edge)) / length(true_edge)
  
  F_score <- 2 * precision * recall / (precision + recall)
  
  TP = length(intersect(pred_edge,true_edge))
  FP = length(intersect(pred_edge,false_edge))
  FN = length(intersect(pred_null_edge,true_edge))
  TN = length(intersect(pred_null_edge,false_edge))
  
  FDR <-  FP / (FP + TP)
  FPR <-  FP / (FP + TN)
  
  # MCC
  # https://www.rdocumentation.org/packages/mltools/versions/0.3.5/topics/mcc
  MCC <- mcc(est_adjacency_matrix_no_diag,true_adjacency_matrix_no_diag)
  MCC <- mcc(TP=TP, FP=FP, TN=TN, FN=FN)   # same
  
  ### SHD: Structrual Hamming Distance
  SHD = shd_input_adj_matrix(true_adjacency_matrix,est_adjacency_matrix)

  return(list(FPR = FPR, FDR = FDR, F_score = F_score, MCC = MCC, SHD = SHD))
  
}

# include diagonals, used for measuring the accuracy of estimating W matrix
output_metrics_with_diag <- function(true_adjacency_matrix,est_adjacency_matrix){
  
  p = dim(true_adjacency_matrix)[1]
  
  true_edge = which(true_adjacency_matrix != 0)
  
  pred_edge = which(est_adjacency_matrix != 0)
  
  false_edge = setdiff(1:p^2,true_edge)
  pred_null_edge = setdiff(1:p^2,pred_edge)
  
  TPR = length(intersect(pred_edge,true_edge))/length(true_edge)
  
  FPR = (length(intersect(pred_edge,false_edge)))/(length(false_edge))
  
  retrieved <- length(pred_edge)
  
  precision <- length(intersect(pred_edge,true_edge)) / retrieved
  
  recall <- length(intersect(pred_edge,true_edge)) / length(true_edge)
  
  F_score <- 2 * precision * recall / (precision + recall)
   
  TP = length(intersect(pred_edge,true_edge))
  FP = length(intersect(pred_edge,false_edge))
  FN = length(intersect(pred_null_edge,true_edge))
  TN = length(intersect(pred_null_edge,false_edge))
  
  FDR <-  FP / (FP + TP)
  FPR <-  FP / (FP + TN)
  
  # MCC
  # https://www.rdocumentation.org/packages/mltools/versions/0.3.5/topics/mcc
  MCC <- mcc(est_adjacency_matrix,true_adjacency_matrix)
  MCC <- mcc(TP=TP, FP=FP, TN=TN, FN=FN)   # same
  
  ### SHD: Structrual Hamming Distance
  SHD = shd_input_adj_matrix(true_adjacency_matrix,est_adjacency_matrix)

  return(list(FPR = FPR, FDR = FDR, F_score = F_score, MCC = MCC, SHD = SHD))
  
}







############################################# DAG estimating functions
###########################################################################################################################
###########################################################################################################################

### First step: Identify ancestral relationships
GAMPI_first_stage.EBIC <- function(X,Y,lambdas=NULL,K_list = 1:(ncol(X)-1),taus = c(0.5,1,2),intercept=TRUE,nlambda=10,lambda.min=0.01,gamma = .5,family="gaussian"){
  
  ### Fit glmtlp
  n = dim(Y)[1]
  p = dim(Y)[2]
  q = dim(X)[2]
  V = matrix(rep(0,q*p),ncol=p)
  family_list = family
  
  for (s in 1:p){
    y = Y[,s]
    if (length(family_list) > 1){
      family = family_list[s]
    }
    
    ### Select tuning parameters
    ebic.fit <- BIC.unconstrained_glmtlp(X,y,lambdas=lambdas,K_list = K_list,taus = taus,intercept = intercept,nlambda=nlambda,lambda.min=lambda.min,gamma = gamma,family=family)
    lambda.optimal = ebic.fit$optpar.lambda
    K.optimal = ebic.fit$optpar.K
    tau.optimal = ebic.fit$optpar.tau
    
    ### Fit model with optimal parameters
    beta_hat = glmtlp(X,y,lambda.optimal,tau = tau.optimal,intercept = intercept,family=family)$beta
    V[,s] = trun_K_glm(X,y,beta_hat,K_list = K.optimal,intercept = intercept,family=family)$beta_matrix
    ### check if V is all 0
    if (sum(abs(V[,s])) == 0){
      beta_hat = glmtlp(X,y,0,tau = tau.optimal,intercept = intercept,family=family)$beta
      V[,s] = trun_K_glm(X,y,beta_hat,K_list = K.optimal,intercept = intercept,family=family)$beta_matrix
    }
  }
  
  ### scale columns of V
  max_each_column = apply(abs(V),2,max)
  col_ind = which(max_each_column > 10)
  scaling_factor =  median(max_each_column)
  for (m in col_ind){
    V[,m] = V[,m] / max_each_column[m] * scaling_factor
  }

  ### Peeling algorithm
  result = peeling(V)
  
  return(list(result=result, V_est = V))
  
}

GAMPI_first_stage.cv <- function(X,Y,lambdas=NULL,fold = 5,type.measure = "deviance",K_list = 1:(ncol(X)-1),taus = c(0.5,1,2),intercept=TRUE,nlambda=10,lambda.min=0.01,cv.rule="1se",family="gaussian"){
  
  ### Fit glmtlp
  n = dim(Y)[1]
  p = dim(Y)[2]
  q = dim(X)[2]
  V = matrix(rep(0,q*p),ncol=p)
  family_list = family
  
  for (s in 1:p){
    y = Y[,s]
    if (length(family_list) > 1){
      family = family_list[s]
    }
    
    ### Select tuning parameters
    cv.fit <- cv.unconstrained_glmtlp_K_tau(X,y,lambdas=lambdas,fold = fold,type.measure=type.measure,K_list = K_list,taus = taus,intercept = intercept,nlambda=nlambda,lambda.min=lambda.min,family=family)
    lambda.optimal = cv.fit$optpar.lambda
    if (cv.rule == "1se"){
      K.optimal = cv.fit$optpar.K.1se
    } else if (cv.rule == "min"){
      K.optimal = cv.fit$optpar.K.min
    }
    tau.optimal = cv.fit$optpar.tau

    ### Fit model with optimal parameters
    beta_hat = glmtlp(X,y,lambda.optimal,tau = tau.optimal,intercept = intercept,family=family)$beta
    V[,s] = trun_K_glm(X,y,beta_hat,K_list = K.optimal,intercept = intercept,family=family)$beta_matrix
    ### check if V is all 0
    if (sum(abs(V[,s])) == 0){
      beta_hat = glmtlp(X,y,0,tau = tau.optimal,intercept = intercept,family=family)$beta
      V[,s] = trun_K_glm(X,y,beta_hat,K_list = K.optimal,intercept = intercept,family=family)$beta_matrix
    }
  }
  
  ### scale columns of V
  max_each_column = apply(abs(V),2,max)
  col_ind = which(max_each_column > 10)
  scaling_factor =  median(max_each_column)
  for (m in col_ind){
    V[,m] = V[,m] / max_each_column[m] * scaling_factor
  }

  ### Peeling algorithm
  result = peeling(V)
  
  return(list(result=result, V_est = V))
  
}


### Second step: Select parent-child relations from ancestors
### No adjustment for confounders with EBIC
select_ancestors.EBIC <- function(X,Y,result,intercept=TRUE,gamma = .5,taus = c(0.5,1,2),family="gaussian"){
  p = ncol(Y)
  family_list = family
  
  ancestral_relation = result$ances_relation
  interv_relation = result$interv_relation
  
  parent_relation_updated = matrix(rep(0,p*p),ncol=p) # parent by child relationship
  
  for (s in 1:p){
    
    if (length(family_list) > 1){
      family = family_list[s]
    }
    
    ### Construct data frame for GLMTLP regression
    Y_anc = Y[,ancestral_relation[,s]!=0]; X_int = X[,interv_relation[,s] != 0];
    X_subset = cbind(Y_anc,X_int)
    y = Y[,s]
    p_anc = sum(ancestral_relation[,s]!=0); p_int = sum(interv_relation[,s] != 0)
    anc_indice = which(ancestral_relation[,s]!=0); int_indice = which(interv_relation[,s] != 0)
    
    is_penalize = c(rep(1,p_anc),rep(0,p_int))
    
    if (p_anc > 0){
      ### Select tuning parameters
      K_list = (p_int): (p_anc + p_int)
      BIC.fit <- BIC.unconstrained_glmtlp_weighted(X_subset,y,K_list = K_list,taus = taus,is_penalize = is_penalize,intercept = intercept,gamma = gamma,family = family)
      lambda.optimal = BIC.fit$optpar.lambda
      K.optimal = BIC.fit$optpar.K
      tau.optimal = BIC.fit$optpar.tau
      
      ### Fit model with optimal parameters
      beta_hat = glmtlp_weighted(X_subset,y,lambda.optimal,tau = tau.optimal,is_penalize = is_penalize,intercept = intercept,family = family)$beta
      beta_glm = trun_K_glm_weighted(X_subset,y,beta_hat,K_list = K.optimal,intercept = intercept,is_penalize = is_penalize,family = family)$beta_matrix
      beta_anc = beta_glm[1:p_anc]; beta_int = beta_glm[(p_anc+1):(p_anc+p_int)]
      
      anc_indice_select = anc_indice[which(beta_anc!=0)]
      parent_relation_updated[anc_indice,s] = beta_anc
    }
  }
  return(parent_relation_updated = parent_relation_updated)
}

### Second step: Select parent-child relations from ancestors
### No adjustment for confounders with cross validation
select_ancestors.cv <- function(X,Y,result,intercept=TRUE,cv.rule="1se",taus = c(0.5,1,2),family="gaussian"){
  p = ncol(Y)
  family_list = family
  
  ancestral_relation = result$ances_relation
  interv_relation = result$interv_relation
  
  parent_relation_updated = matrix(rep(0,p*p),ncol=p) # parent by child relationship
  
  for (s in 1:p){
    
    if (length(family_list) > 1){
      family = family_list[s]
    }
    
    ### Construct data frame for GLMTLP regression
    Y_anc = Y[,ancestral_relation[,s]!=0]; X_int = X[,interv_relation[,s] != 0];
    X_subset = cbind(Y_anc,X_int)
    y = Y[,s]
    p_anc = sum(ancestral_relation[,s]!=0); p_int = sum(interv_relation[,s] != 0)
    anc_indice = which(ancestral_relation[,s]!=0); int_indice = which(interv_relation[,s] != 0)
    
    is_penalize = c(rep(1,p_anc),rep(0,p_int))
    
    if (p_anc > 0){
      ### Select tuning parameters
      K_list = (p_int): (p_anc + p_int)
      cv.fit <- cv.unconstrained_glmtlp_K_tau_weighted(X_subset,y,fold = 5,K_list = K_list,taus = taus,is_penalize = is_penalize,intercept = intercept,family = family)
      lambda.optimal = cv.fit$optpar.lambda
      if (cv.rule == "1se"){
        K.optimal = cv.fit$optpar.K.1se
      } else if (cv.rule == "min"){
        K.optimal = cv.fit$optpar.K.min
      }
      tau.optimal = cv.fit$optpar.tau
      
      ### Fit model with optimal parameters
      beta_hat = glmtlp_weighted(X_subset,y,lambda.optimal,tau = tau.optimal,is_penalize = is_penalize,intercept = intercept,family = family)$beta
      beta_glm = trun_K_glm_weighted(X_subset,y,beta_hat,K_list = K.optimal,intercept = intercept,is_penalize = is_penalize,family = family)$beta_matrix
      beta_anc = beta_glm[1:p_anc]; beta_int = beta_glm[(p_anc+1):(p_anc+p_int)]
      
      anc_indice_select = anc_indice[which(beta_anc!=0)]
      parent_relation_updated[anc_indice,s] = beta_anc
    }
  }
  return(parent_relation_updated = parent_relation_updated)
}



GAMPI_no_deconf.EBIC <- function(X,Y,lambdas=NULL,K_list = 1:(ncol(X)-1),taus = c(0.5,1,2),intercept=TRUE,nlambda=10,lambda.min=0.01,gamma=c(0.5,0.5),family="gaussian"){
  # first stage
  result_first_stage = GAMPI_first_stage.EBIC(X,Y,lambdas=lambdas,K_list = K_list,taus = taus,intercept=intercept,nlambda=nlambda,lambda.min=lambda.min,gamma = gamma[1],family=family)
  # second stage
  result_second_stage = select_ancestors.EBIC(X,Y,result_first_stage$result,intercept=intercept,gamma = gamma[2],taus=taus,family=family)
  return(list(causal_relation = result_second_stage, interv_relation = result_first_stage$result$interv_relation, V_est = result_first_stage$V))
}

GAMPI_no_deconf.cv <- function(X,Y,lambdas=NULL,fold = 5,type.measure = "deviance",K_list = 1:(ncol(X)-1),taus = c(0.5,1,2),intercept=TRUE,nlambda=10,lambda.min=0.01,cv.rule.first_stage="1se",cv.rule.second_stage="1se",family="gaussian"){
  # first stage
  result_first_stage = GAMPI_first_stage.cv(X,Y,lambdas=lambdas,fold = fold,type.measure=type.measure,K_list = K_list,taus = taus,intercept=intercept,nlambda=nlambda,lambda.min=lambda.min,cv.rule=cv.rule.first_stage,family=family)
  # second stage
  result_second_stage = select_ancestors.cv(X,Y,result_first_stage$result,intercept=intercept,cv.rule=cv.rule.second_stage,taus=taus,family=family)
  return(list(causal_relation = result_second_stage, interv_relation = result_first_stage$result$interv_relation, V_est = result_first_stage$V))
}

GAMPI_no_deconf <- function(X,Y,lambdas=NULL,K_list = 1:(ncol(X)-1),taus = c(0.5,1,2),intercept=TRUE,nlambda=10,lambda.min=0.01,method="EBIC",cv.rule.first_stage="1se",cv.rule.second_stage="1se",ebic.gamma=c(0.5,0.5),family="gaussian"){
  
  K_list = 1:min(10,K_list[length(K_list)]) ### save computation
  
  if (method == "CV"){
    result = GAMPI_no_deconf.cv(X,Y,lambdas=lambdas,K_list = K_list,taus = taus,intercept=intercept,nlambda=nlambda,lambda.min=lambda.min,cv.rule.first_stage=cv.rule.first_stage,cv.rule.second_stage=cv.rule.second_stage,family=family)
  } else if (method == "EBIC"){
    result = GAMPI_no_deconf.EBIC(X,Y,lambdas=lambdas,K_list = K_list,taus = taus,intercept=intercept,nlambda=nlambda,lambda.min=lambda.min,gamma=ebic.gamma,family=family)
    if ( ncol(Y) < 100 && sum(result$causal_relation!=0) < 0.1 * ncol(Y)){
      ebic.gamma = c(0,0)
      result = GAMPI_no_deconf.EBIC(X,Y,lambdas=lambdas,K_list = K_list,taus = taus,intercept=intercept,nlambda=nlambda,lambda.min=lambda.min,gamma=ebic.gamma,family=family)
      print("Use BIC criterion")
    }
  }
  return(result)
}



###########################################################################################################################
###########################################################################################################################
### Second step: Select parent-child relations from ancestors
### adjusting for confounders using residual inclusion with EBIC
select_ancestors.EBIC_DRI <- function(X,Y,result,intercept=TRUE,gamma = .5,taus = c(0.5,1,2),family="gaussian"){
  n = nrow(X); p = ncol(Y)
  family_list = family

  ancestral_relation = result$ances_relation
  interv_relation = result$interv_relation
  
  parent_relation_deconfound = matrix(rep(0,p*p),ncol=p) # parent by child relationship
  
  g <- graph_from_adjacency_matrix(abs(ancestral_relation), mode="directed")
  ### Get topological order
  topo_order <- topo_sort(g)
  
  confounding_effect = matrix(rep(0,n*p),nrow = n, ncol = p)  # estimated confounding effect, n by p
  
  ### start peeling
  for (j in topo_order){
    
    if (length(family_list) > 1){
      family = family_list[j]
    }
    
    anc_indice = which(ancestral_relation[,j]!=0); int_indice = which(interv_relation[,j] != 0);
    Y_anc = Y[,anc_indice]; X_int = X[,int_indice]; confounders = confounding_effect[,anc_indice];

    ### Construct data frame for GLMTLP regression
    X_subset = cbind(Y_anc,X_int,confounders)
    y = Y[,j]
    p_anc = sum(ancestral_relation[,j]!=0); p_int = sum(interv_relation[,j] != 0); p_conf = sum(ancestral_relation[,j]!=0);
    is_penalize = c(rep(1,p_anc),rep(0,p_int),rep(1,p_conf))
    
    if (p_anc == 0 & p_int>0){ ### root node
      d = data.frame(y,Y_anc,X_int,confounders)
      if (intercept){
        fit = glm(y ~ X_int, data=d,family = family)
        confounding_effect[,j] = Y[,j] - fit$fitted.values
      } else{
        fit = glm(y ~ X_int - 1, data=d,family = family)
        confounding_effect[,j] = Y[,j] - fit$fitted.values
      }
    } else if (p_anc > 0){
      ### Select tuning parameters
      K_list = (p_int): (p_anc + p_int + p_conf)
      BIC.fit <- BIC.unconstrained_glmtlp_weighted(X_subset,y,K_list = K_list,taus = taus,is_penalize = is_penalize,intercept = intercept,gamma = gamma,family = family)
      lambda.optimal = BIC.fit$optpar.lambda
      K.optimal = BIC.fit$optpar.K
      tau.optimal = BIC.fit$optpar.tau
      
      ### Fit model with optimal parameters
      beta_hat = glmtlp_weighted(X_subset,y,lambda.optimal,tau = tau.optimal,is_penalize=is_penalize,intercept = intercept,family = family)$beta
      result_beta_glm = trun_K_glm_weighted_include_confounder(X_subset,y,beta_hat,K_list = K.optimal,intercept = intercept,is_penalize = is_penalize,p_anc = p_anc,p_int = p_int,p_conf = p_conf,family = family)
      beta_glm = result_beta_glm$beta_matrix
      beta_anc = beta_glm[1:p_anc]; beta_int = beta_glm[(p_anc+1):(p_anc+p_int)]
      
      parent_relation_deconfound[anc_indice,j] = beta_anc;
      
      ### Estimate residuals
      confounding_effect[,j] = Y[,j] - result_beta_glm$fitted_values
    }
  }
  return(parent_relation_updated = parent_relation_deconfound)
}

### Second step: Select parent-child relations from ancestors
### adjusting for confounders using residual inclusion with cross validation
select_ancestors.cv_DRI <- function(X,Y,result,intercept=TRUE,cv.rule="1se",taus = c(0.5,1,2),family="gaussian"){
  n = nrow(X); p = ncol(Y)
  family_list = family

  ancestral_relation = result$ances_relation
  interv_relation = result$interv_relation
  
  parent_relation_deconfound = matrix(rep(0,p*p),ncol=p) # parent by child relationship
  
  g <- graph_from_adjacency_matrix(abs(ancestral_relation), mode="directed")
  ### Get topological order
  topo_order <- topo_sort(g)
  
  confounding_effect = matrix(rep(0,n*p),nrow = n, ncol = p)  # estimated confounding effect, n by p
  
  ### start peeling
  for (j in topo_order){
    
    if (length(family_list) > 1){
      family = family_list[j]
    }
    
    anc_indice = which(ancestral_relation[,j]!=0); int_indice = which(interv_relation[,j] != 0);
    Y_anc = Y[,anc_indice]; X_int = X[,int_indice]; confounders = confounding_effect[,anc_indice];

    ### Construct data frame for GLMTLP regression
    X_subset = cbind(Y_anc,X_int,confounders)
    y = Y[,j]
    p_anc = sum(ancestral_relation[,j]!=0); p_int = sum(interv_relation[,j] != 0); p_conf = sum(ancestral_relation[,j]!=0);
    is_penalize = c(rep(1,p_anc),rep(0,p_int),rep(1,p_conf))
    
    if (p_anc == 0 & p_int>0){ ### root node
      d = data.frame(y,Y_anc,X_int,confounders)
      if (intercept){
        fit = glm(y ~ X_int, data=d,family = family)
        confounding_effect[,j] = Y[,j] - fit$fitted.values
      } else{
        fit = glm(y ~ X_int - 1, data=d,family = family)
        confounding_effect[,j] = Y[,j] - fit$fitted.values
      }
    } else if (p_anc > 0){
      ### Select tuning parameters
      K_list = (p_int): (p_anc + p_int + p_conf)
      cv.fit <- cv.unconstrained_glmtlp_K_tau_weighted(X_subset,y,fold = 5,K_list = K_list,taus = taus,is_penalize = is_penalize,intercept = intercept,family = family)
      lambda.optimal = cv.fit$optpar.lambda
      if (cv.rule == "1se"){
        K.optimal = cv.fit$optpar.K.1se
      } else if (cv.rule == "min"){
        K.optimal = cv.fit$optpar.K.min
      }
      tau.optimal = cv.fit$optpar.tau
      
      ### Fit model with optimal parameters
      beta_hat = glmtlp_weighted(X_subset,y,lambda.optimal,tau = tau.optimal,is_penalize=is_penalize,intercept = intercept,family = family)$beta
      result_beta_glm = trun_K_glm_weighted_include_confounder(X_subset,y,beta_hat,K_list = K.optimal,intercept = intercept,is_penalize = is_penalize,p_anc = p_anc,p_int = p_int,p_conf = p_conf,family = family)
      beta_glm = result_beta_glm$beta_matrix
      beta_anc = beta_glm[1:p_anc]; beta_int = beta_glm[(p_anc+1):(p_anc+p_int)]
      
      parent_relation_deconfound[anc_indice,j] = beta_anc;
      
      ### Estimate residuals
      confounding_effect[,j] = Y[,j] - result_beta_glm$fitted_values
    }
  }
  return(parent_relation_updated = parent_relation_deconfound)
}

### Truncated GLM with confounders
trun_K_glm_weighted_include_confounder <- function(X,y,beta_tlp,K_list = 1:ncol(X),intercept=TRUE,is_penalize = rep(1,ncol(X)),p_anc,p_int,p_conf,family="gaussian"){
  
  n = nrow(X)
  q = ncol(X)
  beta_tlp[is_penalize==0] = max(abs(beta_tlp)) + 1 # make unpenalized coefs the largest ones so that they are always included in the model

  if (!intercept){
    beta_matrix = matrix(rep(0,q*length(K_list)), ncol = length(K_list))
    fitted_values = matrix(rep(0,n*length(K_list)), ncol = length(K_list))
    
    if (sum(abs(beta_tlp)) != 0){
      for (s in 1:length(K_list)){
        K = K_list[s]
        if (K<q){
          indice = order(abs(beta_tlp), decreasing = T)[K+1]
          ### find nonzero index
          act_set = which(abs(beta_tlp) > abs(beta_tlp[indice]))
          if (length(act_set)==0){
            act_set = order(abs(beta_tlp), decreasing = T)[1:K]
          }
        } else if (K==q){
          act_set = which(beta_tlp != 0)
        }
        
        act_set_Y = intersect(act_set,1:p_anc) # parents that are selected
        act_set_conf = p_anc + p_int + act_set_Y # confounders associated with the selected parents
        act_set = union(act_set,act_set_conf)
        
        ### refit
        fit = glm(y~X[,act_set]-1,family = family)
        beta_refit = fit$coefficients
        beta = rep(0,q)
        beta[act_set] = beta_refit
        beta_matrix[,s] = beta
        fitted_values[,s] = fit$fitted.values
      } 
    }
    return(list(beta_matrix=beta_matrix,fitted_values=fitted_values))
  } else { # Fit with intercept
    
    beta_matrix = matrix(rep(0,q*length(K_list)), ncol = length(K_list))
    beta0_vec = rep(0,length(K_list))
    fitted_values = matrix(rep(0,n*length(K_list)), ncol = length(K_list))
    
    if (sum(abs(beta_tlp)) != 0){
      for (s in 1:length(K_list)){
        K = K_list[s]
        if (K<q){
          indice = order(abs(beta_tlp), decreasing = T)[K+1]
          ### find nonzero index
          act_set = which(abs(beta_tlp) > abs(beta_tlp[indice]))
          if (length(act_set)==0){
            act_set = order(abs(beta_tlp), decreasing = T)[1:K]
          }
        } else if (K==q){
          act_set = which(beta_tlp != 0)
        }
        
        act_set_Y = intersect(act_set,1:p_anc) # parents that are selected
        act_set_conf = p_anc + p_int + act_set_Y # confounders associated with the selected parents
        act_set = union(act_set,act_set_conf)
        
        ### refit
        fit = glm(y~X[,act_set],family = family)
        beta_refit = fit$coefficients[-1]
        beta = rep(0,q)
        beta[act_set] = beta_refit
        beta_matrix[,s] = beta
        beta0_vec[s] = fit$coefficients[1]
        fitted_values[,s] = fit$fitted.values
      } 
    }
    return(list(beta_matrix=beta_matrix,beta0_vec=beta0_vec,fitted_values=fitted_values))
  }
}


### Second step: Select parent-child relations from ancestors
### adjusting for confounders using predictor substitution with EBIC
select_ancestors.EBIC_DPS <- function(X,Y,result,intercept=TRUE,gamma = .5,taus = c(0.5,1,2),family="gaussian"){
  n = nrow(X); p = ncol(Y)
  family_list = family

  ancestral_relation = result$ances_relation
  interv_relation = result$interv_relation
  
  parent_relation_deconfound = matrix(rep(0,p*p),ncol=p) # parent by child relationship
  
  g <- graph_from_adjacency_matrix(abs(ancestral_relation), mode="directed")
  ### Get topological order
  topo_order <- topo_sort(g)
  
  Y_imputed = matrix(rep(0,n*p),nrow = n, ncol = p)  # imputed value, n by p
  
  ### start peeling
  for (j in topo_order){
    
    if (length(family_list) > 1){
      family = family_list[j]
    }
    
    anc_indice = which(ancestral_relation[,j]!=0); int_indice = which(interv_relation[,j] != 0);
    Y_anc = Y_imputed[,anc_indice]; X_int = X[,int_indice];
    
    ### Construct data frame for GLMTLP regression
    X_subset = cbind(Y_anc,X_int)
    y = Y[,j]
    p_anc = sum(ancestral_relation[,j]!=0); p_int = sum(interv_relation[,j] != 0);
    is_penalize = c(rep(1,p_anc),rep(0,p_int))
    
    if (p_anc == 0 & p_int>0){ ### root node
      d = data.frame(y,Y_anc,X_int)
      if (intercept){
        fit = glm(y ~ X_int, data=d,family = family)
        Y_imputed[,j] = fit$fitted.values
      } else{
        fit = glm(y ~ X_int - 1, data=d,family = family)
        Y_imputed[,j] = fit$fitted.values
      }
    } else if (p_anc > 0){
      ### Select tuning parameters
      K_list = (p_int): (p_anc + p_int)
      BIC.fit <- BIC.unconstrained_glmtlp_weighted(X_subset,y,K_list = K_list,taus = taus,is_penalize = is_penalize,intercept = intercept,gamma = gamma,family = family)
      lambda.optimal = BIC.fit$optpar.lambda
      K.optimal = BIC.fit$optpar.K
      tau.optimal = BIC.fit$optpar.tau
      
      ### Fit model with optimal parameters
      beta_hat = glmtlp_weighted(X_subset,y,lambda.optimal,tau = tau.optimal,is_penalize = is_penalize,intercept = intercept,family = family)$beta
      result_beta_glm = trun_K_glm_weighted(X_subset,y,beta_hat,K_list = K.optimal,intercept = intercept,is_penalize = is_penalize,family = family)
      beta_glm = result_beta_glm$beta_matrix
      beta_anc = beta_glm[1:p_anc]; beta_int = beta_glm[(p_anc+1):(p_anc+p_int)]
      
      parent_relation_deconfound[anc_indice,j] = beta_anc;
      
      ### Impute Y
      Y_imputed[,j] = result_beta_glm$fitted_values
    }
  }
  return(parent_relation_updated = parent_relation_deconfound)
}

### Second step: Select parent-child relations from ancestors
### adjusting for confounders using predictor substitution with cross validation
select_ancestors.cv_DPS <- function(X,Y,result,intercept=TRUE,cv.rule="1se",taus = c(0.5,1,2),family="gaussian"){
  n = nrow(X); p = ncol(Y)
  family_list = family

  ancestral_relation = result$ances_relation
  interv_relation = result$interv_relation
  
  parent_relation_deconfound = matrix(rep(0,p*p),ncol=p) # parent by child relationship
  
  g <- graph_from_adjacency_matrix(abs(ancestral_relation), mode="directed")
  ### Get topological order
  topo_order <- topo_sort(g)
  
  Y_imputed = matrix(rep(0,n*p),nrow = n, ncol = p)  # imputed value, n by p
  
  ### start peeling
  for (j in topo_order){
    
    if (length(family_list) > 1){
      family = family_list[j]
    }
    
    anc_indice = which(ancestral_relation[,j]!=0); int_indice = which(interv_relation[,j] != 0);
    Y_anc = Y_imputed[,anc_indice]; X_int = X[,int_indice];
    
    ### Construct data frame for GLMTLP regression
    X_subset = cbind(Y_anc,X_int)
    y = Y[,j]
    p_anc = sum(ancestral_relation[,j]!=0); p_int = sum(interv_relation[,j] != 0);
    is_penalize = c(rep(1,p_anc),rep(0,p_int))
    
    if (p_anc == 0 & p_int>0){ ### root node
      d = data.frame(y,Y_anc,X_int)
      if (intercept){
        fit = glm(y ~ X_int, data=d,family = family)
        Y_imputed[,j] = fit$fitted.values
      } else{
        fit = glm(y ~ X_int - 1, data=d,family = family)
        Y_imputed[,j] = fit$fitted.values
      }
    } else if (p_anc > 0){
      ### Select tuning parameters
      K_list = (p_int): (p_anc + p_int)
      cv.fit <- cv.unconstrained_glmtlp_K_tau_weighted(X_subset,y,fold = 5,K_list = K_list,taus = taus,is_penalize = is_penalize,intercept = intercept,family = family)
      lambda.optimal = cv.fit$optpar.lambda
      if (cv.rule == "1se"){
        K.optimal = cv.fit$optpar.K.1se
      } else if (cv.rule == "min"){
        K.optimal = cv.fit$optpar.K.min
      }
      tau.optimal = cv.fit$optpar.tau
      
      ### Fit model with optimal parameters
      beta_hat = glmtlp_weighted(X_subset,y,lambda.optimal,tau = tau.optimal,is_penalize = is_penalize,intercept = intercept,family = family)$beta
      result_beta_glm = trun_K_glm_weighted(X_subset,y,beta_hat,K_list = K.optimal,intercept = intercept,is_penalize = is_penalize,family = family)
      beta_glm = result_beta_glm$beta_matrix
      beta_anc = beta_glm[1:p_anc]; beta_int = beta_glm[(p_anc+1):(p_anc+p_int)]
      
      parent_relation_deconfound[anc_indice,j] = beta_anc;
      
      ### Impute Y
      Y_imputed[,j] = result_beta_glm$fitted_values
    }
  }
  return(parent_relation_updated = parent_relation_deconfound)
}





### Adjust for confounders: residual inclusion with EBIC
GAMPI.EBIC_DRI <- function(X,Y,lambdas=NULL,K_list = 1:(ncol(X)-1),taus = c(0.5,1,2),intercept=TRUE,nlambda=10,lambda.min=0.01,gamma=c(0.5,0.5),family="gaussian"){
  # first stage
  result_first_stage = GAMPI_first_stage.EBIC(X,Y,lambdas=lambdas,K_list = K_list,taus = taus,intercept=intercept,nlambda=nlambda,lambda.min=lambda.min,gamma = gamma[1],family=family)
  # second stage
  result_second_stage = select_ancestors.EBIC_DRI(X,Y,result_first_stage$result,intercept=intercept,gamma = gamma[2],taus=taus,family=family)
  return(list(causal_relation = result_second_stage, interv_relation = result_first_stage$result$interv_relation, V_est = result_first_stage$V))
}

### Adjust for confounders: predictor substitution with EBIC
GAMPI.EBIC_DPS <- function(X,Y,lambdas=NULL,K_list = 1:(ncol(X)-1),taus = c(0.5,1,2),intercept=TRUE,nlambda=10,lambda.min=0.01,gamma=c(0.5,0.5),family="gaussian"){
  # first stage
  result_first_stage = GAMPI_first_stage.EBIC(X,Y,lambdas=lambdas,K_list = K_list,taus = taus,intercept=intercept,nlambda=nlambda,lambda.min=lambda.min,gamma = gamma[1],family=family)
  # second stage
  result_second_stage = select_ancestors.EBIC_DPS(X,Y,result_first_stage$result,intercept=intercept,gamma = gamma[2],taus=taus,family=family)
  return(list(causal_relation = result_second_stage, interv_relation = result_first_stage$result$interv_relation, V_est = result_first_stage$V))
}

### Adjust for confounders: residual inclusion with cross validation
GAMPI.cv_DRI <- function(X,Y,lambdas=NULL,fold = 5,type.measure = "deviance",K_list = 1:(ncol(X)-1),taus = c(0.5,1,2),intercept=TRUE,nlambda=10,lambda.min=0.01,cv.rule.first_stage="1se",cv.rule.second_stage="1se",family="gaussian"){
  # first stage
  result_first_stage = GAMPI_first_stage.cv(X,Y,lambdas=lambdas,fold = fold,type.measure=type.measure,K_list = K_list,taus = taus,intercept=intercept,nlambda=nlambda,lambda.min=lambda.min,cv.rule=cv.rule.first_stage,family=family)
  # second stage
  result_second_stage = select_ancestors.cv_DRI(X,Y,result_first_stage$result,intercept=intercept,cv.rule=cv.rule.second_stage,taus=taus,family=family)
  return(list(causal_relation = result_second_stage, interv_relation = result_first_stage$result$interv_relation, V_est = result_first_stage$V))
}

### Adjust for confounders: predictor substitution with cross validation
GAMPI.cv_DPS <- function(X,Y,lambdas=NULL,fold = 5,type.measure = "deviance",K_list = 1:(ncol(X)-1),taus = c(0.5,1,2),intercept=TRUE,nlambda=10,lambda.min=0.01,cv.rule.first_stage="1se",cv.rule.second_stage="1se",family="gaussian"){
  # first stage
  result_first_stage = GAMPI_first_stage.cv(X,Y,lambdas=lambdas,fold = fold,type.measure=type.measure,K_list = K_list,taus = taus,intercept=intercept,nlambda=nlambda,lambda.min=lambda.min,cv.rule=cv.rule.first_stage,family=family)
  # second stage
  result_second_stage = select_ancestors.cv_DPS(X,Y,result_first_stage$result,intercept=intercept,cv.rule=cv.rule.second_stage,taus=taus,family=family)
  return(list(causal_relation = result_second_stage, interv_relation = result_first_stage$result$interv_relation, V_est = result_first_stage$V))
}





### Main functions
GAMPI_DRI <- function(X,Y,lambdas=NULL,K_list = 1:(ncol(X)-1),taus = c(0.5,1,2),intercept=TRUE,nlambda=10,lambda.min=0.01,method="EBIC",cv.rule.first_stage="1se",cv.rule.second_stage="1se",ebic.gamma=c(0.5,0.5),family="gaussian"){
  
  K_list = 1:min(10,K_list[length(K_list)]) ### save computation
  
  if (method == "CV"){
    result = GAMPI.cv_DRI(X,Y,lambdas=lambdas,K_list = K_list,taus = taus,intercept=intercept,nlambda=nlambda,lambda.min=lambda.min,cv.rule.first_stage=cv.rule.first_stage,cv.rule.second_stage=cv.rule.second_stage,family=family)
  } else if (method == "EBIC"){
    result = GAMPI.EBIC_DRI(X,Y,lambdas=lambdas,K_list = K_list,taus = taus,intercept=intercept,nlambda=nlambda,lambda.min=lambda.min,gamma=ebic.gamma,family=family)
    if ( ncol(Y) < 100 && sum(result$causal_relation!=0) < 0.1 * ncol(Y)){
      ebic.gamma = c(0,0)
      result = GAMPI.EBIC_DRI(X,Y,lambdas=lambdas,K_list = K_list,taus = taus,intercept=intercept,nlambda=nlambda,lambda.min=lambda.min,gamma=ebic.gamma,family=family)
      print("Use BIC criterion")
    }
  }
  return(result)
}

GAMPI_DPS <- function(X,Y,lambdas=NULL,K_list = 1:(ncol(X)-1),taus = c(0.5,1,2),intercept=TRUE,nlambda=10,lambda.min=0.01,method="EBIC",cv.rule.first_stage="1se",cv.rule.second_stage="1se",ebic.gamma=c(0.5,0.5),family="gaussian"){
  
  K_list = 1:min(10,K_list[length(K_list)]) ### save computation
  
  if (method == "CV"){
    result = GAMPI.cv_DPS(X,Y,lambdas=lambdas,K_list = K_list,taus = taus,intercept=intercept,nlambda=nlambda,lambda.min=lambda.min,cv.rule.first_stage=cv.rule.first_stage,cv.rule.second_stage=cv.rule.second_stage,family=family)
  } else if (method == "EBIC"){
    result = GAMPI.EBIC_DPS(X,Y,lambdas=lambdas,K_list = K_list,taus = taus,intercept=intercept,nlambda=nlambda,lambda.min=lambda.min,gamma=ebic.gamma,family=family)
    if ( ncol(Y) < 100 && sum(result$causal_relation!=0) < 0.1 * ncol(Y)){
      ebic.gamma = c(0,0)
      result = GAMPI.EBIC_DPS(X,Y,lambdas=lambdas,K_list = K_list,taus = taus,intercept=intercept,nlambda=nlambda,lambda.min=lambda.min,gamma=ebic.gamma,family=family)
      print("Use BIC criterion")
    }
  }
  return(result)
}



