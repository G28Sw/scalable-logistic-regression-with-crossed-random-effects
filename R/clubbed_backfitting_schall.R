## Helper Functions for scalable logistic regression with crossed random effect


#------------------------------------------------------------------------------------
## Clear the environment ##
#------------------------------------------------------------------------------------
rm(list = ls())


#------------------------------------------------------------------------------------
## Install the required packages into ##
## R if it's not already installed.   ##
#------------------------------------------------------------------------------------
if(!require(dplyr)) install.packages("dplyr")
if(!require(tibble)) install.packages("tibble")

#------------------------------------------------------------------------------------
## Load the packages into the R environment ##
#------------------------------------------------------------------------------------
library(dplyr)
library(tibble) # For generating from the multivariate normal distribution.



#------------------------------------------------------------------------------------
## Functions for clubbed backfitting ##
#------------------------------------------------------------------------------------

preclub=function(X,f,w,lambda){
  ###Precompute the clubbing operator matrix, same size as X
  ## X nxp matrix, including intercept if appropriate
  ## w n vector of weights
  ## f a factor variable of length n with k classes, stored as an integer
  ## lambda -shrinkage parameter
  wX=w*X
  wX=as_tibble(wX)
  d=tibble(f,w,wX)
  pwsum = d %>%
    group_by(f) %>%
    summarize_all(sum,na.rm=TRUE)
  pwsum=data.matrix(pwsum[,-1])
  denom = drop(pwsum[,1]+lambda)
  Xbar = pwsum[,-1]/denom
  Xad = w*(X-Xbar[f,])
  return(list(operator = solve(t(X)%*%Xad,t(Xad)), denom = denom))
}

club = function(X,f,w,y,lambda,Xad=NULL){
  ## This function will do the clubbing, and produce beta and b
  ## typically y is a residual
  ## if Xad is supplied, it has been precomputed

  if(is.null(Xad))Xad=preclub(X,f,w,lambda)$operator
  beta=drop(Xad%*%y)
  wr=as_tibble(w*(y-drop(X%*%beta)))
  d=tibble(f,w,wr)
  pwsum=d %>%
    group_by(f) %>%
    summarize_all(sum,na.rm=TRUE)
  pwsum=data.matrix(pwsum[,-1])
  b=pwsum[,-1]/drop(pwsum[,1]+lambda)
  return(list(beta=beta,coef=b))
}


expit<-function(t){
  return(1/(1+exp(-t)))
}

backfitting_inner_loop<-function(x,z,f1,f2,weight,weight_a,weight_b,beta_hat,a_hat,b_hat,
                                 lambda_a_hat,lambda_b_hat,max_iter=100,epsilon=1e-8,
                                 Xclub_a,Xclub_b){
  #it solves the penalized least square problem using backfitting
  iter = 0
  relative_change = epsilon+2

  predx = x%*%beta_hat
  preda = a_hat[f1]
  predb = b_hat[f2]

  while((relative_change > epsilon)&(iter<=max_iter)){
    iter = iter+1
    deltaf=0

    predx_old = predx
    preda_old = preda
    predb_old = predb

    betaa=club(x,f1,weight,z-b_hat[f2],lambda_a_hat,Xclub_a)
    a_hat=betaa$coef

    betab=club(x,f2,weight,z-a_hat[f1],lambda_b_hat,Xclub_b)
    beta_hat=betab$beta
    b_hat=betab$coef

    predx = x%*%beta_hat
    preda = a_hat[f1]
    predb = b_hat[f2]

    deltaf=mean((predx_old-predx)^2)+mean((preda_old-preda)^2)+mean((predb_old-predb)^2)
    dd = mean((predx_old+preda_old+predb_old)^2)
    relative_change = deltaf/dd
    # cat("inner loop relative change :",relative_change,"\n")

  }
  v1star = sum(1/(weight_a+lambda_a_hat))*lambda_a_hat
  v2star = sum(1/(weight_b+lambda_b_hat))*lambda_b_hat
  sigma_sq_a_hat = sum(a_hat^2)/(length(a_hat)-v1star)
  sigma_sq_b_hat = sum(b_hat^2)/(length(b_hat)-v2star)
  # cat("sigma_sq_a",sigma_sq_a_hat,"\n")
  # cat("sigma_sq_b",sigma_sq_b_hat,"\n")
  # cat("beta",beta_hat,"\n")


  pred = predx+preda+predb

  return(list(betahat = as.vector(beta_hat), ahat = a_hat, bhat = b_hat,
              sigma_a_sqhat = sigma_sq_a_hat, sigma_b_sqhat = sigma_sq_b_hat,
              inner_iter=iter,fitted_val=pred,v1star=v1star,v2star=v2star))
}


backfitting_outer_loop<-function(x,y,f1,f2,betahat = rep(0,ncol(x)),sigma_a_sqhat=1,sigma_b_sqhat=1,
                                 epsilon=1e-8,max_iter=100,over_dispersed=TRUE,trace.it=TRUE,epsilon_backfitting=1e-8){
  #it estimates the parameter using modified schall algorithm
  iter = 0
  nf1 = length(unique(f1))
  nf2 = length(unique(f2))
  f1 = as.factor(f1)
  f2 = as.factor(f2)
  levels(f1) = 1:nf1
  levels(f2) = 1:nf2
  ahat = rep(0,nf1)
  bhat = rep(0,nf2)
  dispersion=1
  rel_change = epsilon+2
  convergence_code = 1
  fitted_val = as.vector(x%*%betahat + ahat[f1] + bhat[f2])

  inner_iteration = c()
  sigma_a_sqhat_vector = c()
  sigma_b_sqhat_vector = c()
  betahat_inner_loop = list()

  while((rel_change>epsilon) & (iter < max_iter)){
    etahat = as.vector(x%*%betahat + ahat[f1] + bhat[f2])
    phat = expit(etahat)
    diag_p_elem = phat*(1-phat)
    weight = as.vector(diag_p_elem/dispersion)

    weight_a = as.vector(tapply(weight,f1,sum))
    weight_b = as.vector(tapply(weight,f2,sum))

    z = as.vector(etahat + ((y-phat)/diag_p_elem))

    lambda_a_hat = 1/sigma_a_sqhat
    lambda_b_hat = 1/sigma_b_sqhat
    sigma_a_sqhat_vector = c(sigma_a_sqhat_vector,sigma_a_sqhat)
    sigma_b_sqhat_vector = c(sigma_b_sqhat_vector,sigma_b_sqhat)

    start_time_factor_a = Sys.time()

    Xclub_a = preclub(x,f1,weight,lambda_a_hat)$operator

    start_time_factor_b = Sys.time()
    Xclub_b = preclub(x,f2,weight,lambda_b_hat)$operator
    scha = backfitting_inner_loop(x,z,f1,f2,weight,weight_a,weight_b,betahat,ahat,bhat,
                                  lambda_a_hat,lambda_b_hat,max_iter,epsilon_backfitting,Xclub_a,Xclub_b)

    betahat_new = scha$betahat
    ahat_new = scha$ahat
    bhat_new = scha$bhat
    sigma_a_sqhat_new = scha$sigma_a_sqhat
    sigma_b_sqhat_new = scha$sigma_b_sqhat
    inner_iteration = c(inner_iteration,scha$inner_iter)
    fitted_val_new = scha$fitted_val
    betahat_inner_loop[[iter+1]] = betahat_new

    change = sum((fitted_val_new-fitted_val)^2)
    rel_change = change/sum(fitted_val^2)

    if(trace.it) cat("Iter",iter,"Rel Change",rel_change,"\n")

    betahat = betahat_new
    sigma_a_sqhat = sigma_a_sqhat_new
    sigma_b_sqhat = sigma_b_sqhat_new
    ahat = ahat_new
    bhat = bhat_new
    fitted_val = fitted_val_new

    # compute dispersion
    etahat = as.vector(x%*%betahat + ahat[f1] + bhat[f2])
    phat = expit(etahat)
    diag_p_elem = phat*(1-phat)
    residual =  y-phat
    denom=nrow(x)-(nf1-scha$v1star)-(nf2-scha$v2star)
    dispersion = sum((residual^2)/diag_p_elem )/denom
    iter = iter + 1
    if(iter==max_iter){convergence_code = 0}
  }
  return(list(beta=betahat,sigma_a=sqrt(sigma_a_sqhat),sigma_b=sqrt(sigma_b_sqhat),dispersion=dispersion,ahat=ahat,bhat=bhat,
              iteration=iter,convergence_code=convergence_code,inner_iteration=inner_iteration,betahat_each_inner_loop_at_convergence = betahat_inner_loop,
              sigma_a_sqhat_vector = sigma_a_sqhat_vector, sigma_b_sqhat_vector = sigma_b_sqhat_vector))
}

level_means<-function(X,f,w,lambda){
  #This function uses tibble to do tapply on multiple columns of a matrix

  wX = w*X
  d= tibble(f,w,as_tibble(wX))
  pwsum=d %>%
    group_by(f) %>%
    summarize_all(sum,na.rm=TRUE)
  pwsum=data.matrix(pwsum[,-1])
  denom = drop(pwsum[,1]+lambda)
  centering_weight = (1/denom)/sum((1/denom))
  num = pwsum[,-1]
  estimate = scale(num, center= colSums(centering_weight*num),scale=F)/denom
  return(estimate)

}

schall_covariance<-function(x,f1,f2,lambda_a_hat,lambda_b_hat,weight,
                            epsilon,max_iter=100,trace.it=TRUE){
  #helper function for covariance
  #it finds covariance by backfitting on each columns of X

  start_time = Sys.time()

  iter = 0
  nf1 = length(unique(f1))
  nf2 = length(unique(f2))

  a_hat = matrix(rep(0,nf1*ncol(x)),nrow=nf1)
  b_hat = matrix(rep(0,nf2*ncol(x)),nrow=nf2)
  preda = a_hat[f1,]
  predb = b_hat[f2,]
  fitted_val = preda + predb

  rel_change = epsilon+2
  convergence_code = 1

  while((rel_change>epsilon) & (iter < max_iter)){
    iter = iter+1
    deltaf=0

    preda_old = preda
    predb_old = predb

    resida = x - predb
    a_hat = level_means(resida,f1,weight,lambda_a_hat)
    preda = a_hat[f1,]

    residb = x - preda
    b_hat = level_means(residb,f2,weight,lambda_b_hat)
    predb = b_hat[f2,]

    deltaf=mean((preda_old-preda)^2)+mean((predb_old-predb)^2)
    dd = mean((preda_old+predb_old)^2)
    rel_change = deltaf/dd
    if(trace.it){cat("relative change while estimating covariance is : ",rel_change,"\n")}
  }

  list(a_hat=a_hat,b_hat=b_hat,fitted_val=preda+predb,iter = iter)
}




var_betahat_glmm <-function(x,f1,f2,trace.it=TRUE,backfit,epsilon=1e-10,max_iter=100){
  #it estimates the cov(\hat{\beta}_{GLMM}) using backfitting

  etahat =  x%*%backfit$beta + backfit$ahat[f1] + backfit$bhat[f2]
  phat = expit(etahat)
  weight_vector = as.vector(phat*(1-phat))#/backfit$dispersion
  sigma_a_sqhat = backfit$sigma_a^2
  sigma_b_sqhat = backfit$sigma_b^2
  sx = schall_covariance(x,f1,f2,1/sigma_a_sqhat,1/sigma_b_sqhat,
                         weight_vector,epsilon,max_iter,trace.it)

  w_resid_sx = weight_vector*(x - sx$fitted_val)

  xt_w_i_minus_s_x = t(x)%*%w_resid_sx
  xt_w_i_minus_s_x_inv = solve(xt_w_i_minus_s_x)

  df1 = tibble(f1,as_tibble(w_resid_sx))
  pwsum1 = df1 %>%group_by(f1) %>%summarize_all(sum)
  za_t_w_resid_sx=data.matrix(pwsum1[,-1])

  df2 = tibble(f2,as_tibble(w_resid_sx))
  pwsum2 = df2 %>%group_by(f2) %>%summarize_all(sum)
  zb_t_w_resid_sx=data.matrix(pwsum2[,-1])

  cov_z = sigma_a_sqhat*t(za_t_w_resid_sx)%*%za_t_w_resid_sx + sigma_b_sqhat*t(zb_t_w_resid_sx)%*%zb_t_w_resid_sx + (t(w_resid_sx)%*%(x - sx$fitted_val))#t(w_resid_sx)%*%(x - sx$fitted_val)
  cov_z_dispersion = sigma_a_sqhat*t(za_t_w_resid_sx)%*%za_t_w_resid_sx + sigma_b_sqhat*t(zb_t_w_resid_sx)%*%zb_t_w_resid_sx + (t(w_resid_sx)%*%(x - sx$fitted_val)*backfit$dispersion)#t(w_resid_sx)%*%(x - sx$fitted_val)

  cov_sandwiched = xt_w_i_minus_s_x_inv%*%cov_z%*%xt_w_i_minus_s_x_inv
  cov_sandwiched_dispersion = xt_w_i_minus_s_x_inv%*%cov_z_dispersion%*%xt_w_i_minus_s_x_inv

  return(list(cov=cov_sandwiched,cov_dispersion=cov_sandwiched_dispersion,iter=sx$iter,phat=phat))
}

var_betahat_logistic_under_glmm <-function(x,f1,f2,trace.it=TRUE,backfit,cov,logistic_fit_w){
  #it estimates cov(\hat{\beta}_{logistic}) under glmm

  sigma_a_sqhat = backfit$sigma_a^2
  sigma_b_sqhat = backfit$sigma_b^2
  w_backfit = cov$phat*(1-cov$phat)

  w_x = logistic_fit_w*x

  df11 = tibble(f1,as_tibble(w_x))
  pwsum11 = df11 %>%group_by(f1) %>%summarize_all(sum)
  za_t_wx=data.matrix(pwsum11[,-1])

  df22 = tibble(f2,as_tibble(w_x))
  pwsum22 = df22 %>%group_by(f2) %>%summarize_all(sum)
  zb_t_wx=data.matrix(pwsum22[,-1])

  cov_zx =  sigma_a_sqhat*t(za_t_wx)%*%za_t_wx + sigma_b_sqhat*t(zb_t_wx)%*%zb_t_wx + t(x)%*%(as.vector(logistic_fit_w^2/w_backfit)*x)
  cov_zx_dispersion = sigma_a_sqhat*t(za_t_wx)%*%za_t_wx + sigma_b_sqhat*t(zb_t_wx)%*%zb_t_wx + t(x)%*%(as.vector(logistic_fit_w^2*backfit$dispersion/w_backfit)*x)

  return(list(cov_zx=cov_zx,cov_zx_dispersion=cov_zx_dispersion))
}
