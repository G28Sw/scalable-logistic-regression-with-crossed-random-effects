## Simulation example for scalable logistic regression with crossed random effect


#------------------------------------------------------------------------------------
## Clear the environment ##
#------------------------------------------------------------------------------------
rm(list = ls())

#------------------------------------------------------------------------------------
## Install the required packages into ##
## R if it's not already installed.   ##
#------------------------------------------------------------------------------------
if(!require(MASS)) install.packages("MASS4")
if(!require(lme4)) install.packages("lme4")
if(!require(tibble)) install.packages("tibble")
if(!require(tidyverse)) install.packages("tidyverse")

#------------------------------------------------------------------------------------
## Load the packages into the R environment ##
#------------------------------------------------------------------------------------
library(lme4)
library(MASS)
library(tibble)
library(tidyverse)

## -----------------------------------------------------------------------------------
# The following avoid unnecessary warnings produced by glmer from lme4
## -----------------------------------------------------------------------------------
options(warn=-1)

#------------------------------------------------------------------------------------
## Source all the functions from "clubbed_backfitting_schall.R" ##
#------------------------------------------------------------------------------------
source("clubbed_backfitting_schall.R")


#------------------------------------------------------------------------------------
## Define the parameters for the simulation setting. ##
#------------------------------------------------------------------------------------


pow = seq(3,5.5,0.25)
pow_glmer = seq(3,4.5,0.25)
S_list = 10^pow
S_list_glmer = 10^pow_glmer
n_list_glmer_len = length(S_list_glmer)
niter = 100
rho = .56
kappa = .56
p = 7
autocorr = 0.5
sigma_a = .8
sigma_b = .4
intercept = -2
Upsilon_upper_lim = 1.2716
autocorr.mat <- function(p = 7, autocorr = 0.5){
  mat <- diag(p)
  return(rho^abs(row(mat)-col(mat)))
}

beta = c(-2,rep(0,7))

#Here we consider the first choice of beta made in page 16. 
#For the second choice change beta to -2+0.5*0:7 

#------------------------------------------------------------------------------------
## The following is a function to simulate the observation pattern $Z$ according to mechanism
#mentioned in Ghosh et al. (2022) and obtain estimates using logistic regression, clubbed
#backfitting and glmer from Bates et al. (2015).
#------------------------------------------------------------------------------------




get_estimates<-function(S_list = 10^seq(3,5.5,0.25), rho = 0.56, kappa = 0.56, beta = c(-2,rep(0,7)), sigma_a = 0.8, sigma_b=0.4, niter = 100, n_list_glmer_len = 7, p = 7, Upsilon_upper_lim = 1.2716, autocorr = 0.5){

  #define lists to store the results

  schall_backfit = list()
  var_schall = list()
  logistic_glm_coef = list()

  glmm_glmer_fixef = list()
  glmm_glmer_nAGQ_0_fixef = list()
  glmm_glmer_sigma = list()
  glmm_glmer_nAGQ_0_sigma = list()

  glmm_glmer_outer_iter = list()
  glmm_glmer_nAGQ_0_outer_iter = list()

  time_schall_backfit = list()
  time_schall_cov= list()
  time_logistic_glm = list()
  time_glmm_glmer = list()
  time_glmm_glmer_nAGQ_0 = list()
  n_list = matrix(rep(0,length(S_list)*niter),ncol=niter)


  for (i in 1:length(S_list)){

    schall_backfit[[i]] = list()
    var_schall[[i]] = list()
    logistic_glm_coef[[i]] = list()

    time_schall_backfit[[i]] = list()
    time_schall_cov[[i]] = list()
    time_logistic_glm[[i]] = list()

    if(i<=n_list_glmer_len){
      glmm_glmer_fixef[[i]] = list()
      glmm_glmer_nAGQ_0_fixef[[i]] = list()

      glmm_glmer_outer_iter[[i]] = list()
      glmm_glmer_nAGQ_0_outer_iter[[i]] = list()

      glmm_glmer_sigma[[i]] = list()
      glmm_glmer_nAGQ_0_sigma[[i]] = list()

      time_glmm_glmer[[i]] = list()
      time_glmm_glmer_nAGQ_0[[i]] = list()
    }


    S = floor(S_list[i])

    nf1 = ceiling(S^rho)
    nf2 = ceiling(S^kappa)
    Z = matrix(rep(0,nf1*nf2),nrow=nf1)



    for(kl in 1:niter){


      set.seed((4*i)+kl+1)

      #generate Z, X, random effects and then the binary response from a crossed random effect model
      for (ik in 1:nf1){
        for(jk in 1:nf2){
          Upsilon=runif(1,1,Upsilon_upper_lim)

          Z[ik,jk] = rbinom(1,1,Upsilon*(S/(nf1*nf2)))

        }
      }

      #The following is a faster implementation, for the purpose of reproducibility I keep the original code used while obtaining the results and use the above

      # Upsilon <- runif(nf1 * nf2, 1, Upsilon_upper_lim)
      # p_seq <-  pmin(Upsilon * (S / (nf1 * nf2)), 1)
      # Z <-  matrix(rbinom(nf1 * nf2, size = 1, prob = p_seq), nf1, nf2)


      n = sum(Z)
      n_list[i,kl] <- n


      Z <- Z[rowSums(Z) != 0, ]
      Znonzero <- Z[, colSums(Z) != 0]
      l <- which(Znonzero!=0,arr.ind = T)
      f1 <- l[,1]
      f2 <- l[,2]


      x <- cbind(1,matrix(mvrnorm(n = n, mu=rep(0,p), Sigma=autocorr.mat(p,autocorr)),n,p))
      nf1 <-  length(unique(f1))
      nf2 <- length(unique(f2))


      a = rnorm(nf1,0,sd=sigma_a)
      b = rnorm(nf2,0,sd=sigma_b)

      eta = as.vector(x%*%beta)+a[f1]+b[f2]
      prob = 1/(1+exp(-eta))
      y = rbinom(n,1,prob)


      #get estimates using plain logistic regression

      start_time_logistic_glm = Sys.time()
      logistic_model = glm(y~-1+x,family="binomial")
      logistic_glm_coef[[i]][[kl]] = coef(logistic_model)
      end_time_logistic_glm = Sys.time()

      time_logistic_glm[[i]][[kl]] = difftime(end_time_logistic_glm,start_time_logistic_glm,units = "secs")



      #get estimates using clubbed backfitting


      start_time_schall_backfit = Sys.time()
      schall_backfit[[i]][[kl]] = backfitting_outer_loop(x,y,f1,f2,betahat = rep(0,(ncol(x))),sigma_a_sqhat=1,sigma_b_sqhat=1,epsilon=1e-8,max_iter=1000,over_dispersed=TRUE,trace.it=FALSE,epsilon_backfitting=1e-8)

      end_time_schall_backfit = Sys.time()

      start_time_schall_cov = Sys.time()
      var_schall[[i]][[kl]] = var_betahat_glmm(x,f1,f2,trace.it=FALSE, schall_backfit[[i]][[kl]],epsilon=1e-8,max_iter=1000)

      end_time_schall_cov = Sys.time()

      time_schall_backfit[[i]][[kl]] = difftime(end_time_schall_backfit,start_time_schall_backfit,units = "secs")

      time_schall_cov[[i]][[kl]] = difftime(end_time_schall_cov,start_time_schall_cov,units = "secs")

      f1_factor = as.factor(f1)
      f2_factor = as.factor(f2)

      #get estimates using glmer from lme4 for small problem size


      if(i<=n_list_glmer_len){

        start_time_lme4 = Sys.time()
        glmer_obj = glmer(y~-1+x+(1|f1_factor)+(1|f2_factor),family=binomial,control = glmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=1000)))
        end_time_lme4 = Sys.time()
        time_glmm_glmer[[i]][[kl]] = difftime(end_time_lme4,start_time_lme4,units = "secs")

        #storing the glmer objects are expensive
        glmm_glmer_fixef[[i]][[kl]] = fixef(glmer_obj)
        glmm_glmer_sigma[[i]][[kl]] = c(sqrt(unlist(summary(glmer_obj)$varcor)[1]),sqrt(unlist(summary(glmer_obj)$varcor)[2]))
        glmm_glmer_outer_iter[[i]][[kl]] = summary(glmer_obj)$optinfo$feval-1

        start_time_lme4_nAGQ_0 = Sys.time()
        glmer_nAGQ_0_obj = glmer(y~-1+x+(1|f1_factor)+(1|f2_factor),family=binomial,nAGQ=0,control = glmerControl(optimizer="Nelder_Mead",optCtrl=list(maxfun=1000)))
        end_time_lme4_nAGQ_0 = Sys.time()
        time_glmm_glmer_nAGQ_0[[i]][[kl]] = difftime(end_time_lme4_nAGQ_0,start_time_lme4_nAGQ_0,units = "secs")

        glmm_glmer_nAGQ_0_fixef[[i]][[kl]]  = fixef(glmer_nAGQ_0_obj)
        glmm_glmer_nAGQ_0_sigma[[i]][[kl]]  = c(sqrt(unlist(summary(glmer_nAGQ_0_obj)$varcor)[1]),sqrt(unlist(summary(glmer_nAGQ_0_obj)$varcor)[2]))
        glmm_glmer_nAGQ_0_outer_iter[[i]][[kl]] = summary(glmer_nAGQ_0_obj)$optinfo$feval-1
      }
    }
  }

  return(list(problem_size = S_list, logistic_coef = logistic_glm_coef, logistic_time = time_logistic_glm,schall_backfit = schall_backfit, schall_cov = var_schall, schall_backfit_time = time_schall_backfit, schall_cov_time = time_schall_cov, glmer_fixef = glmm_glmer_fixef, glmer_sigma = glmm_glmer_sigma, glmer_iter = glmm_glmer_outer_iter, glmer_time = time_glmm_glmer, glmer_nAGQ_0_fixef = glmm_glmer_nAGQ_0_fixef, glmer_nAGQ_0_sigma = glmm_glmer_nAGQ_0_sigma, glmer_iter_nAGQ_0 = glmm_glmer_nAGQ_0_outer_iter, glmer_time_nAGQ_0 = time_glmm_glmer_nAGQ_0, n_list = n_list))
}

#-----------------------------------------------------------------------------------------
## Next, we obtain the estimates and find the MSE based on the estimates and known parameter. ##
#----------------------------------------------------------------------------------------

results_mse = get_estimates(S_list, rho, kappa, beta, sigma_a, sigma_b, niter, n_list_glmer_len, p, Upsilon_upper_lim, autocorr)


logistic_coef_matrix = matrix(unlist(results_mse$logistic_coef),ncol=p+1,byrow=TRUE)
logistic_residual = data.frame(cbind("S"=rep(results_mse$problem_size, each=niter), t(t(logistic_coef_matrix) - beta)))
logistic_mse = logistic_residual %>%
  group_by(S) %>%
  dplyr::summarize(across(everything(), ~sum(.x^2)/n()))


backfit_coef_matrix = matrix(,nrow(logistic_coef_matrix),ncol(logistic_coef_matrix))
backfit_iter = matrix(,nrow(logistic_coef_matrix))
for (i in 1:length(S_list)){
  for(kl in 1:niter){
    backfit_coef_matrix[((i-1)*niter)+kl,] = results_mse$schall_backfit[[i]][[kl]]$beta
    backfit_iter[((i-1)*niter)+kl,] = results_mse$schall_backfit[[i]][[kl]]$iteration
  }
}
backfit_schall_residual = data.frame(cbind("S"=rep(results_mse$problem_size,each=niter),t(t(backfit_coef_matrix) - beta)))
backfit_schall_mse = backfit_schall_residual %>%
  group_by(S) %>%
  dplyr::summarize(across(everything(), ~sum(.x^2)/n()))



glmer_coef_matrix = matrix(unlist(results_mse$glmer_fixef),ncol=p+1,byrow=TRUE)
glmer_residual = data.frame(cbind("S"=rep(results_mse$problem_size[1:n_list_glmer_len],each=niter), t(t(glmer_coef_matrix) - beta)))
glmer_mse = glmer_residual %>%
  group_by(S) %>%
  dplyr::summarize(across(everything(), ~sum(.x^2)/n()))




glmer_nAGQ_0_coef_matrix = matrix(unlist(results_mse$glmer_nAGQ_0_fixef),ncol=p+1,byrow=TRUE)
glmer_nAGQ_0_residual = data.frame(cbind("S"=rep(results_mse$problem_size[1:n_list_glmer_len],each=niter), t(t(glmer_nAGQ_0_coef_matrix) - beta)))
glmer_nAGQ_0_mse = glmer_nAGQ_0_residual %>%
  group_by(S) %>%
  dplyr::summarize(across(everything(), ~sum(.x^2)/n()))

#-----------------------------------------------------------------------------------------
## Now, we plot the MSE. The following generates the two plots in figure 2 of section 5
#----------------------------------------------------------------------------------------

s = -1
s2 = -min(rho,kappa)
avg_n_list = rowMeans(results_mse$n_list)

plot(log10(avg_n_list),log10(rowMeans(backfit_schall_mse[,3:9])),ylab="log10(MSE)",xlab="log10(average sample size)",main="log10(MSE) with log10(average sample size) without intercept",ylim=c(-4.4,-1.5),pch=2)
lines(log10(avg_n_list),s*log10(avg_n_list)-s*log10(avg_n_list)[1]+log10(rowMeans(backfit_schall_mse[,3:9])[1]),col="black",lty=1,lwd=1)
points(log10(avg_n_list[1:n_list_glmer_len]),log10(rowMeans(glmer_mse[,3:9])),pch=16)
points(log10(avg_n_list[1:n_list_glmer_len]),log10(rowMeans(glmer_nAGQ_0_mse[,3:9])))
points(log10(avg_n_list),log10(rowMeans(logistic_mse[,3:9])),pch=4)
legend("bottomleft",c("Schall MSE","GLMER MSE","GLMER nAGQ 0 MSE","Logistic MSE", "Line with Slope = -1"),col=c("black","black","black","black","black"),
       lty=c(NA,NA,NA,NA,1),lwd =c(1,1,1,1,1),pch = c(2,16,1,4,NA))

plot(log10(avg_n_list),log10(pull(backfit_schall_mse, V2)),ylab="log10(MSE)",xlab="log10(average sample size)",main="log10(MSE) with log10(average sample size) only intercept",ylim=c(-4,-1),pch=2)
lines(log10(avg_n_list),s*log10(avg_n_list)-s*log10(avg_n_list)[1]+log10(pull(backfit_schall_mse, V2)[1]),col="black",lty=1,lwd=1)
lines(log10(avg_n_list),s2*log10(avg_n_list)-s2*log10(avg_n_list)[1]+log10(pull(backfit_schall_mse, V2)[1]),col="black",lty=2,lwd=2)
points(log10(avg_n_list[1:n_list_glmer_len]),log10(pull(glmer_mse, V2)),pch=16)
points(log10(avg_n_list[1:n_list_glmer_len]),log10(pull(glmer_nAGQ_0_mse, V2)),pch=1)
points(log10(avg_n_list),log10(pull(logistic_mse, V2)),pch=4)
legend("bottomleft",c("Schall MSE","GLMER MSE","GLMER nAGQ 0 MSE","Logistic MSE", "Line with Slope = -1",paste0("Line with Slope = ",s2)),col=c("black","black","black","black","black","black"),
       lty=c(NA,NA,NA,NA,1,2),lwd =c(1,1,1,1,2,2),pch = c(2,16,1,4,NA,NA))

