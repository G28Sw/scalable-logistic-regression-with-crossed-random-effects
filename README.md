# Overview

The cost of both generalized least squares (GLS) and Gibbs sampling in a crossed random effects model is superlinear in N, the number of observations. Ghosh et al. (2021) developed a backfitting algorithm that reduces the cost to O(N). Here we extend that method to a generalized linear mixed model for logistic regression. This is based on adaptation of penalized quasi-likelihood algorithm of Schall(1991) using clubbed backfitting.

This repository contains codes for the article: "Scalable logistic regression with crossed random effects"

The article link: https://arxiv.org/abs/2105.13747

# Example 
Below is a basic example of using clubbed backfitting for logistic regression with crossed random effect. 


```r
source("clubbed_backfitting_schall.R")
set.seed(100)
p=7
n=40000
X=matrix(rnorm(n*p),n,p)
X=cbind(1,X)
nf1=400 #corresponds to rho ~ 0.56
Fa=sample(seq(nf1),n,replace=T)
cat("R is ",length(unique(Fa)))
cat("min a",min(table(Fa)))#make sure there are good enough number of points for PQL to be good
nf2=400 #corresponds to kappa ~ 0.56
Fb=sample(seq(nf2),n,replace=T)
cat("C is ",length(unique(Fb)))
cat("min b",min(table(Fb)))#make sure there are good enough number of points for PQL to be good

a = rnorm(nf1,0,0.8)
b = rnorm(nf2,0,0.4)
beta = c(-2,rep(0,7))
mu = X%*%beta + a[Fa] + b[Fb]

y=rbinom(n,1,1/(1+exp(-mu)))
```

We use clubbed backfitting with default starting values.
```r
resu=backfitting_outer_loop(X,y,Fa,Fb)
cat("estimate of beta is ",resu$beta,"\n")
cat("estimate of sigma a is ",resu$sigma_a,"\n")
cat("estimate of sigma b is ",resu$sigma_b,"\n")
```


