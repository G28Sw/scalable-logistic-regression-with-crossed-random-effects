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


resu=backfitting_outer_loop(X,y,Fa,Fb,max_iter=1000,epsilon = 1e-8)
cat("estimate of beta is ",resu$beta,"\n")
cat("estimate of sigma a is ",resu$sigma_a,"\n")
cat("estimate of sigma b is ",resu$sigma_b,"\n")

