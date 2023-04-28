library(MASS)
library(survival)
library(SGL)
library(glmnet)
n=200
p=30
q=0
m=5
rho = 0.5 # correlation between adjacent covariates in a block
rho.ar = TRUE # rho^|i-j| block correlation structure; first order autoregressive
library(pec)

# Function to extract survival probability predictions from SGL
predictSurvProb.my_sgl <- function (object, newdata, times, ...) 
{
  require(SGL)
  x <- as.matrix(newdata[,3:(2+object$ngenes+object$nclinical)])
  lp <- predictSGL(object$model, newX=x, lam=object$lambda_index)
  lp
  #
  require(survival)
  survival.coxph <- getFromNamespace("coxph", ns = "survival")
  survival.survfit.coxph <- getFromNamespace("survfit.coxph", ns = "survival")
  survival.summary.survfit <- getFromNamespace("summary.survfit", ns = "survival")
  #     
  f <- survival.coxph(Surv(newdata$time, newdata$status) ~ lp)
  survfit.object <- survival.survfit.coxph(f, newdata=data.frame(lp=lp), 
                                           se.fit=FALSE, conf.int=FALSE)
  inflated.pred <- survival.summary.survfit(survfit.object, times = times)
  p <- t(inflated.pred$surv)
  if ((miss.time <- (length(times) - NCOL(p))) > 0) 
    p <- cbind(p, matrix(rep(NA, miss.time * NROW(p)), nrow = NROW(p)))
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)) 
    stop("Prediction failed")
  p
}
# Function to extract survival probability predictions from elnet/lasso
predictSurvProb.my_elnet <- function (object, newdata, times, ...) 
{
  require(glmnet)
  x <- as.matrix(newdata[,3:(2+object$ngenes+object$nclinical)])
  lp <- predict(object$model, type="link", newx=x, s=object$lambda_index)
  
  require(survival)
  survival.coxph <- getFromNamespace("coxph", ns = "survival")
  survival.survfit.coxph <- getFromNamespace("survfit.coxph", ns = "survival")
  survival.summary.survfit <- getFromNamespace("summary.survfit", ns = "survival")
  #     
  f <- survival.coxph(Surv(newdata$time, newdata$status) ~ lp)
  survfit.object <- survival.survfit.coxph(f, newdata=data.frame(lp=lp), 
                                           se.fit=FALSE, conf.int=FALSE)
  inflated.pred <- survival.summary.survfit(survfit.object, times = times)
  p <- t(inflated.pred$surv)
  if ((miss.time <- (length(times) - NCOL(p))) > 0) 
    p <- cbind(p, matrix(rep(NA, miss.time * NROW(p)), nrow = NROW(p)))
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)) 
    stop("Prediction failed")
  p
}

# find scale parameter for exponential cox
Expo <- function (times, surv) {
  z1 <- -log(surv[1])
  t1 <- times[1]
  lambda <- z1/(t1)
  list(scale =lambda)
}

# effects
#bg = c(0.75, 0.75, 0.5, 0.5, 0.25, -0.25, rep(0, p-6))
#bc = c(-1.0, 1.0, 0.3, 0, -0.3, rep(0,q-5))

bX <- c(0.75, 0.75, 0.5, 0.5, 0.25, -0.25, rep(0, p-6))
# Correlation structure for genomic variables
Sigma <- diag(1, p)
if(rho != 0){
  k <- sum(bX!=0)     			# number of related covariates
  m <- min(m, (p-k) %/% k)  # number of correlated but unrelated vars per related var
  
  for (i in 1:k){ # loop through correlated vars (1:k)
    cov.mat       <- matrix(rho, ncol=m+1, nrow=m+1)  
    diag(cov.mat) <- rep(1, m+1) 
    if(rho.ar){
      # rho^|i-j| auto regressive correlation structure
      for(index in 1:(m+1)){
        cov.mat[index,] <- c(rev(cumprod(rep(rho, index-1))), 1, cumprod(rep(rho, m+1-index)))
      }
    }
    
    # position i in block of related covs (1:k) plus i-th block of correlated covs
    Sigma[c(i, k+(1:m)+(i-1)*m), c(i, k+(1:m)+(i-1)*m)] <- cov.mat
  }
}

# index <- ceiling(1:(p+q) / ((p+q)/m))
index <- c(c(1:6,rep(1:6,rep(4,6))))

################
# seeding data:
set.seed(12345)

# covariates
# genomic
means = rep(0, p)
# Sigma = diag(1, p)
Xg <- mvrnorm(n, means, Sigma)

# clinical
#x1 = rbinom(n=n,size=1, prob=0.7)
#x2 = rbinom(n=n,size=1, prob=0.3)
#x3 = rnorm(n=n, mean=0, sd=1)
#x4 = rnorm(n=n, mean=0, sd=1)
#x5 = rnorm(n=n, mean=0, sd=1)
#Xc <- cbind(x1,x2,x3,x4,x5)
X <- data.frame(Xg)
#names(X) = c(paste("G",1:p, sep=""), paste("C",1:q, sep=""))
names(X) = c(paste("G",1:p, sep=""))
# censoring function
# administrative censoring: uniform data entry.
ACT = 36; FUT = 72
cens.start = FUT; cens.end = ACT+FUT
cens <- runif(n, cens.start, cens.end)
# # baseline hazard
# # median survival time =36 ; S(36) = 1/2 - def of median survival time; S(t)=ho*exp(beta*X); -> calculate h0 
# med.os = 36
# h0 <- round(log(2)/med.os,2)
# #hazard function = h0*exp(bX*X)
# h <- h0 * exp( as.matrix(X) %*% bX )
# # times
# survival distribution (exponential, survival probs 0.5 at 36 months)
surv <- Expo(times=c(36), surv=c(0.5))
dt = -log(runif(n)) * (1/surv$scale) * exp(-as.matrix(X) %*% bX)
status <- ifelse(dt <= cens, 1, 0)
time <- pmin(dt, cens)


################
# Test data:
set.seed(6789)

# covariates
# genomic
Xg. <- mvrnorm(n, means, Sigma)

# clinical
#x1. = rbinom(n=n,size=1, prob=0.7)
#x2. = rbinom(n=n,size=1, prob=0.3)
#x3. = rnorm(n=n, mean=0, sd=1)
#x4. = rnorm(n=n, mean=0, sd=1)
#x5. = rnorm(n=n, mean=0, sd=1)
#Xc. <- cbind(x1,x2,x3,x4,x5)
X. <- data.frame(Xg.)
#names(X.) = c(paste("G",1:p, sep=""), paste("C",1:q, sep=""))
names(X.) = c(paste("G",1:p, sep=""))
# censoring function
# administrative censoring: uniform data entry.
cens. <- runif(n, cens.start, cens.end)
dt. = -log(runif(n)) * (1/surv$scale) * exp(-as.matrix(X.) %*% bX)
status. <- ifelse(dt <= cens., 1, 0)
time. <- pmin(dt., cens.)

################
#Finalise both data sets

#center and scale (all variables, including clinical covariables):
sd.train <- apply(X,2,sd) 
X  <- scale(X,  scale=TRUE)
X. <- scale(X., scale=sd.train)

train <- data.frame("time"=time, "status"=status,  X)
test <- data.frame("time"=time., "status"=status., X.)


################################################
# Fit models

set.seed(101112)

################
# Oracle
oracle.model=coxph(Surv(time, status) ~ G1+G2+G3+G4+G5+G6, data=train)
oracle.summary <- summary(oracle.model)

################
# SGL with default alpha=0.95 ~ almost lasso
sgl_cv_lasso <- SGL::cvSGL(data = list(x =X, time=time, status=status),
                           index = index,
                           type = 'cox',
                           verbose = TRUE,
                           nlam = 50,
                           standardize = FALSE
)

# SGL with default alpha=0.05 ~ almost group lasso
sgl_cv_glasso <- SGL::cvSGL(data = list(x =X, time=time, status=status),
                            index = index,
                            type = 'cox',
                            verbose = TRUE,
                            nlam = 50,
                            standardize = FALSE,
                            alpha=0.05
                            
)
# to see the death times of the surviving candidates (status==1), print sgl_cv_lasso$fit$death.times)

# choose beta corresponding to lambda which has minimum lldiff (maximum likelihood)
# lldiff => (vector of cross validated negative log likelihoods (squared error loss in the linear case, along the regularization path))
index_of_lambda_with_min_lldiff = which(sgl_cv_lasso$lldiff == min(sgl_cv_lasso$lldiff))
lambda_with_min_lldiff = sgl_cv_lasso$lambdas[index_of_lambda_with_min_lldiff]
final_betas = sgl_cv_lasso$fit$beta[,index_of_lambda_with_min_lldiff]
plot(sgl_cv_lasso)
sgl.summary <- summary(sgl_cv_lasso)

index_of_lambda_with_min_lldiff2 = which(sgl_cv_glasso$lldiff == min(sgl_cv_glasso$lldiff))
lambda_with_min_lldiff2 = sgl_cv_glasso$lambdas[index_of_lambda_with_min_lldiff2]
final_betas2 = sgl_cv_glasso$fit$beta[,index_of_lambda_with_min_lldiff2]
plot(sgl_cv_glasso)
sgl.summary <- summary(sgl_cv_glasso)

# NOTE: need to provide a dummy lambda (0.5 here) since the SGL library has a bug if we provide just one lambda in lambdas
sgl_fit = SGL(data=list(x =X, time=time, status=status),
              index=index,
              type = 'cox',
              verbose = TRUE,
              lambdas = c(lambda_with_min_lldiff, 0.5)
)

sgl_fit2 = SGL(data=list(x =X, time=time, status=status),
               index=index,
               type = 'cox',
               verbose = TRUE,
               lambdas = c(lambda_with_min_lldiff2, 0.5)
)

################
# Elastic Net
elnet_cv <- cv.glmnet(x=as.matrix(train[,-(1:2)]), y=Surv(train$time, train$status), family="cox", 
                      alpha=0.5, standardize=FALSE)
final_betas_elnet <- coef(elnet_cv, s=elnet_cv$lambda.min)

################
# Lasso
lasso_cv <- cv.glmnet(x=as.matrix(train[,-(1:2)]), y=Surv(train$time, train$status), family="cox", 
                      alpha=1, standardize=FALSE)
final_betas_lasso <- coef(lasso_cv, s=lasso_cv$lambda.min)


sgl.object <- list(model=sgl_fit, lambda_index=1, ngenes=p, nclinical=q)
class(sgl.object) <- "my_sgl"

sgl2.object <- list(model=sgl_fit2, lambda_index=1, ngenes=p, nclinical=q)
class(sgl2.object) <- "my_sgl"

elnet.object <- list(model=elnet_cv$glmnet.fit, lambda_index=elnet_cv$lambda.min, ngenes=p, nclinical=q)
class(elnet.object) <- "my_elnet"

lasso.object <- list(model=lasso_cv$glmnet.fit, lambda_index=lasso_cv$lambda.min, ngenes=p, nclinical=q)
class(lasso.object) <- "my_elnet"

# finding probability error curves
pe <- pec(list(
  #"Clinical only"=coxph(Surv(time,status)~C1+C2+C3+C4+C5, data=train, x=TRUE),
  "Oracle model"=coxph(Surv(time,status)~G1+G2+G3+G4+G5+G6, data=train, x=TRUE),
  "Frequentist elnet"=elnet.object,
  "Frequentist lasso"=lasso.object,
  "Frequentist sgl_lasso"=sgl.object,
  "Frequentist sgl_glasso"=sgl2.object),
  data=test, formula=Surv(time,status)~1, times=72)

# finding concorance index
cind <- cindex(list(
  #"Clinical only"=coxph(Surv(time,status)~C1+C2+C3+C4+C5, data=train, x=TRUE),
  "Oracle model"=coxph(Surv(time,status)~G1+G2+G3+G4+G5+G6, data=train, x=TRUE),
  "Frequentist elnet"=elnet.object,
  "Frequentist lasso"=lasso.object,
  "Frequentist sgl_lasso"=sgl.object,
  "Frequentist sgl_glasso"=sgl2.object),
  data=test, formula=Surv(time,status)~1, times=72)

# summarizing prediction error curves (aka integrate Brier scores, ibs)
ibs <- crps(pe, times=72)

# plotting IBS
xlim <- c(0,72)
ylim <- c(0,0.35)
cols <- c("black","purple","goldenrod","green3", "chocolate","pink")
plot(pe, xlim=xlim, ylim=ylim, xlab="Time t", legend=FALSE, lwd=2, col=cols, models=c(1:6,7))
legend("topleft", bty="n", lwd=2, cex=0.75, col=cols,
       legend=c(paste("Reference (IBS = ",round(ibs[1],3),")", sep=""),
                #paste("Clinical only (IBS = ",round(ibs[2],3),")", sep=""),
                paste("Oracle model (IBS = ",round(ibs[2],3),")", sep=""),
                paste("Frequentist elnet (IBS = ",round(ibs[3],3),")", sep=""),
                paste("Frequentist lasso (IBS = ",round(ibs[4],3),")", sep=""),
                paste("Frequentist sgl_lasso (IBS = ",round(ibs[5],3),")", sep=""),
                paste("Frequentist sgl_glasso (IBS = ",round(ibs[6],3),")", sep="")
       ))