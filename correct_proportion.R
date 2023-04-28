# helper function to calculate the proportion of correctly identified coefficients
findProportion <- function(model, correct){
  # find number of non-zero coefficients for all lambdas
  num_non_zero_coeff <- apply(model$beta, 2, function(x){sum(x!=0)})
  print(num_non_zero_coeff)
    # find min lambda where  the number of nonzero coefficients chosen in the 
  #     fit match (or >) the true number of nonzero coefficients in the generative model
  chosenPenaltyParamIndex <- min(which(num_non_zero_coeff >= length(correct)))
  print(chosenPenaltyParamIndex)
  # calculate proportion of correctly identified parameters for the provided model
  non_zero_coeff_from_model <- which(model$beta[,chosenPenaltyParamIndex] != 0)
  print(non_zero_coeff_from_model)
  correctly_identified_coeff <- intersect(correct,non_zero_coeff_from_model)
  print(correctly_identified_coeff)
  propCorrectL <- length(correctly_identified_coeff)/num_non_zero_coeff[chosenPenaltyParamIndex]
  print(propCorrectL)
  return(propCorrectL)
}


##### main
# n is number of observations, p is the number of covariates,
# m is number of groups and g is number of generative groups
# 1. Simulation 1: n = 150,m = 100,p = 5000, (g = 1, 2, 3)
# 2. Simulation 2: n = 200, m = 400, p = 40000, (g = 1, 2, 3)



# install.packages("glmnet")
# install.packages("SGL")

# glmnet for elasticNet, SGL for sparse group lasso
library("glmnet")
library("SGL")
library("pls") # for stdize


# initialize a reproducible pseudorandom number generator
set.seed(100)


# simulated values of m, n, p and g
n <- c(60,70,150,200)
p <- c(100,200,1200,2000)
m <- c(10,50,100,200)
group_size <- p/m 
g <- c(1,2,3)

num_simulations = 4

# empty arrays to store final proportions
Elnet_res <- array(0, c(num_simulations,3))
SGL_res <- array(0, c(num_simulations,3))


# two simulations with each simulation for 1,2 and 3 groups
for(i in 1:num_simulations){
  for(num_groups in 1:3){
    index <- ceiling(1:p[i] / group_size[i])   ### needed for SGL
    
    # The columns of X are constructed as iid Gaussian.
    # The columns are created and merged in batches of the 
    # group_size to keep each group normally distributed
    X <- matrix(rnorm(n[i]*group_size[i]), nrow = n[i])
    for(cols in 2:m[i]){
      X <- cbind(X, X <- matrix(rnorm(n[i]*group_size[i]), nrow = n[i]))
    }
    X <- stdize(X)
    
    # beta = (1,2,...,5,0,0...,0) for g = 1
    # beta = (1,2,...,5,0,0...0,1,2...5,0,..,0) for g = 2
    beta <- rep(0,p[i]) 
    correct <- NULL
    for(g_i in 1:num_groups){
      beta[1:5 + (g_i-1)*group_size[i]] = seq(1,5)
      # correct non-zero parameters column indexes
      correct = c(correct, 1:5+(g_i-1)*group_size[i])
    }
    # noise to signal ratio = 0.5
    y <- stdize(X %*% beta) + 0.5*rnorm(n[i])
    
    # fit elasticNet and calculate proportions
    data <- list(x = X, y = y)
    elasticNet <- glmnet(X,y, family = "gaussian", nlam = 100, alpha = 0.5)
    proportionElnet <- findProportion(elasticNet, correct)
    
    
    # fit sparse group lasso and calculate proportions
    sparse_group_lasso <- SGL::SGL(data = data,
                                   index = index,
                                   type = 'linear',
                                   alpha = 0.95,
                                   verbose = TRUE,
                                   thresh = 0.001,
                                   gamma = 0.8,
                                   nlam = 100) #NOTE: change it to 1000. 10 is kept temporarily to speed up the process
    proportionSGL <-  findProportion(sparse_group_lasso, correct)
    
    # store results
    Elnet_res[i,num_groups] <- proportionElnet
    SGL_res[i,num_groups] <- proportionSGL
  }
}

Elnet_res
SGL_res
