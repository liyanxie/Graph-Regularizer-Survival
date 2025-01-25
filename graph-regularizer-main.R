# Implementation of Graph-regularized Cox

rm(list = ls())
library(survival)
library(grpreg)
library(glmnet)
library(MASS)
library(coxed) 
library(Rlab)
library(ncvreg)
library(penalized)
library(lqmm)
library(ggplot2)
library(ggraph)
library(igraph)
library(reshape2)
library(pheatmap)

# set working directory
# setwd("~/Working Directory")

source("graph-lasso-helper.R")

# graph structures
ERDOS = "Erdos"
RING = "Ring"
COMMUNITIES = "Communities" 

########################################################################
# inputs graph

graph_structure = ERDOS # COMMUNITIES #ERDOS #RING

num_iter = 1

EDGE_MAGNITUDE  = 0.5 # off-diagonal connectivity value in B matrix
T_MAX = 100 # max follow up time in survival
CENSOR_RATE = 0.3
SIZE_TRAIN = 5000
SIZE_TEST = 1000
P = 100 # variable dimension

# sparse graph parameters
MAX_DEGREE_VALUE = 10 

PROB_ERDOS = 0.1 # probablity of entry of B being EDGE_MAGNITUDE else 0
SEEDS = 1:1000

ind = 1:4 

# community graph parameters
PROB_INNER_COMMU = 0.7 #0.95 #0.5

PROB_OUTER_COMMU = 0.01 
COMMU1 =30  
COMMU2 = COMMU1+30  
COMMU3 = COMMU2+30 

# ring graph parameters
PROB_RING_ADD_EDGE = 0

##################################################################################################

if (graph_structure == ERDOS) {
  B = create_B_sparse(PROB_ERDOS, EDGE_MAGNITUDE, P)
  
  omega = B
  diag(omega) = 1
  
  Sigma0 = solve(omega)
  Sigma = nearPD(Sigma0)$mat
  
  degree = colSums(omega)
  max_degree = sapply(sort(degree, index.return=TRUE), `[`, length(degree) - ind + 1)
  max_ind = max_degree[, 'ix']
  
  Sigma_xy = rep(0, P)
  Sigma_xy[max_ind] = MAX_DEGREE_VALUE
  beta_true = omega %*% Sigma_xy
} 

if (graph_structure == RING) {
  EDGE_MAGNITUDE = 0.7
  B = matrix(0, P, P)
  for (i in 1:P) {
    for (j in 1:P) {
      B[i,j] = rbern(n = 1, prob = PROB_RING_ADD_EDGE) * EDGE_MAGNITUDE
      if (abs(i - j) < 2) {
        B[i, j] = EDGE_MAGNITUDE
      }
    }
  }
  B[1, P] = EDGE_MAGNITUDE
  B[P, 1] = EDGE_MAGNITUDE
  diag(B) = rep(1,P)
  
  omega = B
  Sigma = solve(omega)
  Sigma = nearPD(Sigma)$mat # make Sigma near positive definite
  
  #Sigma_xy = rep(1, P)
  degree = colSums(omega) # column sum
  max_degree = sapply(sort(degree, index.return=TRUE), `[`, length(degree) - ind + 1) # the maximum four degrees and their index, from 100th-96th
  max_ind = max_degree[, 'ix'] # the index of the maximum four degree
  
  Sigma_xy = rep(0, P)
  Sigma_xy[max_ind] = MAX_DEGREE_VALUE
  
  beta_true = omega %*% Sigma_xy
}

if (graph_structure == COMMUNITIES) {
  B = create_B_commu(PROB_INNER_COMMU, PROB_OUTER_COMMU, COMMU1, COMMU2, COMMU3, EDGE_MAGNITUDE, P)
  
  omega = B
  diag(omega) = 0.1 #1
  
  Sigma = solve(omega)
  Sigma = nearPD(Sigma,posd.tol = 1e-4)$mat # nearPD(Sigma,posd.tol = 1e-4)
  
  degree = colSums(omega) # column sum
  max_degree = sapply(sort(degree, index.return=TRUE), `[`, length(degree) - ind + 1) # the maximum four degrees and their index, from 100th-96th
  max_ind = max_degree[, 'ix'] # the index of the maximum four degree
  
  Sigma_xy = rep(0, P)
  Sigma_xy[max_ind] = MAX_DEGREE_VALUE
  beta_true = omega %*% Sigma_xy ## the true parameter value beta
} 

##################################################################################################
# graph lasso set up 

# N[[i]] is the set of neighbors for variable i, including i itself
N <- vector("list", P)
for (i in 1:P) {
  for (j in 1:P) {
    if (omega[i,j] != 0) {
      N[[i]] <- c(N[[i]], j)
    } 
  }
}

# group indicator for X tilde
# this is the same length as X tilde
# indicates the group that an element of X tilde belongs to
# this is used for identifying groups in the group lasso function call
# total should be sum_i (1 + i's number of neighbors)
X_tilde_group_indicator = c() # a long vector: repeat i for group_size_for_all_vars(i) times
group_size_for_all_vars = c() # a list of P, i-th entry indicates the number of neighbors (including itself) of node-i in the inverse covariance matrix omega;
for (i in 1:P) {
  group_size = length(N[[i]])
  group_size_for_all_vars = c(group_size_for_all_vars, group_size)
  rep_each_group = rep(i, group_size)
  X_tilde_group_indicator = c(X_tilde_group_indicator, rep_each_group)
}

group_size_cumsum = cumsum(group_size_for_all_vars) # cusum of group size

# add zero padding for the non-neighbors
V = list()
for (i in 1:P){
  V[[i]] = rep(0, P)
}

# X_tilde, V_tilde are the expanded form of X and V
# index for X_tilde and V_tilde
index_pool = 1:length(X_tilde_group_indicator)
V_tilde_ind = list(index_pool[1:group_size_cumsum[1]])
for (i in 2:P) {
  V_tilde_ind[[i]] = index_pool[(group_size_cumsum[i - 1] + 1) : group_size_cumsum[i]]
} # the list of index for 100 groups, the index is in the long vector X_tilde_group_indicator 


##################################################################################################
# simulation

#print(paste0("Graph structure is ", graph_structure, ", p = ", P, ", PROB_ERDOS=", PROB_ERDOS))

start_time = Sys.time()

for (i in 1:(num_iter)) {
  print(paste0("Current iteration: ", i, ", current time: ", Sys.time()))
  set.seed(SEEDS[i])
  
  # initialize all values
  tmp_result = c()
  beta_graph = rep(0, P) #matrix(0, nrow = num_iter, ncol = P)
  cmax_graph = 0.5
  
  # data generation
  X_train = mvrnorm(n = SIZE_TRAIN, mu = rep(0, P), Sigma)
  X_test = mvrnorm(n = SIZE_TEST, mu = rep(0, P), Sigma)
  
  
  # covariate-independent censoring
  simdata_train = sim.survdata(N = SIZE_TRAIN, T = T_MAX, type = "none", hazard.fun = NULL,
                               num.data.frames = 1, fixed.hazard = TRUE, knots = 8,
                               spline = TRUE, X = X_train, beta = beta_true, censor = CENSOR_RATE, censor.cond = FALSE)
  data_train = simdata_train$data
  
  simdata_test = sim.survdata(N = SIZE_TEST, T = T_MAX, type = "none", hazard.fun = NULL,
                              num.data.frames = 1, fixed.hazard = T, knots = 8,
                              spline = TRUE, X = X_test, beta = beta_true, censor = CENSOR_RATE, censor.cond = FALSE)
  data_test = simdata_test$data

  fit_cox = coxph(Surv(y, failed) ~ ., data = data_train) 
  y_train = fit_cox$y
  
  #############################################################################################################################
  # graph lasso
  
  time_now = Sys.time()
  X_N_train = list()
  for (j in 1:P){
    X_N_train[[j]] = X_train[,N[[j]]]
  }
  
  # X is the fully expanded model matrix
  X_tilde_train = data.frame(matrix(unlist(X_N_train), nrow = SIZE_TRAIN, byrow = F))
  fit_graph = grpsurv(X_tilde_train, y_train, X_tilde_group_indicator)
  cv_graph <- cv.grpsurv(X_tilde_train, y_train, X_tilde_group_indicator)
  V_tilde = coef(cv_graph) # Beta at minimum CVE
  
  X_N_test = list()
  for (j in 1:P){
    X_N_test[[j]] = X_test[,N[[j]]]
  }
  X_tilde_test = data.frame(matrix(unlist(X_N_test), nrow = SIZE_TEST, byrow = F))
  
  # update
  for (j in 1:P){
    V[[j]][N[[j]]] = V_tilde[V_tilde_ind[[j]]] 
    # V[[j]] is the zero padding vector of length p
    # N is the neighbor indicator
  }
  for (j in 1:P){
    beta_graph[j] = sum(V[[j]])
  }
  cmax_graph = get_cv_cindex_max(fit_graph, cv_graph, as.matrix(X_tilde_test), data_test) # c-index
}

end_time = Sys.time()

elapsed = end_time - start_time

print(elapsed)

