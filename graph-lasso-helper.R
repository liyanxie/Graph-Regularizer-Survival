# Helper functions for graph regularizer
# by Xi He
# ISyE @ GaTech
# 6/30/2020
# xi.isabel.he@gmail.com

# for erods, for every pair of variables, probability of being neighbor = PROB_ERDOS

create_B_sparse = function(probability, edge_mag, dimension) {
  B = matrix(0, dimension, dimension)
  set.seed(0)
  predictor_relation_vector = rbern(n = choose(n = dimension, k = 2), prob = probability) * edge_mag
  B[lower.tri(B)] = predictor_relation_vector
  B = t(B)
  B[lower.tri(B)] = predictor_relation_vector
  return(B)
}

create_B_commu = function(inner_prob, outer_prob, commu1, commu2, commu3, edge_mag, dimension) {
  B = matrix(0, dimension, dimension)
  set.seed(0)
  for (i in 1:dimension) {
    for (j in i:dimension) {
      if(((i<=commu1) && (j<=commu1))
         || (((commu1<i) & (i<=commu2)) && ((commu1<j) & (j<=commu2)))
         || (((commu2<i) & (i<=commu3)) && ((commu2<j) & (j<=commu3)))){
        B[i,j] = rbern(n=1, prob = inner_prob) * edge_mag
      } else {
        B[i,j] = rbern(n=1, prob = outer_prob) * edge_mag
      }
    }
  }
  B = B + t(B)
  diag(B) = rep(1, dimension)
  return(B)
}

# get omega (inverse covariance of predictors)
# B is p*p matrix
get_omega = function(B, p) {
  if ((get_kappa_given_lambda(B, lambda_lower_bound, p) > p) 
      && (get_kappa_given_lambda(B, lambda_upper_bound, p) < p)) {
    while (abs((get_kappa_given_lambda(B, lambda_lower_bound, p) - p)) > kappa_threshold
           && abs((get_kappa_given_lambda(B, lambda_upper_bound, p) - p)) > kappa_threshold) {
      midpoint = 0.5 * (lambda_lower_bound + lambda_upper_bound)
      if (get_kappa_given_lambda(B, midpoint, p) < p) {
        lambda_upper_bound = midpoint
      } else {
        lambda_lower_bound = midpoint
      }
    }
    if (abs(get_kappa_given_lambda(B, lambda_lower_bound, p) - p) <= kappa_threshold) {
      lambda = lambda_lower_bound
    } else {
      lambda = lambda_upper_bound
    } 
    diag(B) = rep(lambda, p)
  } else {
    print('Error: Reset lambda initial bounds.')
  }
  omega = B/lambda
  return(omega)
}

# get condition number of omega given delta as in omega = B + delta*I
get_kappa_given_lambda = function(B_matrix, delta, p_dimension) {
  diag(B_matrix) = rep(delta, p_dimension)
  condition_num = kappa(B_matrix) # condition number
  return(condition_num)
}

# compute predicted probs from beta estimates
get_probs_from_beta = function(X_train, X_test, beta) {
  rM = (colMeans(X_train) %*% as.numeric(beta))[1,1]
  r = as.matrix(X_test) %*% as.numeric(beta)
  probs = r - rM
  return(probs)
}

# find lambda that returns max cindex in cv
get_cv_cindex_max = function(fit_model, cv_model, X_test, data_test) {
  c_tmp = rep(0, length(cv_model$lambda))
  for (j in 1:length(cv_model$lambda)) {
    probs = predict(fit_model, X_test, cv_model$lambda[j], type = "link")
    c_tmp[j] = get_cindex(probs, data_test)
  }
  cmax = max(c_tmp)
  return(cmax)
}

# get_cv_cindex_max_from_beta returns the same result as get_cv_cindex_max
# it is used for models that don't have fit_model and cv_model
get_cv_cindex_max_from_beta = function(betas, X_train, X_test, data_test) {
  c_tmp = rep(0, nrow(betas))
  for (j in 1:nrow(betas)) {
    rM = (colMeans(X_train) %*% as.numeric(betas[j,]))[1,1]
    r = as.matrix(X_test) %*% as.numeric(betas[j,])
    probs = r - rM
    c_tmp[j] = get_cindex(probs, data_test)
  }
  cmax = max(c_tmp)
  return(cmax)
  
}

# alasso max cindex, fused 
# get_cv_cindex_max_alasso should return the same result as get_cv_cindex_max
# it is used for models that don't have fit_model and cv_model
get_cv_cindex_max_alasso = function(newbeta.s, X_train, X_test, data_test) {
  c_tmp = rep(0, nrow(newbeta.s))
  for (j in 1:nrow(newbeta.s)) {
    rM = (colMeans(X_train) %*% as.numeric(newbeta.s[j,]))[1,1]
    r = as.matrix(X_test) %*% as.numeric(newbeta.s[j,])
    probs = r - rM
    c_tmp[j] = get_cindex(probs, data_test)
  }
  cmax = max(c_tmp)
  return(cmax)
  
}



# compute c-index
get_cindex = function(probs, data_test) {
  cindex = as.numeric(1 - rcorr.cens(probs, with(data_test, Surv(y, failed)), outx=FALSE)[1])
  return(cindex)
}

# fit model and get result
get_performance_result = function(Sigma, beta_true, i, censor_rate) {
  tmp_result = c()
  
  beta_scad = beta_cox = beta_lasso = beta_ridge = beta_enet = beta_graph = beta_alasso = rep(0, P)#matrix(0, nrow = num_iter, ncol = P)
  l2_scad = l2_cox = l2_lasso = l2_ridge = l2_enet = l2_graph = l2_alasso = -1
  l1_scad = l1_cox = l1_lasso = l1_ridge = l1_enet = l1_graph = l1_alasso = -1
  linf_scad = linf_cox = linf_lasso = linf_ridge = linf_enet = linf_graph = linf_alasso = -1
  l2_scad_rel = l2_cox_rel = l2_lasso_rel = l2_ridge_rel = l2_enet_rel = l2_graph_rel = l2_alasso_rel = -1
  l1_scad_rel = l1_cox_rel = l1_lasso_rel = l1_ridge_rel = l1_enet_rel = l1_graph_rel = l1_alasso_rel = -1
  linf_scad_rel = linf_cox_rel = linf_lasso_rel = linf_ridge_rel = linf_enet_rel = linf_graph_rel = linf_alasso_rel = -1
  rpe_scad = rpe_cox = rpe_lasso = rpe_ridge = rpe_enet = rpe_graph = rpe_alasso = -1
  rpl_scad = rpl_cox = rpl_lasso = rpl_ridge = rpl_enet = rpl_graph = rpl_alasso = 0
  fpr_scad = fpr_cox = fpr_lasso = fpr_ridge = fpr_enet = fpr_graph = fpr_alasso = 0
  fnr_scad = fnr_cox = fnr_lasso = fnr_ridge = fnr_enet = fnr_graph = fnr_alasso = 0
  cmax_scad = c_cox = cmax_lasso = cmax_ridge = cmax_enet = cmax_graph = cmax_alasso = 0.5
  
  
  X_train = mvrnorm(n = SIZE_TRAIN, mu = rep(0, P), Sigma)
  X_test = mvrnorm(n = SIZE_TEST, mu = rep(0, P), Sigma)
  
  simdata_train = sim.survdata(N = SIZE_TRAIN, T = T_MAX, type = "none", hazard.fun = NULL,
                               num.data.frames = 1, fixed.hazard = T, knots = 8,
                               spline = TRUE, X = X_train, beta = beta_true, censor = censor_rate, censor.cond = FALSE)
  data_train = simdata_train$data
  
  simdata_test = sim.survdata(N = SIZE_TEST, T = T_MAX, type = "none", hazard.fun = NULL,
                              num.data.frames = 1, fixed.hazard = T, knots = 8,
                              spline = TRUE, X = X_test, beta = beta_true, censor = censor_rate, censor.cond = FALSE)
  data_test = simdata_test$data
  
  
  #############################################################################################################################
  # cox model without penalty
  print("cox")
  
  fit_cox = coxph(Surv(y, failed) ~ ., data = data_train) 
  beta_tmp = fit_cox$coefficients
  beta_tmp = beta_tmp[!is.na(beta_tmp)]
  beta_cox[1:length(beta_tmp)] = beta_tmp # estimated beta from Cox
  
  # c-index 
  probs_cox = get_probs_from_beta(X_train, X_test, beta_cox)
  c_cox = get_cindex(probs_cox, data_test)
  
  y_train = fit_cox$y
  
  
  
  #############################################################################################################################
  # classical lasso
  print("lasso")
  
  fit_lasso = glmnet(as.matrix(X_train), Surv(data_train$y, data_train$failed), family = "cox")
  cv_lasso = cv.glmnet(as.matrix(X_train), Surv(data_train$y, data_train$failed), family = "cox")
  beta_tmp = as.vector(coef(fit_lasso, s = cv_lasso$lambda.min))
  beta_tmp = beta_tmp[!is.na(beta_tmp)]
  beta_lasso[1:length(beta_tmp)] = beta_tmp
  
  # c-index
  cmax_lasso = get_cv_cindex_max(fit_lasso, cv_lasso, X_test, data_test)
  
  
  #############################################################################################################################
  # ridge regression
  print("ridge")
  
  fit_ridge = glmnet(as.matrix(X_train), Surv(data_train$y, data_train$failed), family = "cox", alpha = 0)
  cv_ridge = cv.glmnet(as.matrix(X_train), Surv(data_train$y, data_train$failed), family = "cox", alpha = 0)
  beta_tmp = as.vector(coef(fit_ridge, s = cv_ridge$lambda.min))
  beta_tmp = beta_tmp[!is.na(beta_tmp)]
  beta_ridge[1:length(beta_tmp)] = beta_tmp
  
  # c-index
  cmax_ridge = get_cv_cindex_max(fit_ridge, cv_ridge, X_test, data_test)
  
  
  
  #############################################################################################################################
  # elastic net
  print("elastic net")
  
  fit_enet = glmnet(as.matrix(X_train), Surv(data_train$y, data_train$failed), family = "cox", alpha = (1/3))
  cv_enet = cv.glmnet(as.matrix(X_train), Surv(data_train$y, data_train$failed), family = "cox", alpha = (1/3))
  beta_tmp = as.vector(coef(fit_enet, s = cv_enet$lambda.min))
  beta_tmp = beta_tmp[!is.na(beta_tmp)]
  beta_enet[1:length(beta_tmp)] = beta_tmp
  
  # c-index_
  cmax_enet = get_cv_cindex_max(fit_enet, cv_enet, X_test, data_test)
  
  
  
  #############################################################################################################################
  # alasso
  # adapted from code by Hao Helen Zhang & Wenbin Lu
  print("alasso")
  source("alasso_functions.R")
  
  x = X_train
  y = data_train[, (P + 1):(P + 2)]
  y$failed = as.numeric(y$failed) * 2
  
  zz = as.matrix(x)
  n = SIZE_TRAIN
  p = P
  NN = 10
  gd = 10
  
  true.sd = sqrt(apply(zz,2,var)*(n-1))
  
  z = apply(zz,2,normalize)
  time1 = as.vector(y[,1])
  delta = as.numeric(y[,2]==2)
  
  set.seed(10)
  time=numeric(length(time1))
  
  ###breaking the ties
  time[delta==1] = time1[delta==1]+runif(sum(delta==1),0,0.1)
  time[delta==0] = time1[delta==0]+runif(sum(delta==0),0.1,0.5)
  
  iter = 100
  tol = 1.0e-10
  ps = numeric(gd)
  
  newbeta.s = matrix(0,nrow=gd,ncol=p)
  ps = numeric(gd)
  GCV = numeric(gd)
  newGCV = numeric(gd)
  
  delta = delta[order(time)]
  z = z[order(time),]
  sd = sqrt(apply(z,2,var)*(n-1))
  time = time[order(time)]
  
  beta.ini = coxph(Surv(time,delta)~z)$coef
  
  ###Computing Adaptive Lasso estimates
  newbeta.s[1,] = beta_lasso
  newGCV[1] = GCV[1]
  
  beta = newbeta.s[1,]
  sd = sd/abs(beta)
  for(j in 2 : gd){
    lam = 3^((j-1)/2) - 1.0
    beta = beta.ini
    ii = 0
    while(ii < NN){
      fn=loglik(n,delta,z,beta)
      G=dloglik(n,delta,z,beta)
      H=ddloglik(n,delta,z,beta)
      
      Xchol = try(chol(H), silent = T)
      if("try-error" %in% class(Xchol)) {
        ii = NN
        print("chol(H) error, break")
        istop = ii
        dx = 0
      } else {
        X = chol(H)
        vecY = forwardsolve(t(X),H%*%beta-G)
        lsbeta = lm(vecY~-1+X)$coef
        vecX = as.vector(X)
        
        beta1=wshoot(p,p,X,vecY,init=lsbeta,sd,lam,iter,tol)
        dx = max(abs(beta1-beta))
        ii = ii + 1
        istop = ii
        if(dx <= 1.0e-5) ii = NN
        beta = beta1
      }
    }
    newbeta.s[j,] = beta
    fn = loglik(n,delta,z,beta)
    cat("alasso ",istop,dx,fn,"\n")
    w = diag(2*abs(beta))
    ginvw = ginv(w)
    A = H + lam*sd*ginvw
    ps[j] = sum(diag(solve(A)%*%H)) - sum(newbeta.s[j,] == 0)
    newGCV[j] = fn/(n*(1-ps[j]/n)^2)
  }
  
  for(j in 1:gd){
    newbeta.s[j,] = newbeta.s[j,]/true.sd
  }
  
  beta_alasso = newbeta.s[rank(newGCV)==1,]
  
  # cindex
  cmax_alasso = get_cv_cindex_max_alasso(newbeta.s, X_train, X_test, data_test)
  
  ####################################
  # scad
  
  print("scad")
  
  
  fit_scad = ncvsurv(as.matrix(X_train), Surv(data_train$y, data_train$failed), penalty = "SCAD")
  cv_scad = cv.ncvsurv(as.matrix(X_train), Surv(data_train$y, data_train$failed))
  beta_tmp = as.vector(coef(fit_scad, cv_scad$lambda.min))
  beta_tmp = beta_tmp[!is.na(beta_tmp)]
  beta_scad[1:length(beta_tmp)] = beta_tmp
  
  # c-index
  cmax_scad = get_cv_cindex_max(fit_scad, cv_scad, X_test, data_test)
 
  #############################################################################################################################
  # graph lasso
  print("graph")
  
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
  
  # c-index
  cmax_graph = get_cv_cindex_max(fit_graph, cv_graph, as.matrix(X_tilde_test), data_test)
  
  
  
  #############################################################################################################################
  
  # l2_cox = sqrt(sum((exp(beta_cox) - beta_true)^2))
  l2_cox = sqrt(sum((beta_cox - beta_true)^2))
  l2_lasso = sqrt(sum((beta_lasso - beta_true)^2))
  l2_ridge = sqrt(sum((beta_ridge - beta_true)^2))
  l2_enet = sqrt(sum((beta_enet - beta_true)^2))
  l2_scad = sqrt(sum((beta_scad - beta_true)^2))
  l2_graph = sqrt(sum((beta_graph - beta_true)^2))
  l2_alasso = sqrt(sum((beta_alasso - beta_true)^2))
  
  l1_cox = sum(abs(beta_cox - beta_true))
  l1_lasso = sum(abs(beta_lasso - beta_true))
  l1_ridge = sum(abs(beta_ridge - beta_true))
  l1_enet = sum(abs(beta_enet - beta_true))
  l1_scad = sum(abs(beta_scad - beta_true))
  l1_graph = sum(abs(beta_graph - beta_true))
  l1_alasso = sum(abs(beta_alasso - beta_true))
  
  linf_cox = max(abs(beta_cox - beta_true))
  linf_lasso = max(abs(beta_lasso - beta_true))
  linf_ridge = max(abs(beta_ridge - beta_true))
  linf_enet = max(abs(beta_enet - beta_true))
  linf_scad = max(abs(beta_scad - beta_true))
  linf_graph = max(abs(beta_graph - beta_true))
  linf_alasso = max(abs(beta_alasso - beta_true))
  
  rpe_cox = (t(beta_cox - beta_true) %*% t(X_test) %*% X_test %*% (beta_cox - beta_true))/(SIZE_TEST)
  rpe_lasso = (t(beta_lasso - beta_true) %*% t(X_test) %*% X_test %*% (beta_lasso - beta_true))/(SIZE_TEST)
  rpe_enet = (t(beta_enet - beta_true) %*% t(X_test) %*% X_test %*% (beta_enet - beta_true))/(SIZE_TEST)
  rpe_scad = (t(beta_scad - beta_true) %*% t(X_test) %*% X_test %*% (beta_scad - beta_true))/(SIZE_TEST)
  rpe_ridge = (t(beta_ridge - beta_true) %*% t(X_test) %*% X_test %*% (beta_ridge - beta_true))/(SIZE_TEST)
  rpe_alasso = (t(beta_alasso - beta_true) %*% t(X_test) %*% X_test %*% (beta_alasso - beta_true))/(SIZE_TEST)
  rpe_graph = (t(beta_graph - beta_true) %*% t(X_test) %*% X_test %*% (beta_graph - beta_true))/(SIZE_TEST)
  
  # rpl_cox = -1000000 #coxph(Surv(y, failed) ~ ., init = beta_cox, control=list('iter.max'=0,timefix = FALSE), data = data_test)$loglik[2] 
  # rpl_lasso = coxph(Surv(y, failed) ~ ., init = beta_lasso, control=list('iter.max'=0,timefix = FALSE), data = data_test)$loglik[2]
  # rpl_enet = coxph(Surv(y, failed) ~ ., init = beta_enet, control=list('iter.max'=0,timefix = FALSE), data = data_test)$loglik[2]
  # rpl_scad = coxph(Surv(y, failed) ~ ., init = beta_scad, control=list('iter.max'=0,timefix = FALSE), data = data_test)$loglik[2]
  # rpl_ridge = coxph(Surv(y, failed) ~ ., init = beta_ridge, control=list('iter.max'=0,timefix = FALSE), data = data_test)$loglik[2]
  # rpl_alasso = coxph(Surv(y, failed) ~ ., init = beta_alasso, control=list('iter.max'=0,timefix = FALSE), data = data_test)$loglik[2]
  # rpl_graph = coxph(Surv(y, failed) ~ ., init = beta_graph, control=list('iter.max'=0,timefix = FALSE), data = data_test)$loglik[2]

  
  true_positive_index = which(beta_true!=0)
  true_negative_index = which(beta_true==0)
  epsilon = 1e-2 # the threshold for determining zero or non-zero
  
  # false positive rate: false positive/total true negatives
  fpr_cox = length(intersect(which(abs(beta_cox)>epsilon),true_negative_index))/length(true_negative_index)
  fpr_lasso = length(intersect(which(abs(beta_lasso)>epsilon),true_negative_index))/length(true_negative_index)
  fpr_ridge = length(intersect(which(abs(beta_ridge)>epsilon),true_negative_index))/length(true_negative_index)
  fpr_enet = length(intersect(which(abs(beta_enet)>epsilon),true_negative_index))/length(true_negative_index)
  fpr_scad = length(intersect(which(abs(beta_scad)>epsilon),true_negative_index))/length(true_negative_index)
  fpr_graph = length(intersect(which(abs(beta_graph)>epsilon),true_negative_index))/length(true_negative_index)
  fpr_alasso = length(intersect(which(abs(beta_alasso)>epsilon),true_negative_index))/length(true_negative_index)
  
  # false negative rate: false negative/total true positives
  fnr_cox = length(intersect(which(abs(beta_cox)<=epsilon),true_positive_index))/length(true_positive_index)
  fnr_lasso = length(intersect(which(abs(beta_lasso)<=epsilon),true_positive_index))/length(true_positive_index)
  fnr_ridge = length(intersect(which(abs(beta_ridge)<=epsilon),true_positive_index))/length(true_positive_index)
  fnr_enet = length(intersect(which(abs(beta_enet)<=epsilon),true_positive_index))/length(true_positive_index)
  fnr_scad = length(intersect(which(abs(beta_scad)<=epsilon),true_positive_index))/length(true_positive_index)
  fnr_graph = length(intersect(which(abs(beta_graph)<=epsilon),true_positive_index))/length(true_positive_index)
  fnr_alasso = length(intersect(which(abs(beta_alasso)<=epsilon),true_positive_index))/length(true_positive_index)

  
  
  l2_result = c(l2_graph, l2_lasso, l2_ridge, l2_enet, l2_scad, l2_alasso, l2_cox)
  l1_result = c(l1_graph, l1_lasso, l1_ridge, l1_enet, l1_scad, l1_alasso, l1_cox)
  linf_result = c(linf_graph, linf_lasso, linf_ridge, linf_enet, linf_scad, linf_alasso, linf_cox)
  rpe_result = c(rpe_graph, rpe_lasso, rpe_ridge, rpe_enet, rpe_scad, rpe_alasso, rpe_cox)
  rpl_result = c(rpl_graph, rpl_lasso, rpl_ridge, rpl_enet, rpl_scad, rpl_alasso, rpl_cox)
  fpr_result = c(fpr_graph,fpr_lasso, fpr_ridge, fpr_enet, fpr_scad, fpr_alasso, fpr_cox)
  fnr_result = c(fnr_graph,fnr_lasso, fnr_ridge, fnr_enet, fnr_scad, fnr_alasso, fnr_cox)
  cindex_result = c(cmax_graph, cmax_lasso, cmax_ridge, cmax_enet, cmax_scad, cmax_alasso, c_cox)
  tmp_result = c(l2_result, l1_result,linf_result,rpe_result, rpl_result, fpr_result, fnr_result, cindex_result)
  
  returevalue = list(tmp_result,rbind(beta_cox,beta_lasso,beta_ridge,beta_enet,beta_scad,beta_graph,beta_alasso))
  return(returevalue)
}
