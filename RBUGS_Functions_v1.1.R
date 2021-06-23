###############################################################
#   General Hierarchical Bayes Package (Gibbs Sampling)       #
#   3/13/2019 Kevin Lattery                                   #
#                                                             #
#   Run Everything to create environments                     #
#   Any changes to code here need to rerun everything         # 
###############################################################

# Attaching an environment multiple times creates duplicates
if (exists("env_code")) detach(env_code)
if (exists("HB_R")) detach(HB_R)
if (exists("env_GR")) detach(env_GR)
if (exists("env_out")) detach(env_out)
env_code <- new.env(parent = emptyenv())
HB_R <- new.env(parent = emptyenv())
env_GR <- new.env(parent = emptyenv())
env_out <- new.env(parent = emptyenv())

###############  Coding Functions Environment ################ 
env_code$catcode <- function(kdata, kcol, codetype = 3, varout = NULL, reflev = NULL, setna = 0, priorcode = c(0, NA)) {
  #colvec may be vector or nx1 matrix
  #codetype 1 = indicator, 2= dummy, 3 = effects
  #reflev of NULL defaults to last level
  if (is.null(varout)){
    varout <- colnames(kdata)[kcol]
    if (is.null(varout)) varout <- paste0("V", kcol)
  }
  colvec <- as.vector(as.matrix(kdata)[,kcol])
  colvec[colvec == priorcode[1]] <- priorcode[2]
  na_vals <- is.na(colvec)
  kmax <- max(colvec, na.rm = TRUE)
  colvec[na_vals] <- kmax # temp set to max
  labels_in <- sort(unique(colvec))
  newval <- match(colvec,labels_in)
  varnames <- paste(varout,labels_in, sep = "_")
  numlevs <- length(labels_in)
  if (is.null(reflev)) {reflev <- numlevs}
  if (numlevs == 1){
    outcode <- as.matrix(colvec) # no coding for 1 level
    code_matrix <- as.matrix(1)
    colnames(code_matrix) <- varnames
  } 
  if (numlevs >= 2) {
    data_re <- newval
    code_matrix <- diag(numlevs)
    colnames(code_matrix) <- varnames
    if (codetype %in% c(2, 3)) {
      code_matrix <- as.matrix(code_matrix[, -reflev]) # ref lev col dropped  
    }
    if (codetype == 3) { code_matrix[reflev,] <- -1 }
    outcode <- as.matrix(code_matrix[newval,])
  }
  outcode[na_vals,] <- setna
  if (numlevs <= 2) prior <- matrix(1)
  if (numlevs >= 3) {
    if (codetype == 1){ # ind
      prior <- diag(numlevs)
    }
    if (codetype == 2){ # dummy
      prior <- matrix(1, nrow = numlevs-1, ncol = numlevs-1) # off-diagnol
      diag(prior) <- 2
    }
    if (codetype == 3){ #effects
      prior <- matrix(-1/numlevs, nrow = numlevs-1, ncol = numlevs-1) # off-diagnol
      diag(prior) <- (numlevs -1)/numlevs
    }
  }
  return(list(outcode = outcode, code_matrix = code_matrix, vnames = varnames, prior = prior))
}

env_code$usercode <- function(kdata, kcol, varout = NULL){
  if (is.null(varout)){
    varout <- colnames(kdata)[kcol]
    if (is.null(varout)) varout <- "user_"
  }
  outcode <- as.matrix(kdata[,kcol])
  numcol <- ncol(outcode)
  if ((length(varout) == 1) & (numcol > 1)){
    varout <- paste0(varout, 1:numcol)
  }
  if (!(length(varout) == numcol)){
    varout <- paste0("user_", 1:numcol)
  }
  colnames(outcode) <- varout
  return(list(outcode = outcode, code_matrix = diag(ncol(outcode)), vnames = varout, prior = diag(ncol(outcode))))
}

env_code$ordmatrix <- function(num_levels) {
  #num_levels must be >= 3, 2 levels is one variable coded -1,1 or 0,1
  negval <- (num_levels - 1) * -1
  ord_matrix <- matrix(rep(1:(num_levels - 1), each = num_levels), nrow = num_levels)
  ord_matrix[1,] <- (negval:-1)
  for (i in 2:(num_levels - 1)) {
    negval <- negval + 1
    ord_matrix[i,] <- c(1:(i - 1), (negval:-1))
  }
  maxsum <- sum(ord_matrix[num_levels,])
  return(ord_matrix / maxsum)
}

env_code$ordcode <- function(kdata, kcol, cut_pts, varout = NULL, setna = 0) {
  # xvec must be vector
  # cut_pts must be sequential vector from low to high
  # NA and values outside cut_pts are set to "setna", default = 0
  # varout is prefix for varout_cut
  # uses function ordmatrix
  if (is.null(varout)){
    varout <- colnames(kdata)[kcol]
    if (is.null(varout)) varout <- "O"
  }
  xvec <- as.vector(as.matrix(kdata)[,kcol])
  bad <- xvec < min(cut_pts) | xvec > max(cut_pts) | is.na(xvec)
  xvec[bad] <- min(cut_pts) # temp set to min
  high <- sapply(cut_pts, function(x) xvec <= x)
  high_col <- apply(high, 1, function(x) match(TRUE, x))
  low_col <- pmax(high_col - 1,1)
  low_pt <- cut_pts[low_col] 
  high_pt <- cut_pts[high_col] 
  dist1 <- (high_pt - xvec)/(high_pt - low_pt)
  dist1[is.na(dist1)] <- 0
  dist2 <- 1 - dist1
  mat0 <- matrix(0, length(xvec), length(cut_pts))
  mat_low <- mat0
  mat_high <- mat0
  mat_low[cbind(1:length(xvec),low_col)] <- dist1 
  mat_high[cbind(1:length(xvec),high_col)] <- dist2 
  rowcode <- mat_low + mat_high
  code_matrix <- ordmatrix(length(cut_pts))
  ordcode <- round(rowcode %*% code_matrix, 5)
  ordcode[bad] <- setna # reset initial NA (default 0)
  vnames <- unlist(lapply(2:length(cut_pts), function(i) paste0(varout, cut_pts[i-1],"_", cut_pts[i])))
  colnames(ordcode) <- vnames
  varnames <- paste0(varout, "_", cut_pts)
  return(list(outcode = ordcode, code_matrix = code_matrix, vnames = varnames, prior = diag(ncol(code_matrix))))
}

env_code$back_code <- function(indcode_list){
  #att_codes is a list of codes for each attribute
  att_codes <- lapply(indcode_list, function(x) x$code_matrix)
  row_sizes <- sapply(att_codes, nrow)
  col_sizes <- sapply(att_codes, ncol)
  code_to_ind <- matrix(0, nrow = sum(row_sizes), ncol = sum(col_sizes))
  row_end <- cumsum(row_sizes)
  col_end <- cumsum(col_sizes)
  row_start <- c(1, row_end[-length(row_end)] + 1)
  col_start <- c(1, col_end[-length(col_end)] + 1)
  for (i in 1:length(att_codes)){
    code_to_ind[(row_start[i]:row_end[i]), (col_start[i]:col_end[i])] <- att_codes[[i]]   
  }
  rownames(code_to_ind) <- do.call(c, lapply(indcode_list, function(x) x$vnames))
  return(code_to_ind)  
}

env_code$get_prior <- function(indcode_list){
  att_codes <- lapply(indcode_list, function(x) x$prior)
  row_sizes <- sapply(att_codes, nrow)
  col_sizes <- sapply(att_codes, ncol)
  result <- matrix(0, nrow = sum(row_sizes), ncol = sum(col_sizes))
  row_end <- cumsum(row_sizes)
  col_end <- cumsum(col_sizes)
  row_start <- c(1, row_end[-length(row_end)] + 1)
  col_start <- c(1, col_end[-length(col_end)] + 1)
  for (i in 1:length(att_codes)){
    result[(row_start[i]:row_end[i]), (col_start[i]:col_end[i])] <- att_codes[[i]]   
  }
  return(result)  
}

env_code$make_codefiles <- function(indcode_list){
  # Makes global vars:
  # code_master, indcode, indprior
  .GlobalEnv$indprior <- get_prior(indcode_list) # of levels effect
  .GlobalEnv$indcode <- do.call(cbind, lapply(indcode_list, function(x) x$outcode)) # coded variables 
  .GlobalEnv$code_master <- back_code(indcode_list)
  colnames(.GlobalEnv$code_master) <- colnames(.GlobalEnv$indcode) # CHECK code_master
}

################################
#      Gelman Rubin            #
################################
env_GR$split_iter <- function(i, iter_hist = iter_hist) {
  post <- iter_hist[[i]] # drop first column with iteration number
  n_split <- floor(nrow(post) / 2)
  post1 <- post[1:n_split,,drop = FALSE]
  post2 <- post[(n_split + 1):(2 * n_split),,drop = FALSE]
  return(list(post1 = post1, post2 = post2))
}

env_GR$GelRubin <- function(chain, col_num) {
  onevar <- do.call(cbind, lapply(chain, function(x) x[, col_num]))
  krow <- nrow(onevar)
  kcol <- ncol(onevar)
  ch_means <- colMeans(onevar, na.rm = TRUE)
  ch_means_m <- mean(ch_means)
  B <- sum((ch_means - ch_means_m) ^ 2) * krow / (kcol - 1)
  ch_var <- apply(onevar, 2, var, na.rm = TRUE)
  W <- mean(ch_var)
  W_adj <- (W * (krow - 1) + B) / krow
  R <- sqrt(W_adj / W)
  return(R)
}


env_GR$GR <- function(iter_hist, vnames = NULL) {
  # iter_hist is list, each element is matrix/df of iterations from separate chain
  
  # Check dimensions
  kdim <- dim(iter_hist[[1]])
  check <- sapply(iter_hist, function(x) sum(!(dim(x) == kdim)))
  if (max(check) > 0){
    message("FAILURE: Data File dimensions of chains do not match")
  } else {
    
    if (is.null(vnames)){
      vnames <- colnames(iter_hist[[1]])
      if (is.null(vnames)) vnames <- paste0("V", 1:ncol(iter_hist[[1]]))
    }
    chain_split <- lapply(1:length(iter_hist), function(x) split_iter(x, iter_hist))
    chain <- list()
    k <- 1
    for (i in 1:length(chain_split)) {
      for (j in 1:length(chain_split[[i]])) {
        chain[[k]] <- chain_split[[i]][[j]]
        k <- k + 1
      }
    }
    Rhat <- sapply(1:ncol(chain[[1]]), function(x) GelRubin(chain, x))
    Conv_R <- as.data.frame(vnames)
    Conv_R$Rhat <- Rhat
    colnames(Conv_R)[2] <- "Rhat"
    Rhat_sort <- Rhat[order(-Rhat)]
    hline <- max(min(Rhat_sort), 1.2)
    ymax <- max(c(Rhat, 1.3)) # 1.2 + cushion
    ymin <- min(c(Rhat, 1))
    plot(1:length(Rhat_sort), Rhat_sort, pch = 19, col = "seagreen", lab = c(10, 8, 3), las = 1, xlab = "MCMC Variables (Sorted by R)", grid(col = "lightgray", lty = "dotted"),
         ylim = c(ymin,ymax))
    abline(h = hline, col = "red")
    bad <- sum(Rhat > 1.2)
    text(x = 1, y = hline, labels = paste0("     ", bad, " bad"), pos = 3, col = "red")
    print(Conv_R)
    return(list(Rhat = Conv_R, chain = chain))
  }
}

env_GR$Chain_Plot <- function(chain, col_num, colors = c("blue3", "tan", "red", "seagreen")) {
  # chain is list where elements 1:2 are 1st chain, 3:4 2nd chain, etc
  chain_plot <- NULL
  vname <- colnames(chain[[1]])[col_num]
  for (i in 1:(length(chain) / 2)) {
    x <- (2 * i) - 1
    chain_plot <- cbind(chain_plot, c(chain[[x]][, col_num], chain[[x + 1]][, col_num]))
  }
  xvec <- 1:nrow(chain_plot)
  plot(xvec, chain_plot[, 1], type = "p", col = colors[1], cex = .5, xlab = "Post Burn-In Iteration", ylab = paste0(" Col", col_num, ": ", vname), ylim = c(min(chain_plot), max(chain_plot)))
  for (i in 2:(ncol(chain_plot))) {
    lines(xvec, chain_plot[, i], col = colors[i], cex = .5)
  }
}

# Below is graphing function I am working on.  Would like to dsiplay multiple charts in one window.  
env_GR$Chain_Plot_Cols <- function(chain, col_nums) {
  # Works only with 4 columns
  kcol <- floor(sqrt(length(col_nums)))
  krow <- length(col_nums) - kcol
  n_zeros <- krow * kcol - length(col_nums)
  if (n_zeros > 0) col_nums <- c(col_nums,rep(0, n_zeros))
  layout(matrix(col_nums, krow, kcol, byrow = TRUE))
  lapply(col_nums, function(i) Chain_Plot(chain, col_nums[i]))
  layout(1)
}

env_GR$plot_bad <- function(GR_out, threshold = 1.2){
  items_plot <- which(GR_out$Rhat[,2] > threshold)
  if (length(items_plot) == 0){
    cat(paste0("No items to plot above threshold of ", threshold))
    cat("\n")
    cat(paste0("Highest Rhat value is ", round(max(GR_out$Rhat[,2]),4)))
  } else{
    lapply(items_plot, function(x)
      Chain_Plot(GR_out$chain, x))
  } 
}

env_GR$prepdata <- function(kdata, dropcol = NULL, keeplast = 10000){
  # Converts columns to numeric
  # Any row with NA is removed
  keepcol <- TRUE
  if (!is.null(dropcol) & dropcol <= ncol(kdata)) keepcol <- -dropcol
  result <- suppressWarnings(as.data.frame(sapply(kdata[,keepcol], as.numeric)))
  rowkeep <- (rowSums(is.na(result)) == 0)
  result <- result[rowkeep,]
  badrows <- nrow(kdata) - nrow(result)
  if (badrows > 0) {
    cat(paste0(nrow(kdata) - nrow(result) , " rows with strings removed"))
    cat("\n")
  }
  keeplast <- min(keeplast, nrow(result))
  result <- tail(result,keeplast)
  return(result)
}


HB_R$PredProb_MNL <- function(data_list, pred_env){
  # Basic MNL 
  # Assumes matrix of betas
  # data_list has elements ind, match_id, idtask_r
  betamatch <- pred_env$betas[data_list$match_id,] # Map unit betas to rows of data
  V <- rowSums(data_list$ind * betamatch)
  V[V >500] <- 500; V[V < -500] <- -500
  U <- exp(V) #Utility   
  tasksum <- rowsum(U, data_list$idtask_r) # Sum U each task
  esum <- (tasksum[data_list$idtask_r,]) # Map Sum to rows
  predprob <- U/esum
  predprob[is.na(predprob) | (predprob < 9.88e-324)] <- 9.88e-324
  return(predprob)
}

HB_R$LogLike_id <- function(data_list, pred_env) { # Standard likelihood
  predprob <- pred_env$pred_func(data_list)
  LL_id <- rowsum(log(predprob) * data_list$dep * data_list$wts, data_list$match_id)
  return(as.vector(LL_id))
}

HB_R$OLS_LL <- function(pred, dep){
  # Convert OLS to LogLikelihood  
  # (dep - pred) ~ N(0, sigma(dep-pred))
  diff <- (dep - pred) ^2
  sigma <- sqrt(mean(diff))
  prefix <- 1/(sigma * sqrt(2*pi))
  inside <- -diff/(2 * sigma^2) # e^inside is part of normal
  ll_row <- log(prefix) + inside 
  return(ll_row)
}
HB_R$Norm_DataCov <- function(data_list, model_i){
  # alpha ~ N(mean_data, priorcov/n)
  mean_data <- colMeans(get(model_env$model_methods[[model_i]]$prior$data, model_env))
  npar <- length(mean_data)
  ndata <- nrow(get(model_env$model_methods[[model_i]]$prior$data, model_env))
  new_alpha <- mean_data +
    t(chol(get(model_env$model_methods[[model_i]]$prior$cov, model_env)/sqrt(ndata))) %*%
    matrix(rnorm(npar), nrow = npar, ncol = 1) 
  assign(model_env$model_methods[[model_i]]$par, as.vector(new_alpha), envir = model_env)  
}

HB_R$IW_Cov <- function(data_list, model_i){
  # cov ~ IW(DF, priorcov)
  ab_diff <- sweep(get(model_env$model_methods[[model_i]]$prior$data, model_env),
                   2, get(model_env$model_methods[[model_i]]$prior$alpha, model_env), "-")
  npar <- ncol(ab_diff)
  ndata <- nrow(ab_diff)
  df <- model_env$model_methods[[model_i]]$prior$df
  # H <- get(model_env$model_methods[[model_i]]$prior$cov, model_env) * (npar + df) + (t(ab_diff) %*% ab_diff)
  H <- get(model_env$model_methods[[model_i]]$prior$cov, model_env) * (npar + ndata + df) + (t(ab_diff) %*% ab_diff)
  T <- t(chol(solve(H)))
  u <- matrix(rnorm(npar * (ndata + npar + df)), 
              nrow = npar, ncol = ndata + npar + df)
  S_Inv <- solve((T %*% u) %*% t(T %*% u)) # T %*% u are vectors from N(0, Inv(H))
  assign(model_env$model_methods[[model_i]]$par, S_Inv, envir = model_env)
}


HB_R$Prior_MVN <- function(model_i, par_now, par_try){
  # returns log(prior_try/prior_now)

  alpha <- get(model_env$model_methods[[model_i]]$prior_data$alpha, model_env)
  par_tryz <- sweep(par_try, 2, alpha, "-", FALSE)
  par_nowz <- sweep(par_now, 2, alpha, "-", FALSE)
  cov_inv <- solve(get(model_env$model_methods[[model_i]]$prior_data$cov, model_env))
  
  LL_rel_dens <- -0.5 * (colSums(t(par_tryz) * (cov_inv %*% t(par_tryz)))) -
                 -0.5 * (colSums(t(par_nowz) * (cov_inv %*% t(par_nowz))))
    
  return(LL_rel_dens) # log(try) - log(now)
}

HB_R$Proposal_RW <- function(model_i, par_now){
  # Random walk with N(0, prop_data$cov)
  ncols <- ncol(par_now); nrows <- nrow(par_now) # For convenience
  par_try <- par_now + (model_env$model_methods[[model_i]]$prop_data$rho *
                      (matrix(rnorm(nrows * ncols), nrow = nrows, ncol = ncols) %*%
                         chol(get(model_env$model_methods[[model_i]]$prop_data$cov, model_env))))
  # No correction for symmetric
  return(list(par_try = par_try,
              log_correction = matrix(0, nrow = nrows, ncol = 1)
              ))
}

HB_R$Proposal_Bound_old <- function(model_i, par_now){
  # Uses beta distribution - replaced this with truncated normal
   # par_now is 1 item or vector or 1 column
   #model_i must have $prop_data$rho, lb, ub as below
  
  
  rho <- model_env$model_methods[[model_i]]$prop_data$rho
  lb <- model_env$model_methods[[model_i]]$prop_data$bound[1]
  ub <- model_env$model_methods[[model_i]]$prop_data$bound[2]

  par_now_unit <- (par_now - lb)/(ub - lb) # Scale to [0,1]
  par_now_unit[par_now_unit > 1] <- .999
  par_now_unit[par_now_unit < 0] <- .001
  
  jump_inv <- 100/rho # Bigger rho means smaller jump_inv = bigger jump 
  jump_inv <- min(max(jump_inv, 10), 500) # jump_inv of 10 is pretty big jump 
  
  par_try_unit <- sapply(as.vector(par_now_unit), function(x) rbeta(1, x * jump_inv + 1, (1-x) * jump_inv + 1))
  par_try <- lb + (par_try_unit * (ub - lb)) # Scaled proposal = rescaled unit proposal
  
  # Hastings correction factor on unit scale
  try_given_now <- sapply(1:length(par_try_unit), function(i) dbeta(par_try_unit[i],
                                                                    par_now_unit[i] * jump_inv + 1,
                                                                    (1-par_now_unit[i]) * jump_inv + 1
  ))
  now_given_try <- sapply(1:length(par_now_unit), function(i) dbeta(par_now_unit[i],
                                                                    par_try_unit[i] * jump_inv + 1,
                                                                    (1-par_try_unit[i]) * jump_inv + 1
  ))
  
  return(list(par_try = as.matrix(par_try),
              log_correction = as.matrix(log(now_given_try) - log(try_given_now))
  ))
}

###   Truncated Normal Functions ###############
HB_R$rnormtrunc <- function(n, mu, cov, lb, ub){
  # n draws from truncated normal
  sigma <- sqrt(min(cov, 1000))
  FA <- pnorm(((lb - mu)/sigma))
  FB <- pnorm(((ub - mu)/sigma))
  result <- mu + sigma * qnorm(runif(n) * (FB - FA) + FA)
  
  # Should not need these but doing for safety
  result[result > ub] <- ub - .001
  result[result < lb] <- lb + .001 
  
  return(result)  
}


HB_R$dnormtrunc <- function(x, mu, cov, lb, ub){
  # density for truncated normal
  # Not usually needed
  sigma <- sqrt(min(cov, 1000))
  t_scale <- pnorm(ub, mu, sigma) -  pnorm(lb, mu, sigma) 
  result <- dnorm(x, mu, sigma)/t_scale
  #result[x < lb] <- 0 # Return something even if input is bad
  #result[x > ub] <- 0  
  return(result)
}

HB_R$Proposal_RW_NT <- function(model_i, par_now){
  # Normal Truncated Proposal
  # par now is 1 item or vector or 1 column
  # model_i must have $prop_data$rho, lb, ub as below
  
  rho <- model_env$model_methods[[model_i]]$prop_data$rho
  rho <- min(rho, 1000)
  lb <- model_env$model_methods[[model_i]]$prop_data$bound[1]
  ub <- model_env$model_methods[[model_i]]$prop_data$bound[2]
  
  par_try <- sapply(as.vector(par_now), function(x) HB_R$rnormtrunc(1, x, 1*rho, lb, ub))
  
  # Hastings correction factor
  try_given_now <- sapply(1:length(par_try), function(i) HB_R$dnormtrunc(par_try[i],
                                                                    par_now[i],
                                                                    rho, lb, ub
  ))
  now_given_try <- sapply(1:length(par_now), function(i) HB_R$dnormtrunc(par_now[i],
                                                                    par_try[i],
                                                                    rho, lb, ub
  ))
  
  return(list(par_try = as.matrix(par_try),
              log_correction = as.matrix(log(now_given_try) - log(try_given_now))
  ))
}

HB_R$NormTrunc_DataCov <- function(data_list, model_i){
  # alpha ~ Truncated Normal(mean_data, priorcov/n)
  mean_data <- colMeans(get(model_env$model_methods[[model_i]]$prior$data, model_env))
  npar <- length(mean_data)
  ndata <- nrow(get(model_env$model_methods[[model_i]]$prior$data, model_env))
  lb <- model_env$model_methods[[model_i]]$bound[1]
  ub <- model_env$model_methods[[model_i]]$bound[2]
  cov <- get(model_env$model_methods[[model_i]]$prior$cov, model_env)
  new_alpha <- rnormtrunc(1, mean_data, sqrt(cov/ndata), lb, ub)
  assign(model_env$model_methods[[model_i]]$par, as.vector(new_alpha), envir = model_env)  
}
#####  End Truncated Normal ###################

HB_R$MH <- function(data_list, model_i){
  # model_env for prior
  # pred_env for prediction/LL
  
  LL_now <- pred_env$LL_func(data_list, pred_env) # LL function call
  # Need this because of changes to par_env since last loop
  
  par_name <- model_env$model_methods[[model_i]]$par
  par_now <- get(par_name, model_env)
  
  propose <-  model_env$model_methods[[model_i]]$prop_func(model_i, par_now) # Proposal function
  # list(par_try, ln(correction))
  par_try <- propose[[1]]
  
  prior_LL <- model_env$model_methods[[model_i]]$prior_func(model_i, par_now, par_try) # Prior/Density function
  # ln(prior_try/prior_now )
  
  par_pred_copy <- get(par_name, pred_env) # Copy of constrained data
    # Now modify parameter in pred_env with constraints/bounds of par_try
  if(!is.null(model_env$model_methods[[model_i]]$constraints)){
    assign(par_name, con_func(par_try, model_env$model_methods[[model_i]]$constraints), envir = pred_env) 
  } else  assign(par_name, par_try, envir = pred_env) 
 
  LL_try <- pred_env$LL_func(data_list, pred_env) # LL function call
  
  #r_log <-  (LL_try + prior_LL[[2]]) - (LL_now + prior_LL[[1]])
  r_log <- (LL_try - LL_now) + prior_LL + propose[[2]]
           #ln(try/now) + ln(d_try/d_now) + ln(correction)
  accept <- r_log >= log(runif(length(r_log)))
  
  model_env$model_methods[[model_i]]$accept$current <- mean(accept) # Update acceptance rate
  #model_env$model_methods[[model_i]]$LL_now[accept] <- LL_try[accept]
  model_env$LL_now <- sum(LL_now)
  
  eval(parse(text = 
               paste0("model_env$",par_name, "[accept,] <- par_try[accept,]")
  )) # model_env gets modified by accepted values
  eval(parse(text = 
               paste0("pred_env$",par_name, "[!accept,] <- par_pred_copy[!accept,]")
  )) #model_par was modified, make unaccepted original
  
  #update rho
  if (model_env$model_methods[[model_i]]$prop_data$update_rho) {
    accept <- mean(accept)
    target <- model_env$model_methods[[model_i]]$accept$target
    if (abs(target - accept) > .05){
      model_env$model_methods[[model_i]]$prop_data$rho <- model_env$model_methods[[model_i]]$prop_data$rho * qnorm(target/2)/qnorm(accept/2)
      # https://support.sas.com/documentation/cdl/en/statug/63033/HTML/default/viewer.htm#statug_mcmc_sect022.htm
    }
  } 
}

HB_R$prep_file <- function(idtaskdep, ind) {
  sort_order <- order(idtaskdep[, 1], idtaskdep[, 2])
  ind <- as.matrix(ind[sort_order,])
  ind[is.na(ind)] <- 0
  dep <- idtaskdep[sort_order, 3]
  idtask <- idtaskdep[sort_order, 1:2]
  idtask_u <- unique(idtask)
  idtask_r <- (match(data.frame(t(idtask)), data.frame(t(idtask_u)))) # unique tasks
  resp_id <- as.vector(unique(idtask_u[, 1]))
  match_id <- match(idtask[, 1], as.matrix(resp_id))
  # Next 3 lines recodes dep to sum to 1
  depsum <- rowsum(dep, idtask_r) # Sum dep each task
  depsum_match <- as.matrix(depsum[idtask_r,]) # Map Sum to rows
  dep <- dep / depsum_match # sum of dep will add to 1
  dep[is.na(dep)] <- 0
  wts <- rep(1, nrow(dep))
  return(list(idtask = idtask, dep = dep, ind = ind, idtask_r = idtask_r, resp_id = resp_id, match_id = match_id, wts = wts))
}

HB_R$prep_file_MD <- function(idtask, best, worst, ind) {
  sort_order <- order(idtask[, 1], idtask[, 2])
  ind <- as.matrix(ind[sort_order,])
  ind[is.na(ind)] <- 0
  dep1 <- as.matrix(best[sort_order])
  dep2 <- as.matrix(worst[sort_order])
  idtask <- idtask[sort_order, 1:2]
  idtask_u <- unique(idtask)
  idtask_r <- (match(data.frame(t(idtask)), data.frame(t(idtask_u)))) # unique tasks
  resp_id <- as.vector(unique(idtask_u[, 1]))
  match_id <- match(idtask[, 1], as.matrix(resp_id))
  wts <- rep(1, nrow(ind))
  return(list(idtask = idtask, dep1 = dep1, dep2 = dep2, ind = ind, idtask_r = idtask_r, resp_id = resp_id, match_id = match_id, wts = wts))
}

HB_R$con_func <- function(x, constraints){
  # x must be a matrix
  
  # Do relative constraints with tying
  if (!is.null(constraints$rel)){
    nBad <- Inf
    while (nBad > 0) {
      nBad <- 0
      for (i in 1:nrow(constraints$rel)) {
        par1 <- constraints$rel[i,1]; par2 <- constraints$rel[i,3]
        if (constraints$rel[i,2] == 1) {
          bad <- x[, par1] < x[, par2]
          x[bad, par1] <- x[bad, par2]
        } else {
          bad <- x[, par1] > x[, par2]
          x[bad, par1] <- x[bad, par2]
        }
        nBad <- nBad + sum(bad)
      } # end for
    } # end while
  }
  # Do pos/neg constraints
  if (!is.null(constraints$bounds)){
    for (i in 1:length(constraints$bounds)){
      lb <- constraints$bounds[[i]][1]
      ub <- constraints$bounds[[i]][2]
      cols <- constraints$bounds[[i]][-1:-2]
      x[,cols] <- lb + (ub - lb)/(1 + exp(-1* x[,cols]))      
    }
  }
  return(x)
}

# Status Update Reports
HB_R$report_base <- function (data_list){
  if (!is.null(dev.list())) dev.off()
  dev.new()
  tpi <- (Sys.time() - model_env$time_beg) /model_env$iter
  tleft <- tpi * (model_env$tot_iter - model_env$iter)
  units(tleft) <- "mins"
  model_env$LL_Hist[,1] <- model_env$LL_Hist[,1] / 1000
  LLPlot <- model_env$LL_Hist
  LLPlot[,1] <- model_env$LL_Hist[,1]/1000
  LLPlot[,2] <- model_env$LL_Hist[,2]/model_env$scale_y
  ymin <- round(min(LLPlot[,2]) * 1.05, 1)
  ymax <- round(max(LLPlot[,2]) * .95, 1)
  interval <- max(round((ymax - ymin)/4,1),.5)
  ymax <- ymin + (4 * interval)
  plot(x = 0, y = 0, type = "n", main = paste0("MCMC Fit -- Time Left: ", format(tleft, digits = 2)),
       xlim = c(0,model_env$tot_iter/1000), xlab = "Iterations x 1000", 
       ylim = c(ymin, ymax), ylab = paste0("Log Likelihood / ",model_env$scale_y), axes = FALSE, cex = 0.5)
  #mtext(paste0("Time to completion: ", format(tleft, digits = 3)), side = 3)
  segments(hb_control$iter_burn/1000, ymin, hb_control$iter_burn/1000, ymax, col = "red", 
           lty = 2, lwd = 2)
  segments(0, 0, model_env$tot_iter, 0, col = "gray", lty = 1, lwd = 1)
  axis(1, at = seq(from = 0, to = model_env$tot_iter/1000, by = (model_env$tot_iter)/10000))
  axis(2, at = seq(from = ymin, to = ymax, by = interval))
  lines(LL_Hist, type = 'l', col = 'blue', cex = .5)
}

#################################################################
#   MAIN PROGRAM
#################################################################

HB_EST <- function(data_list, pred_func, LL_func, model_pars, model_methods, hb_control, ContinueOld = FALSE){
  # Create environments with and Initial Data 
  if (!ContinueOld) { #Setup model_env
    .GlobalEnv$model_env <- list2env(c(model_pars$par_func, model_pars$par_others), parent = baseenv())
    model_env$model_methods <- model_methods
    model_env$singlepar <- (length(model_pars$par_func) == 1)
    
    .GlobalEnv$pred_env <- list2env(model_pars$par_func, parent = baseenv())
    pred_env$pred_func <- pred_func
    pred_env$LL_func <- LL_func
    environment(pred_env$pred_func) <- pred_env # set environment of function
    for (i in 1:length(model_methods)){ # Make pars in par_func constrained
      if (model_methods[[i]]$par %in% ls(pred_env)){
        if(!is.null(model_methods[[i]]$constraints)){
          par_name <- model_methods[[i]]$par
          par_new <- con_func(get(par_name, pred_env), model_env$model_methods[[i]]$constraints)
          assign(par_name, par_new, envir = pred_env)
        }
      }
    }
  } # End of setup New
  model_env$LL_now <- sum(pred_env$LL_func(data_list, pred_env))
  model_env$scale_y <- round(model_env$LL_now/-90,-1) # scale for y-axis plot
  model_env$LL_Hist <- c(1, model_env$LL_now)
  model_env$tot_iter <- hb_control$iter_burn + hb_control$iter_sample
  
  # Set up storage
  track_detail <- sapply(hb_control$track, function(x){
    if (x %in% names(model_pars$par_func)) {
      result <- paste0("pred_env$",x)
    } else
      if (x %in% names(model_pars$par_others)) {
        result <- paste0("model_env$",x)
      } else result <- x
      return(result)
  })
  rel_iter <- 1
  model_env$track <- vector("list", length(track_detail))
  names(model_env$track) <- hb_control$track
  
  model_env$time_beg <- Sys.time()  
  for (iter in 1:model_env$tot_iter){
    ############ UPDATES  ############################
    model_env$iter <- iter
    for (i in 1:length(model_env$model_methods)){
      do.call(model_env$model_methods[[i]]$method, list(data_list, i))      
    }
    
    #if (iter == hb_control$iter_burn) model_betas$prop_data$update_rho == FALSE # Keep rho constant
    
    ############ Storage  #######################        
    if (iter > hb_control$iter_burn){
      if ((iter - hb_control$iter_burn -1)%%hb_control$thin == 0) {
        for (k in 1:length(hb_control$track)){
          eval(parse(text =
                       paste0("model_env$track[[k]][[rel_iter]] <- ", track_detail[[k]])               
          ))
        }
        rel_iter <- rel_iter + 1 
      }
    }
    ############# Report #########################
    if ((iter%%hb_control$report_interval == 0) | (iter == 30)) {
      model_env$LL_Hist <- rbind(model_env$LL_Hist, c(iter,model_env$LL_now))
      hb_control$report_function(data_list)
    }
  } 
  elapsed <- Sys.time() - model_env$time_beg
  print(elapsed)
  close.screen(all = TRUE)
}


attach(env_code)
attach(HB_R)
attach(env_GR)
attach(env_out)

