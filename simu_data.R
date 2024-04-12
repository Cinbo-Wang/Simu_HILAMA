
# Simu 1 ------------------------------------------------------------------
## generate parameters----------------------------------------------
generate_param <- function(p,q,eta,rho,s,r_p,r_pq,r_h=1){
  Sigma_E <- matrix(0,nrow=p,ncol=p)
  for(i in 1:(p-1)){
    for(j in (i+1):p){
      Sigma_E[i,j] <- rho**abs(i-j)
    }
  }
  Sigma_E <- Sigma_E + t(Sigma_E) + diag(p)
  
  Psi <- matrix(rnorm(n=s*p,mean=0,sd=1),
                nrow=s,byrow = T)
  
  if((p-ceiling(p*r_h))>0){
    for(i in 1:s){
      loc_active <- sample(1:p,size=ceiling(p*r_h),replace = F)
      loc_zero <- setdiff(1:p,loc_active)
      Psi[i,loc_zero] <- 0
    }
  }
  
  Theta <- matrix(0,nrow=p,ncol=q)
  row_active <- sort(sample(1:p,size=ceiling(p*r_p),replace = F))
  col_active_ls <- vector('list',length(row_active))
  for(i in 1:length(row_active)){
    col_active_ls[[i]] <- sample(1:q,size=sample(5:20,size=1))
  }
  col_active <- sort(unique(unlist(col_active_ls)))
  for(i in 1:length(row_active)){
    row_loc <- row_active[i]
    for(col_loc in col_active_ls[[i]]){
      Theta[row_loc,col_loc] <- rnorm(1,0.8,sqrt(0.1))*sample(c(1,-1),1)
    }
  }
  
  Phi <- matrix(rnorm(s*q,mean=eta,sd=1)*sample(c(1,-1),size=s*q,replace = T),
                nrow=s,byrow=T)
  if((q-ceiling(q*r_h))>0){
    for(i in 1:s){
      loc_active <- sample(1:q,size=ceiling(q*r_h),replace = F)
      loc_zero <- setdiff(1:q,loc_active)
      Phi[i,loc_zero] <- 0
    }
  }
  
  col_silence <- setdiff(1:q,col_active)
  s10 <- min(length(col_silence),ceiling(q*r_pq*0.2))
  if(s10>0){
    med_loc_false <- sample(col_silence,size=s10)
  }else{
    med_loc_false <- NULL
  }
  col_active_prob <- table(unlist(col_active_ls)) / length(unlist(col_active_ls))
  med_loc_true <- sample(as.numeric(names(col_active_prob)),size=q * r_pq - s10,prob = col_active_prob)
  beta <- rep(0,q)
  for(i in c(med_loc_true,med_loc_false)){
    beta[i] <-  rnorm(1,0.8,sqrt(0.1))*sample(c(1,-1),1) 
  }
  
  exp_loc_true <- sample(1:p,ceiling(p*r_pq))
  gamma  <- rep(0,p);
  for(i in exp_loc_true){
    gamma[i] <- rnorm(1,0.8,sqrt(0.1))*sample(c(1,-1),1) 
  }
  
  phi <- rnorm(s,eta,1)*sample(c(1,-1),size=s,replace = T)
  
  return(list(
    Psi = Psi, Sigma_E=Sigma_E,Theta=Theta,Phi=Phi,
    gamma = gamma,beta=beta,phi = phi
  ))
  
}

rho_vector <- c(0.4,0.8)
p_vector <- c(100,200,400)
q_vector <- c(50,100,200)

eta_vector <- c(0.5,1.5)
s <- 5; r_p <- r_pq <- 0.1; r_h <- 1

rep_data_comb_list <- list()
comb_num_df <- c()

iter_num <- 1
comb_param_ls <- list()
param_comb_df <- NULL
for(p in p_vector ){
  for(rho in rho_vector){
    for(q in q_vector){
      for(eta in eta_vector){
        param_ls <- generate_param(p,q,eta,rho,s,r_p,r_pq,r_h)
        comb_param_ls[[iter_num]] <- param_ls

        param_comb_df <- rbind(param_comb_df,
                               c(iter_num,rho,eta,p,q,s,r_p,r_pq,r_h))
        iter_num <- iter_num + 1
        if(iter_num %% 10 ==0) cat('iter num:',iter_num,'\n')
      }
    }
  }
}


colnames(param_comb_df) <- c('iter_num','rho','eta','p','q','s','r_p','r_pq','r_h')
param_comb_df <- as.data.frame(param_comb_df)

comb_param_ls$param_comb_df <- param_comb_df

saveRDS(comb_param_ls,file='../simuData/comb_param_ls_simu1.rds')

## generate datasets----------------
generate_dataset <- function(n,param_ls){
  s <- nrow(param_ls[['Psi']]); p <- ncol(param_ls[['Psi']]); q <- ncol(param_ls[['Theta']])
  
  H <- matrix(rnorm(n*s,0,1),nrow=n,byrow=T)
  E <- MASS::mvrnorm(n=n,mu=rep(0,p),Sigma = param_ls[['Sigma_E']])
  X <- H %*% param_ls[['Psi']] + E
  
  E_tilde <- matrix(rnorm(n*q,0,1),nrow=n,byrow = T)
  M <- X %*% param_ls[['Theta']] + H %*% param_ls[['Phi']] + E_tilde
  
  e <- matrix(rnorm(n,0,1),nrow=n)
  Y <- X %*% param_ls[['gamma']] + M %*% param_ls[['beta']] + H %*% param_ls[['phi']]+e
  
  return(list(
    X = X, M = M, Y = Y, Z = H[,1:3],H = H[,4:5],
    Theta=param_ls[['Theta']],beta=param_ls[['beta']],
    gamma =param_ls[['gamma']]
  ))
}
# 产生数据
n_vector <- c(200,400)
comb_param_ls <- readRDS('../simuData/comb_param_ls_simu1.rds')
data_comb_df <- NULL
rep_data_comb_list <- list()
for(rep_num in 1:100){
  data_comb_list <- list()
  data_num <- 1
  for(iter_num in 1:(length(comb_param_ls)-1)){
    param_ls <- comb_param_ls[[iter_num]]
    param_vec <- as.numeric(comb_param_ls$param_comb_df[iter_num,])
    
    for(n in n_vector){
      if(rep_num == 1){
        data_comb_df <- rbind(data_comb_df,c(data_num,n,param_vec))
      }
      data_ls <- generate_dataset(n,param_ls)
      
      data_comb_list[[data_num]] <- data_ls
      data_num <- data_num + 1
    }
    
  }
  if(rep_num %% 10 == 0) cat('rep num:',rep_num,'\n')
  rep_data_comb_list[[rep_num]] <- data_comb_list
}
data_comb_df <- as.data.frame(data_comb_df)
colnames(data_comb_df) <- c('data_num','n',colnames(comb_param_ls$param_comb_df))

rep_data_comb_list[['comb_type']] <- data_comb_df

saveRDS(rep_data_comb_list,file='../simuData/rep_comb_data_ls_simu1.rds')



# Simu 2 --------------------------------------
## generate parameters----------------------------------------------

generate_param <- function(p,q,eta,rho,s,r_p,r_pq,r_h,p_z=0){
  # r_h ：proportion of variables affected by hidden confounder
  
  Sigma_E <- matrix(0,nrow=p,ncol=p)
  for(i in 1:(p-1)){
    for(j in (i+1):p){
      Sigma_E[i,j] <- (rho)**abs(i-j)
    }
  }
  Sigma_E <- Sigma_E + t(Sigma_E) + diag(p)
  
  Psi <- matrix(rnorm(n=s*p,mean=0,sd=1),
                nrow=s,byrow = T)
  
  if((p-ceiling(p*r_h))>0){
    for(i in 1:s){
      loc_active <- sample(1:p,size=ceiling(p*r_h),replace = F)
      loc_zero <- setdiff(1:p,loc_active)
      Psi[i,loc_zero] <- 0
    }
  }
  
  Psi_z <- matrix(rnorm(n=p_z*p,mean=0,sd=1),
                  nrow=p_z,byrow = T)
  
  Theta <- matrix(0,nrow=p,ncol=q)
  row_active <- sort(sample(1:p,size=ceiling(p*r_p),replace = F))
  col_active_ls <- vector('list',length(row_active))
  for(i in 1:length(row_active)){
    col_active_ls[[i]] <- sample(1:q,size=sample(5:50,size=1))
  }
  col_active <- sort(unique(unlist(col_active_ls)))
  for(i in 1:length(row_active)){
    row_loc <- row_active[i]
    for(col_loc in col_active_ls[[i]]){
      Theta[row_loc,col_loc] <- rnorm(1,0.8,sqrt(0.1))*sample(c(1,-1),1)
    }
  }
  
  Phi <- matrix(rnorm(s*q,mean=eta,sd=1)*sample(c(1,-1),size=s*q,replace = T),
                nrow=s,byrow=T)
  if((q-ceiling(q*r_h))>0){
    for(i in 1:s){
      loc_active <- sample(1:q,size=ceiling(q*r_h),replace = F)
      loc_zero <- setdiff(1:q,loc_active)
      Phi[i,loc_zero] <- 0
    }
  }
  
  Phi_z <- matrix(rnorm(p_z*q,mean=eta,sd=1)*sample(c(1,-1),size=p_z*q,replace = T),
                  nrow=p_z,byrow=T)
  
  
  col_silence <- setdiff(1:q,col_active)
  s10 <- min(length(col_silence),ceiling(q*r_pq*0.2))
  if(s10>0){
    med_loc_false <- sample(col_silence,size=s10)
  }else{
    med_loc_false <- NULL
  }
  col_active_prob <- table(unlist(col_active_ls)) / length(unlist(col_active_ls))
  med_loc_true <- sample(as.numeric(names(col_active_prob)),size=q * r_pq - s10,prob = col_active_prob) 
  
  beta <- rep(0,q)
  for(i in c(med_loc_true,med_loc_false)){
    beta[i] <-  rnorm(1,0.8,sqrt(0.1))*sample(c(1,-1),1) 
  }
  
  exp_loc_true <- sample(1:p,ceiling(p*r_pq))
  gamma  <- rep(0,p);
  for(i in exp_loc_true){
    gamma[i] <- rnorm(1,0.8,sqrt(0.1))*sample(c(1,-1),1) 
  }
  
  phi <- rnorm(s,eta,1)*sample(c(1,-1),size=s,replace = T)
  
  phi_z <- rnorm(p_z,eta,1)*sample(c(1,-1),size=p_z,replace = T)
  
  return(list(
    Psi = Psi, Sigma_E=Sigma_E,Theta=Theta,Phi=Phi,
    gamma = gamma,beta=beta,phi = phi,
    Psi_z = Psi_z, Phi_z = Phi_z, phi_z = phi_z
  ))
  
}

s <- 2; p_z <- 3
r_p <- r_pq <- 0.1; 
rho <- 0.6; p <- 100; q <- 100
eta <- 1
r_h_vector <- seq(0,1,0.1)


iter_num <- 1
comb_param_ls <- list()
param_comb_df <- NULL
for(r_h in r_h_vector){
  param_ls <- generate_param(p,q,eta,rho,s,r_p,r_pq,r_h,p_z=p_z)
  comb_param_ls[[iter_num]] <- param_ls
  param_comb_df <- rbind(param_comb_df,
                         c(iter_num,rho,eta,p,q,s,r_p,r_pq,r_h))
  iter_num <- iter_num + 1
  if(iter_num %% 10 ==0) cat('iter num:',iter_num,'\n')
  
}


colnames(param_comb_df) <- c('iter_num','rho','eta','p','q','s','r_p','r_pq','r_h')
param_comb_df <- as.data.frame(param_comb_df)
comb_param_ls$param_comb_df <- param_comb_df

saveRDS(comb_param_ls,file='../simuData/comb_param_ls_simu2.rds')

## generate datasets----------------
generate_dataset <- function(n,param_ls){
  s <- nrow(param_ls[['Psi']]); p <- ncol(param_ls[['Psi']]); q <- ncol(param_ls[['Theta']])
  p_z <- nrow(param_ls[['Psi_z']])
  
  H <- matrix(rnorm(n*s,0,1),nrow=n,byrow=T)
  Z <- matrix(rnorm(n*p_z,0,1),nrow=n,byrow=T)
  E <- MASS::mvrnorm(n=n,mu=rep(0,p),Sigma = param_ls[['Sigma_E']])
  X <- Z %*% param_ls[['Psi_z']] + H %*% param_ls[['Psi']] + E
  
  E_tilde <- matrix(rnorm(n*q,0,1),nrow=n,byrow = T)
  M <- X %*% param_ls[['Theta']] + Z %*% param_ls[['Phi_z']] + H %*% param_ls[['Phi']] + E_tilde
  
  e <- matrix(rnorm(n,0,1),nrow=n)
  Y <- X %*% param_ls[['gamma']] + M %*% param_ls[['beta']] + Z %*% param_ls[['phi_z']]
  H %*% param_ls[['phi']]+e
  
  return(list(
    X = X, M = M, Y = Y,Z = Z, H = H,
    Theta=param_ls[['Theta']],beta=param_ls[['beta']],
    gamma =param_ls[['gamma']]
  ))
}

# 产生数据
n_vector <- 300
comb_param_ls <- readRDS('../simuData/comb_param_ls_simu2.rds')

rep_data_comb_list <- list()
data_comb_df <- NULL
for(rep_num in 1:100){
  data_comb_list <- list()
  data_comb_num <- 1
  for(iter_num in 1:(length(comb_param_ls)-1)){
    param_ls <- comb_param_ls[[iter_num]]
    param_vec <- as.numeric(comb_param_ls$param_comb_df[iter_num,])
    
    for(n in n_vector){
      if(rep_num == 1){
        data_comb_df <- rbind(data_comb_df,c(data_comb_num,n,param_vec))
      }
      data_ls <- generate_dataset(n,param_ls)
      
      data_comb_list[[data_comb_num]] <- data_ls
      data_comb_num <- data_comb_num + 1
    }
    
  }
  if(rep_num %% 10 == 0) cat('rep num:',rep_num,'\n')
  rep_data_comb_list[[rep_num]] <- data_comb_list
}

data_comb_df <- as.data.frame(data_comb_df)
colnames(data_comb_df) <- c('data_num','n',colnames(comb_param_ls$param_comb_df))

rep_data_comb_list[['comb_type']] <- data_comb_df

saveRDS(rep_data_comb_list,file='../simuData/rep_comb_data_ls_simu2.rds')


# Simu3--------------------------------

## generate parameters-----------------------

generate_param <- function(p,q,eta,rho,s,r_p,r_pq,r_h,p_z=0,kappa=0.8){
  Sigma_E <- matrix(0,nrow=p,ncol=p)
  for(i in 1:(p-1)){
    for(j in (i+1):p){
      Sigma_E[i,j] <- (rho)**abs(i-j)
    }
  }
  Sigma_E <- Sigma_E + t(Sigma_E) + diag(p)
  
  Psi <- matrix(rnorm(n=s*p,mean=0,sd=1),
                nrow=s,byrow = T)
  
  if((p-ceiling(p*r_h))>0){
    for(i in 1:s){
      loc_active <- sample(1:p,size=ceiling(p*r_h),replace = F)
      loc_zero <- setdiff(1:p,loc_active)
      Psi[i,loc_zero] <- 0
    }
  }
  Psi_z <- matrix(rnorm(n=p_z*p,mean=0,sd=1),
                  nrow=p_z,byrow = T)
  
  Theta <- matrix(0,nrow=p,ncol=q)
  row_active <- sort(sample(1:p,size=ceiling(p*r_p),replace = F))
  col_active_ls <- vector('list',length(row_active))
  for(i in 1:length(row_active)){
    col_active_ls[[i]] <- sample(1:q,size=sample(5:50,size=1))
  }
  col_active <- sort(unique(unlist(col_active_ls)))
  for(i in 1:length(row_active)){
    row_loc <- row_active[i]
    for(col_loc in col_active_ls[[i]]){
      Theta[row_loc,col_loc] <- rnorm(1,kappa,sqrt(0.1))*sample(c(1,-1),1)
    }
  }
  
  Phi <- matrix(rnorm(s*q,mean=eta,sd=1)*sample(c(1,-1),size=s*q,replace = T),
                nrow=s,byrow=T)
  if((q-ceiling(q*r_h))>0){
    for(i in 1:s){
      loc_active <- sample(1:q,size=ceiling(q*r_h),replace = F)
      loc_zero <- setdiff(1:q,loc_active)
      Phi[i,loc_zero] <- 0
    }
  }
  
  Phi_z <- matrix(rnorm(p_z*q,mean=eta,sd=1)*sample(c(1,-1),size=p_z*q,replace = T),
                  nrow=p_z,byrow=T)
  
  
  col_silence <- setdiff(1:q,col_active)
  s10 <- min(length(col_silence),ceiling(q*r_pq*0.2))
  if(s10>0){
    med_loc_false <- sample(col_silence,size=s10)
  }else{
    med_loc_false <- NULL
  }
  col_active_prob <- table(unlist(col_active_ls)) / length(unlist(col_active_ls))
  med_loc_true <- sample(as.numeric(names(col_active_prob)),size=q * r_pq - s10,prob = col_active_prob) # 被target越多，越有可能成为true mediator
  
  beta <- rep(0,q)
  for(i in c(med_loc_true,med_loc_false)){
    beta[i] <-  rnorm(1,kappa,sqrt(0.1))*sample(c(1,-1),1) 
  }
  
  exp_loc_true <- sample(1:p,ceiling(p*r_pq))
  gamma  <- rep(0,p);
  for(i in exp_loc_true){
    gamma[i] <- rnorm(1,kappa,sqrt(0.1))*sample(c(1,-1),1) 
  }
  
  phi <- rnorm(s,eta,1)*sample(c(1,-1),size=s,replace = T)
  
  phi_z <- rnorm(p_z,eta,1)*sample(c(1,-1),size=p_z,replace = T)
  
  return(list(
    Psi = Psi, Sigma_E=Sigma_E,Theta=Theta,Phi=Phi,
    gamma = gamma,beta=beta,phi = phi,
    Psi_z = Psi_z, Phi_z = Phi_z, phi_z = phi_z
  ))
  
}

s <- 2; p_z <- 3
r_p <- r_pq <- 0.1; 
rho <- 0.6; q <- 100
eta <- 1

p_vec <- c(100,300)
rh_vec <- c(0,0.05,0.1,0.2,0.4,0.8,1)
kappa_vec <- seq(0.1,1.5,0.05)

iter_num <- 1
comb_param_ls <- list()
param_comb_df <- NULL
for(p in p_vec){
  for(r_h in rh_vec){
    for(kappa in kappa_vec){
      param_ls <- generate_param(p,q,eta,rho,s,r_p,r_pq,r_h=r_h,p_z=p_z,kappa=kappa)
      param_ls$p <- p
      param_ls$r_h <- r_h
      param_ls$kappa <- kappa
      
      comb_param_ls[[iter_num]] <- param_ls
      param_comb_df <- rbind(param_comb_df,
                             c(iter_num,rho,eta,p,q,s,r_p,r_pq,r_h,kappa))
      iter_num <- iter_num + 1
      if(iter_num %% 10 ==0) cat('iter num:',iter_num,'\n')
    }
    
  }
  
}


colnames(param_comb_df) <- c('iter_num','rho','eta','p','q','s','r_p','r_pq','r_h','kappa')
param_comb_df <- as.data.frame(param_comb_df)
comb_param_ls$param_comb_df <- param_comb_df
saveRDS(comb_param_ls,file='../simuData/comb_param_ls_simu3.rds')

## generate datasets----------------
generate_dataset <- function(n,param_ls){
  s <- nrow(param_ls[['Psi']]); p <- ncol(param_ls[['Psi']]); q <- ncol(param_ls[['Theta']])
  p_z <- nrow(param_ls[['Psi_z']])
  
  H <- matrix(rnorm(n*s,0,1),nrow=n,byrow=T)
  Z <- matrix(rnorm(n*p_z,0,1),nrow=n,byrow=T)
  E <- MASS::mvrnorm(n=n,mu=rep(0,p),Sigma = param_ls[['Sigma_E']])
  X <- Z %*% param_ls[['Psi_z']] + H %*% param_ls[['Psi']] + E
  
  E_tilde <- matrix(rnorm(n*q,0,1),nrow=n,byrow = T)
  M <- X %*% param_ls[['Theta']] + Z %*% param_ls[['Phi_z']] + H %*% param_ls[['Phi']] + E_tilde
  
  e <- matrix(rnorm(n,0,1),nrow=n)
  Y <- X %*% param_ls[['gamma']] + M %*% param_ls[['beta']] + Z %*% param_ls[['phi_z']]
  H %*% param_ls[['phi']]+e
  
  return(list(
    X = X, M = M, Y = Y,Z = Z, H = H,
    Theta=param_ls[['Theta']],beta=param_ls[['beta']],
    gamma =param_ls[['gamma']]
  ))
}

n <- 300
comb_param_ls <- readRDS('../simuData/comb_param_ls_simu3.rds')

for(rep_num in 1:100){
  for(iter_num in 1:(length(comb_param_ls)-1)){
    param_ls <- comb_param_ls[[iter_num]]
    kappa <- param_ls$kappa; r_h <- param_ls$r_h;p <- param_ls$p
    
    data_comb_list <- generate_dataset(n,param_ls)
    file_path <- paste0("../simuData/simu3/data_ls_p_",p,"_rh_",r_h,"_kappa_",kappa,'_rep_',rep_num,'.rda')
    saveRDS(data_comb_list,file=file_path)
  }
  if(rep_num %% 10 == 0) cat('rep num:',rep_num,'\n')
}


















