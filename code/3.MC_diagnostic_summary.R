library(BayesianTools)
library(dplyr)
library(Hmisc)
library(matrixStats)
library(sn)
library(stringr)

setwd("/home/quan/Documents/NotreDame/Alex/Yellow fever/yf_age_exposure/code/")

#Load data:
sero_data = read.csv("../data/yf_sero_data_with_coverage_mod.csv")
Study_id_seq = rle(as.character(sero_data$STUDY_ID))
one_age_group_data = Study_id_seq$values[Study_id_seq$lengths == 1]
sero_data = filter(sero_data, !(STUDY_ID %in% one_age_group_data))
case.noti.LR.df = read.csv("../data/yf_case_data_by_age_group_with_coverage_LR.csv")
VE = 0.975
Study_time = c((sero_data %>% group_by(STUDY_ID) %>% summarise(YEAR = mean(YEAR)))$YEAR,
               (case.noti.LR.df %>% group_by(study_id) %>% summarise(year = mean(year)))$year)
#All model runs:
all_model_type = c("Constant", "One outbreak", "Two outbreaks")

#Utility functions:
source("util_function.R")

#Generating summary:
par_sel_sum_l = c()
FOI_time_sum_data_l = c()
#pos_age_sum_data_l = c()
rho_sum_data_l = c()
Age_depend_sum_data_l = c()
Rhat_l = c()
WAIC_all = c()
sero_data_fit_l = c()
case_LR_data_fit_l = c()
p_m_l = c()
for(which.scenario in 1:length(all_model_type)){
  print(which.scenario)
  #model runs:
  model_type = all_model_type[which.scenario]
  #Parameter id:
  par_id(model_type)
  #FOI function: outbreak model:
  FOI_function_model_type(model_type)
  
  
  #Par credible interval; 
  par_sel = readRDS(paste0("../output/MCMC_output/MCMC_output_sero_case_noWHO_morethanoneagegroup_", model_type,"_extracted.Rdata"))
  gc()
  
  if(model_type == "Constant"){
    start_prob = 1/20
    end_prob = 1/2
  } else if(model_type == "One outbreak"){
    start_prob = 1/2
    end_prob = 1
  } else {
    start_prob = 1/4
    end_prob = 1
  }

  par_sel_thin = matrix(nrow = 1002, ncol = ncol(par_sel))
  for(ii in 1:ncol(par_sel)){
    post_chains = matrix(par_sel[,ii], ncol = 3, nrow = nrow(par_sel)/3, byrow = T)
    sel_it = round(nrow(post_chains)*start_prob) %>% seq(end_prob*nrow(post_chains), by = round(nrow(post_chains)*(end_prob - start_prob)/334))
    par_sel_thin[, ii] = post_chains[sel_it,] %>% t %>% as.vector()
    #matplot(post_chains[sel_it,])
    #Rhat = c(Rhat, rstan::Rhat(post_chains[sel_it,]))
  }
  #names(Rhat) = par_names
  par_sel = par_sel_thin
  n_samples = nrow(par_sel)
  #par_sel = getSample(par_sel, start = 1e7/3/2, thin = 5e3) %>% tail(1002)
  # par_sel = readRDS(paste0("output/MCMC_output/crc/Twalk_output_sero_case_noWHO_", model_type,".Rdata"))
  # par_sel = getSample(par_sel, start = 2e7/2/2, thin = 1e4) %>% tail(1000)
  
  FOI_par = par_sel[,FOI_par_id]; 
  if(model_type == "One outbreak"){
    FOI_par[,alpha_id] <- 10^(FOI_par[,alpha_id])
  } else if(model_type == "Two outbreaks") {
    FOI_par[,c(alpha1_id, alpha2_id)] <- 10^(FOI_par[,c(alpha1_id, alpha2_id)])
  } else if(model_type == "Constant"){
    FOI_par[,FOI_par_id] <- 10^(FOI_par[,FOI_par_id])
  }
  
  FOI_time_array = array(dim = c(100, N_study, n_samples))
  for(i in 1:n_samples){
    FOI_time_array[,,i] = FOI_func(unlist(FOI_par[i,])) %>% matrix(nrow = 100, ncol = N_study)
  }
  
  Age_depend_par = par_sel[,Age_depend_id];
  age_dep_term_array = array(dim = c(100, N_study, n_samples))
  for(i in 1:n_samples){
    age_dep_term_array[,,i] = age_dep_func(unlist(Age_depend_par[i,]))
  }
  if(model_type == "Constant"){
    saveRDS(age_dep_term_array[,1,], "../output/Age_exp_curve_Constant.rds")
  }
  
  #Include the effect of bg infection and age dependent:
  FOI_array = ((FOI_time_array)*age_dep_term_array)
  
  par_sel[,rho_case_LR_id] = inverse_logit(par_sel[,rho_case_LR_id])
  
  #Sum par:
  FOI_time_sum = apply(FOI_time_array, c(1, 2), quantile95cri)
  FOI_time_age_sum = apply(FOI_time_array*age_dep_term_array, c(1, 2), quantile95cri)
  FOI_time_sum_data = lapply(1:N_study, function(x) rbind(t(FOI_time_sum[,,x]), t(FOI_time_age_sum[,,x]))) %>% 
    Reduce(rbind, .) %>% 
    data.frame(Study_id = rep(names_study, each = 100*2), time = 0:99, par = rep(c("FOI_time", "FOI_time.age"), each = 100),model_type)
  Age_depend_sum = apply(age_dep_term_array, 1, quantile95cri)
  Age_depend_sum_data = data.frame(t(Age_depend_sum), age = 0:99, model_type)
  FOI_time_sum_data_l = c(FOI_time_sum_data_l, list(FOI_time_sum_data))
  Age_depend_sum_data_l = c(Age_depend_sum_data_l, list(Age_depend_sum_data))
  rho_sum_data_l = c(rho_sum_data_l, list(data.frame(Study_id = case_LR_study_id, datatype = "LR", apply(par_sel[,rho_case_LR_id], 2, quantile95cri) %>% t, model_type)))
  
  #Pred pos proportion:
  # pos_age_array = apply(FOI_array, c(2, 3), cumsum)
  # pos_age_array = 1 - exp(-pos_age_array)
  # pos_age_sum = apply(pos_age_array, c(1, 2), quantile95cri)
  # pos_age_sum_data = lapply(1:N_study, function(x) t(pos_age_sum[,,x])) %>% 
  #   Reduce(rbind, .) %>% 
  #   data.frame(Study_id = rep(names_study, each = 100), age = 0:99, model_type, bg, age_dependent)
  # pos_age_sum_data_l = c(pos_age_sum_data_l, list(pos_age_sum_data))
  
  #Datafit:
  all_datafit = lapply(1:n_samples, function(x) gen_func(unlist(par_sel[x,])))
  p_m = sapply(all_datafit, function(x) x$sero)
  p_m_l = c(p_m_l, list(p_m))
  case_LR_fit = sapply(all_datafit, function(x) x$case_LR)
  sero_data_fit = sapply(all_datafit, function(x) x$sero) %>% apply(1, quantile95cri) %>% t %>%
    data.frame(select(sero_data, ISO, YEAR, AGE_LOWER, AGE_UPPER, SAMPLE_SIZE, POSITIVE, STUDY_ID), model_type)
  sero_data_fit_l = c(sero_data_fit_l, list(sero_data_fit))
  case_LR_data_fit = sapply(all_datafit, function(x) x$case_LR) %>% apply(1, quantile95cri) %>% t %>%
    cbind(case.noti.LR.or.df, model_type, .) %>% data.frame()
  case_LR_data_fit_l = c(case_LR_data_fit_l, list(case_LR_data_fit))

  #WAIC
  #Calculate WAIC:
  LL_samples = sapply(1:n_samples, function(x) {
    c(dbinom(sero_data$POSITIVE, sero_data$SAMPLE_SIZE, p_m[,x], log = T),
      dpois(case_by_study_agegroup, lambda = case_LR_fit[,x], log = T))
  })
  lppd_v <- sapply(1:nrow(LL_samples) , function(i) logSumExp(LL_samples[i,]) - log(n_samples))
  pWAIC_v <- sapply(1:nrow(LL_samples) , function(i) var(LL_samples[i,]))
  WAIC = -2*(sum(lppd_v) - sum(pWAIC_v))
  WAIC_all = c(WAIC_all, WAIC)
  # #data specific WAIC:
  # WAIC_v = -2*(lppd_v - pWAIC_v)
  # WAIC_each_data_l = c(WAIC_each_data_l, list(data.frame(Study_id = names_study, WAIC = aggregate(WAIC_v, list(n_row_study), sum)$x, model_type = model_type)))
  # lppd_each_data_l = c(lppd_each_data_l, list(data.frame(Study_id = names_study, WAIC = aggregate(lppd_v, list(n_row_study), sum)$x, model_type = model_type)))
  # pWAIC_each_data_l = c(pWAIC_each_data_l, list(data.frame(Study_id = names_study, WAIC = aggregate(pWAIC_v, list(n_row_study), sum)$x, model_type = model_type)))
  
  #Sum all par:
  par_sel[,FOI_par_id] = FOI_par; 
  #reorder the timing parameter for the Two outbreaks model:
  if(which.scenario == 3){
    for(i in 1:length(time1_id)){
      print(i)
      par_sel[,c(time1_id[i], time2_id[i])] = cbind(par_sel[,time1_id][,i], par_sel[,time2_id][,i]) %>% apply(1, sort) %>% t
    } 
  }
  par_sel[,Age_depend_id[2]] = exp(par_sel[,Age_depend_id[2]])
  if(which.scenario == 2){
    par_sel[,c(time_id)] = matrix(Study_time, ncol = length(Study_time), nrow = nrow(par_sel), byrow = T) - par_sel[, time_id]
  }
  if(which.scenario == 3){
    par_sel[,c(time1_id, time2_id)] = matrix(Study_time, ncol = length(Study_time)*2, nrow = nrow(par_sel), byrow = T) - par_sel[,c(time1_id, time2_id)]
  }
  
  par_sel_sum_l = c(par_sel_sum_l, list(data.frame(par_names = par_names, apply(par_sel, 2, quantile95cri) %>% t, model_type)))

  #Rhat:
  Rhat = c()
  for(ii in 1:ncol(par_sel)){
    post_chains = matrix(par_sel[,ii], ncol = 3, nrow = n_samples/3, byrow = T)
    Rhat = c(Rhat, rstan::Rhat(post_chains))
  }
  names(Rhat) = par_names
  Rhat_l = c(Rhat_l, list(Rhat))
  
}

par_sel_sum = Reduce(rbind, par_sel_sum_l)
FOI_time_sum_data =  Reduce(rbind, FOI_time_sum_data_l)
Age_depend_sum_data =  Reduce(rbind, Age_depend_sum_data_l)
sero_data_fit =  Reduce(rbind, sero_data_fit_l)
case_LR_data_fit =  Reduce(rbind, case_LR_data_fit_l)

#assess convergence:
lapply(Rhat_l, mean)
lapply(Rhat_l, function(x) (which(x > 1.1)))

matplot(matrix(par_sel[,15], ncol = 3, nrow = n_samples/3, byrow = T), type = "l")

#WIAC:
WAIC_all

saveRDS(list(par_sel_sum = par_sel_sum, FOI_time_sum_data = FOI_time_sum_data, Age_depend_sum_data = Age_depend_sum_data,
             sero_data_fit = sero_data_fit, case_LR_data_fit = case_LR_data_fit), "../output/par_sum_fit.rds")
