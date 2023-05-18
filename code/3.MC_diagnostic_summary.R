library(dplyr)
library(matrixStats)
library(stringr)
library(sn)
library(rstan)

setwd("/home/quan/Documents/NotreDame/Alex/Yellow fever/yf_age_exposure/code/")

#Load data:
all_data_subset = c("SA", "AF", "SA_AF")
model_type = 'Constant' #all_model_type[which.scenario]

#Load data:
sero_data = read.csv("../data/yf_sero_data_with_coverage_mod.csv")
sero_data = filter(sero_data, YEAR > 1980) #only get data from 1980
sero_study_id = rle(as.character(sero_data$STUDY_ID))
sero_data = filter(sero_data, STUDY_ID %in% sero_study_id$values[which(sero_study_id$lengths > 1)]) # study with age group > 1
case.noti.LR.df = read.csv("../data/yf_case_data_by_age_group_with_coverage_LR.csv")
case.noti.PAHO.df = read.csv("../data/yf_case_data_by_age_group_with_coverage_PAHO.csv") %>% filter(!(study_id %in% c("ECU_PAHO", "ARG_PAHO")))
readRDS("../data/FOI_rho_prior.rds") %>% list2env(globalenv())
prior_gamma = readRDS("../data/gamma_lnorm_priors_constant.rds")
#Vaccine efficacy: alpha and beta of beta distribution
a_VE = 22.812; b_VE = 1.306
Study_time = c((sero_data %>% group_by(STUDY_ID) %>% summarise(YEAR = mean(YEAR)))$YEAR,
               (case.noti.LR.df %>% group_by(study_id) %>% summarise(year = mean(year)))$year)
# #All model runs:
# all_model_type = c("Constant", "One outbreak", "Two outbreaks")


#Generating summary:
par_sel_sum_l = c()
FOI_time_sum_data_l = c()
#pos_age_sum_data_l = c()
rho_sum_data_l = c()
Age_depend_sum_data_l = c()
Rhat_l = c()
WAIC_all = c()
WAIC_study_l = c()
sero_data_fit_l = c()
case_data_fit_l = c()
p_m_l = c()
Age_dept_mode_l = c()

n_samples = 1e3
for(which.scenario in 1:length(all_data_subset)){
  print(which.scenario)
  # #model runs:
  # model_type = all_model_type[which.scenario]
  data_subset = all_data_subset[which.scenario]
  
  sero_data = read.csv("../data/yf_sero_data_with_coverage_mod.csv")
  sero_data = filter(sero_data, YEAR > 1980) #only get data from 1980
  sero_study_id = rle(as.character(sero_data$STUDY_ID))
  sero_data = filter(sero_data, STUDY_ID %in% sero_study_id$values[which(sero_study_id$lengths > 1)]) # study with age group > 1
  
  source("util_function.R")
  par_id(model_type)
  FOI_function_model_type(model_type)
  if(which.scenario != 2){
    case.noti.df[is.na(case.noti.df$cov),]$pop = 1e-10
    case.noti.df[is.na(case.noti.df$cov),]$cov = 0
  }
  
  #Par credible interval; 
  stan_fit_ob = readRDS(paste0("../output/MCMC_output/MCMC_output_1e4_3pois_stepsize2_6PAHO_sero_case_age_depend_", model_type,"_datasubset_", data_subset, "_rstan.Rdata"))
  #stan_fit_ob = readRDS(paste0("../output/MCMC_output/MCMC_output_1e4_3pois_6PAHO_sero_case_age_depend_Constant_datasubset_AF_rstan.Rdata"))
  
  # summary(stan_fit_ob)$summary
  # traceplot(stan_fit_ob, pars = c("age_dep1", "age_dep2", "age_dep3", "lp__"))
  # traceplot(stan_fit_ob, pars = c("gamma"))
  # traceplot(stan_fit_ob, pars = c("rho_case"))
  # traceplot(stan_fit_ob, pars = c("VE"))
  # pairs(stan_fit_ob, pars =c("VE" ,"age_dep1", "age_dep2", "age_dep3"))
  # stan_dens(stan_fit_ob, pars = c("age_dep1", "age_dep2", "age_dep3"))
  
  par_sel = extract(stan_fit_ob, pars = c("gamma", "age_dep1", "age_dep2", "age_dep3", "rho_case", "VE", "lp__")) %>% lapply(tail, n_samples) %>% Reduce(cbind, .)
  
  # par_sel = extract(stan_fit_ob, pars = c("gamma", "age_dep1", "age_dep2", "age_dep3", "rho_case", "VE")) %>% lapply(tail, n = n_samples) %>%
  #   Reduce(cbind, .)
  
  FOI_par = par_sel[,FOI_par_id]
  
  FOI_time_array = array(dim = c(100, N_study, n_samples))
  for(i in 1:n_samples){
    FOI_time_array[,,i] = FOI_func(unlist(FOI_par[i,])) %>% matrix(nrow = 100, ncol = N_study)
  }
  
  Age_depend_par = par_sel[,Age_depend_id]
  age_dep_term_array = array(dim = c(100, N_study, n_samples))
  for(i in 1:n_samples){
    age_dep_term_array[,,i] = age_dep_func(unlist(Age_depend_par[i,]))
  }
  #matplot(age_dep_term_array[,1,], type = "l")
  if(model_type == "Constant"){
    saveRDS(age_dep_term_array[,1,], paste0("../output/Age_exp_curve_", data_subset, "_Constant.rds"))
  }

  #Mode of age dependence:
  Age_dept_mode = sapply(1:n_samples, function(x){
    optim(50, fn = function(y) -dsn(y, Age_depend_par[x,1], (Age_depend_par[x,2]), Age_depend_par[x,3]),
          lower = 0, upper = 100, method = "Brent")$par
  })
  Age_dept_mode_l = c(Age_dept_mode_l, list(Age_dept_mode))
  #Include the effect of bg infection and age dependent:
  FOI_array = ((FOI_time_array)*age_dep_term_array)
  
  rho_case = inverse_logit(par_sel[,rho_case_id])
  
  #Sum par:
  FOI_time_sum = apply(FOI_time_array, c(1, 2), quantile95cri)
  FOI_time_age_sum = apply(FOI_time_array*age_dep_term_array, c(1, 2), quantile95cri)
  FOI_time_sum_data = lapply(1:N_study, function(x) rbind(t(FOI_time_sum[,,x]), t(FOI_time_age_sum[,,x]))) %>% 
    Reduce(rbind, .) %>% 
    data.frame(Study_id = rep(names_study, each = 100*2), time = 0:99, par = rep(c("FOI_time", "FOI_time.age"), each = 100),model_type, data_subset)
  Age_depend_sum = apply(age_dep_term_array, 1, quantile95cri)
  Age_depend_sum_data = data.frame(t(Age_depend_sum), age = 0:99, model_type, data_subset)
  FOI_time_sum_data_l = c(FOI_time_sum_data_l, list(FOI_time_sum_data))
  Age_depend_sum_data_l = c(Age_depend_sum_data_l, list(Age_depend_sum_data))
  rho_sum_data_l = c(rho_sum_data_l, list(data.frame(Study_id = case_study_id, datatype = "LR", apply(par_sel[,rho_case_id], 2, quantile95cri) %>% t, model_type, data_subset)))
  
  #Datafit:
  all_datafit = lapply(1:n_samples, function(x) gen_func(unlist(par_sel[x,])))
  p_m = sapply(all_datafit, function(x) x$sero)
  p_m_l = c(p_m_l, list(p_m))
  case_fit = sapply(all_datafit, function(x) x$case)
  if(data_subset != "SA"){
    sero_data_fit = sapply(all_datafit, function(x) x$sero) %>% apply(1, quantile95cri) %>% t %>%
      data.frame(select(sero_data, ISO, YEAR, AGE_LOWER, AGE_UPPER, SAMPLE_SIZE, POSITIVE, STUDY_ID), model_type, data_subset)
    sero_data_fit_l = c(sero_data_fit_l, list(sero_data_fit))
  }
  case_data_fit = case_fit %>% apply(1, quantile95cri) %>% t %>%
    cbind(case.noti.or.df, model_type, data_subset, .) %>% data.frame()
  case_data_fit_l = c(case_data_fit_l, list(case_data_fit))
  
  #WAIC:
  if(data_subset != "SA"){
    LL_samples = sapply(1:n_samples, function(x) {
      c(dbinom(sero_data$POSITIVE, sero_data$SAMPLE_SIZE, p_m[,x], log = T),
      dpois(case.noti.or.df$case, case_fit[,x], log = T))
    })
    } else {
      LL_samples = sapply(1:n_samples, function(x) dpois(case.noti.or.df$case, case_fit[,x], log = T))
  }
  lppd_l <- sapply(1:nrow(LL_samples), function(i) logSumExp(LL_samples[i,]) - log(n_samples))# %>% sum
  pWAIC_l <- sapply(1:nrow(LL_samples) , function(i) var(LL_samples[i,]))# %>% sum
  WAIC_study = aggregate(-2*(lppd_l - pWAIC_l), list(n_row_study), sum)$x
  WAIC_study_l = c(WAIC_study_l, WAIC_study)
  WAIC = sum(WAIC_study)
  WAIC_all = c(WAIC_all, WAIC)

  #Sum all par:
  sum_stan = summary(stan_fit_ob)$summary
  sum_par = sum_stan[-which(rownames(sum_stan) == "lp__"), c( "2.5%", "50%", "97.5%")] %>% data.frame
  sum_par$par = row.names(sum_par)
  par_sel_sum_l = c(par_sel_sum_l, list(data.frame(sum_par, model_type, data_subset)))

  #Rhat:
  Rhat = sum_stan[-nrow(sum_stan),"Rhat"]
  Rhat_l = c(Rhat_l, list(Rhat))
}

par_sel_sum = Reduce(rbind, par_sel_sum_l)
FOI_time_sum_data =  Reduce(rbind, FOI_time_sum_data_l)
Age_depend_sum_data =  Reduce(rbind, Age_depend_sum_data_l)
sero_data_fit =  Reduce(rbind, sero_data_fit_l)
case_data_fit =  Reduce(rbind, case_data_fit_l)

#assess convergence:
lapply(Rhat_l, mean)
lapply(Rhat_l, function(x) (which(x > 1.1)))

saveRDS(list(par_sel_sum = par_sel_sum, FOI_time_sum_data = FOI_time_sum_data, Age_depend_sum_data = Age_depend_sum_data, Age_dept_mode = Age_dept_mode_l, WAIC_study_l = WAIC_study_l,
             sero_data_fit = sero_data_fit, case_data_fit = case_data_fit), "../output/par_sum_fit.rds")

