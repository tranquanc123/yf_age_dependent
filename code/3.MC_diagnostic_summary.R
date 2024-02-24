library(dplyr)
library(matrixStats)
library(stringr)
library(sn)
library(rstan)
library(extraDistr)

setwd("/home/quan/Documents/NotreDame/Alex/Yellow fever/yf_age_exposure/code/")

#Load data:
model_type = "Constant" #all_model_type[which.scenario]
data_subset = "SA_AF"

#Load data:
sero_data = read.csv("../data/yf_sero_data_with_coverage_mod.csv")
sero_data = filter(sero_data, YEAR > 1980) #only get data from 1980
sero_study_id = rle(as.character(sero_data$STUDY_ID))
sero_data = filter(sero_data, STUDY_ID %in% sero_study_id$values[which(sero_study_id$lengths > 1)]) # study with age group > 1
case.noti.LR.df = read.csv("../data/yf_case_data_by_age_group_with_coverage_LR.csv")
case.noti.PAHO.df = read.csv("../data/yf_case_data_by_age_group_with_coverage_PAHO.csv") %>% filter(!(study_id %in% c("ECU_PAHO", "ARG_PAHO")))
readRDS("../data/FOI_rho_prior.rds") %>% list2env(globalenv())
prior_gamma = readRDS("../data/gamma_lnorm_priors_constant_nonorm.rds")
#Vaccine efficacy: alpha and beta of beta distribution
a_VE = 22.812; b_VE = 1.306

#different subset of data:
case.noti.df = rbind(case.noti.LR.df, case.noti.PAHO.df)
#case.noti.df = case.noti.df[case.noti.df$case != 0,]
case.noti.df[is.na(case.noti.df$cov),]$pop = 1e-10
case.noti.df[is.na(case.noti.df$cov),]$cov = 0

###Study index:
#Sero 
sero_study_rle = rle(as.character(sero_data$STUDY_ID))
N_sero_study = length(sero_study_rle$lengths)
nrow_sero_study = mapply(rep, 1:N_sero_study, each = sero_study_rle$lengths) %>% unlist()
sero_study_id = sero_study_rle$values

#Case
case_study_rle = rle(as.character(case.noti.df$study_id))
N_case_study = length(case_study_rle$lengths)
nrow_case_study = mapply(rep, 1:N_case_study, each = case_study_rle$lengths) %>% unlist()
case_study_id = case_study_rle$values

N_study = N_sero_study + N_case_study
names_study = c(sero_study_id, case_study_id)

#Speed up calculation the likelihood for sero data:
ind.pop = which(substr(names(sero_data),0,4)=='POP_')
ind.cov = which(substr(names(sero_data),0,4)=='COV_')
sel_0_cov_study = which(sero_data$STUDY_ID %in% c("ETH-10-2014-TSEGAYE", "NGA-8-2008-BABA")) #no coverage for Baba and Tsegaye
sero_data[sel_0_cov_study,ind.cov] = 0
n_FOI_sero = nrow(sero_data)
age_range_l = 
  sapply(1:n_FOI_sero, function(x){
    age_range = (sero_data$AGE_LOWER[x]:sero_data$AGE_UPPER[x])
    age_range[age_range >= 100] <- 99
    return(age_range)
  })
age_range_v = unlist(age_range_l)
nrow_sero_study_age = mapply(rep, nrow_sero_study, lengths(age_range_l)) %>% unlist
pop_l = 
  sapply(1:n_FOI_sero, function(x){
    age_range = age_range_l[[x]] + 1
    pop = sero_data[x,ind.pop[age_range]]
    return(pop)
  })
pop_v = unlist(pop_l)

cov_l = 
  sapply(1:n_FOI_sero, function(x){
    age_range = age_range_l[[x]] + 1
    cov = sero_data[x,ind.cov[age_range]]
    return(cov)
  })
cov_v = unlist(cov_l)
cov_v[cov_v < 0]<-1e-16
cov_v[cov_v > 1]<-1 - 1e-16
n_row_sero = lengths(age_range_l)
n_row_sero_l = unlist(mapply(rep, 1:length(n_row_sero), n_row_sero))
pop_by_row = aggregate(pop_v, list(n_row_sero_l), sum)$x


#Speed up calculation for case LR data:
year_case_max = aggregate(case.noti.df$year, list(nrow_case_study), max)$x[nrow_case_study]
FOI_case_id_row = year_case_max - case.noti.df$year + case.noti.df$age + 1
FOI_case_id_row_max = max(FOI_case_id_row)
study_agegroup = interaction(case.noti.df$age_group_ID, case.noti.df$study_id) %>% as.vector %>% factor(level = rle(.)$values)
case_by_study_agegroup = aggregate(case.noti.df$case, list(study_agegroup), unique)$x
length_case_study = rle(as.character(nrow_case_study))$lengths

n_row_study_case = rle(as.character(study_agegroup))$values %>% str_split("[.]") %>% sapply(function(x) x[2]) %>% rle
n_row_study = c(mapply(rep, 1:N_sero_study, each = sero_study_rle$lengths) %>% unlist, 
                mapply(rep, N_sero_study + 1:N_case_study, each = n_row_study_case$lengths) %>% unlist)

#retreive back the original case LR data:
case.noti.or.df = unique(case.noti.df[,c("age_group_ID", "country", "study_id")])
case.noti.or.df$age_min = aggregate(case.noti.df$age, by = list(case.noti.df$age_group_ID), min)$x
case.noti.or.df$age_max = aggregate(case.noti.df$age, by = list(case.noti.df$age_group_ID), max)$x
case.noti.or.df$case = case_by_study_agegroup

#helper function:
#Age dependent function
age_dep_func <- function(sn_par){
  age_dep = dsn(0:99, sn_par[1], sn_par[2], sn_par[3])  #1: postition (range 0:99) 2: scale (positive) 3:skewness (real number)
  return(age_dep)
}
inverse_logit <- function(x) exp(x)/(exp(x) + 1)
quantile95cri <- function(x) quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))
#function to generate cases and serological data:####
gen_func <- function(par){
  
  FOI_par = unlist(par[names_par == "gamma"]); 
  
  FOI_m = matrix(unlist(FOI_par), nrow = 100, ncol = N, byrow = T) 
  
  if(model_sel != "noage"){
    age_dep_term_array = array(dim = c(100, N_coun))
    for(study in 1:N_coun){
      age_dep1 = par[which(names_par == "age_dep1")[study]]
      age_dep2 = par[which(names_par == "age_dep2")[study]]
      age_dep3 = par[which(names_par == "age_dep3")[study]]
      age_dep_term_array[,study] = age_dep_func(c(age_dep1, age_dep2, age_dep3))
    }
  } else {
    age_dep_term_array = array(1, dim = c(100, N_coun))
  }
  
  
  #Include the effect of bg infection and age dependent:
  FOI_m = ((FOI_m)*age_dep_term_array[,age_dep_coun_id]) 
  
  VE = par[names_par == "VE"]
  
  rho_case = inverse_logit(unlist(par[names_par == "rho_case"]))
  
  if(data_subset != "SA"){
    FOI_sero_m = FOI_m[, 1:N_sero_study]
    FOI_sero_cumsum = sapply(1:length(age_range_v), function(x) sum(FOI_sero_m[1:(age_range_v[x] + 1), nrow_sero_study_age[x]]))
    seropos = 1-exp(-FOI_sero_cumsum)
    sus_v = (cov_v*VE + (1-cov_v*VE) * seropos) * pop_v
    p = aggregate(sus_v, list(n_row_sero_l), sum)$x/pop_by_row
    p[p > 1 - 1e-16] = 1 - 1e-16
    p[p < 1e-16] = 1e-16
  } else {
    p = NULL
  }
  
  ###Case LR data:
  #case data from LR:
  FOI_case_m = FOI_m[, N_sero_study + 1:N_case_study]
  FOI_case_m = rbind(FOI_case_m, matrix(FOI_case_m[100,], nrow = FOI_case_id_row_max - 100, ncol = N_case_study, byrow = T))
  FOI_case_v = sapply(1:length(nrow_case_study), function(i) FOI_case_m[FOI_case_id_row[i], nrow_case_study[i]])
  FOI_case_m_cumsum = apply(FOI_case_m, 2, cumsum)
  FOI_case_int_v = sapply(1:length(nrow_case_study), function(x) {
    FOI_case_id_int = FOI_case_id_row[x] - 1
    if(FOI_case_id_int == 0) return(1)
    return(FOI_case_m_cumsum[FOI_case_id_int,nrow_case_study[x]]) 
  })
  
  l_rho_case = mapply(rep, rho_case, each = length_case_study) %>% unlist #rep(rho_case, each = N_case_study*100)
  case_sus = case.noti.df$pop*(1 - case.noti.df$cov*VE)
  exp_cases = case_sus*exp(-FOI_case_int_v)*(1 - exp(-FOI_case_v))*l_rho_case
  
  exp_case_by_age_group = aggregate(exp_cases, by = list(study_agegroup), FUN = sum)$x
  exp_case_by_age_group[exp_case_by_age_group < 1e-3] = 1e-3
  
  return(list(sero = p, case = exp_case_by_age_group))
}

#prior for gamma:
coun_list_study = sapply(case_study_id, function(x) filter(case.noti.df, study_id == x)$country %>% unique)
coun_list_study_prior = coun_list_study
coun_list_study_prior[!(coun_list_study_prior %in% prior_gamma$ISO)] <- "BRA"
sero_list_study = sapply(sero_study_id, function(x) filter(sero_data, STUDY_ID == x)$ISO %>% unique)
match_prior_gamma = match(c(sero_list_study, coun_list_study_prior), prior_gamma$ISO)
prior_gamma_input = prior_gamma[match_prior_gamma,]
#age-dep by country:
all_coun = unique(c(sero_list_study, coun_list_study))
age_dep_coun_id = match(c(sero_list_study, coun_list_study), all_coun)

n_samples = 1e3

model_sel = "noage"
#Par credible interval; 
#stan_fit_ob = readRDS(paste0("../output/MCMC_output/MCMC_output_1e4_3pois_covmod_6PAHO_sero_case_age_depend_", model_type,"_datasubset_", data_subset, "_rstan.Rdata"))

if(model_sel == "globalcurve"){
  stan_fit_ob = readRDS(paste0("../output/MCMC_output/MCMC_output_1e4_globalcurve_init_rstan.Rdata"))
} else if(model_sel == "counspec") {
  stan_fit_ob = readRDS(paste0("../output/MCMC_output/MCMC_output_1e4_heir_init_rstan.Rdata"))
} else if(model_sel == "noage"){
  stan_fit_ob = readRDS(paste0("../output/MCMC_output/MCMC_output_1e4_noage_init_rstan.Rdata"))
}

sum_stan = summary(stan_fit_ob)$summary
summary(stan_fit_ob)$c_summary["lp__",, ]
#traceplot(stan_fit_ob, pars = c("age_dep1", "age_dep2", "age_dep3", "lp__"))
# traceplot(stan_fit_ob, pars = c("gamma"))
# traceplot(stan_fit_ob, pars = c("rho_case"))
# traceplot(stan_fit_ob, pars = c("VE"))
# pairs(stan_fit_ob, pars =c("VE" ,"age_dep1_mean", "age_dep2_mean", "age_dep3_mean"))
# stan_dens(stan_fit_ob, pars = c("age_dep1", "age_dep2", "age_dep3"))

if(model_sel == "globalcurve"){
  par_sel = extract(stan_fit_ob, pars = c("gamma", "age_dep1", "age_dep2", "age_dep3", 
                                          #"age_dep1_mean", "age_dep2_mean", "age_dep3_mean", "age_dep1_sd", "age_dep2_sd", "age_dep3_sd",
                                          "rho_case", "VE", "lp__"))
} else if(model_sel == "counspec") {
  par_sel = extract(stan_fit_ob, pars = c("gamma", "age_dep1", "age_dep2", "age_dep3", 
                                          "age_dep1_mean", "age_dep2_mean", "age_dep3_mean", "age_dep1_sd", "age_dep2_sd", "age_dep3_sd",
                                          "rho_case", "VE", "lp__"))
} else if(model_sel == "noage") {
  par_sel = extract(stan_fit_ob, pars = c("gamma", #"age_dep1", "age_dep2", "age_dep3", 
                                          #"age_dep1_mean", "age_dep2_mean", "age_dep3_mean", "age_dep1_sd", "age_dep2_sd", "age_dep3_sd",
                                          "rho_case", "VE", "lp__"))
}

times_par = sapply(par_sel, ncol); times_par[is.na(times_par)] = 1; names_par = mapply(rep, names(par_sel), each = times_par) %>% unlist
par_sel = par_sel %>% Reduce(cbind, .)
colnames(par_sel) = names_par
#par_sel = par_sel[par_sel[,"lp__"] > -1000,]
par_sel = tail(par_sel, n_samples)
n_samples = nrow(par_sel)

if(model_sel == "noage"){
  par_sel[,which(names_par == "gamma")] = 10^par_sel[,which(names_par == "gamma")]
}
FOI_par = par_sel[,which(names_par == "gamma")]

N = N_study;        #number studies
if(model_sel %in% c("globalcurve", "noage")){
  N_coun = 1#length(all_coun) #=> change for global curve model
} else if(model_sel == "counspec") {
  N_coun = length(all_coun) #=> change for global curve model
}


FOI_time_array = array(dim = c(100, N, n_samples))
for(i in 1:n_samples){
  FOI_time_array[,,i] = matrix(unlist(FOI_par[i,]), nrow = 100, ncol = N, byrow = T)
}

if(model_sel != "noage"){
  Age_depend_par = par_sel[,grep("age_dep", names_par)]
}

if(model_sel == "globalcurve"){
  age_dep_term_pop = array(dim = c(100, n_samples))
  age_dep1_pop = par_sel[,which(names_par == "age_dep1")]#lapply(1:n_samples, function(x) rnorm(1e3, par_sel[x,which(names_par == "age_dep1_mean")], par_sel[x,which(names_par == "age_dep1_sd")])) %>% unlist
  age_dep2_pop = par_sel[,which(names_par == "age_dep2")]#lapply(1:n_samples, function(x) par_sel[x,which(names_par == "age_dep2_mean")] + rhnorm(1e3, par_sel[x,which(names_par == "age_dep2_sd")])) %>% unlist
  age_dep3_pop = par_sel[,which(names_par == "age_dep3")]#lapply(1:n_samples, function(x) rnorm(1e3, par_sel[x,which(names_par == "age_dep3_mean")], par_sel[x,which(names_par == "age_dep3_sd")])) %>% unlist
  age_dep2_pop = age_dep2_pop[age_dep2_pop > 0]
  for(i in 1:length(age_dep2_pop)){
    age_dep_term_pop[,i] = age_dep_func(c(age_dep1_pop[i], age_dep2_pop[i], age_dep3_pop[i]))
  }
} else if(model_sel == "counspec") {
  age_dep_term_array = array(dim = c(100, N_coun, n_samples))
  for(study in 1:N_coun){
    age_dep1 = par_sel[,which(names_par == "age_dep1")[study]]
    age_dep2 = par_sel[,which(names_par == "age_dep2")[study]]
    age_dep3 = par_sel[,which(names_par == "age_dep3")[study]]
    for(i in 1:n_samples){
      age_dep_term_array[,study,i] = age_dep_func(c(age_dep1[i], age_dep2[i], age_dep3[i]))
    }
  }
  
  age_dep_term_pop = array(dim = c(100, 1e6))
  age_dep1_pop = lapply(1:n_samples, function(x) rnorm(1e3, par_sel[x,which(names_par == "age_dep1_mean")], par_sel[x,which(names_par == "age_dep1_sd")])) %>% unlist
  age_dep2_pop = lapply(1:n_samples, function(x) par_sel[x,which(names_par == "age_dep2_mean")] + rhnorm(1e3, par_sel[x,which(names_par == "age_dep2_sd")])) %>% unlist
  age_dep3_pop = lapply(1:n_samples, function(x) rnorm(1e3, par_sel[x,which(names_par == "age_dep3_mean")], par_sel[x,which(names_par == "age_dep3_sd")])) %>% unlist
  age_dep2_pop = age_dep2_pop[age_dep2_pop > 0]
  for(i in 1:length(age_dep2_pop)){
    age_dep_term_pop[,i] = age_dep_func(c(age_dep1_pop[i], age_dep2_pop[i], age_dep3_pop[i]))
  }
} else if(model_sel == "noage") {
  age_dep_term_array = array(1, dim = c(100, N_coun, n_samples))
}
#matplot(age_dep_term_pop[,1:1e4], type = "l")
if(model_type == "Constant" & model_sel != "noage"){
  saveRDS(age_dep_term_pop, paste0("../output/Age_exp_", model_sel, "_", data_subset, "_Constant.rds"))
}

#Mode of age dependence:
if(model_sel %in% c("globalcurve", "noage")){
  sel_n_samples = n_samples
} else if(model_sel == "counspec") {
  sel_n_samples = 1e4
}

Age_dept_mode = sapply(1:sel_n_samples, function(x){
  optim(50, fn = function(y) -dsn(y, age_dep1_pop[x], (age_dep2_pop[x]), age_dep3_pop[x]),
        lower = 0, upper = 100, method = "Brent")$par
})

# Age_dept_coun_mode = sapply(1:(16*1e3), function(x){
#   age_dep1 = par_sel[,which(names_par == "age_dep1")][x]
#   age_dep2 = par_sel[,which(names_par == "age_dep2")][x]
#   age_dep3 = par_sel[,which(names_par == "age_dep3")][x]
#   optim(50, fn = function(y) -dsn(y, age_dep1, age_dep2, age_dep3),
#         lower = 0, upper = 100, method = "Brent")$par
# })

# Age_dept_coun_mode_sum = data.frame(country = all_coun, Age_dept_coun_mode %>% matrix(nrow = n_samples) %>% apply(2, quantile95cri) %>% t)
# Age_dept_coun_mode_sum[order(Age_dept_coun_mode_sum$X50.),]
# Age_dept_mode_l = c(Age_dept_mode_l, list(Age_dept_mode))
# 
#ratio between high and lowest:
if(model_sel != "noage"){
  Age_dept_low = sapply(1:sel_n_samples, function(x){
    optim(50, fn = function(y) dsn(y, age_dep1_pop[x], (age_dep2_pop[x]), age_dep3_pop[x]),
          lower = 0, upper = 100, method = "Brent")$par
  })
  Age_dept_ratio_hl = dsn(Age_dept_mode, age_dep1_pop[1:sel_n_samples], age_dep2_pop[1:sel_n_samples], age_dep3_pop[1:sel_n_samples])/dsn(Age_dept_low, age_dep1_pop[1:sel_n_samples], age_dep2_pop[1:sel_n_samples], age_dep3_pop[1:sel_n_samples])
  #Age_high_to_low_ratio_l = c(Age_high_to_low_ratio_l, list(Age_dept_ratio_hl))
  #ratio between highest and youngest:
  Age_dept_ratio_hy = (dsn(Age_dept_mode, age_dep1_pop[1:sel_n_samples], age_dep2_pop[1:sel_n_samples], age_dep3_pop[1:sel_n_samples])/dsn(0, age_dep1_pop[1:sel_n_samples], age_dep2_pop[1:sel_n_samples], age_dep3_pop[1:sel_n_samples])) %>% quantile95cri()
  Age_dept_ratio_ho = (dsn(Age_dept_mode, age_dep1_pop[1:sel_n_samples], age_dep2_pop[1:sel_n_samples], age_dep3_pop[1:sel_n_samples])/dsn(100, age_dep1_pop[1:sel_n_samples], age_dep2_pop[1:sel_n_samples], age_dep3_pop[1:sel_n_samples])) %>% quantile95cri()
}

#by country
# Age_dept_coun_ratio_hl = data.frame(coun = all_coun, (apply(age_dep_term_array, c(2, 3), max)/apply(age_dep_term_array, c(2, 3), min)) %>% apply(1, quantile95cri) %>% t %>% round(1))
# Age_dept_coun_ratio_hl[order(Age_dept_coun_ratio_hl$X50.),]
#Include the effect of bg infection and age dependent:
if(model_sel == "globalcurve"){
  age_dep_term_array = array(dim = dim(FOI_time_array))
  age_dep_coun_id = 1:N_study
  for(i in age_dep_coun_id){
    age_dep_term_array[,i,] = age_dep_term_pop
  }
} else if(model_sel == "noage"){
  age_dep_coun_id = rep(1, N_study)
  #FOI_time_array = 10^FOI_time_array
}

FOI_array = ((FOI_time_array)*age_dep_term_array[,age_dep_coun_id,])

rho_case = inverse_logit(par_sel[,which(colnames(par_sel) == "rho_case")])

#Sum par:
FOI_time_sum = apply(FOI_time_array, c(1, 2), quantile95cri)
FOI_time_age_sum = apply(FOI_array, c(1, 2), quantile95cri)
FOI_time_sum_data = lapply(1:N_study, function(x) rbind(t(FOI_time_sum[,,x]), t(FOI_time_age_sum[,,x]))) %>% 
  Reduce(rbind, .) %>% 
  data.frame(Study_id = rep(names_study, each = 100*2), time = 0:99, par = rep(c("FOI_time", "FOI_time.age"), each = 100),model_type, data_subset)
if(model_sel != "noage"){
  age_dep_term_array_scaled = age_dep_term_array/array(matrix(apply(age_dep_term_array, 2, mean), ncol = dim(age_dep_term_array)[2], nrow = dim(age_dep_term_array)[1], byrow = T), dim = dim(age_dep_term_array))
  
  if(model_sel == "counspec"){
    Age_depend_sum = apply(age_dep_term_array_scaled, c(1, 2), quantile95cri) %>% apply(1, rbind)
    Age_depend_sum_data = data.frame(Age_depend_sum, age = 0:99,
                                     #Study_id = rep(names_study, each = 100),
                                     country = rep(all_coun, each = 100),
                                     model_type, data_subset)
  }
  Age_depend_global_sum = apply(age_dep_term_pop[,!apply(age_dep_term_pop, 2, function(x) any(is.na(x)))], 1, quantile95cri) %>% apply(1, rbind)
  Age_depend_global_sum_data = data.frame(Age_depend_global_sum, age = 0:99, 
                                          #Study_id = rep(names_study, each = 100),
                                          #country = rep(all_coun, each = 100), 
                                          model_type, data_subset)
}


rho_sum_data = data.frame(Study_id = case_study_id, datatype = "LR", apply(rho_case, 2, quantile95cri) %>% t, model_type, data_subset)

# ggplot(Age_depend_sum_data, aes(x = age, y = X50.))+
#   geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "gray")+
#   geom_line()+
#   facet_wrap(.~country, scale = "free_y")
# 
# ggplot(Age_depend_sum_pop_data, aes(x = age, y = X50.))+
#   geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = "gray")+
#   geom_line()
# Age_dep_sum_pop = apply(age_dep_term_pop[,which(!is.na(age_dep_term_pop[1,]))], 1, quantile95cri)
# Age_depend_sum_pop_data = data.frame(t(Age_dep_sum_pop), age = 0:99)

#Datafit:
if(model_sel == "globalcurve"){
  age_dep_coun_id = 1
}
all_datafit = lapply(1:n_samples, function(x) gen_func(unlist(par_sel[x,])))
p_m = sapply(all_datafit, function(x) x$sero)
case_fit = sapply(all_datafit, function(x) x$case)
sero_data_fit = sapply(all_datafit, function(x) x$sero) %>% apply(1, quantile95cri) %>% t %>%
    data.frame(select(sero_data, ISO, YEAR, AGE_LOWER, AGE_UPPER, SAMPLE_SIZE, POSITIVE, STUDY_ID), model_type, data_subset)
case_data_fit = case_fit %>% apply(1, quantile95cri) %>% t %>%
  cbind(case.noti.or.df, model_type, data_subset, .) %>% data.frame()

# ggplot(sero_data_fit, aes(x = (AGE_LOWER + AGE_UPPER)/2))+
#   geom_pointrange(aes(y = X50., ymin = X2.5., ymax = X97.5.), color = "red")+
#   geom_point(aes(y = POSITIVE/SAMPLE_SIZE))+
#   facet_wrap(.~study, scale = "free_y")
# 
# ggplot(case_data_fit, aes(x = (age_l + age_u)/2))+
#   geom_pointrange(aes(y = X50., ymin = X2.5., ymax = X97.5.), color = "red")+
#   geom_point(aes(y = case))+
#   facet_wrap(.~study, scale = "free_y")
# 
#WAIC:
LL_samples = rbind(sapply(1:n_samples, function(x) dbinom(sero_data_fit$POSITIVE, sero_data_fit$SAMPLE_SIZE, p_m[,x], log = T)),
             sapply(1:n_samples, function(x) dpois(case.noti.or.df$case, case_fit[,x], log = T)))
lppd_v <- sapply(1:nrow(LL_samples), function(i) logSumExp(LL_samples[i,]) - log(n_samples))
pWAIC_v <- sapply(1:nrow(LL_samples) , function(i) var(LL_samples[i,]))
WAIC_v = -2*(lppd_v - pWAIC_v)
WAIC_data = aggregate(WAIC_v, by = list(c(sero_data_fit$STUDY_ID, case.noti.or.df$study_id)), sum)
WAIC_data = WAIC_data[match(unique(c(sero_data_fit$STUDY_ID, case.noti.or.df$study_id)), WAIC_data$Group.1),]
WAIC = sum(WAIC_v)

#Sum all par:
sum_par = t(apply(par_sel, 2, quantile95cri)) %>% data.frame
sum_par = sum_par[-which(rownames(sum_par) == "lp__"),]
sum_par$par = head(names_par, -1)

#Rhat:
Rhat = sum_stan[-nrow(sum_stan),"Rhat"]
which(Rhat > 1.1)

if(model_sel == "counspec"){
  saveRDS(list(par_sel_sum = sum_par, WAIC_data = WAIC_data, WAIC = WAIC, WAIC_v= WAIC_v,
               Age_depend_global_sum_data = Age_depend_global_sum_data, 
               Age_depend_sum_data = Age_depend_sum_data, 
               Age_dept_mode = quantile95cri(Age_dept_mode), Age_high_to_low_ratio = quantile95cri(Age_dept_ratio_hl), Age_dept_ratio_hy = Age_dept_ratio_hy, Age_dept_ratio_ho = Age_dept_ratio_ho, 
               sero_data_fit = sero_data_fit, case_data_fit = case_data_fit), "../output/par_sum_fit_heir.rds")
} else if(model_sel == "globalcurve"){
  saveRDS(list(par_sel_sum = sum_par, WAIC_data = WAIC_data, WAIC = WAIC, WAIC_v= WAIC_v,
               Age_depend_global_sum_data = Age_depend_global_sum_data, 
               #Age_depend_sum_data = Age_depend_sum_data, 
               Age_dept_mode = quantile95cri(Age_dept_mode), Age_high_to_low_ratio = quantile95cri(Age_dept_ratio_hl), Age_dept_ratio_hy = Age_dept_ratio_hy, Age_dept_ratio_ho = Age_dept_ratio_ho, 
               sero_data_fit = sero_data_fit, case_data_fit = case_data_fit), "../output/par_sum_fit_heir_globalcurve.rds")
} else if(model_sel == "noage"){
  saveRDS(list(par_sel_sum = sum_par, WAIC_data = WAIC_data, WAIC = WAIC, WAIC_v= WAIC_v,
               #Age_depend_global_sum_data = Age_depend_global_sum_data, 
               #Age_depend_sum_data = Age_depend_sum_data, 
               #Age_dept_mode = quantile95cri(Age_dept_mode), Age_high_to_low_ratio = quantile95cri(Age_dept_ratio_hl), Age_dept_ratio_hy = Age_dept_ratio_hy, Age_dept_ratio_ho = Age_dept_ratio_ho, 
               sero_data_fit = sero_data_fit, case_data_fit = case_data_fit), "../output/par_sum_fit_heir_noage.rds")
}



