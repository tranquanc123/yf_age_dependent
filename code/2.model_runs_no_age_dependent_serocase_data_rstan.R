library(dplyr)
library(rstan)
library(extraDistr)
library(matrixStats)

#setwd("/home/quan/Documents/NotreDame/Alex/Yellow fever/yf_age_exposure/code/")
#running parallel in crc, which.scenario is the index:

#all_model_type = c('Constant', "One_outbreak", "Two_outbreaks")
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
#prior_gamma = readRDS("../data/gamma_lnorm_priors_constant_nonorm.rds")
#Vaccine efficacy: alpha and beta of beta distribution
a_VE = 22.812; b_VE = 1.306

#different subset of data:
case.noti.df = rbind(case.noti.LR.df, case.noti.PAHO.df)
#case.noti.df = case.noti.df[case.noti.df$case != 0,]
case.noti.df[is.na(case.noti.df$cov),]$pop = 1e-10
case.noti.df[is.na(case.noti.df$cov),]$cov = 0

###Study index:
#Sero 
sero_study_id = rle(as.character(sero_data$STUDY_ID))
N_sero_study = length(sero_study_id$lengths)
nrow_sero_study = mapply(rep, 1:N_sero_study, each = sero_study_id$lengths) %>% unlist()
sero_study_id = sero_study_id$values

#Case
case_study_id = rle(as.character(case.noti.df$study_id))
N_case_study = length(case_study_id$lengths)
nrow_case_study = mapply(rep, 1:N_case_study, each = case_study_id$lengths) %>% unlist() %>% as.vector()
case_study_id = case_study_id$values

N_study = N_sero_study + N_case_study

#All studies run
all_scen = data.frame(study_id = c(sero_study_id, case_study_id), data_type = c(rep("sero", N_sero_study), rep("case", N_case_study)))

for(which.scenario in 1:N_study){
  study = all_scen$study_id[which.scenario]
  data_type = all_scen$data_type[which.scenario]
  print(study)
  
  if(data_type == "sero"){
    sero_data_study = filter(sero_data, STUDY_ID == study)
    #Speed up calculation the likelihood for sero data:
    ind.pop = which(substr(names(sero_data_study),0,4)=='POP_')
    ind.cov = which(substr(names(sero_data_study),0,4)=='COV_')
    if(study %in% c("ETH-10-2014-TSEGAYE", "NGA-8-2008-BABA")){
      sel_0_cov_study = which(sero_data_study$STUDY_ID %in% c("ETH-10-2014-TSEGAYE", "NGA-8-2008-BABA")) #no coverage for Baba and Tsegaye
      sero_data_study[sel_0_cov_study,ind.cov] = 0 
    }
    n_FOI_sero = nrow(sero_data_study)
    age_range_l = 
      sapply(1:n_FOI_sero, function(x){
        age_range = (sero_data_study$AGE_LOWER[x]:sero_data_study$AGE_UPPER[x])
        age_range[age_range >= 100] <- 99
        return(age_range)
      })
    age_range_v = unlist(age_range_l)
    
    pop_l = 
      sapply(1:n_FOI_sero, function(x){
        age_range = age_range_l[[x]] + 1
        pop = sero_data_study[x,ind.pop[age_range]]
        return(pop)
      })
    pop_v = unlist(pop_l)
    cov_l = 
      sapply(1:n_FOI_sero, function(x){
        age_range = age_range_l[[x]] + 1
        cov = sero_data_study[x,ind.cov[age_range]]
        return(cov)
      })
    cov_v = unlist(cov_l)
    cov_v[cov_v < 0]<-1e-16
    cov_v[cov_v > 1]<-1 - 1e-16
    n_row_sero = lengths(age_range_l)
    n_row_sero_stop = cumsum(n_row_sero)
    n_row_sero_start = n_row_sero_stop - n_row_sero + 1
    Nrow_sero = length(n_row_sero)
    
    #stan code for sero
    stan_code = "
data {
  int<lower=0> N_row_sero_age; //number of rows in every studies for every age
  vector[N_row_sero_age] age_range_v; //age range of every row
  vector[N_row_sero_age] cov_v; //Coverage for every age row
  vector[N_row_sero_age] pop_v; // population for every age row
  int N_row_sero_age_range[N_row_sero_age];
  int<lower=0> Nrow_sero; //number of row in sero dataset 
  int n_row_sero_start[Nrow_sero]; 
  int n_row_sero_stop[Nrow_sero];
  int n_row_sero[Nrow_sero];
  int sero_pos[Nrow_sero];
  int sero_size[Nrow_sero];
  
  //prior
  real a_VE;
  real b_VE;
  real gamma_priors_mean;
  real gamma_priors_sd;
}

parameters {
  real gamma;
  real<lower=0, upper = 1> VE;
}

model {
  vector[N_row_sero_age] gamma_v;
  //matrix[100, N_sero_study] FOI_sero_v;
  //vector[N_row_sero_age] FOI_sero_cumsum;
  vector[N_row_sero_age] sus_pop;
  vector[N_row_sero_age] seropos_prob;
  real p[Nrow_sero];
  real p_ele;
  
  // prior:
  gamma ~ normal(gamma_priors_mean, gamma_priors_sd);
  VE ~ beta(a_VE, b_VE);

  gamma_v = rep_vector(pow(10.0, gamma), N_row_sero_age);
  
  seropos_prob = 1-exp(- gamma_v .* age_range_v);
  sus_pop = (cov_v * VE + (1 - cov_v * VE) .* seropos_prob) .* pop_v;
  
  for(iii in 1:Nrow_sero){
    int n_row_sero_seq[n_row_sero[iii]];
    n_row_sero_seq = N_row_sero_age_range[n_row_sero_start[iii]:n_row_sero_stop[iii]];
    p_ele = sum(sus_pop[n_row_sero_seq])/sum(pop_v[n_row_sero_seq]);
    if (p_ele > 1 - 1e-16)
      p_ele = 1 - 1e-16;
    if (p_ele < 1e-16)
      p_ele = 1e-16; 
    p[iii] = p_ele;
  }
  
  // LL:
  target += binomial_lpmf(sero_pos| sero_size, p);
}
"

#prior for gamma:
sero_list_study = sero_data_study$ISO %>% unique
prior_gamma_input = FOI_prior %>% filter(ISO == sero_list_study)

data_stan = list(
  N_row_sero_age = length(age_range_v), #number of rows in every studies for every age
  age_range_v = age_range_v, #age range of every row
  pop_v = as.numeric(pop_v), # population for every age row
  cov_v = as.numeric(cov_v), #coverage for every age row
  N_row_sero_age_range = 1:length(age_range_v),
  Nrow_sero = Nrow_sero, #number of row in sero dataset 
  n_row_sero_start = n_row_sero_start, 
  n_row_sero_stop = n_row_sero_stop,
  n_row_sero = n_row_sero,
  sero_pos = sero_data_study$POSITIVE,
  sero_size = sero_data_study$SAMPLE_SIZE,
  
  #prior
  a_VE = a_VE, b_VE = b_VE,
  gamma_priors_mean = prior_gamma_input$X1,
  gamma_priors_sd = prior_gamma_input$X2
)

#initial pars:
init_f <- function(n = 1){
  list(gamma = rnorm(n,  prior_gamma_input$X1, prior_gamma_input$X2),
       VE = rbeta(n, a_VE, b_VE)
  )
}
  }

if(data_type == "case"){
  case.noti.df_country = filter(case.noti.df, study_id == study)
  #Speed up calculation for case LR data:
  year_case_max = max(case.noti.df_country$year)
  FOI_case_id_row = year_case_max - case.noti.df_country$year + case.noti.df_country$age + 1
  FOI_case_id_row_max = max(FOI_case_id_row)
  study_agegroup = interaction(case.noti.df_country$age_group_ID, case.noti.df_country$study_id) %>% as.vector()
  case_by_study_agegroup = aggregate(case.noti.df_country$case, list(case.noti.df_country$age_group_ID), unique)$x
  length_case_study = nrow(case.noti.df_country)
  study_agegroup_org = rle(as.character(study_agegroup))$lengths
  study_agegroup_stop = cumsum(study_agegroup_org)
  study_agegroup_start = study_agegroup_stop - study_agegroup_org + 1
  
  #Constant model:####
  stan_code = "
data {
  int length_case_study;
  vector[length_case_study] FOI_case_id_row;
  vector[length_case_study] pop; //Population
  vector[length_case_study] cov; //coverage
  int<lower=0> Nrow_case; 
  int<lower=0> study_agegroup_org[Nrow_case];
  int<lower=0> study_agegroup_start[Nrow_case];
  int<lower=0> study_agegroup_stop[Nrow_case];
  int case_by_study_agegroup[Nrow_case];
  
  //prior
  real a_VE;
  real b_VE;
  vector[2] reporting_prop_prior;
  real gamma_priors_mean;
  real gamma_priors_sd;
}

parameters {
  real gamma;
  real rho_case;
  real<lower=0, upper = 1> VE;
}

model {
  vector[length_case_study] gamma_v;
  real trans_rho_case;
  vector[length_case_study] l_rho_case;
  vector[length_case_study] FOI_case_int_v; 
  vector[length_case_study] exp_cases;
  vector[Nrow_case] exp_case_by_age_group;
  vector[Nrow_case] exp_cases_mod;
  
  // prior:
  gamma ~ normal(gamma_priors_mean, gamma_priors_sd);
  rho_case ~ normal(reporting_prop_prior[1], reporting_prop_prior[2]);
  VE ~ beta(a_VE, b_VE);
  
  gamma_v = rep_vector(pow(10.0, gamma), length_case_study);
  FOI_case_int_v = (FOI_case_id_row - 1.0) .* gamma_v;
  trans_rho_case = inv_logit(rho_case);
  l_rho_case = rep_vector(trans_rho_case, length_case_study);
  exp_cases = exp(-FOI_case_int_v) .* (1.0 - exp(-gamma_v)) .* l_rho_case .* pop .* (1 - cov*VE);
  
  for(xviii in 1:Nrow_case){
    exp_case_by_age_group[xviii] = sum(exp_cases[study_agegroup_start[xviii]:study_agegroup_stop[xviii]]);
  }

  for(xx in 1:(Nrow_case)){
    if(exp_case_by_age_group[xx] < 1e-3)
      exp_cases_mod[xx] = 1e-3;
    else 
      exp_cases_mod[xx] = exp_case_by_age_group[xx];
  }
  
  // LL:
  target += poisson_lpmf(case_by_study_agegroup| exp_cases_mod);
}
"

#prior for gamma:
coun_list_study = case.noti.df_country$country %>% unique
if(coun_list_study %in% unique(case.noti.PAHO.df$country)) coun_list_study = "BRA"
prior_gamma_input = filter(FOI_prior, ISO == coun_list_study)

data_stan = list(
  FOI_case_id_row = FOI_case_id_row,
  length_case_study = length_case_study,
  #length_case_study_start = length_case_study_start,
  #length_case_study_stop = length_case_study_stop,
  Nrow_case = length(study_agegroup_org), 
  pop = case.noti.df_country$pop, 
  cov = case.noti.df_country$cov,
  study_agegroup_org = study_agegroup_org,
  study_agegroup_start = study_agegroup_start,
  study_agegroup_stop = study_agegroup_stop,
  case_by_study_agegroup = case_by_study_agegroup,
  
  #prior
  a_VE = a_VE, b_VE = b_VE,
  reporting_prop_prior = reporting_prop_prior, #c(1.05393, 800.78830), #reporting_prop_prior,
  gamma_priors_mean = prior_gamma_input$X1,
  gamma_priors_sd = prior_gamma_input$X2
)

#initial pars:
init_f <- function(n = 1){
  list(gamma = rnorm(n,  prior_gamma_input$X1, prior_gamma_input$X2),
       rho_case = rnorm(n, reporting_prop_prior[1], reporting_prop_prior[2]),
       VE = rbeta(n, a_VE, b_VE)
  )
}
}

stan_fit <- stan(
  model_code = stan_code,  # Stan program
  data = data_stan,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 2500,          # number of warmup iterations per chain
  iter = 5e3,            # total number of iterations per chain
  cores = 4,  # number of cores (could use one per chain)
  thin = 10,
  control = list(adapt_delta = 0.99),
  init = init_f,
  #refresh = 0             # no progress shown
)
#traceplot(stan_fit)
saveRDS(stan_fit, file = paste0("../output/MCMC_output/MCMC_output_1e4_noage_", data_type, "_", study, "_init_rstan.Rdata"))

}

###Extract the output####
library(dplyr)
library(rstan)
library(extraDistr)

#setwd("/home/quan/Documents/NotreDame/Alex/Yellow fever/yf_age_exposure/code/")
#running parallel in crc, which.scenario is the index:

#all_model_type = c('Constant', "One_outbreak", "Two_outbreaks")
model_type = "Constant" #all_model_type[which.scenario]
data_subset = "SA_AF"

inverse_logit <- function(x) exp(x)/(exp(x) + 1)
quantile95cri <- function(x) quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))

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
sero_study_id = rle(as.character(sero_data$STUDY_ID))
N_sero_study = length(sero_study_id$lengths)
nrow_sero_study = mapply(rep, 1:N_sero_study, each = sero_study_id$lengths) %>% unlist()
sero_study_id = sero_study_id$values

#Case
case_study_id = rle(as.character(case.noti.df$study_id))
N_case_study = length(case_study_id$lengths)
nrow_case_study = mapply(rep, 1:N_case_study, each = case_study_id$lengths) %>% unlist() %>% as.vector()
case_study_id = case_study_id$values

N_study = N_sero_study + N_case_study

ind.pop = which(substr(names(sero_data),0,4)=='POP_')
ind.cov = which(substr(names(sero_data),0,4)=='COV_')

#All studies run
all_scen = data.frame(study_id = c(sero_study_id, case_study_id), data_type = c(rep("sero", N_sero_study), rep("case", N_case_study)))

par_sum_l = c()
par_dist_l = c()
sero_data_fit_l = c()
case_data_fit_l = c()
WAIC_v_l = c()
WAIC_l = c()

#Read output:
for(which.scenario in 1:N_study){
  study = all_scen$study_id[which.scenario]
  data_type = all_scen$data_type[which.scenario]
  print(study)
  
  output = readRDS(paste0("../output/MCMC_output/MCMC_output_1e4_noage_", data_type, "_", study, "_init_rstan.Rdata"))
  
  #summary parameters:
  par_sum = summary(output)$summary[-which(rownames(summary(output)$summary) == "lp__"),c("2.5%", "50%", "97.5%", "Rhat")] %>% data.frame()
  par_sum$Study_id = study
  par_sum$data_type = data_type
  par_sum_l = c(par_sum_l, list(par_sum))
  
  #par extract
  par_sel = extract(output, pars = "lp__", include = FALSE) %>% Reduce(cbind, .)
  par_sel[,1] = 10^par_sel[,1]
  par_dist_l = c(par_dist_l, list(par_sel))
  
  n_samples = nrow(par_sel)
  #data fit and WAIC:
  if(data_type == "sero"){
    
    sel_sero_data = filter(sero_data, STUDY_ID == study)
    if(study %in% c("ETH-10-2014-TSEGAYE", "NGA-8-2008-BABA")){
      sel_0_cov_study = which(sel_sero_data$STUDY_ID %in% c("ETH-10-2014-TSEGAYE", "NGA-8-2008-BABA")) #no coverage for Baba and Tsegaye
      sel_sero_data[sel_0_cov_study,ind.cov] = 0 
    }
    
    gamma = par_sel[,1];
    VE = par_sel[,2]
    p_m  = matrix(nrow = nrow(sel_sero_data), ncol = nrow(par_sel))
    for(study_row in 1:nrow(sel_sero_data)){
      age_range = sel_sero_data$AGE_LOWER[study_row]:sel_sero_data$AGE_UPPER[study_row] + 1
      cov_v = sel_sero_data[study_row, ind.cov][age_range] %>% sapply(function(x) x*VE) %>% t
      pop_v = sel_sero_data[study_row, ind.pop][age_range] %>% as.vector() %>% unlist
      FOI_sero_cumsum = sapply(gamma, function(x) (age_range)*x)
      seropos = 1-exp(-FOI_sero_cumsum)
      sus_v = pop_v * (cov_v + (1-cov_v) * seropos)
      p = colSums(sus_v)/sum(pop_v)
      p[p > 1 - 1e-16] = 1 - 1e-16
      p[p < 1e-16] = 1e-16
      p_m[study_row,] = p
    }
    p_m_sum = apply(p_m, 1, quantile95cri)
    sero_data_fit = cbind(data_subset, model_type, study, data_type, sel_sero_data[,1:8], data.frame(t(p_m_sum)))
    sero_data_fit_l = c(sero_data_fit_l, list(sero_data_fit))
    
    #WAIC:
    LL_samples = sapply(1:n_samples, function(x) dbinom(sel_sero_data$POSITIVE, sel_sero_data$SAMPLE_SIZE, p_m[,x], log = T))
    lppd_v <- sapply(1:nrow(sero_data_fit), function(i) logSumExp(LL_samples[i,]) - log(n_samples))
    pWAIC_v <- sapply(1:nrow(sero_data_fit) , function(i) var(LL_samples[i,]))
    WAIC_v = -2*(lppd_v - pWAIC_v)
    WAIC_v_l = c(WAIC_v_l, WAIC_v)
    WAIC_l = c(WAIC_l, sum(WAIC_v))
    
  } else if(data_type == "case"){
    
    sel_case_data = filter(case.noti.df, study_id == study)
    gamma = par_sel[,1];
    rho = inverse_logit(par_sel[,2]);
    VE = par_sel[,3]
    
    FOI_case_id_row = max(sel_case_data$year) - sel_case_data$year + sel_case_data$age + 1
    FOI_case_int_v = sapply(gamma, function(x) (FOI_case_id_row - 1) * x) %>% t
    pop_m = matrix(sel_case_data$pop, nrow = n_samples, ncol = nrow(sel_case_data))
    cov_m = matrix(sel_case_data$cov, nrow = n_samples, ncol = nrow(sel_case_data))
    exp_cases = exp(-FOI_case_int_v) * (1 - exp(-gamma)) * rho * pop_m * (1 - cov_m*VE)
    exp_case_by_age_group = aggregate(t(exp_cases), by = list(sel_case_data$age_group_ID), sum)[,-1]
    exp_case_by_age_group[exp_case_by_age_group < 1e-3] = 1e-3
    exp_case_by_age_group_sum = apply(exp_case_by_age_group, 1, quantile95cri)
    
    case_by_study_agegroup_i = aggregate(sel_case_data$case, by = list(sel_case_data$age_group_ID), unique)$x
    case_data_fit = data.frame(data_subset, model_type, study, data_type, 
                              age_l = aggregate(sel_case_data$age, by = list(sel_case_data$age_group_ID), min)$x, 
                              age_u = aggregate(sel_case_data$age, by = list(sel_case_data$age_group_ID), max)$x, 
                              t(exp_case_by_age_group_sum), case = case_by_study_agegroup_i)
    case_data_fit_l = c(case_data_fit_l, list(case_data_fit))
    
    #WAIC:
    LL_samples = sapply(1:n_samples, function(x) dpois(case_by_study_agegroup_i, exp_case_by_age_group[,x], log = T))
    lppd_v <- sapply(1:nrow(exp_case_by_age_group), function(i) logSumExp(LL_samples[i,]) - log(n_samples))
    pWAIC_v <- sapply(1:nrow(exp_case_by_age_group) , function(i) var(LL_samples[i,]))
    WAIC_v = -2*(lppd_v - pWAIC_v)
    WAIC_v_l = c(WAIC_v_l, WAIC_v)
    WAIC_l = c(WAIC_l, sum(WAIC_v))
  }
}

par_sum = par_sum_l %>% lapply(function(x) {
  x$par = rownames(x)
  return(x)
  }) %>% Reduce(rbind, .)
sero_data_fit = Reduce(rbind, sero_data_fit_l)
case_data_fit = Reduce(rbind, case_data_fit_l)
all_scen$WAIC = WAIC_l

saveRDS(list(par_sum = par_sum, par_dist = par_dist_l, 
             sero_data_fit = sero_data_fit, case_data_fit = case_data_fit, 
             WAIC = all_scen, WAIC_v = WAIC_v_l),
        "../output/par_sum_fit_noage.rds")
