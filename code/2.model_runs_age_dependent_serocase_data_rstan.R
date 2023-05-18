library(dplyr)
library(rstan)

#setwd("/home/quan/Documents/NotreDame/Alex/Yellow fever/yf_age_exposure/code/")
#running parallel in crc, which.scenario is the index:
args = commandArgs(trailingOnly=TRUE)
which.scenario = as.numeric(args[1])
#which.scenario = 1

#all_model_type = c('Constant', "One_outbreak", "Two_outbreaks")
all_data_subset = c("SA", "AF", "SA_AF")
model_type = "Constant" #all_model_type[which.scenario]
data_subset = all_data_subset[which.scenario]

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

#different subset of data:
if(which.scenario == 1){
  case.noti.df = rbind(filter(case.noti.LR.df, country == "BRA"), case.noti.PAHO.df)
  sero_data = sero_data[NULL,]
} else if(which.scenario == 2){
  case.noti.df = filter(case.noti.LR.df, country != "BRA")
} else {
  case.noti.df = rbind(case.noti.LR.df, case.noti.PAHO.df)
}

case.noti.df = case.noti.df[case.noti.df$case != 0,]
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

if(which.scenario == 1){
  N_study = N_case_study
} else {
  N_study = N_sero_study + N_case_study
}

if(which.scenario != 1){
#Speed up calculation the likelihood for sero data:
ind.pop = which(substr(names(sero_data),0,4)=='POP_')
ind.cov = which(substr(names(sero_data),0,4)=='COV_')
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
n_row_sero_stop = cumsum(n_row_sero)
n_row_sero_start = n_row_sero_stop - n_row_sero + 1
Nrow_sero = length(n_row_sero)
}
#n_row_sero_l = unlist(mapply(rep, 1:length(n_row_sero), n_row_sero))

#Speed up calculation for case LR data:
year_case_max = aggregate(case.noti.df$year, list(nrow_case_study), max)$x[nrow_case_study]
FOI_case_id_row = year_case_max - case.noti.df$year + case.noti.df$age + 1
FOI_case_id_row_max = max(FOI_case_id_row)
study_agegroup = interaction(case.noti.df$age_group_ID, case.noti.df$study_id) %>% as.vector()
case_by_study_agegroup = aggregate(case.noti.df$case, list(study_agegroup), unique)$x
length_case_study = rle(as.character(nrow_case_study))$lengths
length_case_study_stop = cumsum(length_case_study)
length_case_study_start = length_case_study_stop - length_case_study + 1
study_agegroup_org = rle(as.character(study_agegroup))$lengths
study_agegroup_stop = cumsum(study_agegroup_org)
study_agegroup_start = study_agegroup_stop - study_agegroup_org + 1

###Stan code ####

#Constant model:####
if(model_type == "Constant" & data_subset == "SA"){
  
  stan_code = "
data {
  int<lower=0> N;          // number studies
  vector[100] age_seq; // age sequence from 0 to 99

  int N_case_study_id[N]; //N_sero_study + 1:N
  int<lower=0> N_nrow_case_study; //length(nrow_case_study)
  int<lower=0> nrow_case_study[N_nrow_case_study];
  int<lower=0> FOI_case_id_row[N_nrow_case_study];
  int length_case_study[N];
  int length_case_study_start[N];
  int length_case_study_stop[N];
  int nrow_FOI_append;
  vector[N_nrow_case_study] pop; //Population
  vector[N_nrow_case_study] cov; //coverage
  int<lower=0> Nrow_case; 
  int<lower=0> study_agegroup_org[Nrow_case];
  int<lower=0> study_agegroup_start[Nrow_case];
  int<lower=0> study_agegroup_stop[Nrow_case];
  int<lower=0> case_by_study_agegroup[Nrow_case];
  
  //prior
  real a_VE;
  real b_VE;
  vector[2] reporting_prop_prior;
  vector[N] gamma_priors_mean;
  vector[N] gamma_priors_sd;
}

parameters {
  vector<lower=0>[N] gamma;
  real age_dep1; //age dependent: position of the peak age
  real<lower=0> age_dep2; //age dependent: spread of the peak
  real age_dep3; //age dependent: skewness
  vector<lower=-20, upper = 20>[N] rho_case;
  real<lower=0, upper = 1> VE;
}

model {
  //vector[100] age_dep_output_org;
  vector[100] age_dep_output;
  matrix[100, N] age_dep_m;
  matrix[N, 100] gamma_m;
  matrix[100, N] FOI_m;

  vector[N] trans_rho_case;
  vector[N_nrow_case_study] l_rho_case;
  matrix[100, N] FOI_case_m;
  matrix[nrow_FOI_append, N] FOI_case_m_appending;
  matrix[100 + nrow_FOI_append, N] FOI_case_m_appended;
  vector[N_nrow_case_study] FOI_case_v;
  vector[N_nrow_case_study] FOI_case_int_v; 
  vector[N_nrow_case_study] exp_cases_prob;
  vector[N_nrow_case_study] exp_cases;
  vector[Nrow_case] exp_case_by_age_group;
  vector[Nrow_case] exp_cases_mod;
  
  // prior:
  gamma ~ lognormal(gamma_priors_mean, gamma_priors_sd);
  age_dep1 ~ normal(50, 125);
  age_dep2 ~ normal(0, 150);
  age_dep3 ~ student_t(2, 0, 0.5);
  rho_case ~ normal(reporting_prop_prior[1], reporting_prop_prior[2]);
  VE ~ beta(a_VE, b_VE);
  
  //Age dependent component of FOI:
  for(age in 1:100){
    age_dep_output[age] = exp(skew_normal_lpdf(age_seq[age]|age_dep1, age_dep2, age_dep3));
  }
  //age_dep_output = exp(age_dep_output_org - rep_vector(log_sum_exp(age_dep_output_org), 100));

  //Combine to get the FOI matrix:
  age_dep_m = rep_matrix(age_dep_output, N); //row: age; col: study
  gamma_m = rep_matrix(gamma, 100); //row: study; col: age
  FOI_m = age_dep_m .* (gamma_m)'; //row: age; col: study
  
  //LL from case data:
  trans_rho_case = inv_logit(rho_case);
  FOI_case_m_appending = rep_matrix(FOI_m[100,], nrow_FOI_append);
  FOI_case_m_appended = append_row(FOI_m, FOI_case_m_appending);
  for(xi in 1:N_nrow_case_study){
    int FOI_case_id_int = FOI_case_id_row[xi] - 1;
    FOI_case_v[xi] = FOI_case_m_appended[FOI_case_id_row[xi], nrow_case_study[xi]];
    if(FOI_case_id_int == 0)
      FOI_case_int_v[xi] = 0.0;
    else
      FOI_case_int_v[xi] = sum(head(FOI_case_m_appended[, nrow_case_study[xi]], FOI_case_id_int));
  }
  
  l_rho_case = trans_rho_case[nrow_case_study];
  exp_cases_prob = exp(-FOI_case_int_v) .* (1 - exp(-FOI_case_v)) .* l_rho_case;
  exp_cases = exp_cases_prob .* pop .* (1 - cov*VE);

  for(xviii in 1:Nrow_case){
    exp_case_by_age_group[xviii] = sum(exp_cases[study_agegroup_start[xviii]:study_agegroup_stop[xviii]]);
  }

  for(xx in 1:(Nrow_case)){
    if(exp_case_by_age_group[xx] < 1e-3)
      exp_cases_mod[xx] = 1e-3;
    else 
      exp_cases_mod[xx] = exp_case_by_age_group[xx];
  }
  //print(poisson_lpmf(case_by_study_agegroup| exp_cases_mod));
  // LL:
  case_by_study_agegroup ~ poisson(exp_cases_mod);
  //target += poisson_lpmf(case_by_study_agegroup| exp_cases_mod);
}
"

#prior for gamma:
coun_list_study = sapply(case_study_id, function(x) filter(case.noti.df, study_id == x)$country %>% unique)
coun_list_study[!(coun_list_study %in% prior_gamma$ISO)] <- "BRA"
match_prior_gamma = match(c( coun_list_study), prior_gamma$ISO)
prior_gamma_input = prior_gamma[match_prior_gamma,]

if(which.scenario != 2){
  case.noti.df[is.na(case.noti.df$cov),]$pop = 1e-10
  case.noti.df[is.na(case.noti.df$cov),]$cov = 0
}

data_stan = list(
  N = N_study,          #number studies
  age_seq = 0:99, #age sequence from 0 to 99
  #n_samples_marLL = n_samples_marLL, #number of sample to calcualte marginal LL
  #time_FOI_empty = matrix(0, nrow = N_study, ncol = 100),
  
  N_case_study = N_case_study, #number of case LR studies
  N_case_study_id = 1:N_case_study, #N_sero_study + N_case_study + 1:N_case_study
  N_nrow_case_study = length(nrow_case_study), #length(nrow_case_study)
  nrow_case_study = nrow_case_study,
  FOI_case_id_row = FOI_case_id_row,
  nrow_FOI_append = max(FOI_case_id_row) - 100,
  length_case_study = length_case_study,
  length_case_study_start = length_case_study_start,
  length_case_study_stop = length_case_study_stop,
  Nrow_case = length(study_agegroup_org), 
  pop = case.noti.df$pop, 
  cov = case.noti.df$cov,
  study_agegroup_org = study_agegroup_org,
  study_agegroup_start = study_agegroup_start,
  study_agegroup_stop = study_agegroup_stop,
  case_by_study_agegroup = case_by_study_agegroup,
  
  #prior
  a_VE = a_VE, b_VE = b_VE,
  reporting_prop_prior = reporting_prop_prior,
  gamma_priors_mean = prior_gamma_input$X1,
  gamma_priors_sd = prior_gamma_input$X2
)

} 

if(model_type == "Constant" & data_subset != "SA"){
  
  stan_code = "
data {
  int<lower=0> N;          // number studies
  vector[100] age_seq; // age sequence from 0 to 99
  
  int N_sero_study;          // number sero studies
  int N_sero_study_id[N_sero_study]; // id of sero studies
  int<lower=0> N_row_sero_age; //number of rows in every studies for every age
  int age_range_v[N_row_sero_age]; //age range of every row
  int nrow_sero_study_age[N_row_sero_age]; //number of age in every row
  vector[N_row_sero_age] cov_v; //Coverage for every age row
  vector[N_row_sero_age] pop_v; // population for every age row
  int N_row_sero_age_range[N_row_sero_age];
  int<lower=0> Nrow_sero; //number of row in sero dataset 
  int n_row_sero_start[Nrow_sero]; 
  int n_row_sero_stop[Nrow_sero];
  int n_row_sero[Nrow_sero];
  int sero_pos[Nrow_sero];
  int sero_size[Nrow_sero];
  
  int N_case_study; //number of case studies
  int N_case_study_id[N_case_study]; //N_sero_study + 1:N_case_study
  int<lower=0> N_nrow_case_study; //length(nrow_case_study)
  int<lower=0> nrow_case_study[N_nrow_case_study];
  int<lower=0> FOI_case_id_row[N_nrow_case_study];
  int length_case_study[N_case_study];
  int length_case_study_start[N_case_study];
  int length_case_study_stop[N_case_study];
  int nrow_FOI_append;
  vector[N_nrow_case_study] pop; //Population
  vector[N_nrow_case_study] cov; //coverage
  int<lower=0> Nrow_case; 
  int<lower=0> study_agegroup_org[Nrow_case];
  int<lower=0> study_agegroup_start[Nrow_case];
  int<lower=0> study_agegroup_stop[Nrow_case];
  int<lower=0> case_by_study_agegroup[Nrow_case];
  
  //prior
  real a_VE;
  real b_VE;
  vector[2] reporting_prop_prior;
  vector[N] gamma_priors_mean;
  vector[N] gamma_priors_sd;
}

parameters {
  //real<lower=-20, upper = 20> FOI[N];
  vector<lower=0>[N] gamma;
  real age_dep1; //age dependent: position of the peak age
  real<lower=0> age_dep2; //age dependent: spread of the peak
  real age_dep3; //age dependent: skewness
  vector<lower=-20, upper = 20>[N_case_study] rho_case;
  real<lower=0, upper = 1> VE;
}

model {
  //vector[100] age_dep_output_org;
  vector[100] age_dep_output;
  matrix[100, N] age_dep_m;
  matrix[N, 100] gamma_m;
  matrix[100, N] FOI_m;
  matrix[100, N_sero_study] FOI_sero_m;
  vector[N_row_sero_age] FOI_sero_cumsum;
  vector[N_row_sero_age] sus_pop;
  vector[N_row_sero_age] seropos_prob;
  real p[Nrow_sero];
  real p_ele;

  vector[N_case_study] trans_rho_case;
  vector[N_nrow_case_study] l_rho_case;
  matrix[100, N_case_study] FOI_case_m;
  matrix[nrow_FOI_append, N_case_study] FOI_case_m_appending;
  matrix[100 + nrow_FOI_append, N_case_study] FOI_case_m_appended;
  vector[N_nrow_case_study] FOI_case_v;
  vector[N_nrow_case_study] FOI_case_int_v; 
  vector[N_nrow_case_study] exp_cases;
  vector[Nrow_case] exp_case_by_age_group;
  vector[Nrow_case] exp_cases_mod;
  
  // prior:
  gamma ~ lognormal(gamma_priors_mean, gamma_priors_sd);
  age_dep1 ~ normal(50, 125);
  age_dep2 ~ normal(0, 150);
  age_dep3 ~ student_t(2, 0, 0.5);
  rho_case ~ normal(reporting_prop_prior[1], reporting_prop_prior[2]);
  VE ~ beta(a_VE, b_VE);
  
  //Age dependent component of FOI:
  for(age in 1:100){
    age_dep_output[age] = exp(skew_normal_lpdf(age_seq[age]|age_dep1, age_dep2, age_dep3));
  }
  //age_dep_output = exp(age_dep_output_org - rep_vector(log_sum_exp(age_dep_output_org), 100));
  
  //Combine to get the FOI matrix:
  age_dep_m = rep_matrix(age_dep_output, N);
  gamma_m = rep_matrix(gamma, 100);
  FOI_m = age_dep_m .* (gamma_m)';
  
  //LL of sero data:
  FOI_sero_m = FOI_m[, N_sero_study_id];
  for(i in 1:N_row_sero_age){
    FOI_sero_cumsum[i] = sum(head(FOI_sero_m[, nrow_sero_study_age[i]], age_range_v[i] + 1));
  }
  
  seropos_prob = 1-exp(-FOI_sero_cumsum);
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
  
  //LL from case data:
  trans_rho_case = inv_logit(rho_case);
  FOI_case_m = FOI_m[, N_case_study_id];
  FOI_case_m_appending = rep_matrix(FOI_case_m[100,], nrow_FOI_append);
  FOI_case_m_appended = append_row(FOI_case_m, FOI_case_m_appending);
  for(xi in 1:N_nrow_case_study){
    int FOI_case_id_int = FOI_case_id_row[xi] - 1;
    FOI_case_v[xi] = FOI_case_m_appended[FOI_case_id_row[xi], nrow_case_study[xi]];
    if(FOI_case_id_int == 0)
      FOI_case_int_v[xi] = 0.0;
    else
      FOI_case_int_v[xi] = sum(head(FOI_case_m_appended[, nrow_case_study[xi]], FOI_case_id_int));
  }
  
  l_rho_case = trans_rho_case[nrow_case_study];
  exp_cases = exp(-FOI_case_int_v) .* (1 - exp(-FOI_case_v)) .* l_rho_case .* pop .* (1 - cov*VE);

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
  target += binomial_lpmf(sero_pos| sero_size, p); //from sero data
  target += poisson_lpmf(case_by_study_agegroup| exp_cases_mod);
}
"

#prior for gamma:
coun_list_study = sapply(case_study_id, function(x) filter(case.noti.df, study_id == x)$country %>% unique)
coun_list_study[!(coun_list_study %in% prior_gamma$ISO)] <- "BRA"
sero_list_study = sapply(sero_study_id, function(x) filter(sero_data, STUDY_ID == x)$ISO %>% unique)
match_prior_gamma = match(c(sero_list_study, coun_list_study), prior_gamma$ISO)
prior_gamma_input = prior_gamma[match_prior_gamma,]

#
if(which.scenario != 2){
  case.noti.df[is.na(case.noti.df$cov),]$pop = 1e-10
  case.noti.df[is.na(case.noti.df$cov),]$cov = 0
}

data_stan = list(
  N = N_study,          #number studies
  age_seq = 0:99, #age sequence from 0 to 99
  #n_samples_marLL = n_samples_marLL, #number of sample to calcualte marginal LL
  #time_FOI_empty = matrix(0, nrow = N_study, ncol = 100),
  
  N_sero_study = N_sero_study,          #number sero studies
  N_sero_study_id = 1:N_sero_study, #id of sero studies
  N_row_sero_age = length(age_range_v), #number of rows in every studies for every age
  age_range_v = age_range_v, #age range of every row
  nrow_sero_study_age = nrow_sero_study_age, #number of age in every row
  pop_v = as.numeric(pop_v), # population for every age row
  cov_v = as.numeric(cov_v), #coverage for every age row
  N_row_sero_age_range = 1:length(age_range_v),
  Nrow_sero = Nrow_sero, #number of row in sero dataset 
  n_row_sero_start = n_row_sero_start, 
  n_row_sero_stop = n_row_sero_stop,
  n_row_sero = n_row_sero,
  sero_pos = sero_data$POSITIVE,
  sero_size = sero_data$SAMPLE_SIZE,
  
  N_case_study = N_case_study, #number of case LR studies
  N_case_study_id = N_sero_study + 1:N_case_study, #N_sero_study + N_case_study + 1:N_case_study
  N_nrow_case_study = length(nrow_case_study), #length(nrow_case_study)
  nrow_case_study = nrow_case_study,
  FOI_case_id_row = FOI_case_id_row,
  nrow_FOI_append = max(FOI_case_id_row) - 100,
  length_case_study = length_case_study,
  length_case_study_start = length_case_study_start,
  length_case_study_stop = length_case_study_stop,
  Nrow_case = length(study_agegroup_org), 
  pop = case.noti.df$pop, 
  cov = case.noti.df$cov,
  study_agegroup_org = study_agegroup_org,
  study_agegroup_start = study_agegroup_start,
  study_agegroup_stop = study_agegroup_stop,
  case_by_study_agegroup = case_by_study_agegroup,
  
  #prior
  a_VE = a_VE, b_VE = b_VE,
  reporting_prop_prior = reporting_prop_prior,
  gamma_priors_mean = prior_gamma_input$X1,
  gamma_priors_sd = prior_gamma_input$X2
)
}

stan_fit <- stan(
  model_code = stan_code,  # Stan program
  data = data_stan,    # named list of data
  chains = 4,             # number of Markov chains
  warmup = 5e1,          # number of warmup iterations per chain
  iter = 1e2,            # total number of iterations per chain
  cores = 4,  # number of cores (could use one per chain)
  thin = 1,
  control = list(adapt_delta = 0.80, stepsize = 5)
  #refresh = 0             # no progress shown
)
traceplot(stan_fit, pars = c("age_dep1", "age_dep2", "age_dep3", "lp__"))
extract(stan_fit)
saveRDS(stan_fit, file = paste0("../output/MCMC_output/MCMC_output_1e4_unipeakprior_sero_case_age_depend_", model_type, "_datasubset_", data_subset, "_rstan.Rdata"))

#Results from crc ####
library(dplyr)

setwd("/home/quan/Documents/NotreDame/Maggie/Time_varying_FOI/")

all_model_type = c('Constant', "One_outbreak", "Two_outbreaks")

#summary par:

for(model_type in all_model_type){
  a = readRDS(paste0("output/MCMC_output/MCMC_output_sero_case_age_depend_", model_type,"_rstan.Rdata"))
  plot(a$age_dep1)
}


#One outbreak model:####
if(model_type == "One_outbreak"){
  stan_code = "
data {
  int<lower=0> N;          // number studies
  vector[100] age_seq; // age sequence from 0 to 99
  matrix[N, 100] time_FOI_empty;
  
  int<lower=0> N_sero_study;          // number sero studies
  int N_sero_study_id[N_sero_study]; // id of sero studies
  int<lower=0> N_row_sero_age; //number of rows in every studies for every age
  int age_range_v[N_row_sero_age]; //age range of every row
  int nrow_sero_study_age[N_row_sero_age]; //number of age in every row
  vector[N_row_sero_age] pop_v; // population for every age row
  vector[N_row_sero_age] cov_v; // coverage for every age row
  int N_row_sero_age_range[N_row_sero_age];
  int<lower=0> Nrow_sero; //number of row in sero dataset 
  int n_row_sero_start[Nrow_sero]; 
  int n_row_sero_stop[Nrow_sero];
  int n_row_sero[Nrow_sero];
  int sero_pos[Nrow_sero];
  int sero_size[Nrow_sero];
  
  int<lower=0> N_case_study; //number of case LR studies
  int N_case_study_id[N_case_LR_study]; //N_sero_study + N_case_study + 1:N_case_LR_study
  int<lower=0> N_nrow_case_LR_study; //length(nrow_case_LR_study)
  int<lower=0> nrow_case_LR_study[N_nrow_case_LR_study];
  int<lower=0> FOI_case_LR_id_row[N_nrow_case_LR_study];
  int length_case_LR_study[N_case_LR_study];
  int length_case_LR_study_start[N_case_LR_study];
  int length_case_LR_study_stop[N_case_LR_study];
  vector[N_nrow_case_LR_study] pop; //Population
  vector[N_nrow_case_LR_study] cov; //coverage
  int<lower=0> Nrow_case_LR; 
  int<lower=0> study_agegroup_org[Nrow_case_LR];
  int<lower=0> study_agegroup_start[Nrow_case_LR];
  int<lower=0> study_agegroup_stop[Nrow_case_LR];
  int<lower=0> case_by_study_agegroup[Nrow_case_LR];
  
  //prior
  real a_VE;
  real b_VE;
  vector[2] reporting_prop_prior;
}

parameters {
  //real<lower=-20, upper = 20> FOI[N];
  vector<lower=-20, upper = 20>[N] alpha;
  vector<lower=1, upper = 98>[N] time;
  real<lower=1, upper = 98> age_dep1; //age dependent: position of the peak age
  real<lower=0, upper = 100> age_dep2; //age dependent: spread of the peak
  real<lower=-100, upper = 100> age_dep3; //age dependent: skewness
  vector<lower=-20, upper = 20>[N_case_LR_study] rho_case_LR;
  real<lower=0, upper = 1> VE;
}

model {
  vector[N] trans_alpha; //transforming alpha to [0,1]
  matrix[N, 100] time_FOI;
  vector[100] age_dep_output_org;
  vector[100] age_dep_output;
  vector[100] age_dep_output_scale;
  matrix[100, N] age_FOI;
  matrix[100,N] FOI_m;
  matrix[100, N_sero_study] FOI_sero_m;
  vector[N_row_sero_age] FOI_sero_cumsum;
  vector[N_row_sero_age] seropos;
  real p[Nrow_sero];
  real p_ele;

  vector[N_case_LR_study] trans_rho_case_LR;
  vector[N_nrow_case_LR_study] l_rho_case_LR;
  matrix[100, N_case_LR_study] FOI_case_LR_m;
  matrix[6, N_case_LR_study] FOI_case_LR_m_appending;
  matrix[106, N_case_LR_study] FOI_case_LR_m_appended;
  vector[N_nrow_case_LR_study] FOI_case_LR_v;
  vector[N_nrow_case_LR_study] FOI_case_LR_int_v; 
  vector[N_nrow_case_LR_study] exp_cases_LR_prob;
  vector[N_nrow_case_LR_study] exp_cases_LR;
  vector[Nrow_case_LR] exp_case_LR_by_age_group;
  vector[Nrow_case_LR] exp_cases_LR_mod;
  
  //Time dependent component of FOI:
  for (j in 1:N) {
    trans_alpha[j] = pow(alpha[j], 10);
  }
  time_FOI = time_FOI_empty;
  for (study in 1:N){
    int whole_time = 1;
    real part_time;
    while (whole_time < (time[study] - 1)) {
      whole_time += 1;
    }
    part_time = time[study] - whole_time;
    time_FOI[study, whole_time + 1] = (1 - part_time) .* trans_alpha[study];
    time_FOI[study, whole_time + 2] = part_time .* trans_alpha[study];
  }
  
  //Age dependent component of FOI:
  for(age in 1:100){
    age_dep_output_org[age] = exp(skew_normal_lpdf(age_seq[age]|age_dep1, age_dep2, age_dep3)); //exp(-((age_seq[age] - age_dep1)^2) ./ (age_dep2^2));
  }
  //age_dep_output_org = exp(-(pow(rep_vector(2.0,100), age_seq - age_dep1))./exp(age_dep2));
  age_dep_output = age_dep_output_org ./ sum(age_dep_output_org);
  age_FOI = rep_matrix(age_dep_output, N); //getting the age dependent component FOI
  
  //Combine to get the FOI matrix:
  FOI_m = (time_FOI').*age_FOI;
  
  //LL of sero data:
  FOI_sero_m = FOI_m[, N_sero_study_id];
  for(i in 1:N_row_sero_age){
    FOI_sero_cumsum[i] = sum(head(FOI_sero_m[, nrow_sero_study_age[i]], age_range_v[i] + 1));
  }
  seropos = 1-exp(-FOI_sero_cumsum);
  for(iii in 1:Nrow_sero){
    int n_row_sero_seq[n_row_sero[iii]];
    n_row_sero_seq = N_row_sero_age_range[n_row_sero_start[iii]:n_row_sero_stop[iii]];
    p_ele = sum(pop_v[n_row_sero_seq] .* (1 - cov_v[n_row_sero_seq]*VE) .* seropos[n_row_sero_seq])/sum(pop_v[n_row_sero_seq] .* (1 - cov_v[n_row_sero_seq]*VE));
    if (p_ele > 1 - 1e-16)
      p_ele = 1 - 1e-16;
    if (p_ele < 1e-16)
      p_ele = 1e-16; 
    p[iii] = p_ele;
  }
  
  //LL from case LR data:
  trans_rho_case_LR = inv_logit(rho_case_LR);
  FOI_case_LR_m = FOI_m[, N_case_LR_study_id];
  FOI_case_LR_m_appending = rep_matrix(FOI_case_LR_m[100,], 6);
  FOI_case_LR_m_appended = append_row(FOI_case_LR_m, FOI_case_LR_m_appending);
  for(xi in 1:N_nrow_case_LR_study){
    int FOI_case_LR_id_int = FOI_case_LR_id_row[xi] - 1;
    
    FOI_case_LR_v[xi] = FOI_case_LR_m_appended[FOI_case_LR_id_row[xi], nrow_case_LR_study[xi]];
    
    if(FOI_case_LR_id_int == 0)
      FOI_case_LR_int_v[xi] = 0.0;
    else
      FOI_case_LR_int_v[xi] = sum(head(FOI_case_LR_m_appended[, nrow_case_LR_study[xi]], FOI_case_LR_id_int));
  }
  
  l_rho_case_LR = trans_rho_case_LR[nrow_case_LR_study];
  exp_cases_LR_prob = exp(-FOI_case_LR_int_v).*(1 - exp(-FOI_case_LR_v)).*l_rho_case_LR;
  exp_cases_LR = exp_cases_LR_prob .* (1 - cov*VE) .* pop;
  for(xviii in 1:Nrow_case_LR){
    exp_case_LR_by_age_group[xviii] = sum(exp_cases_LR[study_agegroup_start[xviii]:study_agegroup_stop[xviii]]);
  }

  for(xx in 1:(Nrow_case_LR)){
    if(exp_case_LR_by_age_group[xx] < 1e-10)
      exp_cases_LR_mod[xx] = 1e-10;
    else 
      exp_cases_LR_mod[xx] = exp_case_LR_by_age_group[xx];
  }
  
  // prior:
  alpha ~ normal(-3, 1.5);
  time ~ uniform(1, 98);
  age_dep1 ~ uniform(1, 98);
  age_dep2 ~ exponential(1);
  age_dep3 ~ normal(0, 10);
  rho_case_LR ~ normal(reporting_prop_prior[1], reporting_prop_prior[2]);

  // LL:
  target += binomial_lpmf(sero_pos| sero_size, p); //from sero data
  target += poisson_lpmf(case_by_study_agegroup| exp_cases_LR_mod);
}
"
}
#Two outbreaks model:####
if(model_type == "Two_outbreaks"){
  
  stan_code = "
data {
  int<lower=0> N;          // number studies
  vector[100] age_seq; // age sequence from 0 to 99
  matrix[N, 100] time_FOI_empty;
  
  int<lower=0> N_sero_study;          // number sero studies
  int N_sero_study_id[N_sero_study]; // id of sero studies
  int<lower=0> N_row_sero_age; //number of rows in every studies for every age
  int age_range_v[N_row_sero_age]; //age range of every row
  int nrow_sero_study_age[N_row_sero_age]; //number of age in every row
  vector[N_row_sero_age] pop_v; // population for every age row
  vector[N_row_sero_age] cov_v; // coverage for every age row
  int N_row_sero_age_range[N_row_sero_age];
  int<lower=0> Nrow_sero; //number of row in sero dataset 
  int n_row_sero_start[Nrow_sero]; 
  int n_row_sero_stop[Nrow_sero];
  int n_row_sero[Nrow_sero];
  int sero_pos[Nrow_sero];
  int sero_size[Nrow_sero];
  
  int<lower=0> N_case_LR_study; //number of case LR studies
  int N_case_LR_study_id[N_case_LR_study]; //N_sero_study + N_case_study + 1:N_case_LR_study
  int<lower=0> N_nrow_case_LR_study; //length(nrow_case_LR_study)
  int<lower=0> nrow_case_LR_study[N_nrow_case_LR_study];
  int<lower=0> FOI_case_LR_id_row[N_nrow_case_LR_study];
  int length_case_LR_study[N_case_LR_study];
  int length_case_LR_study_start[N_case_LR_study];
  int length_case_LR_study_stop[N_case_LR_study];
  vector[N_nrow_case_LR_study] pop; //Population
  vector[N_nrow_case_LR_study] cov; //coverage
  int<lower=0> Nrow_case_LR; 
  int<lower=0> study_agegroup_org[Nrow_case_LR];
  int<lower=0> study_agegroup_start[Nrow_case_LR];
  int<lower=0> study_agegroup_stop[Nrow_case_LR];
  int<lower=0> case_by_study_agegroup[Nrow_case_LR];
  
  //prior
  real a_VE;
  real b_VE;
  vector[2] reporting_prop_prior;
}

parameters {
  //real<lower=-20, upper = 20> FOI[N];
  vector<lower=-20, upper = 20>[N] alpha1;
  vector<lower=1, upper = 98>[N] time1;
  vector<lower=-20, upper = 20>[N] alpha2;
  vector<lower=1, upper = 98>[N] time2;
  real<lower=1, upper = 98> age_dep1; //age dependent: position of the peak age
  real<lower=0, upper = 100> age_dep2; //age dependent: spread of the peak
  real<lower=-100, upper = 100> age_dep3; //age dependent: skewness
  vector<lower=-20, upper = 20>[N_case_LR_study] rho_case_LR;
  real<lower=0, upper = 1> VE;
}

model {
  vector[N] trans_alpha1; //transforming alpha to [0,1]
  vector[N] trans_alpha2; //transforming alpha to [0,1]
  matrix[N, 100] time_FOI;
  vector[100] age_dep_output_org;
  vector[100] age_dep_output;
  vector[100] age_dep_output_scale;
  matrix[100, N] age_FOI;
  matrix[100,N] FOI_m;
  matrix[100, N_sero_study] FOI_sero_m;
  vector[N_row_sero_age] FOI_sero_cumsum;
  vector[N_row_sero_age] seropos;
  real p[Nrow_sero];
  real p_ele;

  vector[N_case_LR_study] trans_rho_case_LR;
  vector[N_nrow_case_LR_study] l_rho_case_LR;
  matrix[100, N_case_LR_study] FOI_case_LR_m;
  matrix[6, N_case_LR_study] FOI_case_LR_m_appending;
  matrix[106, N_case_LR_study] FOI_case_LR_m_appended;
  vector[N_nrow_case_LR_study] FOI_case_LR_v;
  vector[N_nrow_case_LR_study] FOI_case_LR_int_v; 
  vector[N_nrow_case_LR_study] exp_cases_LR_prob;
  vector[N_nrow_case_LR_study] exp_cases_LR;
  vector[Nrow_case_LR] exp_case_LR_by_age_group;
  vector[Nrow_case_LR] exp_cases_LR_mod;
  
  //Time dependent component of FOI:
  for (j in 1:N) {
    trans_alpha1[j] = pow(10.0, alpha1[j]);
    trans_alpha2[j] = pow(10.0, alpha2[j]);
  }
  time_FOI = time_FOI_empty;
  for (study in 1:N){
    int whole_time1 = 1;
    real part_time1;
    int whole_time2 = 1;
    real part_time2;
    while (whole_time1 < (time1[study] - 1)) {
      whole_time1 += 1;
    }
    while (whole_time2 < (time2[study] - 1)) {
      whole_time2 += 1;
    }
    part_time1 = time1[study] - whole_time1;
    part_time2 = time2[study] - whole_time2;
    
    time_FOI[study, whole_time1 + 1] = (1 - part_time1) .* trans_alpha1[study];
    time_FOI[study, whole_time1 + 2] = part_time1 .* trans_alpha1[study];
    
    time_FOI[study, whole_time2 + 1] = (1 - part_time2) .* trans_alpha2[study];
    time_FOI[study, whole_time2 + 2] = part_time2 .* trans_alpha2[study];
  }

  //Age dependent component of FOI:
  for(age in 1:100){
    age_dep_output_org[age] = exp(skew_normal_lpdf(age_seq[age]|age_dep1, age_dep2, age_dep3)); //exp(-((age_seq[age] - age_dep1)^2) ./ (age_dep2^2));
  }
  //age_dep_output_org = exp(-(pow(rep_vector(2.0,100), age_seq - age_dep1))./exp(age_dep2));
  age_dep_output = age_dep_output_org ./ sum(age_dep_output_org);
  age_FOI = rep_matrix(age_dep_output, N); //getting the age dependent component FOI
  
  //Combine to get the FOI matrix:
  FOI_m = (time_FOI').*age_FOI;
  
  //LL of sero data:
  FOI_sero_m = FOI_m[, N_sero_study_id];
  for(i in 1:N_row_sero_age){
    FOI_sero_cumsum[i] = sum(head(FOI_sero_m[, nrow_sero_study_age[i]], age_range_v[i] + 1));
  }
  seropos = 1-exp(-FOI_sero_cumsum);
  for(iii in 1:Nrow_sero){
    int n_row_sero_seq[n_row_sero[iii]];
    n_row_sero_seq = N_row_sero_age_range[n_row_sero_start[iii]:n_row_sero_stop[iii]];
    p_ele = sum(pop_v[n_row_sero_seq] .* (1 - cov_v[n_row_sero_seq]*VE) .* seropos[n_row_sero_seq])/sum(pop_v[n_row_sero_seq] .* (1 - cov_v[n_row_sero_seq]*VE));
    if (p_ele > 1 - 1e-16)
      p_ele = 1 - 1e-16;
    if (p_ele < 1e-16)
      p_ele = 1e-16; 
    p[iii] = p_ele;
  }
  //LL from case LR data:
  trans_rho_case_LR = inv_logit(rho_case_LR);
  FOI_case_LR_m = FOI_m[, N_case_LR_study_id];
  FOI_case_LR_m_appending = rep_matrix(FOI_case_LR_m[100,], 6);
  FOI_case_LR_m_appended = append_row(FOI_case_LR_m, FOI_case_LR_m_appending);
  for(xi in 1:N_nrow_case_LR_study){
    int FOI_case_LR_id_int = FOI_case_LR_id_row[xi] - 1;
    
    FOI_case_LR_v[xi] = FOI_case_LR_m_appended[FOI_case_LR_id_row[xi], nrow_case_LR_study[xi]];
    
    if(FOI_case_LR_id_int == 0)
      FOI_case_LR_int_v[xi] = 0.0;
    else
      FOI_case_LR_int_v[xi] = sum(head(FOI_case_LR_m_appended[, nrow_case_LR_study[xi]], FOI_case_LR_id_int));
  }
  
  l_rho_case_LR = trans_rho_case_LR[nrow_case_LR_study];
  exp_cases_LR_prob = exp(-FOI_case_LR_int_v).*(1 - exp(-FOI_case_LR_v)).*l_rho_case_LR;
  exp_cases_LR = exp_cases_LR_prob .* pop .* (1 - cov * VE);
  for(xviii in 1:Nrow_case_LR){
    exp_case_LR_by_age_group[xviii] = sum(exp_cases_LR[study_agegroup_start[xviii]:study_agegroup_stop[xviii]]);
  }

  for(xx in 1:(Nrow_case_LR)){
    if(exp_case_LR_by_age_group[xx] < 1e-10)
      exp_cases_LR_mod[xx] = 1e-10;
    else 
      exp_cases_LR_mod[xx] = exp_case_LR_by_age_group[xx];
  }
  
  // prior:
  alpha1 ~ normal(-3, 1.5);
  time1 ~ uniform(1, 98);
  alpha2 ~ normal(-3, 1.5);
  time2 ~ uniform(1, 98);
  age_dep1 ~ uniform(1, 98);
  age_dep2 ~ exponential(1);
  age_dep3 ~ normal(0, 10);
  rho_case_LR ~ normal(reporting_prop_prior[1], reporting_prop_prior[2]);

  // LL:
  target += binomial_lpmf(sero_pos| sero_size, p); //from sero data
  target +=  poisson_lpmf(case_by_study_agegroup| exp_cases_LR_mod);
}
"
}