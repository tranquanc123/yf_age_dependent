library(BayesianTools)
library(dplyr)
library(Hmisc)
library(matrixStats)
library(extraDistr)
library(sn)
library(stringr)

#setwd("/home/quan/Documents/NotreDame/Alex/Yellow fever/yf_age_exposure/code/")

#running parallel in crc, which.scenario is the index:
args = commandArgs(trailingOnly=TRUE)
which.scenario = as.numeric(args[1])

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
  sero_data = sero_data[NULL,]
  case.noti.df = rbind(filter(case.noti.LR.df, country == "BRA"), case.noti.PAHO.df)
} else if(which.scenario == 2){
  case.noti.df = filter(case.noti.LR.df, country != "BRA")
} else {
  case.noti.df = rbind(case.noti.LR.df, case.noti.PAHO.df)
}


#Utility functions:
source("util_function.R")

#model runs:
model_type = "Constant"
#Parameter id:
par_id(model_type)
#FOI function: outbreak model:
FOI_function_model_type(model_type)

if(which.scenario != 2){
  case.noti.df[is.na(case.noti.df$cov),]$pop = 1e-10
  case.noti.df[is.na(case.noti.df$cov),]$cov = 0
}

#LL function:
LL_func <- function(par){
  generated_cases_sero = gen_func(par)
  
  if(data_subset == "SA"){
    LL_sero = 0
  } else {
    #Serological data:
    LL_sero = dbinom(sero_data$POSITIVE, sero_data$SAMPLE_SIZE, unlist(generated_cases_sero$sero), log = T) %>% sum
  }
  
  #Case LR data:
  LL_case = sum(dpois(x = case_by_study_agegroup, lambda = generated_cases_sero$case, log = T))
  
  return(LL_sero + LL_case)
}

###MCMC run
#prior for gamma:
coun_list_study = sapply(case_study_id, function(x) filter(case.noti.df, study_id == x)$country %>% unique)
coun_list_study[!(coun_list_study %in% prior_gamma$ISO)] <- "BRA"
sero_list_study = sapply(sero_study_id, function(x) filter(sero_data, STUDY_ID == x)$ISO %>% unique)
match_prior_gamma = match(c(sero_list_study, coun_list_study), prior_gamma$ISO)
prior_gamma_input = prior_gamma[match_prior_gamma,]


par_prior_dens = function(par){
  prior = dnorm(par[Age_depend_id][1], 50, 125, log = T) + 
    dhnorm(par[Age_depend_id][2], 150, log = T) + 
    dt(par[Age_depend_id][3], 2, 0.5, log = T) +
    dbeta(par[VE_id], a_VE, b_VE, log = T) +
    sum(dnorm(par[c(rho_case_id)], reporting_prop_prior[1], reporting_prop_prior[2], log = T))
  if(model_type == "Constant"){
    prior = prior + sum(dlnorm(par[c(FOI_par_id)], prior_gamma_input$X1, prior_gamma_input$X2, log = T))
  } else if(model_type == "One outbreak"){
    prior = prior + sum(dnorm(par[c(alpha_id)], 0, 10, log = T)) + sum(dunif(par[time_id], 0, 99, log = T)) 
  } else {
    prior = prior + sum(dnorm(par[c(alpha1_id, alpha2_id)], 0, 10, log = T)) + sum(dunif(par[c(time1_id, time2_id)], 1, 99, log = T))
  }
  return(prior)
}

init_values = function(n = 1) {
  all_initial_value = cbind(rnorm(n, 50, 125), rhnorm(n, 150), rt(n, 2, 0.5),
                            matrix(rnorm(n*(N_case_study), reporting_prop_prior[1], reporting_prop_prior[2]), nrow = 1),
                            rbeta(n, a_VE, b_VE))
  if(model_type == "Constant"){
    all_initial_value = cbind(matrix(rlnorm(n*N_study,  prior_gamma_input$X1, prior_gamma_input$X2), nrow = 1), all_initial_value)
  } else if(model_type == "One outbreak"){
    all_initial_value = cbind(matrix(c(rnorm(n*N_study,  0, 10), runif(n*N_study, 0, 99)), nrow = 1), all_initial_value)
  } else {
    all_initial_value = cbind(matrix(c(rnorm(n*N_study,  0, 10), runif(n*N_study, 0, 99), rnorm(n*N_study,  0, 10), runif(n*N_study, 0, 99)), nrow = 1), all_initial_value)
  }
  return(all_initial_value)
}

l_prior = c(c(-Inf, 0, -Inf), rep(-Inf, N_case_study), 0)
u_prior = c(c(Inf, Inf, Inf), rep(Inf, N_case_study), 1) 
m_prior =  c(c(0, 1, 0), rep(0, N_case_study), 0.5)

if(model_type == "Constant"){
  l_prior = c(rep(-Inf, N_study), l_prior)
  u_prior = c(rep(Inf, N_study), u_prior)
  m_prior = c(rep(0, N_study), m_prior)
  
  l_prior[FOI_par_id] = 0
  u_prior[FOI_par_id] = Inf
  m_prior[FOI_par_id] = 1
} else if(model_type == "One outbreak"){
  l_prior = c(rep(c(-Inf, 1), each = N_study), l_prior)
  u_prior = c(rep(c(Inf, 99), each = N_study), u_prior)
  m_prior = c(rep(c(0, 50), each = N_study), m_prior)
} else {
  l_prior = c(rep(rep(c(-Inf, 1), each = N_study), 2), l_prior)
  u_prior = c(rep(rep(c(Inf, 99), each = N_study), 2), u_prior)
  m_prior = c(rep(rep(c(0, 50), each = N_study), 2), m_prior)
}

par_prior = createPrior(density = par_prior_dens, sampler = init_values,
                        lower = l_prior, upper = u_prior, best =  m_prior)

bayesanSetup <- createBayesianSetup(likelihood = LL_func, prior = par_prior)
iter = 1e4
settings = list(iterations = iter)
a <- runMCMC(bayesianSetup = bayesanSetup, sampler = "DREAMzs", settings = settings)
saveRDS(a, file = paste0("../output/MCMC_output/MCMC_PAHO_1e7_prior_SA_vs_AF_output_sero_case_", model_type, "_datasubset_", data_subset, ".Rdata"))