library(BayesianTools)
library(dplyr)
library(Hmisc)
library(matrixStats)
library(sn)
library(stringr)

setwd("/home/quan/Documents/NotreDame/Alex/Yellow fever/yf_age_exposure/code/")

#running parallel in crc, which.scenario is the index:
args = commandArgs(trailingOnly=TRUE)
which.scenario = as.numeric(args[1])

#Load data:
sero_data = read.csv("../data/yf_sero_data_with_coverage_mod.csv")
Study_id_seq = rle(as.character(sero_data$STUDY_ID))
one_age_group_data = Study_id_seq$values[Study_id_seq$lengths == 1]
sero_data = filter(sero_data, !(STUDY_ID %in% one_age_group_data))
case.noti.LR.df = read.csv("../data/yf_case_data_by_age_group_with_coverage_LR.csv")

VE = 0.975

#All model runs:
all_model_type = c("Constant", "One outbreak", "Two outbreaks")

#Utility functions:
source("util_function.R")

#model runs:
model_type = all_model_type[which.scenario]
#Parameter id:
par_id(model_type)
#FOI function: outbreak model:
FOI_function_model_type(model_type)

#LL function:
LL_func <- function(par){
  generated_cases_sero = gen_func(par)
  #Serological data:
  LL_sero = dbinom(sero_data$POSITIVE, sero_data$SAMPLE_SIZE, unlist(generated_cases_sero$sero), log = T) %>% sum
  #Case LR data:
  LL_case_LR =sum(dpois(x = case_by_study_agegroup, lambda = generated_cases_sero$case_LR, log = T)) %>% sum
  
  #Hierarchical terms:
  if(model_type == "One outbreak"){
    LL_hierarchical = sum(dnorm(par[FOI_par_id][alpha_id], par[alpha_avg_id], par[alpha_sd_id], log = T))
  } else if(model_type == "Two outbreaks") {
    LL_hierarchical = sum(dnorm(par[FOI_par_id][c(alpha1_id, alpha2_id)], par[alpha_avg_id], par[alpha_sd_id], log = T))
  } else {
    LL_hierarchical = 0
  }
  return(LL_sero + LL_case_LR + LL_hierarchical)
}

###MCMC run
if(model_type == "Constant"){
  par_prior_dens = function(par){
    sum(dnorm(par[c(FOI_par_id, FOI_avg_id)], -3, 1, log = T)) + dexp(par[FOI_sd_id], 1, log = T) +
      dunif(par[Age_depend_id][1], 0, 99, log = T) + dnorm(par[Age_depend_id][2], 0, 1, log = T) + dnorm(par[Age_depend_id][3], 0, 1, log = T) +
      sum(dnorm(par[c(rho_case_LR_id)], 0, 1.5, log = T))
  }
  init_values = function(n = 1) c(rnorm(n*N_study, -3, 2), rnorm(n, -3, 2), rexp(n, 1), 
                                  c(runif(n, 0, 99), rnorm(n, 0, 1), rnorm(n, 0, 1)),
                                  rnorm(n*(N_case_LR_study), 0, 1.5))
  par_prior = createPrior(density = par_prior_dens, sampler = init_values,
                          lower = c(rep(-Inf, N_study), -Inf, 0, c(0, -Inf, -Inf), rep(-Inf, N_case_LR_study)), 
                          upper = c(rep(Inf, N_study), Inf, 10, c(99, Inf, Inf), rep(Inf, N_case_LR_study)), 
                          best =  c(rep(0, N_study), 0, 1, c(50, 0, 0), rep(0, N_case_LR_study)))
} else if(model_type == "One outbreak"){
  par_prior_dens = function(par){
    sum(dnorm(par[c(alpha_id, alpha_avg_id)], -3, 1, log = T)) + sum(dunif(par[time_id], 0, 99, log = T)) + dexp(par[alpha_sd_id], 1, log = T) +
      dunif(par[Age_depend_id][1], 0, 99, log = T) + dnorm(par[Age_depend_id][2], 0, 1, log = T) + dnorm(par[Age_depend_id][3], 0, 1, log = T) +
      sum(dnorm(par[c(rho_case_LR_id)], 0, 1.5, log = T))
  }
  init_values = function(n = 1) {c(rnorm(n*N_study, -3, 2), runif(n*N_study, 1, 99), rnorm(n, -3, 2), rexp(n, 1),
                                   c(runif(n, 1, 99), rnorm(n, 0, 1), rnorm(n, 0, 1)), rnorm(n*(N_case_LR_study), 0, 1.5))}
  par_prior = createPrior(density = par_prior_dens, sampler = init_values,
                          lower = c(rep(c(-Inf, 1), each = N_study),-Inf, 0, c(0, -Inf, -Inf), rep(-Inf, N_case_LR_study)),
                          upper = c(rep(c(Inf, 99), each = N_study), Inf, 10, c(99, Inf, Inf), rep(Inf, N_case_LR_study)), 
                          best =  c(rep(c(0, 50), each = N_study), 0, 1, c(50, 0, 0), rep(0, N_case_LR_study)))
} else {
  par_prior_dens = function(par){
    sum(dnorm(par[c(alpha1_id, alpha2_id, alpha_avg_id)], -3, 1, log = T)) + sum(dunif(par[c(time1_id, time2_id)], 1, 99, log = T)) + dexp(par[alpha_sd_id], 1, log = T) +
      dunif(par[Age_depend_id][1], 0, 99, log = T) + dnorm(par[Age_depend_id][2], 0, 1, log = T) + dnorm(par[Age_depend_id][3], 0, 1, log = T)+
      sum(dnorm(par[c(rho_case_LR_id)], 0, 1.5, log = T))
  }
  init_values = function(n = 1){
    c(rnorm(n*N_study, -3, 2), runif(n*N_study, 1, 99), rnorm(n*N_study, -3, 2), runif(n*N_study, 1, 99),
      rnorm(n, -3, 2), rexp(n, 1), c(runif(n, 1, 99), rnorm(n, 0, 1), rnorm(n, 0, 1)),
      rnorm(n*(N_case_LR_study), 0, 1.5))
  }
  par_prior = createPrior(density = par_prior_dens, sampler = init_values,
                          lower = c(rep(rep(c(-Inf, 1), each = N_study), 2), -Inf, 0, c(1, -Inf, -Inf),
                                    rep(-Inf, N_case_LR_study)), 
                          upper = c(rep(rep(c(Inf, 99), each = N_study), 2), Inf, Inf, c(99, Inf, Inf),
                                    rep(Inf, N_case_LR_study)), 
                          best =  c(rep(rep(c(0, 50), each = N_study), 2) , 0, 1, c(50, 0, 0),
                                    rep(0, N_case_LR_study)))
}

bayesanSetup <- createBayesianSetup(likelihood = LL_func, prior = par_prior)
iter = 2e7
settings = list(iterations = iter)
a <- runMCMC(bayesianSetup = bayesanSetup, sampler = "DREAMzs", settings = settings)
saveRDS(a, file = paste0("../output/MCMC_output/MCMC_output_sero_case_", model_type, ".Rdata"))