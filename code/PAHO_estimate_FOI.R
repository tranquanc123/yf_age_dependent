setwd("/home/quan/Documents/NotreDame/Alex/Yellow fever/yf_age_exposure/code/")

library(BayesianTools)
library(dplyr)
library(sn)
library(stringr)
library(rstan)
library(matrixStats)

args = commandArgs(trailingOnly=TRUE)
which.scenario = as.numeric(args[1])

#Load data:
sero_data = read.csv("../data/yf_sero_data_with_coverage_mod.csv")
sero_data = filter(sero_data, YEAR > 1980) #only get data from 1980
sero_study_id = rle(as.character(sero_data$STUDY_ID))
sero_data = filter(sero_data, STUDY_ID %in% sero_study_id$values[which(sero_study_id$lengths > 1)]) # study with age group > 1
case.noti.LR.df = read.csv("../data/yf_case_data_by_age_group_with_coverage_LR.csv")
case.noti.PAHO.df = read.csv("../data/yf_case_data_by_age_group_with_coverage_PAHO.csv")# %>% filter(!(study_id %in% c("ECU_PAHO", "ARG_PAHO")))

#Prior:
readRDS("../data/FOI_rho_prior.rds") %>% list2env(globalenv())
prior_gamma = readRDS("../data/gamma_lnorm_priors_constant_nonorm.rds")
a_VE = 22.812; b_VE = 1.306

#all scenarios
SA_study_subset = c(unique(filter(case.noti.LR.df, country == "BRA")$study_id), unique(case.noti.PAHO.df$study_id))
AF_study_subset = c(unique(sero_data$STUDY_ID), unique(filter(case.noti.LR.df, country != "BRA")$study_id))
all_model_type = c("Alternative", "Constant")
all_data_subset  = c("SA", "AF")

all_scen = expand.grid(model_type = all_model_type, data_subset = all_data_subset)
all_scen = all_scen[unlist(mapply(rep, 1:4, c(length(SA_study_subset), length(SA_study_subset), length(AF_study_subset), length(AF_study_subset)))),]
all_scen$study = c(SA_study_subset, SA_study_subset, AF_study_subset, AF_study_subset)
all_scen$all_data_type = c(rep(rep("case", length(SA_study_subset)), 2), rep(c(rep("sero", 4), rep("case", length(AF_study_subset) - 4)), 2))
rownames(all_scen) = 1:nrow(all_scen)

#get a specific scenario:
model_type = all_scen$model_type[which.scenario] # 
data_subset = all_scen$data_subset[which.scenario] #"SA"
study = all_scen$study[which.scenario]

source("util_function.R")
par_id(model_type)
FOI_function_model_type(model_type)

if(data_subset != "AF"){
  case.noti.df[is.na(case.noti.df$cov),]$pop = 1e-10
  case.noti.df[is.na(case.noti.df$cov),]$cov = 0
}

#read in stan object:
data_subset_age = all_data_subset[which(all_data_subset != data_subset)]
if(model_type == "Constant"){
  age_dep = matrix(1/100, nrow = 100, ncol = 1000)
} else if(model_type == "Alternative"){
  age_dep = readRDS(paste0("../output/Age_exp_curve_", data_subset_age, "_Constant.rds"))
}
n_samples = ncol(age_dep)

#LL function
data_type = all_scen$all_data_type[which.scenario]

if(data_type == "sero"){
  sel_study_data = filter(sero_data, STUDY_ID == study)
  study_i_row = which(sero_data$STUDY_ID %in% study) #which(nrow_sero_study_age == which.scenario)
  sel_study_id = which(n_row_sero_l %in% study_i_row)
} else {
  sel_study_data = filter(case.noti.df, study_id == study)
  study_i_row = which(case.noti.df$study_id %in% study)
  case_by_study_agegroup_i = (sel_study_data %>% group_by(age_group_ID) %>% summarise(case = unique(case)))$case
}


LL_validate <- function(par){
  gamma_val = par[1]
  if(data_type == "case") {
    rho_val = inverse_logit(par[2])
    VE_val = par[3]
  } else {
    VE_val = par[2] 
  }
  
  FOI = gamma_val*age_dep
  
  if(data_type == "sero"){
    FOI_sero_cumsum = apply(FOI, 2, cumsum)
    seropos = 1 - exp(-FOI_sero_cumsum[age_range_v[sel_study_id],])
    sus_v = (cov_v[sel_study_id]*VE_val + (1-cov_v[sel_study_id]*VE_val) * seropos) * pop_v[sel_study_id]
    
    p = aggregate(sus_v, list(n_row_sero_l[sel_study_id]), sum)[,-1]/pop_by_row[study_i_row]
    p[p > 1 - 1e-16] = 1 - 1e-16
    p[p < 1e-16] = 1e-16
    
    LL_sero = sapply(1:1e3, function(x) dbinom(sel_study_data$POSITIVE, sel_study_data$SAMPLE_SIZE, p[,x], log = T)) %>% rowMeans() %>% sum
  } else {
    LL_sero = 0
  }
  
  if(data_type == "case"){
    ###Case LR data:
    #case data from LR:
    FOI_case_m = rbind(FOI, matrix(FOI[100,], nrow = FOI_case_id_row_max - 100, ncol = n_samples, byrow = T))
    FOI_case_v = sapply(study_i_row, function(x) FOI_case_m[FOI_case_id_row[x], ]) %>% t
    FOI_case_m_cumsum = apply(FOI_case_m, 2, cumsum)
    FOI_case_int_v = sapply(study_i_row, function(x) {
      FOI_case_id_int = FOI_case_id_row[x] - 1
      if(FOI_case_id_int == 0) return(rep(1, n_samples))
      return(FOI_case_m_cumsum[FOI_case_id_int,]) 
    }) %>% t
    case_sus = sel_study_data$pop*(1 - sel_study_data$cov*VE_val)
    exp_cases = case_sus*exp(-FOI_case_int_v)*(1 - exp(-FOI_case_v))*rho_val
    
    exp_case_by_age_group = aggregate(exp_cases, by = list(study_agegroup[study_i_row]), FUN = sum)[,-1]
    exp_case_by_age_group[exp_case_by_age_group < 1e-10] = 1e-10
    
    LL_case = sapply(1:n_samples, function(x) dpois(case_by_study_agegroup_i, exp_case_by_age_group[,x], log = T)) %>% rowMeans() %>% sum
  } else {
    LL_case = 0
  }
   
  return(LL_sero + LL_case)
}

#model run
if(data_type == 'case'){
  prior_gamma_input = filter(prior_gamma, ISO == unique(sel_study_data$country))
  if(nrow(prior_gamma_input) == 0) prior_gamma_input = filter(prior_gamma, ISO == "BRA")
} else {
  prior_gamma_input = filter(prior_gamma, ISO == unique(sel_study_data$ISO))
}

if(data_type == "case"){
  par_prior_dens = function(par){
    dlnorm(par[1], prior_gamma_input$X1, prior_gamma_input$X2, log = T) + dnorm(par[2], reporting_prop_prior[1], reporting_prop_prior[2], log = T) + dbeta(par[3], a_VE, b_VE, log = T)
  }
  init_values = function(n = 1) c(rlnorm(n, prior_gamma_input$X1, prior_gamma_input$X2), rnorm(n, reporting_prop_prior[1], reporting_prop_prior[2]), rbeta(n, a_VE, b_VE))
  par_prior = createPrior(density = par_prior_dens, sampler = init_values,
                          lower = c(0, -Inf, 0), upper = c(Inf, Inf, 1), best =  c(1, 0, 0.5))
} else {
  par_prior_dens = function(par){
    dlnorm(par[1], prior_gamma_input$X1, prior_gamma_input$X2, log = T) + dbeta(par[2], a_VE, b_VE, log = T)
  }
  init_values = function(n = 1) c(rlnorm(n, prior_gamma_input$X1, prior_gamma_input$X2), rbeta(n, a_VE, b_VE))
  par_prior = createPrior(density = par_prior_dens, sampler = init_values,
                          lower = c(0, 0), upper = c(Inf, 1), best =  c(1, 0.5))
}

bayesanSetup <- createBayesianSetup(likelihood = LL_validate, prior = par_prior)
iter = 2e4
settings = list(iterations = iter)
output <- runMCMC(bayesianSetup = bayesanSetup, sampler = "DREAMzs", settings = settings)
saveRDS(output, paste0("../output/MCMC_output/Validate_", data_subset, "_", model_type,"_", study,"_est.rds"))

#get par from the output:
setwd("/home/quan/Documents/NotreDame/Alex/Yellow fever/yf_age_exposure/code/")

library(BayesianTools)
library(dplyr)
library(sn)
library(stringr)
library(rstan)
library(matrixStats)

#Load data:
sero_data = read.csv("../data/yf_sero_data_with_coverage_mod.csv")
sero_data = filter(sero_data, YEAR > 1980) #only get data from 1980
sero_study_id = rle(as.character(sero_data$STUDY_ID))
sero_data = filter(sero_data, STUDY_ID %in% sero_study_id$values[which(sero_study_id$lengths > 1)]) # study with age group > 1
case.noti.LR.df = read.csv("../data/yf_case_data_by_age_group_with_coverage_LR.csv")
case.noti.PAHO.df = read.csv("../data/yf_case_data_by_age_group_with_coverage_PAHO.csv") #%>% filter(!(study_id %in% c("ECU_PAHO", "ARG_PAHO")))

#all scenarios
SA_study_subset = c(unique(filter(case.noti.LR.df, country == "BRA")$study_id), unique(case.noti.PAHO.df$study_id))
AF_study_subset = c(unique(sero_data$STUDY_ID), unique(filter(case.noti.LR.df, country != "BRA")$study_id))
all_model_type = c("Alternative", "Constant")
all_data_subset  = c("SA", "AF")

all_scen = expand.grid(model_type = all_model_type, data_subset = all_data_subset)
all_scen = all_scen[unlist(mapply(rep, 1:4, c(length(SA_study_subset), length(SA_study_subset), length(AF_study_subset), length(AF_study_subset)))),]
all_scen$study = c(SA_study_subset, SA_study_subset, AF_study_subset, AF_study_subset)
all_scen$all_data_type = c(rep(rep("case", length(SA_study_subset)), 2), rep(c(rep("sero", 4), rep("case", length(AF_study_subset) - 4)), 2))
rownames(all_scen) = 1:nrow(all_scen)
all_scen = all_scen[1:28, ]

gen_validate <- function(par, age_dep_id){
  gamma_val = par[1]
  if(data_type == "case") {
    rho_val = inverse_logit(par[2])
    VE_val = par[3]
  } else {
    VE_val = par[2] 
  }
  
  FOI = gamma_val*age_dep[,age_dep_id]
  
  if(data_type == "sero"){
    FOI_sero_cumsum = cumsum(FOI)
    sel_study_id = which(n_row_sero_l %in% study_i_row)
    seropos = 1 - exp(-FOI_sero_cumsum[age_range_v[sel_study_id]])
    sus_v = (cov_v[sel_study_id]*VE_val + (1-cov_v[sel_study_id]*VE_val) * seropos) * pop_v[sel_study_id]
    
    p = aggregate(sus_v, list(n_row_sero_l[sel_study_id]), sum)$x/pop_by_row[study_i_row]
    p[p > 1 - 1e-16] = 1 - 1e-16
    p[p < 1e-16] = 1e-16
  } else {
    p = NULL
  }
  
  if(data_type == "case"){
    ###Case LR data:
    #case data from LR:
    FOI_case_m = c(FOI, matrix(FOI[100], nrow = FOI_case_id_row_max - 100, byrow = T))
    FOI_case_v = sapply(study_i_row, function(x) FOI_case_m[FOI_case_id_row[x]])
    FOI_case_m_cumsum = cumsum(FOI_case_m)
    FOI_case_int_v = sapply(study_i_row, function(x) {
      FOI_case_id_int = FOI_case_id_row[x] - 1
      if(FOI_case_id_int == 0) return(1)
      return(FOI_case_m_cumsum[FOI_case_id_int]) 
    })
    
    case_sus = sel_study_data$pop*(1 - sel_study_data$cov*VE_val)
    exp_cases = case_sus*exp(-FOI_case_int_v)*(1 - exp(-FOI_case_v))*rho_val
    
    exp_case_by_age_group = aggregate(exp_cases, by = list(study_agegroup[study_i_row]), FUN = sum)[,-1]
    exp_case_by_age_group[exp_case_by_age_group < 1e-10] = 1e-10
  } else {
    exp_case_by_age_group = NULL
  }
  
  return(list(sero = p, case = exp_case_by_age_group))
}

iter = 2e4
n_samples = 1000

study_age_proj_l = c()
study_par_l = c()
WAIC_l = c()
Rhat_l = c()

for(i in 1:nrow(all_scen)){
  print(i)
  data_subset = all_scen$data_subset[i]
  age_exp_subset = setdiff(c("SA", "AF"), data_subset)
  model_type = all_scen$model_type[i]
  data_type = all_scen$all_data_type[i]
  study = all_scen$study[i]
  
  sero_data = read.csv("../data/yf_sero_data_with_coverage_mod.csv")
  sero_data = filter(sero_data, YEAR > 1980) #only get data from 1980
  sero_study_id = rle(as.character(sero_data$STUDY_ID))
  sero_data = filter(sero_data, STUDY_ID %in% sero_study_id$values[which(sero_study_id$lengths > 1)]) # study with age group > 1
  
  source("util_function.R")
  par_id(model_type)
  FOI_function_model_type(model_type)
  
  if(data_subset == "SA"){
    case.noti.df[is.na(case.noti.df$cov),]$pop = 1e-10
    case.noti.df[is.na(case.noti.df$cov),]$cov = 0
  }
  
  output <- readRDS(paste0("../output/MCMC_output/Validate_", data_subset, "_", model_type,"_", study,"_est.rds"))
  par_output = getSample(output, start = round(iter/3/2), thin = 10) %>% tail(n_samples)
  
  #Summary FOI, rho, VE:
  study_gamma = quantile95cri(par_output[,1])
  if(data_type == "case"){
    study_rho = par_output[,2] %>% inverse_logit %>% quantile95cri
    study_VE = par_output[,3] %>% quantile95cri
    study_par = data.frame(data_subset, model_type, study, data_type, 
                           rbind(t(study_gamma), t(study_rho), t(study_VE)), par = c("~italic(gamma)", "~italic(rho)", "~italic(VE)"))
  } else {
    study_VE = par_output[,2] %>% quantile95cri
    study_par = data.frame(data_subset, model_type, study, data_type, 
                           rbind(t(study_gamma), t(study_VE)), par = c("~italic(gamma)", "~italic(VE)"))
  }
  study_par_l = c(study_par_l, list(study_par))
  
  #Project age case pattern and data fit:
  if(model_type == "Constant"){
    age_dep = array(1/100, dim = c(100, n_samples))
  } else {
    age_dep = readRDS(paste0("../output/Age_exp_curve_", age_exp_subset, "_Constant.rds"))
  }
  
  if(data_type == "sero"){
    sel_study_data = filter(sero_data, STUDY_ID == study)
    study_i_row = which(sero_data$STUDY_ID %in% study) #which(nrow_sero_study_age == which.scenario)
  } else {
    sel_study_data = filter(case.noti.df, study_id == study)
    study_i_row = which(case.noti.df$study_id %in% study)
    case_by_study_agegroup_i = (sel_study_data %>% group_by(age_group_ID) %>% summarise(case = unique(case)))$case
  }
  
  #Project case - sero:
  gen_data = lapply(1:n_samples, function(x) gen_validate(unlist(par_output[x,], x)))
  if(data_type == "sero"){
    gen_data_sum = sapply(gen_data, function(x) x$sero) %>% apply(1, quantile95cri)
    study_age_proj = data.frame(data_subset, model_type, study, data_type, sel_study_data, t(gen_data_sum))
  } else {
    gen_data_sum =  sapply(gen_data, function(x) x$case) %>% apply(1, quantile95cri)
    study_age_proj = data.frame(data_subset, model_type, study, data_type, 
                                age_l = aggregate(sel_study_data$age, by = list(sel_study_data$age_group_ID), min)$x, 
                                age_u = aggregate(sel_study_data$age, by = list(sel_study_data$age_group_ID), max)$x, 
                                t(gen_data_sum), case = case_by_study_agegroup_i)
  }
  study_age_proj_l = c(study_age_proj_l, list(study_age_proj))
  
  #WAIC:
  if(data_type == "sero"){
    gen_sero_data = sapply(gen_data, function(x) x$sero)
    LL_samples = sapply(1:n_samples, function(x) dbinom(sel_study_data$POSITIVE, sel_study_data$SAMPLE_SIZE, gen_sero_data[,x], log = T))
    nrow_study = nrow(sel_study_data)
  } else {
    gen_case_data = sapply(gen_data, function(x) x$case)
    LL_samples = sapply(1:n_samples, function(x) dpois(case_by_study_agegroup_i, gen_case_data[,x], log = T))
    nrow_study = length(case_by_study_agegroup_i)
  }
  
  lppd <- sapply(1:nrow_study, function(i) logSumExp(LL_samples[i,]) - log(n_samples)) %>% sum
  pWAIC <- sapply(1:nrow_study , function(i) var(LL_samples[i,])) %>% sum
  WAIC = -2*(lppd - pWAIC)
  WAIC_l = c(WAIC_l, WAIC)
  
  Rhat_cal = sapply(1:ncol(par_output), function(x){
    matrix(par_output[1:(n_samples - 1), x], ncol = 3, byrow = T) %>% Rhat
  })
  Rhat_l = c(Rhat_l, list(Rhat_cal))
}

study_age_sero_proj = study_age_proj_l[which(all_scen$all_data_type == "sero")] %>% Reduce(rbind, .)
study_age_case_proj = study_age_proj_l[which(all_scen$all_data_type == "case")] %>% Reduce(rbind, .)
study_par = Reduce(rbind, study_par_l)
all_scen$WAIC = WAIC_l
Rhat_l[which(sapply(Rhat_l, function(x) any(x > 1.05)))]

write.csv(study_par, "../output/Validate_par.csv", row.names = F)
saveRDS(list(study_age_sero_proj = study_age_sero_proj, study_age_case_proj = study_age_case_proj), "../output/Validate_datafit.rds")
write.csv(all_scen, "../output/WAIC_validate.csv", row.names = F)

#Plotting
library(dplyr)
library(ggplot2)
library(Hmisc)
setwd("/home/quan/Documents/NotreDame/Alex/Yellow fever/yf_age_exposure/code/")

model_order = c("Constant", "Alternative")
model_colors = c("#8073ac", "#e08214")

#WAIC:
WAIC_validate = read.csv("../output/WAIC_validate.csv")
WAIC_validate = (WAIC_validate %>% group_by(study) %>% summarise(model_type = model_type, data_subset = data_subset, study = study, delta_WAIC = WAIC - min(WAIC)))
ggplot(WAIC_validate, aes(y = study, x = log10(delta_WAIC), fill = model_type))+
  geom_bar(stat = "identity", position = position_dodge(width = 0.5))+
  facet_wrap(data_subset~., scale = "free")

#SA only:
WAIC_validate = read.csv("../output/WAIC_validate.csv")
WAIC_validate_SA = WAIC_validate %>% filter(data_subset == "SA")
WAIC_validate_SA$WAIC = round(WAIC_validate_SA$WAIC)
cbind(WAIC_validate_SA[1:14,c(3,5)], WAIC_No_age = WAIC_validate_SA[15:28,5])
#cbind(WAIC_validate_SA[1:12,c(5)], WAIC_No_age = WAIC_validate_SA[13:24,5]) %>% colSums()

# source_names = unique(WAIC_validate$study)
# sero_names = source_names[13:16] %>% strsplit("-") %>% sapply(function(x) paste0(str_to_title(x[4]), " el al ", x[3]))
# case_LR_names = source_names[-c(7:16)] %>% strsplit("_"); author_LR = sapply(case_LR_names, function(x) x[1]) ; author_LR[15] = "De Cock"; author_LR[6] = "Ribeiro"; author_LR[2] = "Saad"; year_LR = sapply(case_LR_names, function(x) tail(x, 1)); year_LR[18:19] <- paste0("1998-", year_LR[18:19])
# case_LR_names = paste0(author_LR, " et al ", year_LR)
# case_PAHO_names = paste0("PAHO - ", c("Bolivia", "Brazil", "Colombia", "Paraguay", "Peru", "Venezuela"))
# source_names = c(case_LR_names[1:6], case_PAHO_names, sero_names, case_LR_names[7:20])

source_names = unique(WAIC_validate$study)
case_LR_names = source_names[1:6] %>% strsplit("_"); author_LR = sapply(case_LR_names, function(x) x[1]) ; author_LR[6] = "Ribeiro"; author_LR[2] = "Saad"; year_LR = sapply(case_LR_names, function(x) tail(x, 1)); 
case_LR_names = paste0(author_LR, " et al ", year_LR)
case_PAHO_names = paste0("PAHO - ", c("Argentina", "Bolivia", "Brazil", "Colombia", "Ecuador", "Paraguay", "Peru", "Venezuela"))
source_names = c(case_LR_names, case_PAHO_names)

#Summary parameters:
Val_par <-read.csv("../output/Validate_par.csv")
Val_par$model_type = factor(Val_par$model_type, levels = model_order)
# Val_par_plot <- ggplot(Val_par, aes(x = X50., y = study, color = model_type))+
#   geom_pointrange(aes(xmin = X2.5., xmax = X97.5.), position = position_dodge(width = 0.7))+
#   facet_grid(data_subset~par, scale = "free", labeller = labeller(par = label_parsed))+
#   scale_color_manual(values = model_colors)+
#   labs(x = "", y = "Study", color = "Models")+
#   theme_bw()+
#   theme(text = element_text(size = 19), 
#         axis.text.y = element_text(size = 14), 
#         axis.text.x = element_text(size = 10),
#         strip.text = element_text(size = 12))
# Val_par_plot
# ggsave("../plots/PAHO_par_plot.png", PAHO_par_plot, width = 8, height = 7)

#SA only:
Val_par_SA = Val_par %>% filter(data_subset == "SA") 
Val_par_SA$study = factor(rep(rep(source_names[1:14], each = 3), 2), levels = c(source_names[1:14]))
Val_par_SA[which(Val_par_SA$par == "~italic(rho)"), c("X2.5.", "X50.", "X97.5.")] <- log10(Val_par_SA[which(Val_par_SA$par == "~italic(rho)"), c("X2.5.", "X50.", "X97.5.")])
Val_par_SA = Val_par_SA %>% mutate(model_type = case_when(model_type == "Alternative" ~ "Africa", TRUE ~ as.character(model_type)),
                                   par = case_when(par == "~italic(rho)" ~ "~italic(rho)~(log10)", TRUE ~ as.character(par)))

Val_par_plot <- ggplot(Val_par_SA, aes(x = X50., y = study, color = model_type))+
  geom_pointrange(aes(xmin = X2.5., xmax = X97.5.), position = position_dodge(width = 0.7))+
  facet_grid(.~par, scale = "free", labeller = labeller(par = label_parsed))+
  scale_color_manual(values = model_colors)+
  labs(x = "", y = "Study", color = "Age exposure")+
  theme_bw()+
  theme(text = element_text(size = 19), 
        axis.text.y = element_text(size = 14), 
        axis.text.x = element_text(size = 10, angle = -90),
        strip.text = element_text(size = 12))
Val_par_plot
ggsave("../plots/Val_par_SA_plot.png", Val_par_plot, width = 8, height = 6)

# Val_age_proj <-readRDS("../output/Validate_datafit.rds")
# Val_age_sero_proj = Val_age_proj$study_age_sero_proj %>% 
#   cbind(CI = binconf(Val_age_proj$study_age_sero_proj$POSITIVE, Val_age_proj$study_age_sero_proj$SAMPLE_SIZE)) %>%
#   mutate(AGE_RANGE = paste0(formatC(AGE_LOWER, width = 2, flag = 0), "-", formatC(AGE_UPPER, width = 2, flag = 0)),
#          AGE_MID = (AGE_LOWER + AGE_UPPER)/2) %>% 
#   mutate(model_type = case_when(model_type == "Constant" ~ "Africa", model_type == "No Age" ~ "Constant")) %>%
#   mutate(study = rep((mapply(rep, source_names[13:16], c(4,5,5,10)) %>% unlist), 2))
# 
# sero_fit_plot = ggplot(Val_age_sero_proj, aes(x = AGE_MID))+
#   geom_pointrange(aes(ymin = CI.Lower, ymax = CI.Upper, y = CI.PointEst))+
#   geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = model_type), alpha = 0.7)+
#   geom_line(aes(y = X50., group = model_type), color = "black")+
#   scale_color_manual(values = c("#82ad58","#4e7c6f"))+
#   scale_fill_manual(values = c("#82ad58","#4e7c6f"))+
#   labs(x = "Age", y = "Positive proportion", color = "Age exposure", fill = "Age exposure")+
#   facet_wrap(.~study, scale = "free_y")+
#   theme_bw()+
#   theme(text = element_text(size = 18), 
#         axis.text.x = element_text(size = 12))
# sero_fit_plot
# 
# Val_age_case_proj = Val_age_proj$study_age_case_proj %>% 
#   mutate(age_m = (age_l + age_u)/2)
# Val_age_case_proj[,7:9] = round(Val_age_case_proj[,7:9])
# ggplot(Val_age_case_proj, aes(x = age_m, fill = model_type))+
#   geom_point(aes(y = case), color = "black")+
#   geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), alpha = 0.6)+
#   geom_line(aes(y = X50., color = model_type))+
#   scale_fill_manual(values = c("#82ad58","#4e7c6f"))+
#   scale_color_manual(values = c("#82ad58","#4e7c6f"))+
#   scale_y_log10()+
#   labs(x = "Age", y = "Cases", color = "Models", fill = "Models")+
#   facet_wrap(study~., scale = "free_y")+
#   theme_bw()+
#   theme(text = element_text(size = 20), 
#         axis.text = element_text(size = 14))


#ggsave("../plots/PAHO_age_case_proj_plot.png", PAHO_age_case_proj_plot, width = 9, height = 6)

#Combining with the original analysis:
# #WAIC:
# model_order = c("Same age exposure", "Alternative age exposure", "Constant age exposure")
# WAIC_same_place_age = readRDS("../output/par_sum_fit.rds")$WAIC_study_l
# WAIC_validate_full = rbind(WAIC_validate, data.frame(model_type = "Same age exposure", filter(WAIC_validate, model_type == "Constant")[,2:4], WAIC = WAIC_same_place_age[1:30])) %>%
#   mutate(model_type = case_when(model_type == "Alternative" ~ "Alternative age exposure",
#                                 model_type == "Constant" ~ "Constant age exposure", 
#                                 TRUE ~ model_type))
# 
# WAIC_validate_full$model_type = factor(WAIC_validate_full$model_type, levels = model_order)
# WAIC_validate_full$study = factor(WAIC_validate_full$study, levels = unique(WAIC_validate_full$study))
# WAIC_validate_full =  (WAIC_validate_full %>% group_by(study) %>% summarise(model_type = model_type, data_subset = data_subset, study = study, delta_WAIC = WAIC - min(WAIC)))
# ggplot(WAIC_validate_full, aes(y = study, x = log10(delta_WAIC), fill = model_type))+
#   geom_bar(stat = "identity", position = position_dodge(width = 1))+
#   facet_wrap(data_subset~., scale = "free")

Val_age_proj <-readRDS("../output/Validate_datafit.rds")
# #sero data fit:
# sero_data_fit_comb = bind_rows(Val_age_sero_proj %>% mutate(model_type = case_when(model_type == "Alternative" ~ "Alternative age exposure", model_type == "Constant" ~ "Constant age exposure")), 
#           sero_data_fit_plot %>% filter(data_subset == "Africa") %>% mutate(model_type = case_when(model_type == "Constant" ~ "Same age exposure"), study = STUDY_ID))
# sero_data_fit_comb$model_type  = factor(sero_data_fit_comb$model_type, levels = model_order)
# ggplot(sero_data_fit_comb, aes(x = AGE_MID))+
#   geom_pointrange(aes(ymin = CI.Lower, ymax = CI.Upper, y = CI.PointEst))+
#   geom_ribbon(aes(ymin = X2.5., ymax = X97.5., color = model_type, fill = model_type), alpha = 0.5)+
#   geom_line(aes(y = X50., group = model_type), color = "black")+
#   #scale_color_manual(values = model_colors)+
#   #scale_fill_manual(values = model_colors)+
#   labs(x = "Age", y = "Positive proportion", color = "Data subset", fill = "Data subset")+
#   facet_wrap(.~STUDY_ID, scale = "free_y")+
#   theme_bw()+
#   theme(text = element_text(size = 18), 
#         axis.text.x = element_text(size = 12))


#case data fit:
# case_data_fit_comb = bind_rows(Val_age_proj$study_age_case_proj %>% mutate(model_type = case_when(model_type == "Alternative" ~ "Alternative age exposure", model_type == "Constant" ~ "Constant age exposure")), 
#                                case_data_fit_data %>% filter(data_subset != "Complete data") %>% mutate(model_type = case_when(model_type == "Constant" ~ "Same age exposure"), 
#                                                                                                  study = study_id, age_l = age_min, age_u = age_max, age_m = AGE_MID))
case_data_fit_comb = bind_rows(Val_age_proj$study_age_case_proj %>% mutate(age_m = (age_l + age_u)/2))
case_data_fit_comb$study = factor(mapply(rep, source_names, each = rle(case_data_fit_comb$study)$lengths) %>% unlist, levels = source_names)

case_data_fit_comb$model_type  = factor(case_data_fit_comb$model_type, levels = model_order) 
case_data_fit_comb = case_data_fit_comb %>%
  mutate(model_type = case_when(model_type == "Alternative" ~ "Africa", TRUE ~ as.character(model_type)))

case_data_fit_comb_plot = ggplot(case_data_fit_comb, aes(x = age_m, fill = model_type))+
  geom_point(aes(y = case), color = "black")+
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), alpha = 0.6)+
  geom_line(aes(y = X50.), color = "black")+
  scale_fill_manual(values = model_colors)+
  #scale_color_manual(values = model_colors)+
  #scale_y_log10()+
  labs(x = "Age", y = "Cases", color = "Age exposure", fill = "Age exposure")+
  facet_wrap(study~., scale = "free_y")+
  theme_bw()+
  theme(text = element_text(size = 20), 
        axis.text = element_text(size = 14))
case_data_fit_comb_plot
ggsave("../plots/case_data_fit_comb_plot.png", case_data_fit_comb_plot, width = 15, height = 8)

# WAIC_table <- read.csv("../output/WAIC_validate_PAHO.csv")
# reshape2::dcast(WAIC_table, Country~moodel_type)
# ggplot(WAIC_table, aes(x = Country, y = WAIC, fill = moodel_type))+
#   geom_bar(stat = "identity", position = position_dodge())+
#   facet_wrap(Country~. , scales = "free" )
# 
# WAIC_prop = WAIC_table %>% group_by(Country) %>% summarise(model_prob = exp(-(WAIC - min(WAIC)))/sum(exp(-(WAIC - min(WAIC)))), model_type = moodel_type)
# ggplot(WAIC_prop, aes(x = Country, y = model_prob, fill = model_type))+
#   geom_bar(stat = "identity")
# 
# WAIC_table %>% group_by(moodel_type) %>% summarise(WAIC = sum(WAIC))
# 
# #Comparing with the FOI from our model:
# load("../../yf_ensemble-main/output/foi_wtd_raw.RData")
# load("../../yf_ensemble-main/output/foi_wtd_with_Continent.RData")
# select.foi.raw = foi.raw.adm1[,,3]
# select.foi.raw = select.foi.raw[(row.names(select.foi.raw) %>% substr(1,3)) %in% mod_coun_list,]
# select.foi.raw = data.frame(select.foi.raw, ISO = substr(row.names(select.foi.raw), 1, 3))
# select.foi.raw.sum = sapply(mod_coun_list, function(x) filter(select.foi.raw,ISO == x) %>% select(starts_with("X")) %>% unlist %>% log10 %>% quantile95cri)
# 
# select.foi.regress = foi.all.wtd[,,3]
# select.foi.regress = select.foi.regress[(row.names(select.foi.regress) %>% substr(1,3)) %in% mod_coun_list,]
# select.foi.regress = data.frame(select.foi.regress, ISO = substr(row.names(select.foi.regress), 1, 3))
# select.foi.regress.sum = sapply(mod_coun_list, function(x) filter(select.foi.regress,ISO == x) %>% select(starts_with("X")) %>% unlist %>% quantile95cri)
# 
# #plot:
# PAHO_vs_model_FOI = data.frame(Country = mod_coun_list, comp = rep(c("No regression", "With regression"), each = length(mod_coun_list)), t(PAHO_FOI), rbind(t(select.foi.raw.sum), t(select.foi.regress.sum)))
# ggplot(PAHO_vs_model_FOI, aes(x = X50., y = X50..1))+
#   geom_pointrange(aes(xmin = X2.5., xmax = X97.5.))+
#   geom_linerange(aes(ymin = X2.5..1, ymax = X97.5..1))
# cor.test(PAHO_FOI[2,], select.foi.raw.sum[2,], method = "spearman")
# cor.test(PAHO_FOI[2,], select.foi.regress.sum[2,], method = "spearman")
