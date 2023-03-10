setwd("/home/quan/Documents/NotreDame/Alex/Yellow fever/yf_age_exposure/code/")

library(BayesianTools)
library(dplyr)
library(sn)
library(stringr)
library(matrixStats)

#Load data:
sero_data = read.csv("../data/yf_sero_data_with_coverage_mod.csv")
Study_id_seq = rle(as.character(sero_data$STUDY_ID))
one_age_group_data = Study_id_seq$values[Study_id_seq$lengths == 1]
sero_data = filter(sero_data, !(STUDY_ID %in% one_age_group_data))
case.noti.LR.df = read.csv("../data/yf_case_data_by_age_group_with_coverage_LR.csv")
#PAHO data:
PAHO_age_case = read.csv('../data/Age_case_PAHO.csv')
year_range = 2000:2014
pop.cov = read.csv("/home/quan/Documents/NotreDame/Alex/Yellow fever/yf_ensemble-main/data/adm_1_pop_and_cov_VIMC.csv")
coun_list = unique(PAHO_age_case$country)
#susceptible pop:
pop.cov.sel = filter(pop.cov, YEAR %in% year_range, ISO %in% c(coun_list, "PRY"))
sus.pop = (select(pop.cov.sel, starts_with("POP"))*(1 - VE*select(pop.cov.sel, starts_with("COV")))) %>%
  aggregate(by = list(ISO = pop.cov.sel$ISO), sum)

#helper function:
inverse_logit <- function(x) exp(x)/(exp(x) + 1)
quantile95cri <- function(x) quantile(x, probs = c(0.025, 0.5, 0.975))
#Vectorizing case data:
agg_id = (unique(PAHO_age_case$age_u) - unique(PAHO_age_case$age_l) + 1) %>% mapply(rep, 1:length(unique(PAHO_age_case$age_u)), each = .) %>% as.vector()

#LL function:
age_dep_term_array = age_dep_term_array[,1:2]
case_gen_fun <- function(par){
  FOI = 10^par[1]; rho = inverse_logit(par[2])
  FOI_time_age = age_dep_term_array*FOI
  pos_age = (1 - exp(-FOI_time_age))*t(exp(-apply(rbind(0, head(FOI_time_age, -1)), 1, cumsum)))
  case_age = pos_age*sus.pop.coun
  case_age_group = aggregate(case_age, by = list(agg_id), sum)[,-1]
  case_age_group_rho = case_age_group*rho
  return(case_age_group_rho)
}
LL_function <- function(par){
  case_age_group_rho = case_gen_fun(par)
  LL = sapply(1:ncol(age_dep_term_array), function(x) sum(dpois(case.coun, case_age_group_rho[,x], log = T))) %>% mean
  return(LL)
}

VE = 0.975
source("util_function.R")

#model runs:
all_model_type = c("Constant", "One outbreak", "Two outbreaks", "No Age")
for(model_type in all_model_type){
  if(model_type == "No Age"){
    age_dep_term_array = array(1/100, dim = c(100, 1002))
  } else {
    #Parameter id:
    par_id(model_type)
    #Age dependent function:
    age_dep_func <- function(sn_par){
      age_dep = dsn(0:99, sn_par[1], exp(sn_par[2]), sn_par[3])  #1: postition (range 0:99) 2: scale (positive) 3:skewness (real number)
      return((age_dep/sum(age_dep)))
    }
    
    #Par credible interval; 
    par_sel = readRDS(paste0("../output/MCMC_output/MCMC_output_sero_case_noWHO_morethanoneagegroup_", model_type,"_estracted.Rdata"))
    n_samples = nrow(par_sel)
    
    Age_depend_par = par_sel[,Age_depend_id];
    #Age_depend_par[,2] <- 10^(Age_depend_par[,2])
    age_dep_term_array = array(dim = c(100, n_samples))
    for(i in 1:n_samples){
      age_dep_term_array[,i] = age_dep_func(unlist(Age_depend_par[i,]))
    }
  }
  
  #model run
  par_prior_dens = function(par){
    dnorm(par[1], -3, 1, log = T) + dnorm(par[2], 0, 1.5, log = T)
  }
  init_values = function(n = 1) c(rnorm(n, -3, 1), rnorm(n, 0, 1.5))
  par_prior = createPrior(density = par_prior_dens, sampler = init_values,
                          lower = c(-10, -10), upper = c(10, 10), best =  c(0, 0))
  bayesanSetup <- createBayesianSetup(likelihood = LL_function, prior = par_prior)
  iter = 2e4
  settings = list(iterations = iter)
  
  output = c()
  for(coun in coun_list){
    if(coun == "PAR"){
      sus.pop.coun = filter(sus.pop, ISO == "PRY") %>% select(starts_with("POP")) %>% as.numeric()
    } else {
      sus.pop.coun = filter(sus.pop, ISO == coun) %>% select(starts_with("POP")) %>% as.numeric()
    }
    
    case.coun = filter(PAHO_age_case, country == coun)$cases; #case.coun[case.coun < 1e-10] = 1e-10
    
    a <- runMCMC(bayesianSetup = bayesanSetup, sampler = "DREAMzs", settings = settings)
    output = c(output, list(a))
  }
  
  saveRDS(output, paste0("../output/MCMC_output/PAHO_", model_type,"_FOI_est.rds"))
}

#get par from the output:
iter = 2e4
n_samples = 1002
PAHO_age_case_proj_l = c()
WAIC_table_l = c()
PAHO_par_l = c()
WAIC_all_coun_l = c()
for(model_type in all_model_type){
  print(model_type)
  output <- readRDS(paste0("../output/MCMC_output/PAHO_", model_type,"_FOI_est.rds"))
  par_output = lapply(output, function(x){
    getSample(x, start = round(iter/3/2), thin = 5) %>% tail(n_samples)
  })
  mod_coun_list = coun_list; mod_coun_list[6] = "PRY"
  names(par_output) = mod_coun_list
  
  #Summary FOI, rho:
  PAHO_FOI = sapply(par_output, function(x) quantile95cri(x[,1]))
  PAHO_rho = sapply(par_output, function(x) x[,2] %>% inverse_logit %>% quantile95cri)
  PAHO_par_l = c(PAHO_par_l, list(data.frame(model_type = model_type, coun = coun_list, rbind(t(PAHO_FOI), t(PAHO_rho)), par = rep(c("~gamma", "~rho"), each = length(coun_list)))))
  
  #Project age case pattern and data fit:
  if(model_type == "No Age"){
    age_dep_term_array = array(1/100, dim = c(100, 1002))
  } else {
    #Parameter id:
    par_id(model_type)
    #Age dependent function:
    age_dep_func <- function(sn_par){
      age_dep = dsn(0:99, sn_par[1], exp(sn_par[2]), sn_par[3])  #1: postition (range 0:99) 2: scale (positive) 3:skewness (real number)
      return((age_dep/sum(age_dep)))
    }
    
    #Par credible interval; 
    par_sel = readRDS(paste0("../output/MCMC_output/MCMC_output_sero_case_noWHO_morethanoneagegroup_", model_type,"_estracted.Rdata"))
    n_samples = nrow(par_sel)
    
    Age_depend_par = par_sel[,Age_depend_id];
    #Age_depend_par[,2] <- 10^(Age_depend_par[,2])
    age_dep_term_array = array(dim = c(100, n_samples))
    for(i in 1:n_samples){
      age_dep_term_array[,i] = age_dep_func(unlist(Age_depend_par[i,]))
    }
  }
  
  project_age_case_l = c()
  WAIC_coun = c()
  LL_model_all_coun_l = c()
  for(i in 1:length(coun_list)){
    coun = coun_list[i]
    if(coun == "PAR"){
      sus.pop.coun = filter(sus.pop, ISO == "PRY") %>% select(starts_with("POP")) %>% as.numeric()
    } else {
      sus.pop.coun = filter(sus.pop, ISO == coun) %>% select(starts_with("POP")) %>% as.numeric()
    }
    
    x = par_output[[i]]
    proj = sapply(1:n_samples, function(y) rowMeans(case_gen_fun(x[y,])))
    project_age_case_l = c(project_age_case_l, list((t(apply(proj, 1, quantile95cri)))))
    
    #LL:
    case.coun = filter(PAHO_age_case, country == coun)$cases;# case.coun[case.coun < 1e-10] = 1e-10
    LL_samples = sapply(1:n_samples, function(x) dpois(case.coun, proj[,x], log = T))
    LL_model_all_coun_l = c(LL_model_all_coun_l, list(LL_samples))
    lppd <- sapply(1:length(case.coun), function(i) logSumExp(LL_samples[i,]) - log(n_samples)) %>% sum
    pWAIC <- sapply(1:length(case.coun) , function(i) var(LL_samples[i,])) %>% sum
    WAIC = -2*(lppd - pWAIC)
    WAIC_coun = c(WAIC_coun, WAIC)
  }
  project_age_case =  Reduce(rbind, project_age_case_l)
  LL_model_all_coun = Reduce(rbind, LL_model_all_coun_l)
  lppd_all_coun = sapply(1:nrow(PAHO_age_case), function(i) logSumExp(LL_model_all_coun[i,]) - log(n_samples)) %>% sum
  pWAIC_all_coun <- sapply(1:nrow(PAHO_age_case) , function(i) var(LL_model_all_coun[i,])) %>% sum
  WAIC_all_coun_l = c(WAIC_all_coun_l, -2*(lppd_all_coun - pWAIC_all_coun))
  
  PAHO_age_case_proj_l = c(PAHO_age_case_proj_l, list(data.frame(cbind(filter(PAHO_age_case, country %in% coun_list), project_age_case), model_type = model_type) %>%
    mutate(age_m = (age_u + age_l)/2)))
  WAIC_table_l = c(WAIC_table_l, list(data.frame(Country = coun_list, WAIC = WAIC_coun, moodel_type = model_type)))
}
PAHO_par = Reduce(rbind, PAHO_par_l)
PAHO_age_case_proj = Reduce(rbind, PAHO_age_case_proj_l)
WAIC_table_l = c(WAIC_table_l, list(data.frame(Country = "Global", WAIC = WAIC_all_coun_l, moodel_type = all_model_type)))
WAIC_table = Reduce(rbind, WAIC_table_l)

write.csv(PAHO_par, "../output/PAHO_par.csv", row.names = F)
write.csv(PAHO_age_case_proj, "../output/PAHO_age_case_proj.csv", row.names = F)
write.csv(WAIC_table, "../output/WAIC_validate_PAHO.csv", row.names = F)

#Plotting
library(ggplot2)
model_order = c("No Age", "One outbreak", "Two outbreaks", "Constant")
model_colors = c("#f36201","#82ad58","#4e7c6f", "#fbbb09")
#Summary parameters:
PAHO_par <-read.csv("../output/PAHO_par.csv")
PAHO_par$model_type = factor(PAHO_par$model_type, levels = model_order)
PAHO_par_plot <- ggplot(PAHO_par, aes(x = X50., y = coun, color = model_type))+
  geom_pointrange(aes(xmin = X2.5., xmax = X97.5.), position = position_dodge(width = 0.7))+
  facet_wrap(~par, ncol = 2, nrow = 1, scale = "free_x", labeller = labeller(par = label_parsed))+
  scale_color_manual(values = model_colors)+
  labs(x = "", y = "Country", color = "Models")+
  theme_bw()+
  theme(text = element_text(size = 19), 
        axis.text.y = element_text(size = 14), 
        axis.text.x = element_text(size = 10),
        strip.text = element_text(size = 12))
ggsave("../plots/PAHO_par_plot.png", PAHO_par_plot, width = 8, height = 7)

PAHO_age_case_proj <-read.csv("../output/PAHO_age_case_proj.csv")
PAHO_age_case_proj$model_type = factor(PAHO_age_case_proj$model_type, levels = model_order)
PAHO_age_case_proj_plot = ggplot(PAHO_age_case_proj, aes(x = age_m, fill = model_type))+
  geom_point(aes(y = cases), color = "black")+
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), alpha = 0.6)+
  geom_line(aes(y = X50.), color = "black")+
  scale_fill_manual(values = model_colors)+
  labs(x = "Age", y = "Cases", fill = "Models")+
  facet_wrap(country~., scale = "free_y")+
  theme_bw()+
  theme(text = element_text(size = 20), 
        axis.text = element_text(size = 14))
ggsave("../plots/PAHO_age_case_proj_plot.png", PAHO_age_case_proj_plot, width = 9, height = 6)

WAIC_table <- read.csv("../output/WAIC_validate_PAHO.csv")
reshape2::dcast(WAIC_table, Country~moodel_type)
ggplot(WAIC_table, aes(x = Country, y = WAIC, fill = moodel_type))+
  geom_bar(stat = "identity", position = position_dodge())+
  facet_wrap(Country~. , scales = "free" )

WAIC_prop = WAIC_table %>% group_by(Country) %>% summarise(model_prob = exp(-(WAIC - min(WAIC)))/sum(exp(-(WAIC - min(WAIC)))), model_type = moodel_type)
ggplot(WAIC_prop, aes(x = Country, y = model_prob, fill = model_type))+
  geom_bar(stat = "identity")

WAIC_table %>% group_by(moodel_type) %>% summarise(WAIC = sum(WAIC))

#Comparing with the FOI from our model:
load("../../yf_ensemble-main/output/foi_wtd_raw.RData")
load("../../yf_ensemble-main/output/foi_wtd_with_Continent.RData")
select.foi.raw = foi.raw.adm1[,,3]
select.foi.raw = select.foi.raw[(row.names(select.foi.raw) %>% substr(1,3)) %in% mod_coun_list,]
select.foi.raw = data.frame(select.foi.raw, ISO = substr(row.names(select.foi.raw), 1, 3))
select.foi.raw.sum = sapply(mod_coun_list, function(x) filter(select.foi.raw,ISO == x) %>% select(starts_with("X")) %>% unlist %>% log10 %>% quantile95cri)

select.foi.regress = foi.all.wtd[,,3]
select.foi.regress = select.foi.regress[(row.names(select.foi.regress) %>% substr(1,3)) %in% mod_coun_list,]
select.foi.regress = data.frame(select.foi.regress, ISO = substr(row.names(select.foi.regress), 1, 3))
select.foi.regress.sum = sapply(mod_coun_list, function(x) filter(select.foi.regress,ISO == x) %>% select(starts_with("X")) %>% unlist %>% quantile95cri)

#plot:
PAHO_vs_model_FOI = data.frame(Country = mod_coun_list, comp = rep(c("No regression", "With regression"), each = length(mod_coun_list)), t(PAHO_FOI), rbind(t(select.foi.raw.sum), t(select.foi.regress.sum)))
ggplot(PAHO_vs_model_FOI, aes(x = X50., y = X50..1))+
  geom_pointrange(aes(xmin = X2.5., xmax = X97.5.))+
  geom_linerange(aes(ymin = X2.5..1, ymax = X97.5..1))
cor.test(PAHO_FOI[2,], select.foi.raw.sum[2,], method = "spearman")
cor.test(PAHO_FOI[2,], select.foi.regress.sum[2,], method = "spearman")
