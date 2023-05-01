#Prior for vac mis model:
library(extraDistr)
library(dplyr)

setwd("/home/quan/Documents/NotreDame/Alex/Yellow fever/yf_age_exposure/code/")

#Reporting proportion:
rep_prop = c()
for(which.scenario in 1:8){
  load(paste0("/home/quan/Documents/NotreDame/Alex/Yellow fever/yf_ensemble-main/output/in_paper/proportion_by_type_", which.scenario,".RData"))
  rep_prop = c(rep_prop, 1 - rdirichlet(1e3, prop.FCU)[,3])
}

inverse_logit <- function(x) exp(x)/(exp(x) + 1)
logit <- function(x) log(x) - log(1-x)
cost_func <- function(par){
  -sum(dnorm(logit(rep_prop), mean = par[1], sd = par[2], log = T))
}

# cost_func <- function(par){
#   -sum(dbeta(rep_prop, par[1], par[2], log = T))
# }

output = optim(par = c(0.1, 0.1), cost_func, method = "BFGS")
output
hist(rep_prop)
hist(inverse_logit(rnorm(8e3,output$par[1], output$par[2] )), breaks = 1e3, xlim = c(0, 0.004))
#hist(rbeta(8e3,output$par[1], output$par[2]), breaks = 1e3, xlim = c(0, 0.004))
reporting_prop = output$par

#FOI
case.df <- read.csv("/home/quan/Documents/NotreDame/Alex/Yellow fever/yf_age_exposure/data/yf_case_data_by_age_group_with_coverage_LR.csv")
sero.df <- read.csv("/home/quan/Documents/NotreDame/Alex/Yellow fever/yf_age_exposure/data/yf_sero_data_with_coverage_mod.csv")
coun_list = unique(c(case.df$country, sero.df$ISO))

load(paste0("/home/quan/Documents/NotreDame/Alex/Yellow fever/yf_ensemble-main/output/in_paper/foi_wtd.RData"))
load(paste0("/home/quan/Documents/NotreDame/Alex/Yellow fever/yf_ensemble-main/output/in_paper/model_foi_1_11.RData"))

cost_func <- function(par){
  -sum(dnorm(sel_FOI, mean = par[1], sd = par[2], log = T))
}

par_all = c()
for(coun in coun_list){
  print(coun)
  sel_FOI = filter(c, ISO == coun) %>% select(starts_with("foi")) %>% unlist %>% as.numeric()
  if(length(sel_FOI) == 0){
    sel_FOI = c %>% select(starts_with("foi")) %>% unlist
  }
  output = optim(par = c(0.1, 0.1), cost_func, method = "BFGS")
  par_all = c(par_all, list(output$par))
}
FOI_data = data.frame(ISO = coun_list, Reduce(rbind, par_all))

saveRDS(list(reporting_prop_prior = reporting_prop, FOI_prior = FOI_data), "../data/FOI_rho_prior.rds")
