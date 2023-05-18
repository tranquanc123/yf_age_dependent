library(sn)
library(extraDistr)
library(matrixStats)
library(ggplot2)
library(dplyr)

####Prior predictive distribution for skew-normal distribution:
n = 1e5
sn_par = data.frame(
  xi = runif(n, -1, 97),
  omega = rhnorm(n, 13),
  alpha = rt(n,2,1/2))
sn_par = data.frame(
  xi = rnorm(n, 50, 125),
  omega = rhnorm(n, 150),
  alpha = rt(n,2,1/2))
sn_gen = sapply(1:nrow(sn_par), function(x) {
  age_dep = dsn(0:99, sn_par$xi[x], sn_par$omega[x], sn_par$alpha[x], log = T)
  return(exp(age_dep))
  })

matplot(sn_gen[, sample(1:n, 1e3)]^0.2, type = "l", lty = 1)

# fold difference between exposure of highest and lowest age
scale_diff = apply(sn_gen,2,function(x){log(diff(log(range(x))),10)})
hist(scale_diff,50,col='gray',
     xlab='fold difference between exposure of high and low age on log10 scale',main='')

# peak age
max_peak = apply(sn_gen,2,which.max)
hist(max_peak[scale_diff > 0],50,col='gray',xlab='peak age',main='')
hist(max_peak[max_peak != 1 & max_peak != 100],50,col='gray',xlab='peak age',main='') #Removing peak age at 1 and 100

# difference between mode and mean
hist(apply(sn_gen,2,function(x){which.max(x)-sum((0:99)*x/sum(x))}),50,col='gray',
     xlab='difference between mode and mean',main='')

####Prior for gamma
#Coun from the studies
setwd("/home/quan/Documents/NotreDame/Alex/Yellow fever/yf_age_exposure/code/")
sero_data = read.csv("../data/yf_sero_data_with_coverage_mod.csv")
case.noti.LR.df = read.csv("../data/yf_case_data_by_age_group_with_coverage_LR.csv")
coun_list = unique(c(sero_data$ISO, case.noti.LR.df$country))
#FOI estimates from SA
load("/home/quan/Documents/NotreDame/Alex/Yellow fever/yf_burden_mod/output/foi_wtd.RData")
load("/home/quan/Documents/NotreDame/Alex/Yellow fever/yf_ensemble-main/output/in_paper/model_foi_1_11.RData")

#Constant g function
g_func <<- function(x) rep(x, each = 100)

#match cumulative sum:
n_cost_func = 1e4
cost_func <- function(par){
  g_par_mat = data.frame(rlnorm(n_cost_func, par[1], par[2]))
  #Cummulative sum from the age dependent model:
  lambda_vec = sapply(1:n_cost_func,function(ii) cumsum(g_par_mat[ii,]*(sn_gen[,ii])))
  return(sum((lambda_vec - cumsum_coun_FOI_sampled[,1:n_cost_func])^2))
}

output_l = c()
for(coun in coun_list){
  print(coun)
  if(coun == "BRA"){
    coun_FOI = foi.all.wtd
  } else {
    coun_FOI = foi.all.wtd[c$ISO == coun,,] 
  }
  #Cummulative sum from the country
  cumsum_coun_FOI_sampled = sapply(coun_FOI[sample(1:prod(dim(coun_FOI)), n_cost_func)], function(x) (0:99)*10^x)
  output <- optim(par = c(1, 1), cost_func)
  output_l = c(output_l, list(output))
}

#
gamma_priors = data.frame(ISO = coun_list, Reduce(rbind, lapply(output_l, function(x) x$par)))

saveRDS(gamma_priors, "../data/gamma_lnorm_priors_constant_nonorm.rds")
#Comparing the FOI:
FOI_vary = data.frame(ISO = rep(coun_list, each = 100), lapply(1:nrow(gamma_priors), function(x) (rlnorm(n_cost_func, gamma_priors$X1[x], gamma_priors$X2[x])*sn_gen[,1:n_cost_func]) %>%
                                            apply(1, quantile, probs = c(0.025, 0.5, 0.975)) %>% t) %>% Reduce(rbind, .),
           FOI_used = "FOI estimated from optimized gamma", age = 0:99)
FOI_SA = data.frame(ISO = rep(coun_list, each = 100), lapply(coun_list, function(x) {
  if(coun == "BRA"){
    coun_FOI = foi.all.wtd
  } else {
    coun_FOI = foi.all.wtd[c$ISO == coun,,] 
  }
  coun_FOI_sum = (10^coun_FOI[sample(1:prod(dim(coun_FOI)), n_cost_func)]) %>% quantile(probs = c(0.025, 0.5, 0.975))
  return(coun_FOI_sum)
}) %>% Reduce(rbind, .), FOI_used = "constant FOI from SAdvance paper", age = 0:99)

comp_FOI_data = rbind(FOI_vary, FOI_SA)

ggplot(comp_FOI_data, aes(x = age, y = X50., color = FOI_used, fill = FOI_used))+
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), alpha = 0.4)+ 
  labs(y = "FOI")+
  facet_wrap(.~ISO)+
  geom_line()

# #Predicted positive proportion:
# pos_FOI_vary = apply((1 - exp(-apply(FOI_vary, 2, cumsum))), 1, quantile, probs = c(0.025, 0.5, 0.975))
# pos_FOI_SA = apply((1 - exp(-cumsum_AGO_FOI_sampled[,1:n_cost_func])) , 1, quantile, probs = c(0.025, 0.5, 0.975))
# comp_pos_FOI_data = data.frame(FOI_used = rep(c("FOI estimated from optimized gamma", "constant FOI from SAdvance paper"), each = ncol(pos_FOI_SA)), 
#                            rbind(t(pos_FOI_vary), t(pos_FOI_SA)), age = 0:99)
# 
# ggplot(comp_pos_FOI_data, aes(x = age, y = X50., color = FOI_used, fill = FOI_used))+
#   geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), alpha = 0.4)+ 
#   geom_line()
#                            
