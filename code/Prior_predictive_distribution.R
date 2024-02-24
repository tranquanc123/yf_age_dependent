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

#
n_cluster = 1e1
n_pop = 1e0
n_coun = 1e4
rm(age_dep_data_coun_sum, age_dep_data_feature, age_dep_term_pop);gc()
age_dep_term_pop = array(dim = c(100, n_cluster, n_pop, n_coun)) #age: 100; number of clusters: 10; 
for(ii in 1:n_pop){
  age_dep1_pop_mean = rnorm(n_cluster, 50, 125)
  age_dep2_pop_mean = rhnorm(n_cluster, 150)
  age_dep3_pop_mean = rt(n_cluster,2,1/2)
  age_dep1_pop_sd = rexp(n_cluster, 1/200)
  age_dep2_pop_sd = rexp(n_cluster, 1/200)
  age_dep3_pop_sd = rexp(n_cluster, 1/10)
  age_dep1_coun = sapply(1:n_cluster, function(x) rnorm(n_coun, age_dep1_pop_mean[x], age_dep1_pop_sd[x]))
  age_dep2_coun = sapply(1:n_cluster, function(x) rnorm(n_coun, age_dep2_pop_mean[x], age_dep2_pop_sd[x]))
  age_dep3_coun = sapply(1:n_cluster, function(x) rnorm(n_coun, age_dep3_pop_mean[x], age_dep3_pop_sd[x]))
 
  #(aov(value ~ cluster, data.frame(cluster = rep(1:10, each = 1e2), value = as.vector(age_dep1_coun))) %>% summary)[[1]]$`Pr(>F)`[1]
  for(i in 1:length(age_dep1_coun)){
    sel_id_i = ((i - 1) %/% n_coun) + 1
    sel_id_ii = i - (sel_id_i - 1)*n_coun
    age_dep_term_pop[, sel_id_i, ii, sel_id_ii] = dsn(0:99, age_dep1_coun[i], age_dep2_coun[i], age_dep3_coun[i], log = T)
  }
}

age_dep_data = cbind(expand.grid(age = 0:99, cluster = 1:n_cluster, pop = 1:n_pop, coun = 1:n_coun), value = as.vector(age_dep_term_pop))
age_dep_data_feature = age_dep_data %>% group_by(cluster, pop, coun) %>% summarise(peak_age = age[which.max(value)],
                                                          diff_high_low = log(diff(range(value)),10),
                                                          diff_mod_mean = (which.max(value) - 1) - sum((0:99)*exp(value - logSumExp(value))))
rm(age_dep_data);gc()
filter(age_dep_data_feature, peak_age != 0, peak_age != 99)

age_dep_data_coun_sum = age_dep_data_feature %>% group_by(cluster, pop) %>% summarise(sd_peak_age = sd(peak_age), sd_diff_high_low = sd(diff_high_low), sd_diff_mod_mean = sd(diff_mod_mean))
age_dep_data_pop_sum = age_dep_data_feature %>% group_by(pop) %>% summarise(sd_peak_age = sd(peak_age), sd_diff_high_low = sd(diff_high_low, na.rm = T), sd_diff_mod_mean = sd(diff_mod_mean, na.rm = T))

hist(filter(age_dep_data_coun_sum, sd_peak_age != 0)$sd_peak_age,50,col='gray',
     xlab='sd of peak age within a country',main='')
hist(age_dep_data_coun_sum$sd_diff_high_low,50,col='gray',
     xlab='sd of fold difference between exposure of high and low age on log10 scale within a country',main='')
hist(age_dep_data_coun_sum$sd_diff_mod_mean,50,col='gray',
     xlab='sd of difference between mode and mean within a country',main='')

hist(age_dep_data_pop_sum$sd_peak_age,50,col='gray',
     xlab='sd of peak age of population draw',main='')
hist(age_dep_data_pop_sum$sd_diff_high_low,50,col='gray',
     xlab='sd of fold difference between exposure of high and low age on log10 scale of population draw',main='')
hist(age_dep_data_pop_sum$sd_diff_mod_mean,50,col='gray',
     xlab='sd of difference between mode and mean of population drawr',main='')

par(mfrow=c(3,4))
for(plot_i in 1:10){
  matplot(exp(age_dep_term_pop[, plot_i,2, ])^0.2, type = "l", lty = 1)
}

par(mfrow=c(3,4))
for(plot_i in 1:10){
  # fold difference between exposure of highest and lowest age
  hist(filter(age_dep_data_feature, cluster == plot_i, pop == 1)$diff_high_low,50,col='gray',
       xlab='fold difference between exposure of high and low age on log10 scale',main='')
}
#
par(mfrow=c(3,4))
for(plot_i in 1:10){
  # peak age
  hist(filter(age_dep_data_feature, cluster == plot_i, pop == 1, peak_age != 0, peak_age != 99)$peak_age,50,col='gray',xlab='peak age',main='') #Removing peak age at 1 and 100
}
#
par(mfrow=c(3,4))
for(plot_i in 1:10){
  # difference between mode and mean
  hist(filter(age_dep_data_feature, cluster == plot_i, pop == 1)$diff_mod_mean,50,col='gray',
       xlab='difference between mode and mean',main='')
}
