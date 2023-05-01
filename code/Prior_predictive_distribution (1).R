library(sn)
library(extraDistr)

####Prior predictive distribution for skew-normal distribution:
sn_par = expand.grid(xi = rnorm(1e2, 0, 1000), omega = rhmorm(1e2, 1000), alpha = rt(1e2,  2, 1/2))
n = 1e3
sn_par = data.frame(
  xi = rnorm(n, 0, 1000),
  omega = rhnorm(n, 1000),
  alpha = rt(n, 2, 1/2))
sn_gen = sapply(1:nrow(sn_par), function(x) dsn(0:99, sn_par$xi[x], sn_par$omega[x], sn_par$alpha[x]))
matplot(sn_gen^0.2, type = "l", lty = 1)
# peak age
hist(apply(sn_gen,2,which.max),50,col='gray',xlab='peak age',main='')
# fold difference between exposure of highest and lowest age on log10 scale
hist(apply(sn_gen,2,function(x){log(diff(log(range(x))),10)}),50,col='gray',
     xlab='fold difference between exposure of high and low age on log10 scale',main='')
# difference between mode and mean
hist(apply(sn_gen,2,function(x){which.max(x)-sum((0:99)*x/sum(x))}),50,col='gray',
     xlab='difference between mode and mean',main='')


#FOI estimates, Population and coverage data:
output <-readRDS("AGO_foi_pop_cov.rds") 
AGO_FOI = output$AGO_FOI #FOI estimates in AGO from science advance paper
AGO_pop = output$AGO_pop #Population data in AGO from POLICI
AGO_cov = output$AGO_cov #Coverage data in AGO from POLICI
AGO_cov[is.na(AGO_cov)] = 0

#Time dependent g function:
g_function_model_type <- function(model_type){
  if(model_type == "Constant"){
    g_func <<- function(x) rep(x, each = 100)
  } else if(model_type == "One outbreak"){
    g_func <<- function(x){
      gamma = x[1];
      t = x[2]
      g_100 = numeric(100)
      whole_time = floor(t); part_time = t - whole_time
      g_100[whole_time + 1] <- (1 - part_time)*gamma
      g_100[whole_time + 2] <- part_time*gamma
      return(g_100)
    }
  } else if(model_type == "Two outbreaks"){
    g_func <<- function(x){
      gamma1 = x[1]; gamma2 = x[2]
      t1 = x[3]; t2 = x[4]
      t_m = sort(c(t1, t2))
      g_100 = numeric(100)
      whole_time_v = floor(t_m); part_time_v = t_m - whole_time_v
      g_100[whole_time_v + 1] <- (1 - part_time_v)*c(gamma1, gamma2)
      g_100[whole_time_v + 2] <- part_time_v*c(gamma1, gamma2)
      return(g_100)
    }
  }
}
#Type in the model type in the argument so it will give you the time depedent g function:
g_function_model_type("Two outbreaks")
g_func

#Age dependent function
age_dep_func <- function(par){
  age_dep = dsn(0:99, par[1], par[2], par[3])  #1: postition  (real number) 2: scale (positive) 3:skewness (real number)
  return((age_dep/sum(age_dep)))
}

#Vaccine efficacy:
VE = 0.975

#Time range of data that was used to estimate FOI from the science advance paper:
sel_year_range = 1980:2014

# example parameters for g
g_par = c(1e-3,1e-2,20,50)

# matrix of parameters for g from possible prior
g_par_mat = cbind(
  exp(rnorm(100,12,2)),
  exp(rnorm(100,12,2)),
  sample(1:100,100,replace=T),
  sample(1:100,100,replace=T))

# function to calculate mean FOI as a function of f and g parameters
lambda_mean = function(g_par,sn_par_rep){
  mean(
    matrix(rep(age_dep_func(as.numeric(sn_par[sn_par_rep,])),times=100),100,100) * # f(a)
    matrix(rep(g_func(g_par),each=100),100,100) * # g(y)
    as.matrix(AGO_pop[1:100,1:100]) * # population
    as.matrix(1 - AGO_cov[1:100,1:100] * VE) / # susceptibility among vaccinated
    sum( # normalize by susceptible population
      as.matrix(AGO_pop[1:100,1:100]) *
      as.matrix(1 - AGO_cov[1:100,1:100] * VE)))
}

# calculate mean FOI for 100 prior samples
lambda_vec = sapply(1:100,function(ii)lambda_mean(g_par_mat[ii,],ii))
hist(log(lambda_vec),10)
