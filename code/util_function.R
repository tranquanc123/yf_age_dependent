#helper function ####
inverse_logit <- function(x) exp(x)/(exp(x) + 1)
quantile95cri <- function(x) quantile(x, probs = c(0.025, 0.5, 0.975))

#Speed up calculation of likelihood of sero and case data ####
#Sero 
sero_study_rle = rle(as.character(sero_data$STUDY_ID))
N_sero_study = length(sero_study_rle$lengths)
nrow_sero_study = mapply(rep, 1:N_sero_study, each = sero_study_rle$lengths) %>% unlist()
sero_study_id = sero_study_rle$values
#Case LR
case_LR_study_rle = rle(as.character(case.noti.LR.df$study_id))
N_case_LR_study = length(case_LR_study_rle$lengths)
nrow_case_LR_study = mapply(rep, 1:N_case_LR_study, each = case_LR_study_rle$lengths) %>% unlist()
case_LR_study_id = case_LR_study_rle$values

N_study = N_sero_study + N_case_LR_study
names_study = c(sero_study_id, case_LR_study_id)

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
n_row_sero_l = unlist(mapply(rep, 1:length(n_row_sero), n_row_sero))

#Speed up calculation for case LR data:
year_case_LR_max = aggregate(case.noti.LR.df$year, list(nrow_case_LR_study), max)$x[nrow_case_LR_study]
FOI_case_LR_id_row = year_case_LR_max - case.noti.LR.df$year + case.noti.LR.df$age + 1
FOI_case_LR_id_row_max = max(FOI_case_LR_id_row)
study_agegroup = interaction(case.noti.LR.df$age_group_ID, case.noti.LR.df$study_id) %>% as.vector %>% factor(level = rle(.)$values)
case_by_study_agegroup = aggregate(case.noti.LR.df$case, list(study_agegroup), unique)$x
length_case_LR_study = rle(as.character(nrow_case_LR_study))$lengths

n_row_study_case_LR = rle(as.character(study_agegroup))$values %>% str_split("[.]") %>% sapply(function(x) x[2]) %>% rle
n_row_study = c(mapply(rep, 1:N_sero_study, each = sero_study_rle$lengths) %>% unlist, 
                mapply(rep, N_sero_study + 1:N_case_LR_study, each = n_row_study_case_LR$lengths) %>% unlist)
# mapply(rep, 1:N_study, each = c(sero_study_rle$lengths, levels(study_WHO_agegroup), n_row_study_case_LR$lengths)) %>% unlist()
#retreive back the original case LR data:
case.noti.LR.or.df = unique(case.noti.LR.df[,c("age_group_ID", "country", "study_id")])
case.noti.LR.or.df$age_min = aggregate(case.noti.LR.df$age, by = list(case.noti.LR.df$age_group_ID), min)$x
case.noti.LR.or.df$age_max = aggregate(case.noti.LR.df$age, by = list(case.noti.LR.df$age_group_ID), max)$x
case.noti.LR.or.df$case = case_by_study_agegroup

#Age dependent function####
age_dep_func <- function(sn_par){
  age_dep = dsn(0:99, sn_par[1], sn_par[2], sn_par[3])  #1: postition (range 0:99) 2: scale (positive) 3:skewness (real number)
  return((age_dep/sum(age_dep)))
}

#parameter id ####
par_id <- function(model_type){
  par_names <<- c()
  if(model_type == "Constant"){
    FOI_par_id <<- 1:N_study
    #FOI_avg_id <<- N_study + 1; FOI_sd_id <<- N_study + 2
    Age_depend_id <<- max(FOI_par_id) + 1:3
    par_names <<- c(par_names, paste0("FOI_", 1:N_study), paste0("Age_dep_", 1:3))
  } else if(model_type == "One outbreak"){
    FOI_par_id <<- 1:(N_study*2)
    alpha_id <<- 1:N_study;
    time_id <<- max(alpha_id) + 1:N_study
    #alpha_avg_id <<- max(time_id) + 1; alpha_sd_id <<- alpha_avg_id + 1
    Age_depend_id <<- max(time_id) + 1:3
    par_names <<- c(par_names,  paste0("FOI_alpha_", 1:N_study),  paste0("FOI_time_", 1:N_study), 
                   paste0("Age_dep_", 1:3))
  } else {
    FOI_par_id <<- 1:(N_study*4)
    alpha1_id <<- 1:N_study;
    time1_id <<- max(alpha1_id) + 1:N_study
    alpha2_id <<- max(time1_id) + 1:N_study;
    time2_id <<- max(alpha2_id) + 1:N_study
    #alpha_avg_id <<- max(time2_id) + 1; alpha_sd_id <<- alpha_avg_id + 1
    Age_depend_id <<- max(time2_id) + 1:3
    par_names <<- c(par_names,  paste0("FOI_alpha1_", 1:N_study),  paste0("FOI_time1_", 1:N_study), paste0("FOI_alpha2_", 1:N_study),  paste0("FOI_time2_", 1:N_study), 
                   paste0("Age_dep_", 1:3))
  }
  
  #VE:
  VE_id <<- max(Age_depend_id) + 1
  par_names <<- c(par_names, "VE")
  
  #reporting proportion
  rho_case_LR_id <<- VE_id + 1:N_case_LR_study
  par_names <<- c(par_names, paste0("rho_case_LR_", 1:N_case_LR_study))
}

#FOI function####
FOI_function_model_type <- function(model_type){
  if(model_type == "Constant"){
    FOI_func <<- function(x) rep(x, each = 100)
  } else if(model_type == "One outbreak"){
    FOI_func <<- function(x){
      alpha = x[alpha_id]; #beta = x[beta_id]; 
      t = x[time_id]
      #if(bg == "N"){
      return(lapply(1:N_study, function(i) {
        FOI_100 = numeric(100)
        whole_time = floor(t[i]); part_time = t[i] - whole_time
        FOI_100[whole_time + 1] <- (1 - part_time)*alpha[i]
        FOI_100[whole_time + 2] <- part_time*alpha[i]
        return(FOI_100)
      }) %>% unlist) 
    }
  } else {
    FOI_func <<- function(x){
      alpha1 = x[alpha1_id]; alpha2 = x[alpha2_id]
      #beta1 = x[beta1_id]; beta2 = x[beta2_id]; 
      t1 = x[time1_id]; t2 = x[time2_id]
      t_m = cbind(t1, t2) %>% apply(1, sort) %>% t
      #if(bg == "N"){
      return(lapply(1:N_study, function(i) {
        FOI_100 = numeric(100)
        whole_time_v = floor(t_m[i, ]); part_time_v = t_m[i, ] - whole_time_v
        FOI_100[whole_time_v + 1] <- (1 - part_time_v)*c(alpha1[i], alpha2[i])
        FOI_100[whole_time_v + 2] <- part_time_v*c(alpha1[i], alpha2[i])
        return(FOI_100)
      }) %>% unlist)
    }
  }
}

#function to generate cases and serological data:####
gen_func <- function(par){
  FOI_par = unlist(par[FOI_par_id]); 
  if(model_type == "One outbreak"){
    FOI_par[alpha_id] <- 10^(FOI_par[alpha_id])
  } else if(model_type == "Two outbreaks") {
    FOI_par[c(alpha1_id, alpha2_id)] <- 10^(FOI_par[c(alpha1_id, alpha2_id)])
  } else if(model_type == "Constant"){
    FOI_par[FOI_par_id] <- 10^(FOI_par[FOI_par_id])
  }
  
  FOI_m = FOI_func(FOI_par) %>% matrix(nrow = 100, ncol = N_study)
  
  Age_depend_par = par[Age_depend_id];
  #Age_depend_par[2] <- 10^(Age_depend_par[2])
  age_dep_term = age_dep_func(Age_depend_par)
  
  #Include the effect of bg infection and age dependent:
  FOI_m = ((FOI_m)*age_dep_term) %>% matrix(nrow = 100, ncol = N_study)
  
  VE = par[VE_id]
  
  rho_case_LR = inverse_logit(unlist(par[rho_case_LR_id]))
  
  FOI_sero_m = FOI_m[, 1:N_sero_study]
  FOI_sero_cumsum = sapply(1:length(age_range_v), function(x) sum(FOI_sero_m[1:(age_range_v[x] + 1), nrow_sero_study_age[x]]))
  seropos = 1-exp(-FOI_sero_cumsum)
  sus_v = (1-cov_v*VE) * pop_v
  p = aggregate(sus_v*seropos, list(n_row_sero_l), sum)$x/aggregate(sus_v, list(n_row_sero_l), sum)$x
  p[p > 1 - 1e-16] = 1 - 1e-16
  p[p < 1e-16] = 1e-16
  
  ###Case LR data:
  #case data from LR:
  FOI_case_LR_m = FOI_m[, N_sero_study + 1:N_case_LR_study]
  FOI_case_LR_m = rbind(FOI_case_LR_m, matrix(FOI_case_LR_m[100,], nrow = FOI_case_LR_id_row_max - 100, ncol = N_case_LR_study, byrow = T))
  FOI_case_LR_v = sapply(1:length(nrow_case_LR_study), function(i) FOI_case_LR_m[FOI_case_LR_id_row[i], nrow_case_LR_study[i]])
  FOI_case_LR_int_v = sapply(1:length(nrow_case_LR_study), function(i) {
    FOI_case_LR_id_int = FOI_case_LR_id_row[i] - 1
    if(FOI_case_LR_id_int == 0) return(1)
    return(sum(FOI_case_LR_m[1:FOI_case_LR_id_int, nrow_case_LR_study[i]]))
  })
  l_rho_case_LR = mapply(rep, rho_case_LR, each = length_case_LR_study) %>% unlist #rep(rho_case_LR, each = N_case_LR_study*100)
  case_LR_sus = case.noti.LR.df$pop*(1 - case.noti.LR.df$cov*VE)
  exp_cases_LR = case_LR_sus*exp(-FOI_case_LR_int_v)*(1 - exp(-FOI_case_LR_v))*l_rho_case_LR
  
  exp_case_LR_by_age_group = aggregate(exp_cases_LR, by = list(study_agegroup), FUN = sum)$x
  exp_case_LR_by_age_group[exp_case_LR_by_age_group < 1e-10] = 1e-10
  
  return(list(sero = p, case_LR = exp_case_LR_by_age_group))
}
