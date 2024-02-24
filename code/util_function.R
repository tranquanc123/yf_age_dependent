#helper function ####
inverse_logit <- function(x) exp(x)/(exp(x) + 1)
quantile95cri <- function(x) quantile(x, probs = c(0.025, 0.5, 0.975))

#different subset of data:
if(data_subset == "SA"){
  case.noti.df = rbind(filter(case.noti.LR.df, country == "BRA"), case.noti.PAHO.df)
  sero_data = sero_data[NULL,]
} else if(data_subset == "AF"){
  case.noti.df = filter(case.noti.LR.df, country != "BRA")
} else {
  case.noti.df = rbind(case.noti.LR.df, case.noti.PAHO.df)
}

#Speed up calculation of likelihood of sero and case data ####
#Sero 
if(data_subset != "SA"){
  sero_study_rle = rle(as.character(sero_data$STUDY_ID))
  N_sero_study = length(sero_study_rle$lengths)
  nrow_sero_study = mapply(rep, 1:N_sero_study, each = sero_study_rle$lengths) %>% unlist()
  sero_study_id = sero_study_rle$values
} else {
  sero_study_rle = rle(as.character(sero_data$STUDY_ID))
  sero_study_id = sero_study_rle$values
  N_sero_study = 0
}

#Case
case_study_rle = rle(as.character(case.noti.df$study_id))
N_case_study = length(case_study_rle$lengths)
nrow_case_study = mapply(rep, 1:N_case_study, each = case_study_rle$lengths) %>% unlist()
case_study_id = case_study_rle$values

N_study = N_sero_study + N_case_study
names_study = c(sero_study_id, case_study_id)

if(data_subset != "SA"){
#Speed up calculation the likelihood for sero data:
ind.pop = which(substr(names(sero_data),0,4)=='POP_')
ind.cov = which(substr(names(sero_data),0,4)=='COV_')
sel_0_cov_study = which(sero_data$STUDY_ID %in% c("ETH-10-2014-TSEGAYE", "NGA-8-2008-BABA")) #no coverage for Baba and Tsegaye
sero_data[sel_0_cov_study,ind.cov] = 0
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
pop_by_row = aggregate(pop_v, list(n_row_sero_l), sum)$x

}

#Speed up calculation for case LR data:
year_case_max = aggregate(case.noti.df$year, list(nrow_case_study), max)$x[nrow_case_study]
FOI_case_id_row = year_case_max - case.noti.df$year + case.noti.df$age + 1
FOI_case_id_row_max = max(FOI_case_id_row)
study_agegroup = interaction(case.noti.df$age_group_ID, case.noti.df$study_id) %>% as.vector %>% factor(level = rle(.)$values)
case_by_study_agegroup = aggregate(case.noti.df$case, list(study_agegroup), unique)$x
length_case_study = rle(as.character(nrow_case_study))$lengths

n_row_study_case = rle(as.character(study_agegroup))$values %>% str_split("[.]") %>% sapply(function(x) x[2]) %>% rle
n_row_study = c(mapply(rep, 1:N_sero_study, each = sero_study_rle$lengths) %>% unlist, 
                mapply(rep, N_sero_study + 1:N_case_study, each = n_row_study_case$lengths) %>% unlist)
# mapply(rep, 1:N_study, each = c(sero_study_rle$lengths, levels(study_WHO_agegroup), n_row_study_case$lengths)) %>% unlist()
#retreive back the original case LR data:
case.noti.or.df = unique(case.noti.df[,c("age_group_ID", "country", "study_id")])
case.noti.or.df$age_min = aggregate(case.noti.df$age, by = list(case.noti.df$age_group_ID), min)$x
case.noti.or.df$age_max = aggregate(case.noti.df$age, by = list(case.noti.df$age_group_ID), max)$x
case.noti.or.df$case = case_by_study_agegroup

#Age dependent function####
age_dep_func <- function(sn_par){
  age_dep = dsn(0:99, sn_par[1], sn_par[2], sn_par[3])  #1: postition (range 0:99) 2: scale (positive) 3:skewness (real number)
  return(age_dep)
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
  
  
  #reporting proportion
  rho_case_id <<- max(Age_depend_id) + 1:N_case_study
  par_names <<- c(par_names, paste0("rho_case_", 1:N_case_study))
  
  #VE:
  VE_id <<- max(rho_case_id) + 1
  par_names <<- c(par_names, "VE")
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
  
  FOI_par = unlist(par[names_par == "gamma"]); 
  
  FOI_m = FOI_func(FOI_par) %>% matrix(nrow = 100, ncol = N_study)
  
  age_dep_term_array = array(dim = c(100, N_study))
  for(study in 1:N_study){
    age_dep1 = par[which(names_par == "age_dep1")[study]]
    age_dep2 = par[which(names_par == "age_dep2")[study]]
    age_dep3 = par[which(names_par == "age_dep3")[study]]
    age_dep_term_array[,study] = age_dep_func(c(age_dep1, age_dep2, age_dep3))
  }
  
  #Include the effect of bg infection and age dependent:
  FOI_m = ((FOI_m)*age_dep_term_array) %>% matrix(nrow = 100, ncol = N_study)
  
  VE = par[names_par == "VE"]
  
  rho_case = inverse_logit(unlist(par[names_par == "rho_case"]))
  
  if(data_subset != "SA"){
    FOI_sero_m = FOI_m[, 1:N_sero_study]
    FOI_sero_cumsum = sapply(1:length(age_range_v), function(x) sum(FOI_sero_m[1:(age_range_v[x] + 1), nrow_sero_study_age[x]]))
    seropos = 1-exp(-FOI_sero_cumsum)
    sus_v = (cov_v*VE + (1-cov_v*VE) * seropos) * pop_v
    p = aggregate(sus_v, list(n_row_sero_l), sum)$x/pop_by_row
    p[p > 1 - 1e-16] = 1 - 1e-16
    p[p < 1e-16] = 1e-16
  } else {
    p = NULL
  }
  
  ###Case LR data:
  #case data from LR:
  FOI_case_m = FOI_m[, N_sero_study + 1:N_case_study]
  FOI_case_m = rbind(FOI_case_m, matrix(FOI_case_m[100,], nrow = FOI_case_id_row_max - 100, ncol = N_case_study, byrow = T))
  FOI_case_v = sapply(1:length(nrow_case_study), function(i) FOI_case_m[FOI_case_id_row[i], nrow_case_study[i]])
  FOI_case_m_cumsum = apply(FOI_case_m, 2, cumsum)
  FOI_case_int_v = sapply(1:length(nrow_case_study), function(x) {
    FOI_case_id_int = FOI_case_id_row[x] - 1
    if(FOI_case_id_int == 0) return(1)
    return(FOI_case_m_cumsum[FOI_case_id_int,nrow_case_study[x]]) 
  })
  
  l_rho_case = mapply(rep, rho_case, each = length_case_study) %>% unlist #rep(rho_case, each = N_case_study*100)
  case_sus = case.noti.df$pop*(1 - case.noti.df$cov*VE)
  exp_cases = case_sus*exp(-FOI_case_int_v)*(1 - exp(-FOI_case_v))*l_rho_case
  
  exp_case_by_age_group = aggregate(exp_cases, by = list(study_agegroup), FUN = sum)$x
  exp_case_by_age_group[exp_case_by_age_group < 1e-3] = 1e-3
  
  return(list(sero = p, case = exp_case_by_age_group))
}
