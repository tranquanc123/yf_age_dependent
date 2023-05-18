#library load:
library(dplyr)
library(ggplot2)
library(stringr)
library(cowplot)
library(Hmisc)

setwd("/home/quan/Documents/NotreDame/Alex/Yellow fever/yf_age_exposure/code/")
list2env(readRDS("../output/par_sum_fit.rds"),globalenv())

#load data:
case.noti.LR.df = read.csv("../data/yf_case_data_by_age_group_with_coverage_LR.csv")

#Study name:
source_names = rle(FOI_time_sum_data$Study_id)$values %>% unique()
sero_names = source_names[13:16] %>% strsplit("-") %>% sapply(function(x) paste0(str_to_title(x[4]), " el al ", x[3]))
case_LR_names = source_names[-c(7:16)] %>% strsplit("_"); author_LR = sapply(case_LR_names, function(x) x[1]) ; author_LR[15] = "De Cock"; author_LR[6] = "Ribeiro"; author_LR[2] = "Saad"; year_LR = sapply(case_LR_names, function(x) tail(x, 1)); year_LR[18:19] <- paste0("1998-", year_LR[18:19])
case_LR_names = paste0(author_LR, " et al ", year_LR)
case_PAHO_names = paste0("PAHO - ", c("Bolivia", "Brazil", "Colombia", "Paraguay", "Peru", "Venezuela"))
source_names = c(case_LR_names[1:6], case_PAHO_names, sero_names, case_LR_names[7:20])
N_study = length(source_names)

#color scheme:
three_subset_color = c("#82ad58","#4e7c6f", "#fbbb09")
Subset_order = c("Africa", "South America", "Complete data")

#Age dependent plot:
Age_depend_sum_plot_data = Age_depend_sum_data %>% 
  mutate(data_subset = case_when(
    data_subset == "SA" ~ Subset_order[2],
    data_subset == "AF" ~ Subset_order[1],
    data_subset == "SA_AF" ~ Subset_order[3]
))
Age_depend_sum_plot_data$data_subset = factor(Age_depend_sum_plot_data$data_subset, levels = Subset_order)
Age_dep_plot = ggplot(Age_depend_sum_plot_data, aes(x = age, fill = data_subset))+
  geom_ribbon(aes(ymin =  X2.5., ymax = X97.5.), alpha = 1)+
  geom_line(aes(y = X50.), color = "black")+
  scale_fill_manual(values = three_subset_color)+
  facet_wrap(data_subset ~ ., scale = "free_y")+
  labs(x = "Age", y = "Exposure risk", fill = "Data subset")+
  theme_bw()+
  theme(text = element_text(size = 20), 
        axis.text.x = element_text(size = 14))

ggsave("../plots/Age_dep_plot.png", Age_dep_plot, width = 10, height = 4)

# #FOI plot:
# FOI_time_sum_data$model_type = factor(FOI_time_sum_data$model_type, levels = Subset_order)
# FOI_time_sum_data = FOI_time_sum_data %>% mutate(par = case_when(par == "FOI_time" ~ "FOI (time component)", par == "FOI_time.age" ~ "FOI (time and\nage component)"))
# FOI_time_plot = ggplot(filter(FOI_time_sum_data), aes(x = time, color = model_type, fill = model_type, group = interaction(model_type, Study_id)))+
#   #geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), alpha = 0.1)+
#   geom_line(aes(y = X50.), size = 1)+
#   scale_y_log10()+
#   scale_color_manual(values = three_subset_color)+
#   labs(x = "Year", y = "FOI", color = "Model")+
#   facet_grid(par~model_type)+
#   theme_bw()+
#   theme(text = element_text(size = 16), 
#         axis.text.x = element_text(size = 12))
# FOI_time_plot
# ggsave("../plots/FOI_plot.png", FOI_time_plot, width = 10, height = 6)

#parameter plot:
par_sel_sum_data_plot =  par_sel_sum %>% 
  mutate(data_subset = case_when(
    data_subset == "SA" ~ Subset_order[2],
    data_subset == "AF" ~ Subset_order[1],
    data_subset == "SA_AF" ~ Subset_order[3]
  )) 
par_sel_sum_data_plot$par = strsplit(par_sel_sum_data_plot$par, "[", fixed=TRUE) %>% sapply(function(x) x[1])
par_sel_sum_data_plot = mutate(par_sel_sum_data_plot, par = case_when(
  par == "gamma" ~ "~italic(gamma)",
  par == "rho_case" ~ "~italic(rho)",
  par == "VE" ~ "~italic(VE)",
  par == "age_dep1" ~ "~italic(xi)",
  par == "age_dep2" ~ "~italic(omega)",
  par == "age_dep3" ~ "~italic(alpha)"
))
par_sel_sum_data_plot$Study_names = c(source_names[1:12], rep(NA, 3), source_names[1:12], NA, #SA
                                      source_names[13:30], rep(NA, 3), source_names[17:30], NA, #AF
                                      source_names[c(13:16, 1:6, 17:30, 7:12)], rep(NA, 3), source_names[c(1:6, 17:30, 7:12)], NA) #SA AF
par_sel_sum_data_plot$Study_names = factor(par_sel_sum_data_plot$Study_names, levels = rev(source_names[c(13:30, 1:12)]))
# par_sel_sum$par_names_group = c(rep("~gamma", N_study), "FOI mean", "FOI sd", paste("Age dependent ", 1:3), rep("~rho", N_study - 8),
#                           rep("~theta", N_study), rep("time", N_study), "FOI mean", "FOI sd", paste("Age dependent ", 1:3), rep("~rho", N_study - 8),
#                           rep(c(rep("~theta", N_study), rep("~time", N_study)), 2), "FOI mean", "FOI sd", paste("Age dependent ", 1:3), rep("~rho", N_study - 8))
# par_sel_sum$par_names_group = factor(par_sel_sum$par_names_group, levels = unique(par_sel_sum$par_names_group))
par_sel_sum_data_plot[par_sel_sum_data_plot$par == "~italic(gamma)",1:3]  = log10(par_sel_sum_data_plot[par_sel_sum_data_plot$par == "~italic(gamma)",1:3])
par_sel_sum_data_plot$data_subset = factor(par_sel_sum_data_plot$data_subset, levels = Subset_order)
Constant_model_par_plot = ggplot(par_sel_sum_data_plot[!(is.na(par_sel_sum_data_plot$Study_names)),], aes(y = Study_names, color = data_subset))+
  geom_pointrange(aes(x = X50., xmin = X2.5., xmax = X97.5.), size = 0.75, position = position_dodge(width = 0.5))+
  scale_color_manual(values = three_subset_color)+
  facet_grid(.~par, scale = "free", labeller = labeller(par = label_parsed))+
  labs(x = "Parameters in log10 scale", y = "Study", color = "Data subset")+
  theme_bw()+
  theme(text = element_text(size = 18),
        axis.text.x = element_text(size = 12))
Constant_model_par_plot
ggsave("../plots/Constant_subset_par_plot.png", Constant_model_par_plot, width = 10, height = 8)

#Age dependent (summary table):
Age_dep_par = par_sel_sum_data_plot[row.names(par_sel_sum_data_plot) %>% grep("age_dep", .),]
Age_dep_par[,1:3] = round(Age_dep_par[,1:3], 3)
Age_dep_par[c(4:6, 1:3, 7:9),]

#VE:
VE_dep_par = par_sel_sum_data_plot %>% filter(par == "~italic(VE)")
VE_dep_par[,1:3] = round(VE_dep_par[,1:3], 3)
VE_dep_par[c(2,1,3),]

# # Constant_model_par_plot = ggplot(par_sel_sum[!(is.na(par_sel_sum$Study_names)),] %>% filter(model_type == "Constant"), aes(y = Study_names))+
# #   geom_pointrange(aes(x = X50., xmin = X2.5., xmax = X97.5.), size = 0.75)+
# #   facet_wrap(par_names_group~., scale = "free_x", labeller = labeller(par_names_group = label_parsed))+
# #   scale_x_log10()+
# #   labs(x = "", y = "Study")+
# #   theme_bw()+
# #   theme(text = element_text(size = 18), 
# #         axis.text.x = element_text(size = 12))
# 
# Constant_model_gamma_plot = ggplot(par_sel_sum[!(is.na(par_sel_sum$Study_names)),] %>% filter(model_type == "Constant", par_names_group == "~gamma"), aes(y = Study_names))+
#   geom_pointrange(aes(x = X50., xmin = X2.5., xmax = X97.5.), size = 0.75, color = three_subset_color[3])+
#   facet_wrap(par_names_group~., scale = "free_x", labeller = labeller(par_names_group = label_parsed))+
#   scale_x_log10()+
#   labs(x = "", y = "")+
#   theme_bw()+
#   theme(text = element_text(size = 18), 
#         axis.text.x = element_text(size = 12))
# Constant_model_rho_plot = ggplot(par_sel_sum[!(is.na(par_sel_sum$Study_names)),] %>% filter(model_type == "Constant", par_names_group == "~rho"), aes(y = Study_names))+
#   geom_pointrange(aes(x = X50., xmin = X2.5., xmax = X97.5.), size = 0.75, color = three_subset_color[3])+
#   facet_wrap(par_names_group~., scale = "free_x", labeller = labeller(par_names_group = label_parsed))+
#   labs(x = "", y = "")+
#   theme_bw()+
#   theme(text = element_text(size = 18), 
#         axis.text.x = element_text(size = 12), 
#         axis.text.y = element_blank(), 
#         axis.ticks.y = element_blank())
# Constant_model_par_plot = plot_grid(Constant_model_gamma_plot, Constant_model_rho_plot, rel_widths = c(1.5,1))
# 
# One_outbreak_model_par_plot = ggplot(par_sel_sum[!(is.na(par_sel_sum$Study_names)),] %>% filter(model_type == "One outbreak"), aes(y = Study_names))+
#   geom_pointrange(aes(x = X50., xmin = X2.5., xmax = X97.5.), size = 0.75, color = three_subset_color[1])+
#   facet_wrap(par_names_group~., scale = "free_x", labeller = labeller(par_names_group = label_parsed))+
#   labs(x = "", y = "")+
#   theme_bw()+
#   theme(text = element_text(size = 18), 
#         axis.text.x = element_text(size = 12))
# 
# One_outbreak_model_alpha_plot = ggplot(par_sel_sum[!(is.na(par_sel_sum$Study_names)),] %>% filter(model_type == "One outbreak", par_names_group == "~alpha"), aes(y = Study_names))+
#   geom_pointrange(aes(x = X50., xmin = X2.5., xmax = X97.5.), size = 0.75, color = three_subset_color[1])+
#   facet_wrap(par_names_group~., scale = "free_x", labeller = labeller(par_names_group = label_parsed))+
#   scale_x_log10(breaks = c(1e-3, 1e-2, 1e-1, 1, 10), labels =  c("1e-3", "1e-2", "1e-1", 1, 10))+
#   labs(x = "", y = "")+
#   theme_bw()+
#   theme(text = element_text(size = 18), 
#         axis.text.x = element_text(size = 12))
# 
# One_outbreak_model_data = par_sel_sum[!(is.na(par_sel_sum$Study_names)),] %>% filter(model_type == "One outbreak", par_names_group != "~alpha")
# One_outbreak_model_data$par_names_group = factor(One_outbreak_model_data$par_names_group, levels = c("time", "~rho"))
# One_outbreak_model_rho_plot = ggplot(filter(One_outbreak_model_data, par_names_group == "~rho"), aes(y = Study_names))+
#   geom_pointrange(aes(x = X50., xmin = X2.5., xmax = X97.5.), size = 0.75, color = three_subset_color[1])+
#   facet_wrap(par_names_group~., scale = "free_x", labeller = labeller(par_names_group = label_parsed))+
#   labs(x = "", y = "")+
#   theme_bw()+
#   theme(text = element_text(size = 18), 
#         axis.text.x = element_text(size = 12),
#         axis.text.y = element_blank(), 
#         axis.ticks.y = element_blank(), 
#         strip.text.x = element_text(margin = margin(0.1, 0, 0.1, 0, "cm")))
# One_outbreak_model_time_plot = ggplot(filter(One_outbreak_model_data, par_names_group == "time"), aes(y = Study_names))+
#   geom_pointrange(aes(x = X50., xmin = X2.5., xmax = X97.5.), size = 0.75, color = three_subset_color[1])+
#   facet_wrap(par_names_group~., scale = "free_x", labeller = labeller(par_names_group = label_parsed))+
#   labs(x = "", y = "")+
#   theme_bw()+
#   theme(text = element_text(size = 18), 
#         axis.text.x = element_text(size = 12),
#         axis.text.y = element_blank(), 
#         axis.ticks.y = element_blank(), 
#         strip.text.x = element_text(margin = margin(0.1, 0, 0.1, 0, "cm")))
# 
# One_outbreak_model_par_plot = plot_grid(One_outbreak_model_alpha_plot, One_outbreak_model_time_plot, One_outbreak_model_rho_plot, rel_widths = c(1.75, 1,  1), nrow = 1)
# 
# Two_outbreaks_model_par = par_sel_sum[!(is.na(par_sel_sum$Study_names)),] %>% filter(model_type == "Two outbreaks")
# Two_outbreaks_model_par$Outbreaks = c(rep(c("second", "first"), each = N_study*2), rep(NA, N_study - 8))
# Two_outbreaks_model_par_plot = ggplot(Two_outbreaks_model_par, aes(y = Study_names, color = Outbreaks))+
#   geom_pointrange(aes(x = X50., xmin = X2.5., xmax = X97.5.), position = position_jitter(width = 0, height = 0.15), size = 0.75)+
#   scale_color_manual(breaks = c("first","second"), values = three_subset_color[1:2], na.value = three_subset_color[2])+
#   facet_wrap(par_names_group~., scale = "free_x", labeller = labeller(par_names_group = label_parsed), nrow = 1)+
#   labs(x = "", y = "")+
#   theme_bw()+
#   theme(text = element_text(size = 18), 
#         axis.text.x = element_text(size = 12))
# 
# Two_outbreaks_model_rho_plot = ggplot(filter(Two_outbreaks_model_par, par_names_group == "~rho"), aes(y = Study_names, color = Outbreaks))+
#   geom_pointrange(aes(x = X50., xmin = X2.5., xmax = X97.5.), position = position_jitter(width = 0, height = 0.15), size = 0.75)+
#   scale_color_manual(breaks = c("first","second"), values = three_subset_color[1:2], na.value = three_subset_color[2])+
#   scale_y_discrete(limits = levels(Two_outbreaks_model_par$Study_names))+
#   facet_wrap(par_names_group~., labeller = labeller(par_names_group = label_parsed), nrow = 1)+
#   labs(x = "", y = "")+
#   theme_bw()+
#   theme(text = element_text(size = 18), 
#         axis.text.x = element_text(size = 12))
# 
# Two_outbreaks_model_alpha_plot = ggplot(filter(Two_outbreaks_model_par, par_names_group == "~alpha"), aes(y = Study_names, color = Outbreaks))+
#   geom_pointrange(aes(x = X50., xmin = X2.5., xmax = X97.5.), position = position_jitter(width = 0, height = 0.15), size = 0.75)+
#   scale_color_manual(breaks = c("first","second"), values = three_subset_color[1:2], na.value = three_subset_color[2])+
#   facet_wrap(par_names_group~., scale = "free_x", labeller = labeller(par_names_group = label_parsed), nrow = 1)+
#   labs(x = "", y = "")+
#   scale_x_log10(breaks = c(1e-7, 1e-5, 1e-3, 1e-1, 1, 10), labels =  c("1e-7", "1e-5", "1e-3", "1e-1", 1, 10))+
#   theme_bw()+
#   theme(text = element_text(size = 18), 
#         axis.text.x = element_text(size = 12),
#         axis.text.y = element_blank(), 
#         axis.ticks.y = element_blank(),
#         legend.position = "none")
# 
# Two_outbreaks_model_time_plot = ggplot(filter(Two_outbreaks_model_par, par_names_group == "~time"), aes(y = Study_names, color = Outbreaks))+
#   geom_pointrange(aes(x = X50., xmin = X2.5., xmax = X97.5.), position = position_jitter(width = 0, height = 0.15), size = 0.75)+
#   scale_color_manual(breaks = c("first","second"), values = three_subset_color[1:2], na.value = three_subset_color[2])+
#   facet_wrap(par_names_group~., scale = "free_x", labeller = labeller(par_names_group = label_parsed), nrow = 1)+
#   labs(x = "", y = "")+
#   theme_bw()+
#   theme(text = element_text(size = 18), 
#         axis.text.x = element_text(size = 12),
#         axis.text.y = element_blank(), 
#         axis.ticks.y = element_blank())
# 
# Two_outbreaks_model_par_plot = plot_grid(Two_outbreaks_model_rho_plot, Two_outbreaks_model_alpha_plot, Two_outbreaks_model_time_plot, rel_widths = c(1.75, 1,  1.5), nrow = 1)
# 
# ggsave("../plots/Constant_par_plot.png", Constant_model_par_plot, width = 10, height = 8)
# ggsave("../plots/1ob_par_plot.png", One_outbreak_model_par_plot, width = 10, height = 8)
# ggsave("../plots/2obs_par_plot.png", Two_outbreaks_model_par_plot, width = 12, height = 8)

#Datafit:
sero_data_fit_plot = sero_data_fit %>%
  cbind(CI = binconf(sero_data_fit$POSITIVE, sero_data_fit$SAMPLE_SIZE)) %>%
  mutate(AGE_RANGE = paste0(formatC(AGE_LOWER, width = 2, flag = 0), "-", formatC(AGE_UPPER, width = 2, flag = 0)),
         AGE_MID = (AGE_LOWER + AGE_UPPER)/2) %>%
  mutate(mod_STUDY_ID = mapply(rep, rep(source_names[13:16], 2), rle(sero_data_fit$STUDY_ID)$lengths) %>% unlist)

sero_data_fit_plot = mutate(sero_data_fit_plot, data_subset = case_when(
  data_subset == "AF" ~ Subset_order[1],
  data_subset == "SA_AF" ~ Subset_order[3]
)) 
sero_data_fit_plot$data_subset = factor(sero_data_fit_plot$data_subset, levels = Subset_order)
sero_fit_plot = ggplot(sero_data_fit_plot, aes(x = AGE_MID))+
  geom_pointrange(aes(ymin = CI.Lower, ymax = CI.Upper, y = CI.PointEst))+
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = data_subset), alpha = 0.5)+
  geom_line(aes(y = X50., group = data_subset), color = "black")+
  scale_color_manual(values = three_subset_color[c(1,3)])+
  scale_fill_manual(values = three_subset_color[c(1,3)])+
  labs(x = "Age", y = "Positive proportion", color = "Data subset", fill = "Data subset")+
  facet_wrap(.~mod_STUDY_ID, scale = "free_y")+
  theme_bw()+
  theme(text = element_text(size = 18), 
        axis.text.x = element_text(size = 12))
sero_fit_plot

case_data_fit_data = case_data_fit %>%
  mutate(AGE_RANGE = paste0(formatC(age_min, width = 2, flag = 0), "-", formatC(age_max, width = 2, flag = 0)),
         AGE_MID = (age_min + age_max)/2) %>% 
  mutate(data_subset = case_when(
    data_subset == "AF" ~ Subset_order[1],
    data_subset == "SA" ~ Subset_order[2],
    data_subset == "SA_AF" ~ Subset_order[3]
  )) %>%
  mutate(mod_STUDY_ID = mapply(rep, c(source_names[-c(13:16)], source_names[c(1:6, 17:30, 7:12)]), rle(case_data_fit$study_id)$lengths) %>% unlist)

case_data_fit_data$data_subset = factor(case_data_fit_data$data_subset, levels = Subset_order)
case_data_fit_data$mod_STUDY_ID = factor(case_data_fit_data$mod_STUDY_ID, levels = source_names[c(13:30, 1:12)])
case_fit_plot = ggplot(case_data_fit_data %>% filter(data_subset != "Complete data"), aes(x = AGE_MID))+
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = data_subset), alpha = 0.5)+
  geom_line(aes(y = X50., group = data_subset), color = "black")+
  geom_point(aes(y = case))+
  scale_color_manual(values = three_subset_color)+
  scale_fill_manual(values = three_subset_color)+
  labs(x = "Age", y = "Positive proportion", color = "Data subset", fill = "Data subset")+
  facet_wrap(.~mod_STUDY_ID, scale = "free_y", ncol = 4)+
  theme_bw()+
  theme(text = element_text(size = 18), 
        axis.text.x = element_text(size = 12))
case_fit_plot

ggsave("../plots/Sero_fit_plot.png", sero_fit_plot, width = 12, height = 8)
ggsave("../plots/Case_fit_plot.png", case_fit_plot, width = 15, height = 8)

# #Project case age pattern given the age exposure ####
# #model runs:
# model_type = "Constant"
# #Parameter id:
# Age_depend_id = 31:33
# #Age dependent function:
# age_dep_func <- function(sn_par){
#   age_dep = dsn(0:99, sn_par[1], exp(sn_par[2]), sn_par[3])  #1: postition (range 0:99) 2: scale (positive) 3:skewness (real number)
#   return((age_dep/sum(age_dep)))
# }
# 
# #Par credible interval; 
# par_sel = readRDS(paste0("../output/MCMC_output/MCMC_output_sero_case_noWHO_morethanoneagegroup_", model_type,"_extracted.Rdata"))
# n_samples = nrow(par_sel)
# 
# Age_depend_par = par_sel[,Age_depend_id];
# #Age_depend_par[,2] <- 10^(Age_depend_par[,2])
# age_dep_term_array = array(dim = c(100, 3, n_samples))
# for(i in 1:n_samples){
#   age_dep_term_array[,,i] = age_dep_func(unlist(Age_depend_par[i,]))
# }
# 
# #FOI and age dependent:
# FOI_const = 0.1
# FOI_array = FOI_const*age_dep_term_array
# 
# #Gen cases
# pop.cov = read.csv("/home/quan/Documents/NotreDame/Alex/Yellow fever/yf_ensemble-main/data/adm_1_pop_and_cov_VIMC.csv")
# all_coun = c("NGA", "BRA")#unique(pop.cov$ISO)
# year_range = 2020
# FOI_range = c(1e-4, 1e-3, 1e-2)
# inf_prop = sapply(FOI_range, function(x) exp(-cumsum(0:99)*x)*(1 - exp(-x)))
# 
# cases_proj_all_coun = array(dim = c(length(all_coun), 100, length(year_range), length(FOI_range), 2)) #age, country, year, scen
# 
# for(scen in 1:2){
#   for(FOI_id in 1:length(FOI_range)){
#     for(year_id in 1:length(year_range)){
#       year = year_range[year_id]
#       sel_pop_cov = sapply(all_coun, function(x) {
#         age_pop = filter(pop.cov, ISO == x, YEAR == year) %>% select(starts_with("POP"))
#         age_pop_sum = colSums(age_pop)
#         age_cov_sum = ((filter(pop.cov, ISO == x, YEAR == year) %>% select(starts_with("COV"))*age_pop) %>% colSums())/age_pop_sum
#         age_cov_sum[is.na(age_cov_sum)] <- 0
#         return(list(age_pop_sum, age_cov_sum))
#       })
#       
#       pop_mat = Reduce(rbind, sel_pop_cov[1, ])
#       if(scen == 1){
#         inf_scen = matrix(inf_prop[,FOI_id], nrow = length(all_coun), ncol = 100, byrow = T)*pop_mat/rowSums(pop_mat)
#       } else {
#         inf_scen = t(sapply(all_coun, function(x) inf_prop[,FOI_id]*(1 - unlist(sel_pop_cov[2,x]))))*pop_mat/rowSums(pop_mat)
#       }
#       
#       cases_proj_all_coun[,, year_id, FOI_id, scen] = inf_scen/apply(inf_scen, 1, sum)
#     }
#   }
# }
# 
# cases_gen_data_plot = data.frame(expand.grid(ISO = all_coun, age = 0:99, year = year_range, FOI = FOI_range, scen = c("No vaccination", "With vaccination")), 
#                                  IR = as.vector(cases_proj_all_coun))
# 
# ggplot(filter(cases_gen_data_plot), aes(x = age, y = IR, color = as.factor(FOI), group = interaction(FOI, scen), linetype = scen))+
#   geom_line()+
#   facet_grid(interaction(FOI, ISO) ~ year, scale = "free")+
#   theme_bw()
# 
# #compare with PAHO data for Brazil:
# BRA_PAHO_data = data.frame(age_l = seq(0, 90, 5), age_u = seq(4, 94, 5), 
#                            cases = c(3,2,11,33,41,31,42,45,33,31,22,12,8,8,3,0,1,0,1)) %>%
#   mutate(age_m = (age_l + age_u)/2)
# 
# plot(BRA_PAHO_data$age_m, BRA_PAHO_data$cases, type = "l")

#