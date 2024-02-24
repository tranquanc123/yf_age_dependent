#library load:
library(dplyr)
library(ggplot2)
library(stringr)
library(cowplot)
library(Hmisc)
library(gridExtra)
library(patchwork)
library(grid)
library(gtable)

setwd("/Users/tranm/Documents/NotreDame/Alex/Yellow fever/yf_age_exposure/code/")
list2env(readRDS("../output/par_sum_fit_heir.rds"),globalenv())

#load data:
sero_data = read.csv("../data/yf_sero_data_with_coverage_mod.csv")
sero_data = filter(sero_data, YEAR > 1980) #only get data from 1980
sero_data = filter(sero_data, SOURCE %in% rle(sero_data$SOURCE)$values[which(rle(sero_data$SOURCE)$lengths > 1)]) # study with age group > 1
case.noti.LR.df = read.csv("../data/yf_case_data_by_age_group_with_coverage_LR.csv")
case.noti.PAHO.df = read.csv("../data/yf_case_data_by_age_group_with_coverage_PAHO.csv") %>% filter(!(study_id %in% c("ECU_PAHO", "ARG_PAHO")))

#different subset of data:
case.noti.df = rbind(case.noti.LR.df, case.noti.PAHO.df)
#case.noti.df = case.noti.df[case.noti.df$case != 0,]
case.noti.df[is.na(case.noti.df$cov),]$pop = 1e-10
case.noti.df[is.na(case.noti.df$cov),]$cov = 0

sero_study_id = rle(as.character(sero_data$STUDY_ID))$values
case_study_id = rle(as.character(case.noti.df$study_id))$values

names_study = c(sero_study_id, case_study_id)

coun_list_study = unique(Age_depend_sum_data$country)
#load prior:
prior_gamma = readRDS("../data/gamma_lnorm_priors_constant_nonorm.rds")
readRDS("../data/FOI_rho_prior.rds") %>% list2env(globalenv())

#Study name:
source_names = names_study
sero_names = source_names[1:4] %>% strsplit("-") %>% sapply(function(x) paste0(str_to_title(x[4]), " el al ", x[3]))
case_LR_names = source_names[c(5:24)] %>% strsplit("_"); author_LR = sapply(case_LR_names, function(x) x[1]) ; author_LR[15] = "De Cock"; author_LR[6] = "Ribeiro"; author_LR[2] = "Saad"; year_LR = sapply(case_LR_names, function(x) tail(x, 1)); year_LR[18:19] <- paste0("1998-", year_LR[18:19])
case_LR_names = paste0(author_LR, " et al ", year_LR)
case_PAHO_names = paste0("PAHO - ", c("Bolivia", "Brazil", "Colombia", "Paraguay", "Peru", "Venezuela"))
source_names = c(sero_names, case_LR_names, case_PAHO_names)#c(case_LR_names[1:6], case_PAHO_names, sero_names, case_LR_names[7:20])
N_study = length(source_names)

quantile95cri <- function(x) quantile(x, probs = c(0.025, 0.25, 0.5, 0.75, 0.975))

#color scheme:
three_subset_color = c("#5C0029","#4e7c6f", "#fbbb09")
Subset_order = c("Africa", "South America", "Complete data")

#Age dependent plot:
par_dist = readRDS("../output/Age_exp_curve_SA_AF_Constant.rds")
# n_rep = 1e2
# rand10_curves = par_dist[,sample(1:ncol(par_dist), n_rep)]
# IQR_curves = apply(par_dist, 1, quantile, c(0.25, 0.75))
# #Age_depend_global_sum_data = cbind(Age_depend_global_sum_data, data.frame(t(IQR_curves)))
# Age_depend_global_sum_rand10_data = rep(list(Age_depend_global_sum_data), 1 + n_rep) %>% Reduce(rbind, .)
# Age_depend_global_sum_rand10_data$X50. = c(Age_depend_global_sum_data$X50., as.vector(rand10_curves))
# #Age_depend_global_sum_rand10_data$X2.5. = c(Age_depend_global_sum_data$X2.5., rep(NA, nrow(Age_depend_global_sum_data)*10))
# #Age_depend_global_sum_rand10_data$X97.5. = c(Age_depend_global_sum_data$X97.5., rep(NA, nrow(Age_depend_global_sum_data)*10))
# Age_depend_global_sum_rand10_data$replicate = rep(c("or", 1:n_rep), each = nrow(Age_depend_global_sum_data))
# #rm(par_dist)
# 
# Age_dep_plot = ggplot(Age_depend_global_sum_data, aes(x = age))+
#   geom_ribbon(aes(ymin =  X2.5., ymax = X97.5.), alpha = 0.5, fill = three_subset_color[3])+
#   geom_ribbon(aes(ymin =  X25., ymax = X75.), alpha = 0.7, fill = three_subset_color[3])+#, color =  "#f03b20")+
#   geom_line(aes(y = X50.), color = "black", size = 1)+
#   scale_y_continuous(limits = c(0, 0.03))+
#   #labs(x = "Age", y = "Exposure risk")+
#   theme_bw()+
#   theme(text = element_text(size = 20), 
#         axis.text.x = element_text(size = 14), axis.title = element_blank())
# 
# Age_dep_plot_replicate = ggplot(Age_depend_global_sum_rand10_data, aes(x = age))+
#   geom_ribbon(aes(ymin =  X2.5., ymax = X97.5.), alpha = 0.5, fill = three_subset_color[3])+
#   geom_line(aes(y = X50., color = replicate), size = 1)+
#   #scale_color_manual(values = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a', "black"))+
#   #labs(x = "Age", y = "Exposure risk")+
#   #scale_y_continuous(limits = c(0, 0.03))+
#   theme_bw()+
#   theme(text = element_text(size = 20), 
#         axis.text.x = element_text(size = 14), axis.title = element_blank(),
#         legend.position = "None")
# 
# Age_dep_plot_comb = ((Age_dep_plot + Age_dep_plot_replicate) + plot_layout() + plot_annotation(tag_levels = "A")) %>% 
#   cowplot::as_gtable() %>%
#     arrangeGrob(left = textGrob("Exposure risk",
#                                 rot = 90,
#                                 vjust = 1,
#                                 gp = gpar(fontsize = 20)),
#                 bottom = textGrob("Age",
#                                   vjust = 0.1,
#                                   hjust = 0.2,
#                                   gp = gpar(fontsize = 20))
#                 )
# plot(Age_dep_plot_comb)
# ggsave("../plots/Age_dep_plot.png", Age_dep_plot_comb, width = 10, height = 5)
# 
# Age_dep_coun_plot = ggplot(Age_depend_sum_data, aes(x = age))+
#   geom_ribbon(aes(ymin =  X2.5., ymax = X97.5.), alpha = 0.5, fill = three_subset_color[3])+
#   geom_ribbon(aes(ymin =  X25., ymax = X75.), alpha = 1, fill = three_subset_color[3])+
#   geom_line(aes(y = X50.), color = "black")+
#   labs(x = "Age", y = "Exposure risk")+
#   facet_wrap(.~country, scales = "free_y")+
#   theme_bw()+
#   theme(text = element_text(size = 20), 
#         axis.text.x = element_text(size = 14))
# Age_dep_coun_plot
# ggsave("../plots/Age_dep_coun_plot.png", Age_dep_coun_plot, width = 12, height = 8)

#scaled risk:
par_dist = par_dist[,1:1e4]
avg_risk_age_global = mean(par_dist)#apply(par_dist, 2, mean)
re_avg_risk_age_global = exp(log(par_dist) - log(avg_risk_age_global))
re_avg_risk_age_global_sum = apply(re_avg_risk_age_global, 1, quantile95cri)
re_avg_risk_age_global_sum_data = data.frame(age = 0:99, t(re_avg_risk_age_global_sum))

Age_dep_plot = ggplot(re_avg_risk_age_global_sum_data, aes(x = age))+
  geom_ribbon(aes(ymin =  X2.5., ymax = X97.5.), alpha = 0.5, fill = three_subset_color[3])+
  geom_ribbon(aes(ymin =  X25., ymax = X75.), alpha = 0.7, fill = three_subset_color[3])+#, color =  "#f03b20")+
  geom_line(aes(y = X50.), color = "black", size = 1)+
  scale_y_continuous(limits = c(0, 3.5))+
  labs(title = "Country-specific")+
  theme_bw()+
  theme(text = element_text(size = 20), 
        axis.text.x = element_text(size = 14), axis.title = element_blank())

n_rand = 1e1
par_dist_10rand = par_dist[,sample(1:1e4, n_rand)]
avg_risk_age_global_10rand = mean(par_dist)
re_avg_risk_age_global_10rand = exp(log(par_dist_10rand) - log(avg_risk_age_global_10rand))
re_avg_risk_age_global_10rand_sum_data = rep(list(re_avg_risk_age_global_sum_data), 1 + n_rand) %>% Reduce(rbind, .)
re_avg_risk_age_global_10rand_sum_data$X50. = c(re_avg_risk_age_global_sum_data$X50., as.vector(re_avg_risk_age_global_10rand))
re_avg_risk_age_global_10rand_sum_data$replicate = rep(c("or", 1:n_rand), each = nrow(re_avg_risk_age_global_sum_data))

Age_dep_plot_replicate = ggplot(re_avg_risk_age_global_10rand_sum_data, aes(x = age))+
  geom_ribbon(aes(ymin =  X2.5., ymax = X97.5.), alpha = 0.5, fill = three_subset_color[3])+
  geom_line(aes(y = X50., color = replicate), size = 1)+
  scale_color_manual(values = c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a', "black"))+
  labs(title = "Country-specific")+
  scale_y_continuous(limits = c(0, 3.5))+
  theme_bw()+
  theme(text = element_text(size = 20), 
        axis.text.x = element_text(size = 14), axis.title = element_blank(),
        legend.position = "None")

#add the same age dept curve model as well:
par_dist = readRDS("../output/Age_exp_globalcurve_SA_AF_Constant.rds")
avg_risk_age_global = mean(par_dist)#apply(par_dist, 2, mean)
re_avg_risk_age_global = exp(log(par_dist) - log(avg_risk_age_global))
re_avg_risk_age_global_sum = apply(re_avg_risk_age_global, 1, quantile95cri)
re_avg_risk_age_global_sum_data_globalcurve = data.frame(age = 0:99, t(re_avg_risk_age_global_sum))
Age_dep_plot_globalcurve = ggplot(re_avg_risk_age_global_sum_data_globalcurve, aes(x = age))+
  geom_ribbon(aes(ymin =  X2.5., ymax = X97.5.), alpha = 0.5, fill = three_subset_color[2])+
  geom_ribbon(aes(ymin =  X25., ymax = X75.), alpha = 0.7, fill = three_subset_color[2])+#, color =  "#f03b20")+
  geom_line(aes(y = X50.), color = "black", size = 1)+
  scale_y_continuous(limits = c(0, 3.5))+
  labs(title = "Geographically invariant")+
  theme_bw()+
  theme(text = element_text(size = 20), 
        axis.text.x = element_text(size = 14), axis.title = element_blank())
Age_dep_plot_comb = ((Age_dep_plot + Age_dep_plot_replicate + Age_dep_plot_globalcurve) + plot_layout() + plot_annotation(tag_levels = "A")) %>% 
  cowplot::as_gtable() %>%
  arrangeGrob(left = textGrob("Exposure risk\n(fold difference compared to the average)",
                              rot = 90,
                              vjust = 1,
                              gp = gpar(fontsize = 15)),
              bottom = textGrob("Age",
                                vjust = 0.1,
                                hjust = 0.2,
                                gp = gpar(fontsize = 20))
  )
plot(Age_dep_plot_comb)
ggsave("../plots/Age_dep_plot.png", Age_dep_plot_comb, width = 15, height = 6)

Age_dep_coun_plot = ggplot(Age_depend_sum_data, aes(x = age))+
  geom_ribbon(aes(ymin =  X2.5., ymax = X97.5.), alpha = 0.5, fill = three_subset_color[3])+
  geom_ribbon(aes(ymin =  X25., ymax = X75.), alpha = 1, fill = three_subset_color[3])+
  geom_line(aes(y = X50.), color = "black")+
  labs(x = "Age", y = "Exposure risk\n(fold difference compared to the average)")+
  facet_wrap(.~country, scales = "free_y")+
  theme_bw()+
  theme(text = element_text(size = 20), 
        axis.text.x = element_text(size = 14))
Age_dep_coun_plot
ggsave("../plots/Age_dep_coun_plot.png", Age_dep_coun_plot, width = 12, height = 8)
# 
# #AF only:
# Age_dep_plot = ggplot(Age_depend_sum_plot_data %>% filter(data_subset == "Africa"), aes(x = age))+
#   geom_ribbon(aes(ymin =  X2.5., ymax = X97.5.), alpha = 1, fill = three_subset_color[1])+
#   geom_line(aes(y = X50.), color = "black")+
#   labs(x = "Age", y = "Exposure risk")+
#   theme_bw()+
#   theme(text = element_text(size = 20), 
#         axis.text.x = element_text(size = 14))
# ggsave("../plots/Age_dep_plot.png", Age_dep_plot, width = 6, height = 4)

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
# par_sel_sum_data_plot =  par_sel_sum %>% 
#   mutate(data_subset = case_when(
#     data_subset == "SA" ~ Subset_order[2],
#     data_subset == "AF" ~ Subset_order[1],
#     data_subset == "SA_AF" ~ Subset_order[3]
#   )) 
# par_sel_sum_data_plot$par = strsplit(par_sel_sum_data_plot$par, "[", fixed=TRUE) %>% sapply(function(x) x[1])
# par_sel_sum_data_plot = mutate(par_sel_sum_data_plot, par = case_when(
#   par == "gamma" ~ "~italic(gamma)",
#   par == "rho_case" ~ "~italic(rho)",
#   par == "VE" ~ "~italic(VE)",
#   par == "age_dep1" ~ "~italic(xi)",
#   par == "age_dep2" ~ "~italic(omega)",
#   par == "age_dep3" ~ "~italic(alpha)"
# ))
# par_sel_sum_data_plot$Study_names = c(#source_names[1:12], rep(NA, 3), source_names[1:12], NA, #SA
#                                       source_names[13:30], rep(NA, 3), source_names[17:30], NA) #AF
#                                       #source_names[c(13:16, 1:6, 17:30, 7:12)], rep(NA, 3), source_names[c(1:6, 17:30, 7:12)], NA) #SA AF
# par_sel_sum_data_plot$Study_names = factor(par_sel_sum_data_plot$Study_names, levels = rev(source_names[c(13:30, 1:12)]))
# #par_sel_sum_data_plot$Study_names = factor(par_sel_sum_data_plot$Study_names, levels = rev(source_names[c(13:30, 1:12)]))
# # par_sel_sum$par_names_group = c(rep("~gamma", N_study), "FOI mean", "FOI sd", paste("Age dependent ", 1:3), rep("~rho", N_study - 8),
# #                           rep("~theta", N_study), rep("time", N_study), "FOI mean", "FOI sd", paste("Age dependent ", 1:3), rep("~rho", N_study - 8),
# #                           rep(c(rep("~theta", N_study), rep("~time", N_study)), 2), "FOI mean", "FOI sd", paste("Age dependent ", 1:3), rep("~rho", N_study - 8))
# # par_sel_sum$par_names_group = factor(par_sel_sum$par_names_group, levels = unique(par_sel_sum$par_names_group))
# #par_sel_sum_data_plot[par_sel_sum_data_plot$par == "~italic(gamma)",1:3]  = log10(par_sel_sum_data_plot[par_sel_sum_data_plot$par == "~italic(gamma)",1:3])
# par_sel_sum_data_plot$data_subset = factor(par_sel_sum_data_plot$data_subset, levels = Subset_order)
# Constant_model_par_plot = ggplot(par_sel_sum_data_plot[!(is.na(par_sel_sum_data_plot$Study_names)),], aes(y = Study_names, color = data_subset))+
#   geom_pointrange(aes(x = X50., xmin = X2.5., xmax = X97.5.), size = 0.75, position = position_dodge(width = 0.5))+
#   scale_color_manual(values = three_subset_color)+
#   facet_grid(.~par, scale = "free", labeller = labeller(par = label_parsed))+
#   labs(x = "Parameters in log10 scale", y = "Study", color = "Data subset")+
#   theme_bw()+
#   theme(text = element_text(size = 18),
#         axis.text.x = element_text(size = 12))
# Constant_model_par_plot
# ggsave("../plots/Constant_subset_par_plot.png", Constant_model_par_plot, width = 10, height = 8)

#SA_AF only:
par_sel_sum_data_plot = mutate(par_sel_sum, par = case_when(
  par == "gamma" ~ "~italic(gamma)",
  par == "rho_case" ~ "~italic(rho)",
  par == "VE" ~ "~italic(VE)",
  par == "age_dep1" ~ "~italic(xi)",
  par == "age_dep2" ~ "~italic(omega)",
  par == "age_dep3" ~ "~italic(alpha)"
))
par_sel_sum_data_plot$Study_names = c(source_names, rep(NA, length(coun_list_study)*3 + 6), source_names[5:30], NA) #SA AF
par_sel_sum_data_plot$Country = c(rep(NA, length(source_names)), rep(coun_list_study, 3), rep(NA, 6 + 26 + 1))
SA_AF_Constant_model_par_plot = par_sel_sum_data_plot[!(is.na(par_sel_sum_data_plot$Study_names)),]
SA_AF_Constant_model_par_plot$Study_names = factor(SA_AF_Constant_model_par_plot$Study_names, 
                                                levels = rev(c(sort(SA_AF_Constant_model_par_plot$Study_names[1:4]), sort(SA_AF_Constant_model_par_plot$Study_names[5:30]))))
SA_AF_Constant_model_gamma =((SA_AF_Constant_model_par_plot %>% filter(par == "~italic(gamma)"))[,1:3])
(SA_AF_Constant_model_gamma[order(SA_AF_Constant_model_gamma$X50.),][c(1,nrow(SA_AF_Constant_model_gamma)),]) %>% round(3)

SA_AF_Constant_model_rho = exp((SA_AF_Constant_model_par_plot %>% filter(par == "~italic(rho)"))[,1:3])/(exp((SA_AF_Constant_model_par_plot %>% filter(par == "~italic(rho)"))[,1:3]) + 1)
(SA_AF_Constant_model_rho[order(SA_AF_Constant_model_rho$X50.),][c(1,nrow(SA_AF_Constant_model_rho)),])

sero_ISO_source = sero_data[,c(1,8)] %>% unique
case_ISO_source = ((case.noti.df)[,2:3]) %>% unique; case_ISO_source$country[21:26] = "BRA"
matched_SA_AF_study_prior = c(match(sero_ISO_source$ISO, prior_gamma$ISO), match(case_ISO_source$country, prior_gamma$ISO))
rho_prior_SA_AF = data.frame(matrix(reporting_prop_prior, nrow = 26, ncol = 2, byrow = T))
prior_gamma_SA_AF = rbind(prior_gamma[matched_SA_AF_study_prior,-1], rho_prior_SA_AF)
prior_gamma_SA_AF_CI = data.frame(prior_UI_l = qnorm(0.025, prior_gamma_SA_AF$X1, prior_gamma_SA_AF$X2), prior_UI_u = qnorm(0.975, prior_gamma_SA_AF$X1, prior_gamma_SA_AF$X2))

SA_AF_prior_Constant_model = SA_AF_Constant_model_par_plot;
SA_AF_prior_Constant_model$X2.5. = prior_gamma_SA_AF_CI$prior_UI_l; SA_AF_prior_Constant_model$X97.5. = prior_gamma_SA_AF_CI$prior_UI_u; SA_AF_prior_Constant_model$X50. = NULL
Constant_model_par_plot = ggplot(data = NULL, aes(y = Study_names))+
  geom_linerange(data = SA_AF_prior_Constant_model,  aes(xmin = X2.5., xmax = X97.5.), linewidth = 2, color = "red", alpha = 0.1)+
  geom_pointrange(data = SA_AF_Constant_model_par_plot, aes(x = X50., xmin = X2.5., xmax = X97.5.), size = 0.25, color = three_subset_color[1])+
  facet_grid(.~par, scale = "free", labeller = labeller(par = label_parsed))+
  labs(x = "Parameters in log10 scale", y = "Study")+
  theme_bw()+
  theme(text = element_text(size = 14),
        axis.text.x = element_text(size = 12))
Constant_model_par_plot
ggsave("../plots/Constant_par_plot.png", Constant_model_par_plot, width = 6, height = 6)

#Age dependent (summary table):
Age_dep_par = par_sel_sum_data_plot[row.names(par_sel_sum_data_plot) %>% grep("age_dep", .),]
Age_dep_par[,1:3] = round(Age_dep_par[,1:3], 3)
Age_dep_par_global = filter(Age_dep_par, is.na(Country))
round(Age_dep_par_global[,1:3], 1)

Age_dept_mode
Age_dept_ratio_ho
Age_dept_ratio_hy
Age_high_to_low_ratio
#Age_dep_par[4:6,] #AF only
Age_dep_par_coun = filter(Age_dep_par, !is.na(Country))
Age_dep_par_coun[rev((Age_dep_par_coun %>% group_by(par) %>% summarise(which.min(X50.)))$`which.min(X50.)`) + (0:2)*length(coun_list_study),]
Age_dep_par_coun[rev((Age_dep_par_coun %>% group_by(par) %>% summarise(which.max(X50.)))$`which.max(X50.)`) + (0:2)*length(coun_list_study),]

Age_dep_par_coun_plot = ggplot(Age_dep_par_coun, aes(y = Country, x = X50.))+
  geom_pointrange(aes(xmin = X2.5., xmax = X97.5.), color = three_subset_color[1])+
  facet_wrap(.~par, scale = "free_x", labeller = labeller(par = label_parsed))+
  labs(x = "Parameters", y = "Country")+
  theme_bw()+
  theme(text = element_text(size = 14),
        axis.text.x = element_text(size = 12))
ggsave("../plots/Age_dep_par_coun_plot.png", Age_dep_par_coun_plot, width = 7, height = 6)

#VE:
VE_dep_par = par_sel_sum_data_plot %>% filter(par == "~italic(VE)")
VE_dep_par[,1:3] = round(VE_dep_par[,1:3], 3)
VE_dep_par
#VE_dep_par[c(2,1,3),]
#VE_dep_par[2,] #AF only

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
  cbind(CI = sapply(c(0.025, 0.5, 0.975), function(x) qbeta(x, sero_data_fit$POSITIVE + 1, sero_data_fit$SAMPLE_SIZE  - sero_data_fit$POSITIVE + 1))) %>%
  mutate(AGE_RANGE = paste0(formatC(AGE_LOWER, width = 2, flag = 0), "-", formatC(AGE_UPPER, width = 2, flag = 0)),
         AGE_MID = (AGE_LOWER + AGE_UPPER)/2) %>%
  mutate(mod_STUDY_ID = mapply(rep, source_names[1:4], rle(sero_data_fit$STUDY_ID)$lengths) %>% unlist)

# sero_data_fit_plot = mutate(sero_data_fit_plot, data_subset = case_when(
#   data_subset == "AF" ~ Subset_order[1],
#   data_subset == "SA_AF" ~ Subset_order[3]
# )) 
#sero_data_fit_plot$data_subset = factor(sero_data_fit_plot$data_subset, levels = Subset_order)
sero_fit_plot = ggplot(sero_data_fit_plot, aes(x = AGE_MID))+
  geom_pointrange(aes(ymin = CI.1, ymax = CI.3, y = CI.2))+
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = data_subset), alpha = 0.5, fill = three_subset_color[3])+
  geom_line(aes(y = X50., group = data_subset), color = "black")+
  #scale_color_manual(values = three_subset_color[c(1,3)])+
  #scale_fill_manual(values = three_subset_color[c(1,3)])+
  labs(x = "Age", y = "Positive proportion", color = "Data subset", fill = "Data subset")+
  facet_wrap(.~mod_STUDY_ID, scale = "free_y")+
  theme_bw()+
  theme(text = element_text(size = 18), 
        axis.text.x = element_text(size = 12))
sero_fit_plot

case_data_fit_data = case_data_fit %>%
  mutate(AGE_RANGE = paste0(formatC(age_min, width = 2, flag = 0), "-", formatC(age_max, width = 2, flag = 0)),
         AGE_MID = (age_min + age_max)/2) %>% 
  # mutate(data_subset = case_when(
  #   data_subset == "AF" ~ Subset_order[1],
  #   data_subset == "SA" ~ Subset_order[2],
  #   data_subset == "SA_AF" ~ Subset_order[3]
  # )) %>%
  mutate(mod_STUDY_ID = mapply(rep, c(source_names[5:30]), rle(case_data_fit$study_id)$lengths) %>% unlist)
  #mutate(mod_STUDY_ID = mapply(rep, c(source_names[-c(13:16)], source_names[c(1:6, 17:30, 7:12)]), rle(case_data_fit$study_id)$lengths) %>% unlist)

#case_data_fit_data$data_subset = factor(case_data_fit_data$data_subset, levels = Subset_order)
case_data_fit_data$mod_STUDY_ID = factor(case_data_fit_data$mod_STUDY_ID, levels = source_names[5:30])
case_fit_plot = ggplot(case_data_fit_data %>% filter(data_subset != "Complete data"), aes(x = AGE_MID))+
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = data_subset), alpha = 0.5, fill = three_subset_color[3])+
  geom_line(aes(y = X50., group = data_subset), color = "black")+
  geom_point(aes(y = case))+
  #scale_color_manual(values = three_subset_color)+
  #scale_fill_manual(values = three_subset_color)+
  labs(x = "Age", y = "Positive proportion", color = "Data subset", fill = "Data subset")+
    facet_wrap(.~mod_STUDY_ID, scale = "free_y", ncol = 5)+
  theme_bw()+
  theme(text = element_text(size = 18), 
        axis.text.x = element_text(size = 12))
case_fit_plot

# #AF only:
# sero_fit_plot = ggplot(sero_data_fit_plot %>% filter(data_subset == "Africa"), aes(x = AGE_MID))+
#   geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), alpha = 1, fill = three_subset_color[1])+
#   geom_pointrange(aes(ymin = CI.Lower, ymax = CI.Upper, y = CI.PointEst))+
#   geom_line(aes(y = X50.), color = "black")+
#   labs(x = "Age", y = "Positive proportion")+
#   facet_wrap(.~mod_STUDY_ID, scale = "free_y")+
#   theme_bw()+
#   theme(text = element_text(size = 18), 
#         axis.text.x = element_text(size = 12))
# case_fit_plot = ggplot(case_data_fit_data %>% filter(!(data_subset %in% c("Complete data", "South America"))), aes(x = AGE_MID))+
#   geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), alpha = 1, fill = three_subset_color[1])+
#   geom_line(aes(y = X50., group = data_subset), color = "black")+
#   geom_point(aes(y = case))+
#   labs(x = "Age", y = "Positive proportion")+
#   facet_wrap(.~mod_STUDY_ID, scale = "free_y", ncol = 4)+
#   theme_bw()+
#   theme(text = element_text(size = 18), 
#         axis.text.x = element_text(size = 12))

ggsave("../plots/Sero_fit_plot.png", sero_fit_plot, width = 8, height = 6)
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
#Alternative models####
heir_country_model = readRDS("../output/par_sum_fit_heir.rds")
heir_globalcurve_model = readRDS("../output/par_sum_fit_heir_globalcurve.rds")
noage_model = readRDS("../output/par_sum_fit_heir_noage.rds")

three_models_name = c("None", "Country-specific", "Geographically invariant")

#Par summary:
par_sel_sum_data_globalcurve_plot = heir_globalcurve_model$par_sel_sum %>% mutate(par = case_when(
  par == "gamma" ~ "~italic(gamma)",
  par == "rho_case" ~ "~italic(rho)",
  par == "VE" ~ "~italic(VE)"#,
  # par == "age_dep1" ~ "~italic(xi)",
  # par == "age_dep2" ~ "~italic(omega)",
  # par == "age_dep3" ~ "~italic(alpha)"
))
par_sel_sum_data_globalcurve_plot$Study_names = c(source_names, rep(NA, 3), source_names[5:30], NA) 
par_sel_sum_data_globalcurve_plot = par_sel_sum_data_globalcurve_plot[!is.na(par_sel_sum_data_globalcurve_plot$Study_names),]
par_sel_sum_data_globalcurve_plot$Model = three_models_name[3]

par_sel_sum_data_noage_plot = noage_model$par_sel_sum %>% mutate(par = case_when(
  par == "gamma" ~ "~italic(gamma)",
  par == "rho_case" ~ "~italic(rho)",
  par == "VE" ~ "~italic(VE)"#,
  # par == "age_dep1" ~ "~italic(xi)",
  # par == "age_dep2" ~ "~italic(omega)",
  # par == "age_dep3" ~ "~italic(alpha)"
))
par_sel_sum_data_noage_plot$Study_names = c(source_names, source_names[5:30], NA)
par_sel_sum_data_noage_plot = par_sel_sum_data_noage_plot[!is.na(par_sel_sum_data_noage_plot$Study_names),]
par_sel_sum_data_noage_plot$Model = three_models_name[1]

par_sel_sum_data_coun_spec_plot = mutate(heir_country_model$par_sel_sum, par = case_when(
  par == "gamma" ~ "~italic(gamma)",
  par == "rho_case" ~ "~italic(rho)"
))
par_sel_sum_data_coun_spec_plot = par_sel_sum_data_coun_spec_plot[!(is.na(par_sel_sum_data_coun_spec_plot$par)),]
par_sel_sum_data_coun_spec_plot$Study_names = c(source_names, source_names[5:30]) #SA AF
par_sel_sum_data_coun_spec_plot$Model = three_models_name[2]

par_sel_sum_data_plot = rbind(par_sel_sum_data_noage_plot, par_sel_sum_data_coun_spec_plot, par_sel_sum_data_globalcurve_plot)
#par_sel_sum_data_plot$Country = c(rep(NA, length(source_names)), rep(coun_list_study, 3), rep(NA, 6 + 26 + 1))

par_sel_sum_data_plot$Study_names = factor(par_sel_sum_data_plot$Study_names, 
                                                   levels = rev(c(sort(source_names[1:4]), sort(source_names[5:30]))))
par_sel_sum_data_plot[which(par_sel_sum_data_plot$par == "~italic(gamma)"),1:5] = log10(par_sel_sum_data_plot[which(par_sel_sum_data_plot$par == "~italic(gamma)"),1:5])
#SA_AF_Constant_model_gamma =((par_sel_sum_data_plot %>% filter(par == "~italic(gamma)"))[,1:3])
#(SA_AF_Constant_model_gamma[order(SA_AF_Constant_model_gamma$X50.),][c(1,nrow(SA_AF_Constant_model_gamma)),]) %>% round(3)

#SA_AF_Constant_model_rho = exp((par_sel_sum_data_plot %>% filter(par == "~italic(rho)"))[,1:3])/(exp((par_sel_sum_data_plot %>% filter(par == "~italic(rho)"))[,1:3]) + 1)
#(SA_AF_Constant_model_rho[order(SA_AF_Constant_model_rho$X50.),][c(1,nrow(SA_AF_Constant_model_rho)),])

sero_ISO_source = sero_data[,c(1,8)] %>% unique
case_ISO_source = ((case.noti.df)[,2:3]) %>% unique; case_ISO_source$country[21:26] = "BRA"
matched_SA_AF_study_prior = c(match(sero_ISO_source$ISO, prior_gamma$ISO), match(case_ISO_source$country, prior_gamma$ISO))
rho_prior_SA_AF = data.frame(matrix(reporting_prop_prior, nrow = 26, ncol = 2, byrow = T))
prior_gamma_heir = rbind(prior_gamma[matched_SA_AF_study_prior,-1], rho_prior_SA_AF)
matched_SA_AF_study_prior_noage = c(match(sero_ISO_source$ISO, FOI_prior$ISO), match(case_ISO_source$country, FOI_prior$ISO))
prior_gamma_noage = rbind(FOI_prior[matched_SA_AF_study_prior_noage,-1], rho_prior_SA_AF)
prior_gamma_SA_AF = rbind(prior_gamma_noage, prior_gamma_heir, prior_gamma_heir)
prior_gamma_SA_AF_CI = data.frame(prior_UI_l = qnorm(0.025, prior_gamma_SA_AF$X1, prior_gamma_SA_AF$X2), prior_UI_u = qnorm(0.975, prior_gamma_SA_AF$X1, prior_gamma_SA_AF$X2))

SA_AF_prior_Constant_model = par_sel_sum_data_plot;
SA_AF_prior_Constant_model$X2.5. = prior_gamma_SA_AF_CI$prior_UI_l; SA_AF_prior_Constant_model$X97.5. = prior_gamma_SA_AF_CI$prior_UI_u; SA_AF_prior_Constant_model$X50. = NULL
Constant_model_par_plot = ggplot(data = NULL, aes(y = Study_names))+
  geom_linerange(data = SA_AF_prior_Constant_model,  aes(xmin = X2.5., xmax = X97.5.), linewidth = 2, color = "red", alpha = 0.1)+
  geom_pointrange(data = par_sel_sum_data_plot, aes(x = X50., xmin = X2.5., xmax = X97.5.,color = Model), size = 0.25)+
  facet_grid(.~Model+ par, scale = "free", labeller = labeller(par = label_parsed))+
  labs(x = "Parameters in log10 scale", y = "Study")+
  scale_color_manual(values = three_subset_color[c(3,1,2)])+
  theme_bw()+
  theme(text = element_text(size = 14),
        axis.text = element_text(size = 10), axis.text.x = element_text(angle = 90), axis.title.y = element_text(margin = margin(0,20,0,0, "pt")))
selgrob <- ggplotGrob(Constant_model_par_plot) 
locations <- grep("strip-t", selgrob$layout$name);locations
strip <- gtable_filter(selgrob, "strip", trim = FALSE);strip

mat   <- matrix(vector("list", length = 6), nrow = 2)
mat[] <- list(zeroGrob())
# The separator for the facets has zero width
res_t <- gtable_matrix("toprow", mat, unit(c(1, 0, 1), "null"), unit(c(1, 1), "null"))
# Adding the first layer
seladd_strip1 <- res_t %>%
  gtable_add_grob(selgrob$grobs[[22]]$grobs[[1]], 1, 1, 1, 3) %>%
  gtable_add_grob(selgrob, ., t = 7,  l = 7,  b = 5,  r = 5, name = c("add-strip"))

seladd_strip2 <- gtable_add_grob(res_t, selgrob$grobs[[24]]$grobs[[1]], 1, 1, 1, 3) %>%
  gtable_add_grob(seladd_strip1, ., t = 7,  l = 9,  b = 7,  r = 11, name = c("add-strip"))

seladd_strip3 <- gtable_add_grob(res_t, selgrob$grobs[[26]]$grobs[[1]], 1, 1, 1, 3) %>%
  gtable_add_grob(seladd_strip2, ., t = 7,  l = 13,  b = 7,  r = 15, name = c("add-strip"))
ggsave("../plots/Constant_par_plot.png", seladd_strip3, width = 12, height = 6)
#WAIC_comp2[,c(3:5)] = round(WAIC_comp2[,c(3:5)])
#WAIC_comp2[,c(3:5)] %>% colSums()

#Age curve of the global:
Age_depend_global_sum_data = heir_globalcurve_model$Age_depend_global_sum_data
Age_dep_globalcurve_plot = ggplot(Age_depend_global_sum_data, aes(x = age))+
  geom_ribbon(aes(ymin =  X2.5., ymax = X97.5.), alpha = 0.5, fill = three_subset_color[2])+
  geom_ribbon(aes(ymin =  X25., ymax = X75.), alpha = 0.7, fill = three_subset_color[2])+#, color =  "#f03b20")+
  geom_line(aes(y = X50.), color = "black", size = 1)+
  #scale_y_continuous(limits = c(0, 0.03))+
  labs(x = "Age", y = "Exposure risk")+
  theme_bw()+
  theme(text = element_text(size = 20), 
        axis.text.x = element_text(size = 14))
Age_dep_globalcurve_plot
ggsave("../plots/Age_dep_globalcurve_plot.png", Age_dep_globalcurve_plot, width = 7, height = 5)

#Age sum
heir_globalcurve_model$Age_dept_mode
heir_globalcurve_model$Age_dept_ratio_hy
heir_globalcurve_model$Age_dept_ratio_ho

#datafit:
# noage_model_sero_data_fit = noage_model$sero_data_fit
# noage_model_sero_data_fit$STUDY_ID = noage_model_sero_data_fit$study
sero_data_fit_all = bind_rows(noage_model$sero_data_fit, heir_country_model$sero_data_fit, heir_globalcurve_model$sero_data_fit)
sero_data_fit_all_plot = sero_data_fit_all %>%
  cbind(CI = sapply(c(0.025, 0.5, 0.975), function(x) qbeta(x, sero_data_fit_all$POSITIVE + 1, sero_data_fit_all$SAMPLE_SIZE  - sero_data_fit_all$POSITIVE + 1))) %>%
  mutate(AGE_RANGE = paste0(formatC(AGE_LOWER, width = 2, flag = 0), "-", formatC(AGE_UPPER, width = 2, flag = 0)),
         AGE_MID = (AGE_LOWER + AGE_UPPER)/2) %>%
  mutate(mod_STUDY_ID = mapply(rep, source_names[1:4], rle(sero_data_fit_all$STUDY_ID)$lengths) %>% unlist)
sero_data_fit_all_plot$Model = factor(rep(three_models_name, each = nrow(noage_model$sero_data_fit)), levels = three_models_name)

sero_fit_all_plot = ggplot(sero_data_fit_all_plot, aes(x = AGE_MID, group = Model))+
  geom_pointrange(aes(ymin = CI.1, ymax = CI.3, y = CI.2))+
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = Model), alpha = 0.7)+
  geom_line(aes(y = X50.), color = "black")+
  scale_color_manual(values = three_subset_color[c(1,3,2)])+
  scale_fill_manual(values = three_subset_color[c(1,3,2)])+
  labs(x = "Age", y = "Positive proportion", color = "Model", fill = "Model")+
  facet_wrap(.~mod_STUDY_ID, scale = "free_y")+
  theme_bw()+
  theme(text = element_text(size = 18), 
        axis.text.x = element_text(size = 12))
sero_fit_all_plot
ggsave("../plots/Sero_fit_plot_alt_models.png", sero_fit_all_plot, width = 12, height = 6)

#noage_model_case_data_fit = noage_model$case_data_fit
#colnames(noage_model_case_data_fit) = c(colnames(noage_model_case_data_fit)[1:2], "study_id", "data_type", "age_min", "age_max", colnames(noage_model_case_data_fit)[7:12])
case_data_fit_all = bind_rows(noage_model$case_data_fit, heir_country_model$case_data_fit, heir_globalcurve_model$case_data_fit)
case_data_fit_all_data = case_data_fit_all %>%
  mutate(AGE_RANGE = paste0(formatC(age_min, width = 2, flag = 0), "-", formatC(age_max, width = 2, flag = 0)),
         AGE_MID = (age_min + age_max)/2) %>% 
  # mutate(data_subset = case_when(
  #   data_subset == "AF" ~ Subset_order[1],
  #   data_subset == "SA" ~ Subset_order[2],
  #   data_subset == "SA_AF" ~ Subset_order[3]
  # )) %>%
  mutate(mod_STUDY_ID = mapply(rep, c(source_names[5:30]), rle(case_data_fit_all$study_id)$lengths) %>% unlist)
#mutate(mod_STUDY_ID = mapply(rep, c(source_names[-c(13:16)], source_names[c(1:6, 17:30, 7:12)]), rle(case_data_fit$study_id)$lengths) %>% unlist)

#case_data_fit_data$data_subset = factor(case_data_fit_data$data_subset, levels = Subset_order)
case_data_fit_all_data$mod_STUDY_ID = factor(case_data_fit_all_data$mod_STUDY_ID, levels = source_names[5:30])
case_data_fit_all_data$Model = factor(rep(three_models_name, each = nrow( noage_model$case_data_fit)), levels = three_models_name)
case_fit_all_plot = ggplot(case_data_fit_all_data, aes(x = AGE_MID))+
  geom_ribbon(aes(ymin = X2.5., ymax = X97.5., fill = Model), alpha = 0.7)+
  geom_line(aes(y = X50., group = Model), color = "black")+
  geom_point(aes(y = case))+
  scale_color_manual(values = three_subset_color[c(1,3,2)])+
  scale_fill_manual(values = three_subset_color[c(1,3,2)])+
  #scale_y_log10()+
  labs(x = "Age", y = "Cases", color = "Model", fill = "Model")+
  facet_wrap(.~mod_STUDY_ID, scale = "free_y", ncol = 5)+
  theme_bw()+
  theme(text = element_text(size = 18), 
        axis.text.x = element_text(size = 12))
case_fit_all_plot
ggsave("../plots/Case_fit_plot_alt_models.png", case_fit_all_plot, width = 19, height = 10)

#WAIC:
WAIC_comp = noage_model$WAIC_data %>% list %>% rep(3) %>% Reduce(rbind, .)
WAIC_comp$x = round(c(noage_model$WAIC_data$x, heir_country_model$WAIC_data$x, heir_globalcurve_model$WAIC_data$x))
WAIC_comp$model = rep(three_models_name, each = nrow(WAIC_comp)/3)
colnames(WAIC_comp) = c("study_id", "WAIC", "model")

ggplot(WAIC_comp, aes(x = study_id, y = log10(WAIC), fill = model))+geom_bar(stat = "identity", position = position_dodge(width = 1))

WAIC_comp2 = noage_model$WAIC_data ;colnames(WAIC_comp2) = c("study_id", three_models_name[1])
WAIC_comp2$Coun_spec = heir_country_model$WAIC_data$x
WAIC_comp2$Global = heir_globalcurve_model$WAIC_data$x
WAIC_comp2[,c(2:4)] = WAIC_comp2[,c(2:4)] %>% apply(1, function(x) x - min(x)) %>% t %>% round(1)
WAIC_comp2$N_data_point = (c(noage_model$sero_data_fit$STUDY_ID, noage_model$case_data_fit$study_id) %>% rle)$lengths
samples_size_all = aggregate(c(noage_model$sero_data_fit$SAMPLE_SIZE, noage_model$case_data_fit$case), by = list(c(noage_model$sero_data_fit$STUDY_ID, noage_model$case_data_fit$study_id)), sum)
WAIC_comp2$N_samples = samples_size_all[match(unique(c(noage_model$sero_data_fit$STUDY_ID, noage_model$case_data_fit$study_id)), samples_size_all$Group.1),]$x
View(WAIC_comp2)

WAIC_comp2$study_id = factor(source_names, levels = source_names)
all_scenarios_data_plot = data.table::melt(WAIC_comp2, id.vars = "study_id", measure.vars = c("None", "Coun_spec", "Global"))
all_scenarios_data_plot = all_scenarios_data_plot %>% mutate(variable = case_when(variable == "None" ~ three_models_name[1],
                                                        variable == "Coun_spec" ~ three_models_name[2],
                                                        variable == "Global" ~ three_models_name[3]))
all_scenarios_data_plot$variable = factor(all_scenarios_data_plot$variable, levels = three_models_name[c(2,3,1)])
WAIC_plot = (ggplot(filter(all_scenarios_data_plot), aes(x = 1, y = 1, fill = round(value)))+
               geom_tile()+
               geom_text(aes(label = round(value)), size = 7)+
               facet_grid(study_id ~ variable, scale = "free", switch = "y")+
               scale_x_continuous(expand = c(0,0))+
               scale_y_continuous(expand = c(0,0))+
               binned_scale(aesthetics = "fill",
                            scale_name = "stepsn", 
                            palette = function(x) c('#5aae61', '#bfd3e6','#9ebcda','#8c96c6','#8c6bb1','#88419d','#810f7c','#4d004b'),
                            breaks =  c(0, 1, 5, 50, 100, 500, 1000),
                            limits = c(0, 1000),
                            show.limits = TRUE, 
                            right = F,
                            guide = "colorsteps")+
               labs(fill = expression(Delta~WAIC), x = "", y = "")+
               theme_bw()+
               theme(axis.text = element_blank(), text = element_text(size = 20),
                     axis.ticks = element_blank(),
                     panel.spacing = unit(0, "mm"),
                     strip.text.y.left = element_text(angle = 0)))
WAIC_plot
ggsave("../plots/WAIC_plot.png", WAIC_plot, width = 12.5, height = 10)

best_model_id = apply(WAIC_comp2[c(2:4)], 1, function(x) which(x == 0))

dif_WAIC_v_noage_global = noage_model$WAIC_v - heir_globalcurve_model$WAIC_v
dif_WAIC_v_noage_coun = noage_model$WAIC_v - heir_country_model$WAIC_v
dif_WAIC_v_global_coun = heir_globalcurve_model$WAIC_v - heir_country_model$WAIC_v

nrow_all_data = mapply(1:N_study, FUN = rep, each = rle(c(noage_model$sero_data_fit$STUDY_ID, noage_model$case_data_fit$study_id))$lengths) %>% unlist
for(study in 1:N_study){
  best_model_id_study = best_model_id[study]
  if(best_model_id_study == 3){
    nrow_data_study_id = which(nrow_all_data == study)
    CI95_diff_1 = as.numeric(WAIC_comp2[study,2]) + c(-1.96, 1.96)*sqrt(length(nrow_data_study_id)*var(dif_WAIC_v_noage_global[nrow_data_study_id]))
    CI95_diff_2 = as.numeric(WAIC_comp2[study,3]) + c(-1.96, 1.96)*sqrt(length(nrow_data_study_id)*var(-dif_WAIC_v_global_coun[nrow_data_study_id]))
    WAIC_comp2[study,c(2:4)] = c(paste0(round(as.numeric(WAIC_comp2[study,2]), 1), " (", paste0(round(CI95_diff_1, 1), collapse = " - "), ")"),
                            paste0(round(as.numeric(WAIC_comp2[study,3]), 1), " (", paste0(round(CI95_diff_2, 1), collapse = " - "), ")"), 
                            0)
  }
  if(best_model_id_study == 2){
    nrow_data_study_id = which(nrow_all_data == study)
    CI95_diff_1 = as.numeric(WAIC_comp2[study,2]) + c(-1.96, 1.96)*sqrt(length(nrow_data_study_id)*var(dif_WAIC_v_noage_coun[nrow_data_study_id]))
    CI95_diff_2 = as.numeric(WAIC_comp2[study,4]) + c(-1.96, 1.96)*sqrt(length(nrow_data_study_id)*var(dif_WAIC_v_global_coun[nrow_data_study_id]))
    WAIC_comp2[study,c(2:4)] = c(paste0(round(as.numeric(WAIC_comp2[study,2]), 1), " (", paste0(round(CI95_diff_1, 1), collapse = " - "), ")"),0,
                                 paste0(round(as.numeric(WAIC_comp2[study,4]), 1), " (", paste0(round(CI95_diff_2, 1), collapse = " - "), ")"))
  }
  if(best_model_id_study == 1){
    nrow_data_study_id = which(nrow_all_data == study)
    CI95_diff_1 = as.numeric(WAIC_comp2[study,3]) + c(-1.96, 1.96)*sqrt(length(nrow_data_study_id)*var(-dif_WAIC_v_noage_coun[nrow_data_study_id]))
    CI95_diff_2 = as.numeric(WAIC_comp2[study,4]) + c(-1.96, 1.96)*sqrt(length(nrow_data_study_id)*var(-dif_WAIC_v_noage_global[nrow_data_study_id]))
    WAIC_comp2[study,c(2:4)] = c(paste0(round(as.numeric(WAIC_comp2[study,3]), 1), " (", paste0(round(CI95_diff_1, 1), collapse = " - "), ")"),0,
                                 paste0(round(as.numeric(WAIC_comp2[study,4]), 1), " (", paste0(round(CI95_diff_2, 1), collapse = " - "), ")"))
  }
}
#all data:
c(sum(dif_WAIC_v_noage_coun), sum(dif_WAIC_v_noage_coun) + c(-1.96, 1.96)*sqrt(length(dif_WAIC_v_noage_coun)*var(dif_WAIC_v_noage_coun)))
c(sum(dif_WAIC_v_global_coun), sum(dif_WAIC_v_global_coun) + c(-1.96, 1.96)*sqrt(length(dif_WAIC_v_noage_coun)*var(dif_WAIC_v_global_coun)))


