############################################## Baboon example
library(STRAND)
library(stringr)
library(ggplot2)
library(psych)

# install_github('ctross/PlvsVltra')
 library(PlvsVltra) # For colors


############################################## 
# Load data
 data(Baboon_Longitudinal_Data)
 d = Baboon_Longitudinal_Data

############################################## 
# Process data into longitudinal form
dat_long = NULL
# Loop over y days of data, making a day-specific network data-set
# Some days are missing outcome data. This is dealt with via the mask layer
# if outcome data is missing, mask[i,j]=1
# currently, missings in the predictors aren't supported in STRAND, but will be eventually
for(y in 1:14){
 # Merge data
 dat_long[[y]] = make_strand_data(
  outcome = list("Affiliative" = d$Affiliative[[y]]),
  exposure = list("Affiliative" = d$Exposure[[y]]),
  mask = list("Affiliative" = d$Mask[[y]]),
  block_covariates = NULL, 
  individual_covariates = d$Individual, 
  dyadic_covariates = list("Presenting" = t(d$Presenting[[y]])),
  longitudinal = TRUE,
  outcome_mode="binomial",
  link_mode="logit"
  )
 }
names(dat_long) = paste("Time", c(1:14))

############################################## 
# Fit model with time-varying slopes, cholesky style
fit_2a = fit_longitudinal_model(
 long_data=dat_long,
 block_regression = ~ 1,
 focal_regression = ~ Age + Sex,
 target_regression = ~ Age + Sex,
 dyad_regression = ~ Presenting,
 coefficient_mode="varying",
 random_effects_mode="fixed",
 bandage_penalty = -1,
  mode="mcmc",
 stan_mcmc_parameters = list(seed = 1, chains = 1, init=0,
    parallel_chains = 1, refresh = 10, iter_warmup = 1000,
    iter_sampling = 1000, max_treedepth = 12)
  )

res_2a = summarize_longitudinal_bsrm_results(fit_2a)

############################################## 
# Fit model with time-varying slopes, l2norm style
fit_2b = fit_longitudinal_model(
 long_data=dat_long,
 block_regression = ~ 1,
 focal_regression = ~ Age + Sex,
 target_regression = ~ Age + Sex,
 dyad_regression = ~ Presenting,
 coefficient_mode="varying",
 random_effects_mode="fixed",
 bandage_penalty = 0.01,
  mode="mcmc",
 stan_mcmc_parameters = list(seed = 1, chains = 1, init=0,
    parallel_chains = 1, refresh = 10, iter_warmup = 1000,
    iter_sampling = 1000, max_treedepth = 12)
  )

res_2b = summarize_longitudinal_bsrm_results(fit_2b)

fit_2a$fit$time()$total
fit_2b$fit$time()$total

############################################## 
# Fit model with time-invariant slopes, cholesky style
fit_2c = fit_longitudinal_model(
 long_data=dat_long,
 block_regression = ~ 1,
 focal_regression = ~ Age + Sex,
 target_regression = ~ Age + Sex,
 dyad_regression = ~ Presenting,
 coefficient_mode="fixed",
 random_effects_mode="fixed",
 bandage_penalty = -1,
  mode="mcmc",
 stan_mcmc_parameters = list(seed = 1, chains = 1, init=0,
    parallel_chains = 1, refresh = 10, iter_warmup = 1000,
    iter_sampling = 1000, max_treedepth = 12),
 priors=NULL
  )

res_2c = summarize_longitudinal_bsrm_results(fit_2c)

############################################## 
# Fit model with time-invariant slopes, l2norm style
fit_2d = fit_longitudinal_model(
 long_data=dat_long,
 block_regression = ~ 1,
 focal_regression = ~ Age + Sex,
 target_regression = ~ Age + Sex,
 dyad_regression = ~ Presenting,
 coefficient_mode="fixed",
 random_effects_mode="fixed",
 bandage_penalty = 0.01,
  mode="mcmc",
 stan_mcmc_parameters = list(seed = 1, chains = 1, init=0,
    parallel_chains = 1, refresh = 10, iter_warmup = 1000,
    iter_sampling = 1000, max_treedepth = 12),
 priors=NULL
  )

res_2d = summarize_longitudinal_bsrm_results(fit_2d)

############################################## 
# Visualize results

# Correlation matrix plots
colors = plvs_vltra("dust_storm", rev=FALSE, elements=c(2,4))
colors = c(colors[1], "grey90", colors[2])
multiplex_plot(fit_2a, type="dyadic", HPDI=0.9, plot = TRUE, export_as_table = FALSE, save_plot = "Baboon_dyadic_cholesky.pdf", height=6, width=7, palette=colors)
multiplex_plot(fit_2a, type="generalized", HPDI=0.9, plot = TRUE, export_as_table = FALSE, save_plot = "Baboon_generalized_cholesky.pdf", height=6, width=7, palette=colors)

multiplex_plot(fit_2b, type="dyadic", HPDI=0.9, plot = TRUE, export_as_table = FALSE, save_plot = "Baboon_dyadic_l2norm.pdf", height=6, width=7, palette=colors)
multiplex_plot(fit_2b, type="generalized", HPDI=0.9, plot = TRUE, export_as_table = FALSE, save_plot = "Baboon_generalized_l2norm.pdf", height=6, width=7, palette=colors)

# Long plots
pal = plvs_vltra("mystic_mausoleum", elements=c(1,9,2))
longitudinal_plot(fit_2a, type="dyadic", save_plot="Baboon_dyadic_long_cholesky.pdf", palette = pal, height=6, width=6.5)
longitudinal_plot(fit_2b, type="dyadic", save_plot="Baboon_dyadic_long_l2norm.pdf", palette = pal, height=6, width=6.5)

pal = plvs_vltra("mystic_mausoleum", elements=c(1,2,9,3,4))
longitudinal_plot(fit_2a, type="generalized", save_plot="Baboon_generalized_long_cholesky.pdf", palette = pal, height=6, width=6.5)
longitudinal_plot(fit_2b, type="generalized", save_plot="Baboon_generalized_long_l2norm.pdf", palette = pal, height=6, width=6.5)

pal = plvs_vltra("mystic_mausoleum", elements=c(1,9,7,5,10,8,6))
longitudinal_plot(fit_2a,type="coefficient", 
    parameter_set = list(
    focal="Age", target="Age", 
    focal="SexMale", target="SexMale",
    dyadic="Presenting"),
    palette=pal,
    normalized=TRUE,
    height=4, width=9,
    save_plot="Slopes_Baboon_long_cholesky.pdf")

longitudinal_plot(fit_2b,type="coefficient", 
    parameter_set = list(
    focal="Age", target="Age", 
    focal="SexMale", target="SexMale",
    dyadic="Presenting"),
    palette=pal,
    normalized=TRUE,
    height=4, width=9,
    save_plot="Slopes_Baboon_long_l2norm.pdf")

############################################## 
# Visualize results
# Make a merged plot

################# Cholesky-style
bab2c_set = c("focal effects coeffs (out-degree), Time 1 - Age", 
"focal effects coeffs (out-degree), Time 1 - SexMale", 
"target effects coeffs (in-degree), Time 1 - Age", 
"target effects coeffs (in-degree), Time 1 - SexMale",
"dyadic effects coeffs, Time 1 - Presenting")

to_add = res_2c$summary[which(res_2c$summary$Variable %in% bab2c_set),]

    parameter_set = list(
    focal="Age", target="Age", 
    focal="SexMale", target="SexMale",
    dyadic="Presenting")

base_set = longitudinal_plot_c(fit_2a, parameters=as.vector(unlist(parameter_set)), type=names(parameter_set),
    plot=FALSE,
    normalized=FALSE,
    export_as_table=TRUE)

to_add$short_names = c("Focal - Time 1 - Age", "Focal - Time 1 - SexMale",
                "Target - Time 1 - Age", "Target - Time 1 - SexMale",
                 "Dyadic - Time 1 - Presenting")
to_add$time_point = rep("Time 0", 5)
to_add$time_point_int = rep(0, 5)
to_add$extra_short_names = c("Focal - Age", "Focal - SexMale",
                "Target - Age", "Target - SexMale",
                 "Dyadic - Presenting")

to_add$type_set = c("Focal", "Focal", "Target", "Target", "Dyadic")

to_add$Model = "Fixed"
base_set$Model = "Varying"

full_set = rbind(base_set,to_add)


full_set2 = full_set[,c("Variable", "time_point", "extra_short_names", "extra_short_names", "Median", "HPDI:0.05", "HPDI:0.95", "Mean", "SD", "time_point_int", "type_set")]
colnames(full_set2) = c("Variable", "Layer", "Target" , "Base" , "Median", "L", "H", "Mean","SD","LayerNumeric","Type")
      
     Diff = as.numeric(full_set2$H)-as.numeric(full_set2$L)   
     full_set2$Median = as.numeric(full_set2$Median)/Diff
     full_set2$L = as.numeric(full_set2$L)/Diff
     full_set2$H =  as.numeric(full_set2$H)/Diff

     full_set2$Model = ifelse(full_set2$"LayerNumeric"==0, "Fixed", "Varying")
     full_set2$Model2 = ifelse(full_set2$"Model"=="Fixed", "Fixed", full_set2$Type)

     full_set2$Model2 = factor(full_set2$Model2)
     full_set2$Model2 = factor(full_set2$Model2, levels=c("Fixed", "Dyadic", "Focal", "Target"))

     full_set2_choleksy = full_set2
    full_set2_choleksy$method = "choleksy"

################# l2norm-style
bab2d_set = c("focal effects coeffs (out-degree), Time 1 - Age", 
"focal effects coeffs (out-degree), Time 1 - SexMale", 
"target effects coeffs (in-degree), Time 1 - Age", 
"target effects coeffs (in-degree), Time 1 - SexMale",
"dyadic effects coeffs, Time 1 - Presenting")

to_add = res_2d$summary[which(res_2d$summary$Variable %in% bab2d_set),]

    parameter_set = list(
    focal="Age", target="Age", 
    focal="SexMale", target="SexMale",
    dyadic="Presenting")

base_set = longitudinal_plot_c(fit_2b, parameters=as.vector(unlist(parameter_set)), type=names(parameter_set),
    plot=FALSE,
    normalized=FALSE,
    export_as_table=TRUE)

to_add$short_names = c("Focal - Time 1 - Age", "Focal - Time 1 - SexMale",
                "Target - Time 1 - Age", "Target - Time 1 - SexMale",
                 "Dyadic - Time 1 - Presenting")
to_add$time_point = rep("Time 0", 5)
to_add$time_point_int = rep(0, 5)
to_add$extra_short_names = c("Focal - Age", "Focal - SexMale",
                "Target - Age", "Target - SexMale",
                 "Dyadic - Presenting")

to_add$type_set = c("Focal", "Focal", "Target", "Target", "Dyadic")

to_add$Model = "Fixed"
base_set$Model = "Varying"

full_set = rbind(base_set,to_add)


full_set2 = full_set[,c("Variable", "time_point", "extra_short_names", "extra_short_names", "Median", "HPDI:0.05", "HPDI:0.95", "Mean", "SD", "time_point_int", "type_set")]
colnames(full_set2) = c("Variable", "Layer", "Target" , "Base" , "Median", "L", "H", "Mean","SD","LayerNumeric","Type")
      
     Diff = as.numeric(full_set2$H)-as.numeric(full_set2$L)   
     full_set2$Median = as.numeric(full_set2$Median)/Diff
     full_set2$L = as.numeric(full_set2$L)/Diff
     full_set2$H =  as.numeric(full_set2$H)/Diff

     full_set2$Model = ifelse(full_set2$"LayerNumeric"==0, "Fixed", "Varying")
     full_set2$Model2 = ifelse(full_set2$"Model"=="Fixed", "Fixed", full_set2$Type)

     full_set2$Model2 = factor(full_set2$Model2)
     full_set2$Model2 = factor(full_set2$Model2, levels=c("Fixed", "Dyadic", "Focal", "Target"))

     full_set2_l2norm = full_set2
     full_set2_l2norm$method = "l2norm"


     full_set_all = rbind(full_set2_choleksy, full_set2_l2norm)

p = ggplot(full_set_all, aes(x=LayerNumeric, y=as.numeric(Median), ymin=as.numeric(L), ymax=as.numeric(H), group=Target, color=Target))+ 
     geom_linerange(size=1, position = position_dodge(width = 0.3)) + facet_grid(method~Model2, scales="free", space="free") + 
     geom_point(size=2, position = position_dodge(width = 0.3))+
     geom_hline(aes(yintercept=0),color="black",linetype="dashed")+
     labs(y="Effect size", x="Time step") + 
     theme(strip.text.x = element_text(size=12,face="bold"), 
      strip.text.y = element_text(size=12,face="bold"),
      axis.text = element_text(size=12),
      axis.title = element_text(size=14, face="bold"))+
     theme(strip.text.y = element_text(angle = 360)) + 
    # coord_flip() + 
     theme(panel.spacing = grid::unit(1, "lines")) + scale_color_manual(values = pal) + 
     theme(legend.position="bottom") + theme(legend.title = element_blank()) + scale_x_continuous(breaks=1:14,expand = c(0, 0.95))

p

ggsave("Slopes_Baboon_merged.pdf", p, height=6, width=13.5)

##################
# Fit model with time-varying slopes, cholesky style
fit_2a$fit$time()$total
mean(as.data.frame(fit_2a$fit$summary(variables ="D_corr"))$ess_bulk, na.rm=TRUE)
mean(as.data.frame(fit_2a$fit$summary(variables ="D_corr"))$ess_tail, na.rm=TRUE)

mean(as.data.frame(fit_2a$fit$summary(variables ="D_corr"))$ess_bulk, na.rm=TRUE)/fit_2a$fit$time()$total
mean(as.data.frame(fit_2a$fit$summary(variables ="D_corr"))$ess_tail, na.rm=TRUE)/fit_2a$fit$time()$total


# Fit model with time-varying slopes, l2norm style
fit_2b$fit$time()$total
mean(as.data.frame(fit_2b$fit$summary(variables ="D_corr"))$ess_bulk, na.rm=TRUE)
mean(as.data.frame(fit_2b$fit$summary(variables ="D_corr"))$ess_tail, na.rm=TRUE)

mean(as.data.frame(fit_2b$fit$summary(variables ="D_corr"))$ess_bulk, na.rm=TRUE)/fit_2b$fit$time()$total
mean(as.data.frame(fit_2b$fit$summary(variables ="D_corr"))$ess_tail, na.rm=TRUE)/fit_2b$fit$time()$total


# Fit model with fixed slopes, cholesky style
fit_2c$fit$time()$total
mean(as.data.frame(fit_2c$fit$summary(variables ="D_corr"))$ess_bulk, na.rm=TRUE)
mean(as.data.frame(fit_2c$fit$summary(variables ="D_corr"))$ess_tail, na.rm=TRUE)

mean(as.data.frame(fit_2c$fit$summary(variables ="D_corr"))$ess_bulk, na.rm=TRUE)/fit_2c$fit$time()$total
mean(as.data.frame(fit_2c$fit$summary(variables ="D_corr"))$ess_tail, na.rm=TRUE)/fit_2c$fit$time()$total


# Fit model with fixed slopes, l2norm style
fit_2d$fit$time()$total
mean(as.data.frame(fit_2d$fit$summary(variables ="D_corr"))$ess_bulk, na.rm=TRUE)
mean(as.data.frame(fit_2d$fit$summary(variables ="D_corr"))$ess_tail, na.rm=TRUE)

mean(as.data.frame(fit_2d$fit$summary(variables ="D_corr"))$ess_bulk, na.rm=TRUE)/fit_2d$fit$time()$total
mean(as.data.frame(fit_2d$fit$summary(variables ="D_corr"))$ess_tail, na.rm=TRUE)/fit_2d$fit$time()$total





