########################################
#
#   Multiplex Bernoulli Analyses - Full  
#
########################################
# Load data
data(RICH_Data)
RICH = RICH_Data

# Structure for sample-size sweep
dat_list_cholesky = NULL
dat_list_l2norm = NULL
dat_list_direct = NULL

samp_seq = seq(23, 93, by=5)
ticker = 0

for(s in samp_seq){
 ticker = 1 + ticker
 Q = s

# Outcomes stored as a labeled list
outcome = list(
 Give = RICH$Give[1:Q,1:Q], 
 Take = RICH$Take[1:Q,1:Q], 
 Reduce = RICH$Reduce[1:Q,1:Q]
)

# Dyadic data as a labeled list
dyad = list(
 Relatedness = RICH$Relatedness[1:Q,1:Q], 
 Friends = RICH$Friends[1:Q,1:Q],
 Marriage = RICH$Marriage[1:Q,1:Q]
)

# Individual data in data-frame
ind = RICH$Individual[1:Q,]

# Individual blocking measures
groups = data.frame(
 Ethnicity = as.factor(ind$Ethnicity), 
 Sex = as.factor(ind$Sex)
)
rownames(groups) = rownames(ind)

# Merge data
dat = make_strand_data(
 outcome = outcome,
 block_covariates = groups, 
 individual_covariates = ind, 
 dyadic_covariates = dyad,
 outcome_mode="bernoulli",
 link_mode="logit",
 multiplex = TRUE
)

dat$node_count = s

dat$model_to_deploy = 1
dat_list_cholesky[[ticker]] = dat

dat$model_to_deploy = 2
dat_list_l2norm[[ticker]] = dat

dat$model_to_deploy = 3
dat_list_direct[[ticker]] = dat
}

dat_list = c(dat_list_cholesky, dat_list_l2norm, dat_list_direct)

####################################################################################### Coeff plots
deploy_on_cores2 = function(dat){
if(dat$model_to_deploy == 1){
############################################### Fast mode
# Model 1 - Full model, all controls
fit = fit_multiplex_model(
 data=dat,
 block_regression = ~ Ethnicity + Sex,
 focal_regression = ~ Age + Wealth + FoodInsecure + Depressed,
 target_regression = ~ Age + Wealth + FoodInsecure + Depressed,
 dyad_regression = ~ Relatedness + Friends + Marriage,
 mode="mcmc",
 bandage_penalty = -1,         # This sets Cholesky style
 stan_mcmc_parameters = list(
   chains = 1, 
   parallel_chains = 1, 
   refresh = 1, 
   iter_warmup = 1000, 
   iter_sampling = 1000, 
   max_treedepth = 12, 
   adapt_delta = 0.98,
   init=0)
)

}

if(dat$model_to_deploy == 2){
############################################### Lasso mode
# Model 2 - Full model, all controls
fit = fit_multiplex_model(
 data=dat,
 block_regression = ~ Ethnicity + Sex,
 focal_regression = ~ Age + Wealth + FoodInsecure + Depressed,
 target_regression = ~ Age + Wealth + FoodInsecure + Depressed,
 dyad_regression = ~ Relatedness + Friends + Marriage,
 mode="mcmc",
 bandage_penalty = 0.01,       # This sets l2norm style
 stan_mcmc_parameters = list(
   chains = 1, 
   parallel_chains = 1, 
   refresh = 1, 
   iter_warmup = 1000, 
   iter_sampling = 1000, 
   max_treedepth = 12, 
   adapt_delta = 0.98,
   init=0)
)

}

if(dat$model_to_deploy == 3){
############################################### Lasso mode
# Model 2 - Full model, all controls
fit = fit_multiplex_model(
 data=dat,
 block_regression = ~ Ethnicity + Sex,
 focal_regression = ~ Age + Wealth + FoodInsecure + Depressed,
 target_regression = ~ Age + Wealth + FoodInsecure + Depressed,
 dyad_regression = ~ Relatedness + Friends + Marriage,
 mode="mcmc",
 bandage_penalty = 0,       # This sets direct style
 stan_mcmc_parameters = list(
   chains = 1, 
   parallel_chains = 1, 
   refresh = 1, 
   iter_warmup = 1000, 
   iter_sampling = 1000, 
   max_treedepth = 12, 
   adapt_delta = 0.98,
   init=0)
)

}

#################################################### Return
 return(fit)
}

############################################################### Fit on server
fit_full = mclapply(dat_list[c(15,30,45)], function(z){
                       deploy_on_cores2(z)
                       }, mc.cores = 3)

########################################################### Covariate effects
res_1 = summarize_strand_results(fit_full[[1]])
res_2 = summarize_strand_results(fit_full[[2]])
res_3 = summarize_strand_results(fit_full[[3]])

vis_1 = strand_caterpillar_plot(res_1, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects"), export_as_table = TRUE, normalized=FALSE)
vis_2 = strand_caterpillar_plot(res_2, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects"), export_as_table = TRUE, normalized=FALSE)
vis_3 = strand_caterpillar_plot(res_3, submodels=c("Focal effects: Out-degree","Target effects: In-degree","Dyadic effects"), export_as_table = TRUE, normalized=FALSE)

vis_1$Site = "cholesky"
vis_2$Site = "l2norm"
vis_3$Site = "direct"

vis_1$Variable = gsub("focal effects coeffs \\(out-degree\\), ", "", vis_1$Variable)
vis_1$Variable = gsub("target effects coeffs \\(in-degree\\), ", "", vis_1$Variable)
vis_1$Variable = gsub("dyadic effects coeffs, ", "", vis_1$Variable)
vis_1$Variable = gsub("dyadic effects ", "", vis_1$Variable)
vis_1$Variable = gsub("focal effects ", "", vis_1$Variable)
vis_1$Variable = gsub("target effects ", "", vis_1$Variable)


vis_2$Variable = gsub("focal effects coeffs \\(out-degree\\), ", "", vis_2$Variable)
vis_2$Variable = gsub("target effects coeffs \\(in-degree\\), ", "", vis_2$Variable)
vis_2$Variable = gsub("dyadic effects coeffs, ", "", vis_2$Variable)
vis_2$Variable = gsub("dyadic effects ", "", vis_2$Variable)
vis_2$Variable = gsub("focal effects ", "", vis_2$Variable)
vis_2$Variable = gsub("target effects ", "", vis_2$Variable)

vis_3$Variable = gsub("focal effects coeffs \\(out-degree\\), ", "", vis_3$Variable)
vis_3$Variable = gsub("target effects coeffs \\(in-degree\\), ", "", vis_3$Variable)
vis_3$Variable = gsub("dyadic effects coeffs, ", "", vis_3$Variable)
vis_3$Variable = gsub("dyadic effects ", "", vis_3$Variable)
vis_3$Variable = gsub("focal effects ", "", vis_3$Variable)
vis_3$Variable = gsub("target effects ", "", vis_3$Variable)

df = rbind(vis_1, vis_2, vis_3)

df$Outcome = ifelse(str_detect(df$Variable, "Take"), "Take",
             ifelse(str_detect(df$Variable, "Give"), "Give",
             ifelse(str_detect(df$Variable, "Reduce"), "Reduce",
                    NA)))

df$Outcome = factor(df$Outcome)
df$Outcome = factor(df$Outcome, levels=c("Give", "Take", "Reduce"))

df$Variable = gsub("Give - ", "", df$Variable)
df$Variable = gsub("Take - ", "", df$Variable)
df$Variable = gsub("Reduce - ", "", df$Variable)

df$Variable = gsub(" - Give", "", df$Variable)
df$Variable = gsub(" - Take", "", df$Variable)
df$Variable = gsub(" - Reduce", "", df$Variable)

df$Variable = gsub("sd", "SD", df$Variable)

df$Variable = gsub("FoodInsecure", "Food Insecure", df$Variable)
df$Variable = gsub("LogWealth", "Log Wealth", df$Variable)

df$Variable = factor(df$Variable)
df$Variable = factor(df$Variable, levels=c("SD", "Age", "Wealth", "Food Insecure", "Depressed", "Relatedness", "Friends", "Marriage"))

df$Method = df$Site

p = ggplot2::ggplot(df, ggplot2::aes(x = Variable, y = Median, group = Method, color=Method,
        ymin = LI, ymax = HI)) + ggplot2::geom_linerange(size = 1,, position = position_dodge(width = 0.7)) + 
        ggplot2::geom_point(size = 2,, position = position_dodge(width = 0.7)) + ggplot2::facet_grid(Submodel ~ 
        Outcome, scales = "free", space = "free") + ggplot2::geom_hline(ggplot2::aes(yintercept = 0), 
        color = "black", linetype = "dashed") + ggplot2::labs(y = "Regression parameters", 
        x = "") + ggplot2::theme(strip.text.x = ggplot2::element_text(size = 12, 
        face = "bold"), strip.text.y = ggplot2::element_text(size = 12, 
        face = "bold"), axis.text = ggplot2::element_text(size = 12), 
        axis.title.y = ggplot2::element_text(size = 14, face = "bold"), 
        axis.title.x = ggplot2::element_blank()) + ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 360)) + 
        ggplot2::coord_flip() + ggplot2::theme(panel.spacing = grid::unit(1, 
        "lines")) + scale_color_manual(values=c("l2norm" = colors2[5], "cholesky" = colors2[7], "direct" = colors2[3])) + theme(legend.position="bottom")

 ggsave("rich_res.pdf",p, width=9, height=5.5)


########################################################### VPCs
VPCs_1 = strand_VPCs(fit_full[[1]], n_partitions = 4)
VPCs_2 = strand_VPCs(fit_full[[2]], n_partitions = 4)
VPCs_3 = strand_VPCs(fit_full[[3]], n_partitions = 4)

df1 = data.frame(do.call(rbind, VPCs_1[[2]]))
colnames(df1) = c("Variable", "Median", "L", "H", "Mean", "SD")
df1$Site = "cholesky"
df1$Submodel = rep(c("Give","Take","Reduce"),each=4)
df1$Variable2 = rep(c("Focal","Target","Dyadic","Error"),3)

df2 = data.frame(do.call(rbind, VPCs_2[[2]]))
colnames(df2) = c("Variable", "Median", "L", "H", "Mean", "SD")
df2$Site = "l2norm"
df2$Submodel = rep(c("Give","Take","Reduce"),each=4)
df2$Variable2 = rep(c("Focal","Target","Dyadic","Error"),3)

df3 = data.frame(do.call(rbind, VPCs_3[[2]]))
colnames(df3) = c("Variable", "Median", "L", "H", "Mean", "SD")
df3$Site = "direct"
df3$Submodel = rep(c("Give","Take","Reduce"),each=4)
df3$Variable2 = rep(c("Focal","Target","Dyadic","Error"),3)

df = rbind(df1, df2, df3)
df$Median = as.numeric(df$Median)
df$L = as.numeric(df$L)
df$H = as.numeric(df$H)

df$Submodel = factor(df$Submodel)
df$Submodel = factor(df$Submodel, levels=c("Give", "Take", "Reduce"))

df$Variable2 = factor(df$Variable2)
df$Variable2 = factor(df$Variable2, levels=rev(c("Focal","Target","Dyadic","Error")))

df$Method = df$Site

p = ggplot2::ggplot(df, ggplot2::aes(x = Variable2, y = Median, group = Method, color=Method,
        ymin = L, ymax = H)) + ggplot2::geom_linerange(size = 1,, position = position_dodge(width = 0.4)) + 
        ggplot2::geom_point(size = 2,, position = position_dodge(width = 0.4)) + ggplot2::facet_grid(Submodel ~ ., scales = "free", space = "free") +
         ggplot2::geom_hline(ggplot2::aes(yintercept = 0), 
        color = "black", linetype = "dashed") + ggplot2::labs(y = "Regression parameters", 
        x = "") + ggplot2::theme(strip.text.x = ggplot2::element_text(size = 12, 
        face = "bold"), strip.text.y = ggplot2::element_text(size = 12, 
        face = "bold"), axis.text = ggplot2::element_text(size = 12), 
        axis.title.y = ggplot2::element_text(size = 14, face = "bold"), 
        axis.title.x = ggplot2::element_blank()) + ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 360)) + 
        ggplot2::coord_flip() + ggplot2::theme(panel.spacing = grid::unit(1, 
        "lines")) + scale_color_manual(values=c("l2norm" = colors2[5], "cholesky" = colors2[7], "direct" = colors2[3])) + theme(legend.position="bottom")

 ggsave("rich_vpc.pdf",p, width=4.5, height=9)


########################################################### Reciprocity
VPCs_1 = strand_VPCs(fit_full[[1]], n_partitions = 4, include_reciprocity = TRUE)
VPCs_2 = strand_VPCs(fit_full[[2]], n_partitions = 4, include_reciprocity = TRUE)
VPCs_3 = strand_VPCs(fit_full[[3]], n_partitions = 4, include_reciprocity = TRUE)

df1 = data.frame(VPCs_1[[3]])
colnames(df1) = c("Variable", "Median", "L", "H", "Mean", "SD")
df1$Site = "cholesky"
df1$Submodel = rep(c("Generalized","Dyadic"),each=15)

df2 = data.frame(VPCs_2[[3]])
colnames(df2) = c("Variable", "Median", "L", "H", "Mean", "SD")
df2$Site = "l2norm"
df2$Submodel = rep(c("Generalized","Dyadic"),each=15)

df3 = data.frame(VPCs_3[[3]])
colnames(df3) = c("Variable", "Median", "L", "H", "Mean", "SD")
df3$Site = "direct"
df3$Submodel = rep(c("Generalized","Dyadic"),each=15)

df = rbind(df1, df2, df3)
df$Median = as.numeric(df$Median)
df$L = as.numeric(df$L)
df$H = as.numeric(df$H)

df$Submodel = factor(df$Submodel)

df$Method = df$Site

p = ggplot2::ggplot(df, ggplot2::aes(x = Variable, y = Median, group = Method, color=Method,
        ymin = L, ymax = H)) + ggplot2::geom_linerange(size = 1,, position = position_dodge(width = 0.6)) + 
        ggplot2::geom_point(size = 2,, position = position_dodge(width = 0.6)) + ggplot2::facet_grid(. ~Submodel, scales = "free", space = "free") +
         ggplot2::geom_hline(ggplot2::aes(yintercept = 0), 
        color = "black", linetype = "dashed") + ggplot2::labs(y = "Regression parameters", 
        x = "") + ggplot2::theme(strip.text.x = ggplot2::element_text(size = 12, 
        face = "bold"), strip.text.y = ggplot2::element_text(size = 12, 
        face = "bold"), axis.text = ggplot2::element_text(size = 12), 
        axis.title.y = ggplot2::element_text(size = 14, face = "bold"), 
        axis.title.x = ggplot2::element_blank()) + ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 360)) + 
        ggplot2::coord_flip() + ggplot2::theme(panel.spacing = grid::unit(1, 
        "lines")) + scale_color_manual(values=c("l2norm" = colors2[5], "cholesky" = colors2[7], "direct" = colors2[3])) + theme(legend.position="bottom")


 ggsave("rich_recip.pdf",p, width=10, height=9)


################ Details for text
# Cholesky
 time_c = fit_full[[1]]$fit$time()$total   # 22279.08
 bulk_c = mean(as.data.frame(fit_full[[1]]$fit$summary(variables ="D_corr"))$ess_bulk,na.rm=TRUE) # 86.50093
 tail_c = mean(as.data.frame(fit_full[[1]]$fit$summary(variables ="D_corr"))$ess_tail,na.rm=TRUE) # 223.436

# L2norm
 time_l = fit_full[[2]]$fit$time()$total   # 17502.75
 bulk_l = mean(as.data.frame(fit_full[[2]]$fit$summary(variables ="D_corr"))$ess_bulk,na.rm=TRUE) # 65.93305
 tail_l = mean(as.data.frame(fit_full[[2]]$fit$summary(variables ="D_corr"))$ess_tail,na.rm=TRUE) # 154.852

# Direct
 time_d = fit_full[[3]]$fit$time()$total   # 23744.78
 bulk_d = mean(as.data.frame(fit_full[[3]]$fit$summary(variables ="D_corr"))$ess_bulk,na.rm=TRUE) # 42.67249
 tail_d = mean(as.data.frame(fit_full[[3]]$fit$summary(variables ="D_corr"))$ess_tail,na.rm=TRUE) # 127.374

eff_bulk_c = bulk_c/time_c
eff_tail_c = tail_c/time_c

eff_bulk_l = bulk_l/time_l
eff_tail_l = tail_l/time_l

eff_bulk_d = bulk_d/time_d
eff_tail_d = tail_d/time_d


eff_bulk_c/eff_bulk_d # 2.160448
eff_bulk_l/eff_bulk_d # 2.096125

eff_tail_c/eff_tail_d # 1.869577
eff_tail_l/eff_tail_d # 1.649293

