##############################################
#
#   Multiplex Bernoulli Analyses - Sweep 
#
#############################################

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

##################################### Now run the simulations
deploy_on_cores = function(dat){
if(dat$model_to_deploy == 1){
############################################### Fast mode
# Model 1 - Full model, all controls
fit_Fast = fit_multiplex_model(
 data=dat,
 block_regression = ~ Ethnicity + Sex,
 focal_regression = ~ Age + Wealth + FoodInsecure + Depressed,
 target_regression = ~ Age + Wealth + FoodInsecure + Depressed,
 dyad_regression = ~ Relatedness + Friends + Marriage,
 mode="mcmc",
 bandage_penalty = -1,
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

res_temp = as.data.frame(fit_Fast$fit$summary(variables ="D_corr"))
res_temp$fit_time = fit_Fast$fit$time()$total
res_temp$node_count = dat$node_count
res_temp$method = "cholesky"
}

if(dat$model_to_deploy == 2){
############################################### Lasso mode
# Model 2 - Full model, all controls
fit_Lasso = fit_multiplex_model(
 data=dat,
 block_regression = ~ Ethnicity + Sex,
 focal_regression = ~ Age + Wealth + FoodInsecure + Depressed,
 target_regression = ~ Age + Wealth + FoodInsecure + Depressed,
 dyad_regression = ~ Relatedness + Friends + Marriage,
 mode="mcmc",
 bandage_penalty = 0.01,
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

res_temp = as.data.frame(fit_Lasso$fit$summary(variables ="D_corr"))
res_temp$fit_time = fit_Lasso$fit$time()$total
res_temp$node_count = dat$node_count
res_temp$method = "l2norm"
}

if(dat$model_to_deploy == 3){
############################################### Direct mode
# Model 2 - Full model, all controls
fit_Direct = fit_multiplex_model(
 data=dat,
 block_regression = ~ Ethnicity + Sex,
 focal_regression = ~ Age + Wealth + FoodInsecure + Depressed,
 target_regression = ~ Age + Wealth + FoodInsecure + Depressed,
 dyad_regression = ~ Relatedness + Friends + Marriage,
 mode="mcmc",
 bandage_penalty = 0,
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

res_temp = as.data.frame(fit_Direct$fit$summary(variables ="D_corr"))
res_temp$fit_time = fit_Direct$fit$time()$total
res_temp$node_count = dat$node_count
res_temp$method = "direct"
}

#################################################### Return
 return(res_temp)
}

############################################################### Fit on server
fit = mclapply(dat_list, function(z){
                       deploy_on_cores(z)
                       }, mc.cores = 45)

res = do.call(rbind, fit)
save(res, file = "res_compare.RData")

################################################################### Process
sum_stats = function(x){return(c(mean(x), HPDI(x, prob=0.9)))}

res_thin = res[which(res$mean != 1),]
res_thin$group = paste0(res_thin$method,"_", res_thin$node_count)
res_thin$method_numeric = ifelse(res_thin$method=="l2norm",1,ifelse(res_thin$method=="cholesky",2,3))

l = list( 
 aggregate(res_thin$method_numeric~res_thin$group, FUN=median),
 aggregate(res_thin$node_count~res_thin$group, FUN=median),
 aggregate(res_thin$fit_time~res_thin$group, FUN=median),
 aggregate(res_thin$ess_bulk~res_thin$group, FUN=sum_stats),
 aggregate(res_thin$ess_tail~res_thin$group, FUN=sum_stats))

res_merged = purrr::reduce(.x = l, merge, by = c('res_thin$group'), all = TRUE)
res_merged = data.frame(as.matrix(data.frame(res_merged)))

colnames(res_merged) = c("group", "method_numeric", "node_count", "fit_time", "ess_bulk_mean", "ess_bulk_05", "ess_bulk_95", "ess_tail_mean", "ess_tail_05", "ess_tail_95")
res_merged$method = ifelse(res_merged$method_numeric==1,"l2norm",ifelse(res_merged$method_numeric==2,"cholesky","direct"))

for(i in 2:10){
 res_merged[,i] = as.numeric(res_merged[,i])
}

res_long = NULL
legal_set = c("D_corr[1,2]", "D_corr[1,3]", "D_corr[2,3]",  "D_corr[1,4]", "D_corr[2,5]", "D_corr[3,6]",  "D_corr[1,5]", "D_corr[1,6]", "D_corr[2,6]")

for(i in 1:length(legal_set)){
res_thin = res[which(res$variable == legal_set[i]),]
res_thin$group = paste0(res_thin$method,"_", res_thin$node_count)
res_thin$method_numeric = ifelse(res_thin$method=="l2norm",1,ifelse(res_thin$method=="cholesky",2,3))
res_long[[i]] = res_thin
}

res_thin = do.call(rbind, res_long)
res_thin$grouping = paste0(res_thin$method, res_thin$variable)

################################################################### Plot
# ggplot(data=res_thin, aes(x=node_count, y=ess_bulk, group=grouping, color=method, fill = method)) + geom_path() 
# ggplot(data=res_thin, aes(x=node_count, y=ess_tail, group=grouping, color=method, fill = method)) + geom_path() 

################################################################### Plot
res_merged$Method = res_merged$method

p1 = ggplot(data=res_merged, aes(x=node_count, y=fit_time, group=Method, color=Method, fill = Method)) +
  geom_line(linewidth=1.5) + ylab("Computation time (seconds)") + xlab("Nodes in network") + 
  scale_color_manual(values=c("l2norm" = colors2[5], "cholesky" = colors2[7], "direct" = colors2[3])) +
  theme(legend.position="bottom")

ggsave("Compare_time.pdf", p1, height=4, width=4)

p2 = ggplot(data=res_merged, aes(x=node_count, y=ess_bulk_mean, group=Method, color=Method, fill = Method)) +
  geom_line(linewidth=1.5) + ylab("Effective samples (bulk)") + xlab("Nodes in network") + ylim(0,1000) +
  scale_color_manual(values=c("l2norm" = colors2[5], "cholesky" = colors2[7], "direct" = colors2[3])) +
  theme(legend.position="bottom")

ggsave("Compare_Bulk_ESS.pdf", p2, height=4, width=4)

p3 = ggplot(data=res_merged, aes(x=node_count, y=ess_tail_mean, group=Method, color=Method, fill = Method)) +
  geom_line(linewidth=1.5) + ylab("Effective samples (tail)") + xlab("Nodes in network") + ylim(0,1000) +
  scale_color_manual(values=c("l2norm" = colors2[5], "cholesky" = colors2[7], "direct" = colors2[3])) +
  theme(legend.position="bottom")

ggsave("Compare_Tail_ESS.pdf", p3, height=4, width=4)
