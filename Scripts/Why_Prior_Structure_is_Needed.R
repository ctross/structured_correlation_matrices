##############################################
#
#   Multiplex Simulation Analyses - Effect of prior structture
#
#############################################
# Make data
N_test = 10
Samp = round(exp(seq(2.5, 5.192, length.out=N_test)),0)
N_switches = 5
N_layers_top = 4
N_datasets = N_test*N_switches

set.seed(42+7)

G_list = vector("list", N_test)
D_list = vector("list", N_test)
N_list = vector("list", N_test)

# Just some random labels for the layers. Need to add names if you want more than 4 layers.
jungian_archetypes =  c("Explorer",
                        "Rebel",
                        "Hero",
                        "Lover")

# Storage for key quantities
dat_list_NoPrior = vector("list", N_datasets)
dat_list_BlockPrior = vector("list", N_datasets)

########################################### Making DR matrix is a bit complicated
make_Dr = function(K){
  if(K>1){
  dr_Rho_A = runif(K, -0.7, 0.7)
  dr_Rho_B = rlkjcorr( 1 , K , eta=1.5 ) 
  dr_Rho_C = rlkjcorr( 1 , K , eta=1.5 ) 

  diag(dr_Rho_B) = dr_Rho_A                    

  dr_Rho = as.matrix(rbind(cbind(dr_Rho_C, dr_Rho_B),       
                 cbind(t(dr_Rho_B), dr_Rho_C)) )  

  dr_Rho = nearPD(dr_Rho, corr = TRUE, keepDiag = TRUE)           
  return(dr_Rho$mat)
  }
  else{
    dr_Rho = rlkjcorr( 1 , K+1 , eta=1.5 ) 
    return(dr_Rho)
  }
}


dr = as.matrix(make_Dr(N_layers_top))  
chol(dr)

 # Make sure all calls to chol succeed

########################################### Set basic params
# Set effect sizes
sr_mu = rep(0, N_layers_top*2)  
sr_sigma = rgamma(N_layers_top*2, 2.5, 2)+1.1
sr_Rho = rlkjcorr( 1 , N_layers_top*2 , eta=1.5 )
dr_mu = rep(0, N_layers_top) 
dr_sigma = rgamma(N_layers_top, 2.5, 2)+2.5
dr_Rho = dr

############################### Now loop over different sample sizes, and simulate network datasets
for(samp in 1:N_test){
  N_id = Samp[samp]

########### Simple intercepts
groups = data.frame(group = factor(rep(1,N_id)))

B = list()

for(i in 1:N_layers_top){
  B[[i]] = list(matrix(runif(1,-10,-2),nrow=1,ncol=1))
}

# Save generative parameters
G_list[[samp]] = sr_Rho
D_list[[samp]] = dr_Rho
N_list[[samp]] = Samp[samp]


########## Simulate network
G = simulate_multiplex_network(
  N_id = N_id,            
  N_layers = N_layers_top,                   
  B = B,                       
  V = 1,       
  groups = groups,                     
  sr_mu = sr_mu,            
  sr_sigma = sr_sigma,                        
  sr_Rho = sr_Rho,                     
  dr_mu = dr_mu,                            
  dr_sigma = dr_sigma,                         
  dr_Rho = dr_Rho,                          
  outcome_mode = "bernoulli",
  link_mode = "logit",                
  individual_predictors = NULL,    
  dyadic_predictors = NULL,        
  individual_effects = NULL,        
  dyadic_effects = NULL           
 )

  for(iter in 1:N_switches){
    new_order_of_individuals = sample(1:N_id, replace=FALSE)

# Create the STRAND data object
outcome = vector("list", N_layers_top)
exposure = vector("list", N_layers_top)

for(i in 1:N_layers_top){
  outcome[[i]] = G$network[i,new_order_of_individuals,new_order_of_individuals]
  exposure[[i]] = G$exposure[i,new_order_of_individuals,new_order_of_individuals]

  rownames(exposure[[i]]) = colnames(exposure[[i]]) = 1:N_id
  rownames(outcome[[i]]) = colnames(outcome[[i]]) = 1:N_id
}

names(outcome) = jungian_archetypes[1:N_layers_top]
names(exposure) = jungian_archetypes[1:N_layers_top]

dat = make_strand_data(outcome = outcome,
                       block_covariates = NULL, 
                       individual_covariates = NULL, 
                       dyadic_covariates = NULL,
                       exposure = exposure,
                       outcome_mode="bernoulli",
                       link_mode="logit",
                       multiplex = TRUE)

dat$node_count = N_id
dat$layer_count = N_layers_top

# Store the STRAND data for the no prior model
dat_list_NoPrior[[(iter-1)*N_test + samp]] = dat
dat_list_NoPrior[[(iter-1)*N_test + samp]]$sr_Rho = sr_Rho
dat_list_NoPrior[[(iter-1)*N_test + samp]]$dr_Rho = dr_Rho
dat_list_NoPrior[[(iter-1)*N_test + samp]]$model_to_deploy = 1

# Store the STRAND data for the block prior model
dat_list_BlockPrior[[(iter-1)*N_test + samp]] = dat
dat_list_BlockPrior[[(iter-1)*N_test + samp]]$sr_Rho = sr_Rho
dat_list_BlockPrior[[(iter-1)*N_test + samp]]$dr_Rho = dr_Rho
dat_list_BlockPrior[[(iter-1)*N_test + samp]]$model_to_deploy = 2
}}

dat_list = c(dat_list_NoPrior, rev(dat_list_BlockPrior))

####################################################
# Run models using a server

# Model function
run_and_parse = function(XXX){
  if(XXX$model_to_deploy == 1){
  fit = fit_multiplex_model(data=XXX,
                          block_regression = ~ 1,
                          focal_regression = ~ 1,
                          target_regression = ~ 1,
                          dyad_regression = ~ 1,
                          mode="mcmc",
                          bandage_penalty = 100,
                          stan_mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 1,
                                                        iter_warmup = 500, iter_sampling = 500,
                                                        max_treedepth = 12, adapt_delta = 0.98, init=0)
 )

  res_temp = as.data.frame(fit$fit$summary(variables ="D_corr"))
  res_temp$fit_time = fit$fit$time()$total
  res_temp$node_count = XXX$node_count
  res_temp$method = "no_prior"
}

  if(XXX$model_to_deploy == 2){
  fit = fit_multiplex_model(data=XXX,
                          block_regression = ~ 1,
                          focal_regression = ~ 1,
                          target_regression = ~ 1,
                          dyad_regression = ~ 1,
                          mode="mcmc",
                          bandage_penalty = 0.01,
                          stan_mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 1,
                                                        iter_warmup = 500, iter_sampling = 500,
                                                        max_treedepth = 12, adapt_delta = 0.98, init=0)
)

  res_temp = as.data.frame(fit$fit$summary(variables ="D_corr"))
  res_temp$fit_time = fit$fit$time()$total
  res_temp$node_count = XXX$node_count
  res_temp$method = "block_prior"
}

res = summarize_strand_results(fit)

G_corr = apply(res$samples$srm_model_samples$G_corr, 2:3, median)
D_corr = apply(res$samples$srm_model_samples$D_corr, 2:3, median)


return(list(G_corr=G_corr, D_corr=D_corr, res_temp=res_temp))
}

# Deploy models on a server in parallel
fit = mclapply(1:(2*N_datasets), function(z){
                run_and_parse(dat_list[[z]])
                }, 
  mc.cores = 2*N_datasets)

######################################################### Test visuals

Res_dr = array(NA,c(N_datasets*2, N_layers_top*2, N_layers_top*2))
Res_model = array(NA,c(N_datasets*2, N_layers_top*2, N_layers_top*2))
Res_IDs = array(NA,c(N_datasets*2, N_layers_top*2, N_layers_top*2))
Res_IDs2 = array(NA,c(N_datasets*2, N_layers_top*2, N_layers_top*2))
Res_N = array(NA,c(N_datasets*2, N_layers_top*2, N_layers_top*2))
Res_Replicate = array(NA,c(N_datasets*2, N_layers_top*2, N_layers_top*2))

IDs2 = IDs = matrix(NA, nrow=N_layers_top*2, ncol=N_layers_top*2)

K = N_layers_top

IDs[1,2] = IDs[K + 1, 2 + K] = "Within 1"
IDs[1,3] = IDs[K + 1, 3 + K] = "Within 2"
IDs[1,4] = IDs[K + 1, 4 + K] = "Within 3"
IDs[2,3] = IDs[K + 2, 3 + K] = "Within 4"
IDs[2,4] = IDs[K + 2, 4 + K] = "Within 5"
IDs[3,4] = IDs[K + 3, 4 + K] = "Within 6"

IDs[1,2 + K] = IDs[2, K + 1] = "Between 1"
IDs[1,3 + K] = IDs[3, K + 1] = "Between 2"
IDs[1,4 + K] = IDs[4, K + 1] = "Between 3"
IDs[2,3 + K] = IDs[3, K + 2] = "Between 4"
IDs[2,4 + K] = IDs[4, K + 2] = "Between 5"
IDs[3,4 + K] = IDs[4, K + 3] = "Between 6"

Order = rep(rep(1:5, each=10),2)

IDs2[1,2] = "Top"
IDs2[1,3] = "Top"
IDs2[1,4] = "Top"
IDs2[2,3] = "Top"
IDs2[2,4] = "Top"
IDs2[3,4] = "Top"

IDs2[K + 1, 2 + K] = "Bottom"
IDs2[K + 1, 3 + K] = "Bottom"
IDs2[K + 1, 4 + K] = "Bottom"
IDs2[K + 2, 3 + K] = "Bottom"
IDs2[K + 2, 4 + K] = "Bottom"
IDs2[K + 3, 4 + K] = "Bottom"

IDs2[1,2 + K] = "Top"
IDs2[1,3 + K] = "Top"
IDs2[1,4 + K] = "Top"
IDs2[2,3 + K] = "Top"
IDs2[2,4 + K] = "Top"
IDs2[3,4 + K] = "Top"

IDs2[2, K + 1] = "Bottom"
IDs2[3, K + 1] = "Bottom"
IDs2[4, K + 1] = "Bottom"
IDs2[3, K + 2] = "Bottom"
IDs2[4, K + 2] = "Bottom"
IDs2[4, K + 3] = "Bottom"


for(i in 1:(N_datasets*2)){
  Res_dr[i,,] = fit[[i]]$D_corr
  Res_model[i,,] = matrix(ifelse(i>50,"With prior block-structure","Without prior block-structure"), nrow=N_layers_top*2, ncol=N_layers_top*2)
  Res_IDs[i,,] = IDs
  Res_IDs2[i,,] = IDs2
  Res_N[i,,] = matrix(fit[[i]]$res_temp$node_count[1], nrow=N_layers_top*2, ncol=N_layers_top*2)
  Res_Replicate[i,,] = matrix(Order[i], nrow=N_layers_top*2, ncol=N_layers_top*2)
}

df_long = data.frame(Value = c(Res_dr),
                     ElementID = c(Res_IDs),
                     ElementID2 = c(Res_IDs2),
                     SampleSize = c(Res_N),
                     Model = c(Res_model),
                     Replicate = c(Res_Replicate)
                 )

df_long = df_long[complete.cases(df_long),]

df_thin = df_long[,c("ElementID", "SampleSize", "Model", "Replicate")]
df_thin = df_thin[!duplicated(df_thin), ]

df_thin$Element1 = NA 
df_thin$Element2 = NA 

for(i in 1:length(df_thin$SampleSize)){
  df_thin$Element1[i] = df_long$Value[which(df_long$SampleSize == df_thin$SampleSize[i] & 
                                            df_long$Model == df_thin$Model[i] &
                                            df_long$ElementID == df_thin$ElementID[i] &
                                            df_long$Replicate == df_thin$Replicate[i] &
                                            df_long$ElementID2 == "Top"
                                            )]

  df_thin$Element2[i] = df_long$Value[which(df_long$SampleSize == df_thin$SampleSize[i] & 
                                            df_long$Model == df_thin$Model[i] &
                                            df_long$ElementID == df_thin$ElementID[i] &
                                            df_long$Replicate == df_thin$Replicate[i] &
                                            df_long$ElementID2 == "Bottom"
                                            )]
}


df_thin$InvSamp = log(1/df_thin$SampleSize) - log(1/50)

df_thin$ElementID2 = ifelse(df_thin$ElementID %in% paste("Within", 1:6), "Within", "Between")
df_thin$Replicate = factor(df_thin$Replicate)
df_thin$SampleSize = factor(df_thin$SampleSize)


df_thin2 = df_thin[which(df_thin$ElementID=="Between 4"),]

colors2 = plvs_vltra("mystic_mausoleum", rev=FALSE, elements=NULL, show=FALSE)

ggplot(df_thin2, aes(x=Element1, y=Element2, group=SampleSize, color=SampleSize, fill=SampleSize)) +
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_point(x=dr[2,7], y=dr[2,7], color="black", size=8, show.legend = FALSE) +
  coord_cartesian(xlim =c(-0.85, 0.85), ylim = c(-0.85, 0.85)) +
  geom_line(alpha=0.9) + 
  scale_color_manual(values = colorRampPalette(c(colors2[c(1,5,6,8,7,3)]))(10), name = "Sample Size") +
  scale_fill_manual(values = colorRampPalette(c(colors2[c(1,5,6,8,7,3)]))(10), name = "Sample Size") +
  geom_point(aes(x=Element1, y=Element2, group=SampleSize, color=SampleSize, fill=SampleSize), shape=21, alpha=0.9) +
  facet_wrap(. ~ Model) 




df_thin2 = df_thin[which(df_thin$ElementID=="Within 1"),]

colors2 = plvs_vltra("mystic_mausoleum", rev=FALSE, elements=NULL, show=FALSE)

p1 = ggplot(df_thin2, aes(x=Element1, y=Element2, group=SampleSize, color=SampleSize, fill=SampleSize)) +
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_point(x=dr[1,2], y=dr[1,2], color="black", size=8, show.legend = FALSE) +
  coord_cartesian(xlim =c(-0.85, 0.85), ylim = c(-0.85, 0.85)) +
  geom_line(alpha=0.9) + 
  scale_color_manual(values = colorRampPalette(c(colors2[c(1,5,6,8,7,3)]))(10), name = "Sample Size") +
  scale_fill_manual(values = colorRampPalette(c(colors2[c(1,5,6,8,7,3)]))(10), name = "Sample Size") +
  geom_point(aes(x=Element1, y=Element2, group=SampleSize, color=SampleSize, fill=SampleSize), shape=21, alpha=0.9) +
  facet_wrap(. ~ Model) 

df_thin2 = df_thin[which(df_thin$ElementID=="Within 2"),]

colors2 = plvs_vltra("mystic_mausoleum", rev=FALSE, elements=NULL, show=FALSE)

p2 = ggplot(df_thin2, aes(x=Element1, y=Element2, group=SampleSize, color=SampleSize, fill=SampleSize)) +
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_point(x=dr[1,3], y=dr[1,3], color="black", size=8, show.legend = FALSE) +
  coord_cartesian(xlim =c(-0.85, 0.85), ylim = c(-0.85, 0.85)) +
  geom_line(alpha=0.9) + 
  scale_color_manual(values = colorRampPalette(c(colors2[c(1,5,6,8,7,3)]))(10), name = "Sample Size") +
  scale_fill_manual(values = colorRampPalette(c(colors2[c(1,5,6,8,7,3)]))(10), name = "Sample Size") +
  geom_point(aes(x=Element1, y=Element2, group=SampleSize, color=SampleSize, fill=SampleSize), shape=21, alpha=0.9) +
  facet_wrap(. ~ Model) 


df_thin2 = df_thin[which(df_thin$ElementID=="Within 3"),]

colors2 = plvs_vltra("mystic_mausoleum", rev=FALSE, elements=NULL, show=FALSE)

p3 = ggplot(df_thin2, aes(x=Element1, y=Element2, group=SampleSize, color=SampleSize, fill=SampleSize)) +
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_point(x=dr[1,4], y=dr[1,4], color="black", size=8, show.legend = FALSE) +
  coord_cartesian(xlim =c(-0.85, 0.85), ylim = c(-0.85, 0.85)) +
  geom_line(alpha=0.9) + 
  scale_color_manual(values = colorRampPalette(c(colors2[c(1,5,6,8,7,3)]))(10), name = "Sample Size") +
  scale_fill_manual(values = colorRampPalette(c(colors2[c(1,5,6,8,7,3)]))(10), name = "Sample Size") +
  geom_point(aes(x=Element1, y=Element2, group=SampleSize, color=SampleSize, fill=SampleSize), shape=21, alpha=0.9) +
  facet_wrap(. ~ Model) 


df_thin2 = df_thin[which(df_thin$ElementID=="Within 4"),]

colors2 = plvs_vltra("mystic_mausoleum", rev=FALSE, elements=NULL, show=FALSE)

p4 = ggplot(df_thin2, aes(x=Element1, y=Element2, group=SampleSize, color=SampleSize, fill=SampleSize)) +
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_point(x=dr[2,3], y=dr[2,3], color="black", size=8, show.legend = FALSE) +
  coord_cartesian(xlim =c(-0.85, 0.85), ylim = c(-0.85, 0.85)) +
  geom_line(alpha=0.9) + 
  scale_color_manual(values = colorRampPalette(c(colors2[c(1,5,6,8,7,3)]))(10), name = "Sample Size") +
  scale_fill_manual(values = colorRampPalette(c(colors2[c(1,5,6,8,7,3)]))(10), name = "Sample Size") +
  geom_point(aes(x=Element1, y=Element2, group=SampleSize, color=SampleSize, fill=SampleSize), shape=21, alpha=0.9) +
  facet_wrap(. ~ Model) 


df_thin2 = df_thin[which(df_thin$ElementID=="Within 5"),]

colors2 = plvs_vltra("mystic_mausoleum", rev=FALSE, elements=NULL, show=FALSE)

p5 = ggplot(df_thin2, aes(x=Element1, y=Element2, group=SampleSize, color=SampleSize, fill=SampleSize)) +
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_point(x=dr[2,4], y=dr[2,4], color="black", size=8, show.legend = FALSE) +
  coord_cartesian(xlim =c(-0.85, 0.85), ylim = c(-0.85, 0.85)) +
  geom_line(alpha=0.9) + 
  scale_color_manual(values = colorRampPalette(c(colors2[c(1,5,6,8,7,3)]))(10), name = "Sample Size") +
  scale_fill_manual(values = colorRampPalette(c(colors2[c(1,5,6,8,7,3)]))(10), name = "Sample Size") +
  geom_point(aes(x=Element1, y=Element2, group=SampleSize, color=SampleSize, fill=SampleSize), shape=21, alpha=0.9) +
  facet_wrap(. ~ Model) 


df_thin2 = df_thin[which(df_thin$ElementID=="Within 6"),]

colors2 = plvs_vltra("mystic_mausoleum", rev=FALSE, elements=NULL, show=FALSE)

p6 = ggplot(df_thin2, aes(x=Element1, y=Element2, group=SampleSize, color=SampleSize, fill=SampleSize)) +
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_point(x=dr[3,4], y=dr[3,4], color="black", size=8, show.legend = FALSE) +
  coord_cartesian(xlim =c(-0.85, 0.85), ylim = c(-0.85, 0.85)) +
  geom_line(alpha=0.9) + 
  scale_color_manual(values = colorRampPalette(c(colors2[c(1,5,6,8,7,3)]))(10), name = "Sample Size") +
  scale_fill_manual(values = colorRampPalette(c(colors2[c(1,5,6,8,7,3)]))(10), name = "Sample Size") +
  geom_point(aes(x=Element1, y=Element2, group=SampleSize, color=SampleSize, fill=SampleSize), shape=21, alpha=0.9) +
  facet_wrap(. ~ Model) 









df_thin2 = df_thin[which(df_thin$ElementID=="Between 1"),]

colors2 = plvs_vltra("mystic_mausoleum", rev=FALSE, elements=NULL, show=FALSE)

b1 = ggplot(df_thin2, aes(x=Element1, y=Element2, group=SampleSize, color=SampleSize, fill=SampleSize)) +
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_point(x=dr[1,6], y=dr[1,6], color="black", size=8, show.legend = FALSE) +
  coord_cartesian(xlim =c(-0.85, 0.85), ylim = c(-0.85, 0.85)) +
  geom_line(alpha=0.9) + 
  scale_color_manual(values = colorRampPalette(c(colors2[c(1,5,6,8,7,3)]))(10), name = "Sample Size") +
  scale_fill_manual(values = colorRampPalette(c(colors2[c(1,5,6,8,7,3)]))(10), name = "Sample Size") +
  geom_point(aes(x=Element1, y=Element2, group=SampleSize, color=SampleSize, fill=SampleSize), shape=21, alpha=0.9) +
  facet_wrap(. ~ Model) 

df_thin2 = df_thin[which(df_thin$ElementID=="Between 2"),]

colors2 = plvs_vltra("mystic_mausoleum", rev=FALSE, elements=NULL, show=FALSE)

b2 = ggplot(df_thin2, aes(x=Element1, y=Element2, group=SampleSize, color=SampleSize, fill=SampleSize)) +
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_point(x=dr[1,7], y=dr[1,7], color="black", size=8, show.legend = FALSE) +
  coord_cartesian(xlim =c(-0.85, 0.85), ylim = c(-0.85, 0.85)) +
  geom_line(alpha=0.9) + 
  scale_color_manual(values = colorRampPalette(c(colors2[c(1,5,6,8,7,3)]))(10), name = "Sample Size") +
  scale_fill_manual(values = colorRampPalette(c(colors2[c(1,5,6,8,7,3)]))(10), name = "Sample Size") +
  geom_point(aes(x=Element1, y=Element2, group=SampleSize, color=SampleSize, fill=SampleSize), shape=21, alpha=0.9) +
  facet_wrap(. ~ Model) 


df_thin2 = df_thin[which(df_thin$ElementID=="Between 3"),]

colors2 = plvs_vltra("mystic_mausoleum", rev=FALSE, elements=NULL, show=FALSE)

b3 = ggplot(df_thin2, aes(x=Element1, y=Element2, group=SampleSize, color=SampleSize, fill=SampleSize)) +
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_point(x=dr[1,8], y=dr[1,8], color="black", size=8, show.legend = FALSE) +
  coord_cartesian(xlim =c(-0.85, 0.85), ylim = c(-0.85, 0.85)) +
  geom_line(alpha=0.9) + 
  scale_color_manual(values = colorRampPalette(c(colors2[c(1,5,6,8,7,3)]))(10), name = "Sample Size") +
  scale_fill_manual(values = colorRampPalette(c(colors2[c(1,5,6,8,7,3)]))(10), name = "Sample Size") +
  geom_point(aes(x=Element1, y=Element2, group=SampleSize, color=SampleSize, fill=SampleSize), shape=21, alpha=0.9) +
  facet_wrap(. ~ Model) 


df_thin2 = df_thin[which(df_thin$ElementID=="Between 4"),]

colors2 = plvs_vltra("mystic_mausoleum", rev=FALSE, elements=NULL, show=FALSE)

b4 = ggplot(df_thin2, aes(x=Element1, y=Element2, group=SampleSize, color=SampleSize, fill=SampleSize)) +
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_point(x=dr[2,7], y=dr[2,7], color="black", size=8, show.legend = FALSE) +
  coord_cartesian(xlim =c(-0.85, 0.85), ylim = c(-0.85, 0.85)) +
  geom_line(alpha=0.9) + 
  scale_color_manual(values = colorRampPalette(c(colors2[c(1,5,6,8,7,3)]))(10), name = "Sample Size") +
  scale_fill_manual(values = colorRampPalette(c(colors2[c(1,5,6,8,7,3)]))(10), name = "Sample Size") +
  geom_point(aes(x=Element1, y=Element2, group=SampleSize, color=SampleSize, fill=SampleSize), shape=21, alpha=0.9) +
  facet_wrap(. ~ Model) 


df_thin2 = df_thin[which(df_thin$ElementID=="Between 5"),]

colors2 = plvs_vltra("mystic_mausoleum", rev=FALSE, elements=NULL, show=FALSE)

b5 = ggplot(df_thin2, aes(x=Element1, y=Element2, group=SampleSize, color=SampleSize, fill=SampleSize)) +
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_point(x=dr[2,8], y=dr[2,8], color="black", size=8, show.legend = FALSE) +
  coord_cartesian(xlim =c(-0.85, 0.85), ylim = c(-0.85, 0.85)) +
  geom_line(alpha=0.9) + 
  scale_color_manual(values = colorRampPalette(c(colors2[c(1,5,6,8,7,3)]))(10), name = "Sample Size") +
  scale_fill_manual(values = colorRampPalette(c(colors2[c(1,5,6,8,7,3)]))(10), name = "Sample Size") +
  geom_point(aes(x=Element1, y=Element2, group=SampleSize, color=SampleSize, fill=SampleSize), shape=21, alpha=0.9) +
  facet_wrap(. ~ Model) 


df_thin2 = df_thin[which(df_thin$ElementID=="Between 6"),]

colors2 = plvs_vltra("mystic_mausoleum", rev=FALSE, elements=NULL, show=FALSE)

b6 = ggplot(df_thin2, aes(x=Element1, y=Element2, group=SampleSize, color=SampleSize, fill=SampleSize)) +
  geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_point(x=dr[3,8], y=dr[3,8], color="black", size=8, show.legend = FALSE) +
  coord_cartesian(xlim =c(-0.85, 0.85), ylim = c(-0.85, 0.85)) +
  geom_line(alpha=0.9) + 
  scale_color_manual(values = colorRampPalette(c(colors2[c(1,5,6,8,7,3)]))(10), name = "Sample Size") +
  scale_fill_manual(values = colorRampPalette(c(colors2[c(1,5,6,8,7,3)]))(10), name = "Sample Size") +
  geom_point(aes(x=Element1, y=Element2, group=SampleSize, color=SampleSize, fill=SampleSize), shape=21, alpha=0.9) +
  facet_wrap(. ~ Model) 







df_thin2a = df_thin[which(df_thin$ElementID=="Between 4"),]
colors2 = plvs_vltra("mystic_mausoleum", rev=FALSE, elements=NULL, show=FALSE)

b4 = ggplot(df_thin2a, aes(x=Element1, y=Element2, group=SampleSize, color=SampleSize, fill=SampleSize)) +
  geom_abline(intercept = 0, slope = 1, linetype="dashed", color="darkred") +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_point(x=dr[2,7], y=dr[2,7], color="black", size=6, show.legend = FALSE) +
  #coord_cartesian(xlim =c(-0.85, 0.85), ylim = c(-0.85, 0.85)) +
  geom_line(alpha=0.9) + 
  scale_color_manual(values = colorRampPalette(c(colors2[c(1,5,6,8,7,3)]))(10), name = "Sample Size") +
  scale_fill_manual(values = colorRampPalette(c(colors2[c(1,5,6,8,7,3)]))(10), name = "Sample Size") +
  geom_point(aes(x=Element1, y=Element2, group=SampleSize, color=SampleSize, fill=SampleSize), shape=21, alpha=0.9) +
  theme(strip.text.x = element_text(size = 12, face = "bold"), 
              strip.text.y = element_text(size = 12, face = "bold"), 
              axis.text = element_text(size = 12), 
              axis.title = element_text(size = 14, face = "bold")) + 
  facet_wrap(. ~ Model) + xlab("Rho(2,7)") + ylab("Rho(3,6)")



df_thin2b = df_thin[which(df_thin$ElementID=="Within 6"),]

colors2 = plvs_vltra("mystic_mausoleum", rev=FALSE, elements=NULL, show=FALSE)

p6 = ggplot(df_thin2b, aes(x=Element1, y=Element2, group=SampleSize, color=SampleSize, fill=SampleSize)) +
  geom_abline(intercept = 0, slope = 1, linetype="dashed", color="darkred") +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  geom_point(x=dr[3,4], y=dr[3,4], color="black", size=6, show.legend = FALSE) +
  coord_cartesian(xlim =c(-0.75, 0.1), ylim = c(-0.75, 0.1)) +
  geom_line(alpha=0.9) + 
  scale_color_manual(values = colorRampPalette(c(colors2[c(1,5,6,8,7,3)]))(10), name = "Sample Size") +
  scale_fill_manual(values = colorRampPalette(c(colors2[c(1,5,6,8,7,3)]))(10), name = "Sample Size") +
  geom_point(aes(x=Element1, y=Element2, group=SampleSize, color=SampleSize, fill=SampleSize), shape=21, alpha=0.9) +
  theme(strip.text.x = element_text(size = 12, face = "bold"), 
              strip.text.y = element_text(size = 12, face = "bold"), 
              axis.text = element_text(size = 12), 
              axis.title = element_text(size = 14, face = "bold")) + 
  facet_wrap(. ~ Model) + xlab("Rho(3,4)") + ylab("Rho(7,8)")


df_thin2a$Model2 = "Weak effect"
df_thin2b$Model2 = "Strong effect"

df_merged = rbind(df_thin2a,df_thin2b)

df_merged$Model2 = factor(df_merged$Model2)
df_merged$Model2 = factor(df_merged$Model2, levels=rev(levels(df_merged$Model2)))

df_merged$SampleSize

pM = ggplot(data=df_merged, aes(x=Element1, y=Element2, group=SampleSize, color=SampleSize, fill=SampleSize)) +
  geom_abline(intercept = 0, slope = 1, linetype="dashed", color="darkred") +
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_vline(xintercept = 0, linetype="dashed") +
  coord_cartesian(xlim =c(-0.75, 0.1), ylim = c(-0.75, 0.1)) +
  geom_point(data = data.frame(Element1=c(dr[2,7], dr[3,4]), Element2=c(dr[2,7], dr[3,4]), Model2 = c("Weak effect", "Strong effect"), SampleSize=NA), color="black", size=6, show.legend = FALSE) +
  geom_line(alpha=0.9) + 
  geom_point(aes(x=Element1, y=Element2, group=SampleSize, color=SampleSize, fill=SampleSize), shape=21, alpha=0.9) +
  scale_color_manual(values = colorRampPalette(c(colors2[c(1,5,6,8,7,3)]))(10), name = "Sample Size", na.translate = FALSE) +
  scale_fill_manual(values = colorRampPalette(c(colors2[c(1,5,6,8,7,3)]))(10), name = "Sample Size", na.translate = FALSE) +
  theme(strip.text.x = element_text(size = 16, face = "bold"), 
              strip.text.y = element_text(size = 16, face = "bold"), 
              axis.text = element_text(size = 16), 
              axis.title = element_text(size = 18, face = "bold")) + 
  facet_grid(Model2 ~ Model) + xlab("Rho(a,b)") + ylab("Rho(a+K, b+K)") + 
  theme(legend.position="bottom") + 
  guides(color=guide_legend(nrow=1), fill=guide_legend(nrow=1)) + theme(legend.text=element_text(size=14))



ggsave("CorrelationEstimateIssue.pdf", pM, width=8.5, height=8.5)
