########################### Parameter sweeps
# Make data
N_time_top = 10
N_nodes_top = 7
N_datasets = N_time_top*N_nodes_top
set.seed(420)

# Storage for key quantities
dat_list_cholesky = vector("list", N_datasets)
dat_list_l2norm = vector("list", N_datasets)

dat_list = vector("list", N_datasets)
dat_list_Rho = vector("list", N_datasets)
G_list = vector("list", N_datasets)
D_list = vector("list", N_datasets)
N_list = vector("list", N_datasets)
L_list = vector("list", N_datasets)
M_list = vector("list", N_datasets)

# Sample sizes and layer number to test
Samp = round(exp(seq(2.7, 4.55,length.out=N_nodes_top)),0)
Time_steps = seq(1, 10, length.out=N_time_top)*2
Time_steps2 = seq(1, 10, length.out=N_time_top)

# Just some random labels for the layers
layer_labels =  paste0("Time ", Time_steps2)

########################################### Making DR matrix is a bit complicated
dr_list = vector("list", N_time_top)

make_Dr = function(N_timesteps){

  dr_seq = seq(from=0.8, to=0.1, length.out=N_timesteps)

  D = matrix(NA, nrow=(N_timesteps*2), ncol=(N_timesteps*2))

  D[1,] = c(c(1, dr_seq[1:(N_timesteps-1)]), 0.6*dr_seq)  
    
  for(k in 1:(N_timesteps-1)){
    K =  N_timesteps - k
    for(m in 1:K){
     D[m, m+k] = D[1, k+1]
     D[m, N_timesteps + m + k] = D[1, N_timesteps + k + 1]
    }
   }
     
  for(m in 1:N_timesteps){
     D[m, m+N_timesteps] = D[1, N_timesteps + 1]
    }

  for(m in 1:(N_timesteps-1)){
   for(n in (m+1):N_timesteps){
     D[m+N_timesteps, n+N_timesteps] = D[m, n] 
    }}

  for(m in 1:(N_timesteps-1)){
   for(n in (m+1):N_timesteps){
     D[n, m+N_timesteps] = D[m, n+N_timesteps]
    }}

  for(i in 1:((2*N_timesteps)-1)){
   for(j in (i+1):(N_timesteps*2)){
     D[j, i] = D[i, j]
    }}

  for(i in 1:((2*N_timesteps))){
     D[i, i] = 1
    }

    D2 = as.matrix(round(D,3))

    D3 = nearPD(D2, corr=TRUE)

    dr_Rho = as.matrix(D3$mat)  

    return(dr_Rho)
}

for(i in 1:N_time_top){
dr_list[[i]] = as.matrix(make_Dr(Time_steps[i]))  
chol(dr_list[[i]])
}
 # Make sure all calls to chol succeed

########################################### Making GR matrix is a bit more complicated
gr_list = vector("list", N_time_top)

make_Gr = function(N_timesteps){

  dr_seq = seq(from=0.8, to=0.1, length.out=N_timesteps)

 G = matrix(NA, nrow=(N_timesteps*2), ncol=(N_timesteps*2))
  
    Vals = seq(from=0.5, to=0.15, length.out=N_timesteps)              # Sender
    Vals2 = seq(from=0.35, to=0.01, length.out=N_timesteps)            # Cross
    Vals3 = rep(0.8, N_timesteps)                                      # GR 
    Vals4 = seq(from=0.7, to=0.35, length.out=N_timesteps)             # Reciever
    Vals5 = seq(from=0.25, to=0.1, length.out=N_timesteps)             # Cross

    target = 0
    for(m in 1:(N_timesteps-1)){
    for(n in (m+1):N_timesteps){
     target = target + 1
     G[m, n] = Vals[target]
     G[m+N_timesteps, n+N_timesteps] = Vals4[target]

     G[n, m+N_timesteps] = Vals2[target]
     G[m, n+N_timesteps] = Vals5[target]
    }}

  for(m in 1:N_timesteps){
     G[m, m+N_timesteps] = Vals3[m]
    }

  for(k in 1:(N_timesteps-1)){
    for(m in 1:(N_timesteps - k)){
     G[m, m+k] = G[1, k+1]
     G[m, N_timesteps + m + k] = G[1, N_timesteps + k + 1]
     G[N_timesteps + m, N_timesteps + m+k] = G[N_timesteps + 1, N_timesteps + k+1]
     G[m + k, N_timesteps + m] = G[k + 1, N_timesteps + 1]
    }
   }

  for(m in 1:N_timesteps){
     G[m, m+N_timesteps] = G[1, N_timesteps + 1]
    }

  for(i in 1:((2*N_timesteps)-1)){
   for(j in (i+1):(N_timesteps*2)){
     G[j, i] = G[i, j]
    }}

  for(i in 1:((2*N_timesteps))){
     G[i, i] = 1
    }

    G2 = as.matrix(round(G,3))

    G3 = nearPD(G2, corr=TRUE)

    sr_Rho = as.matrix(G3$mat)  

    return(sr_Rho)
}

for(i in 1:N_time_top){
gr_list[[i]] = as.matrix(make_Gr(Time_steps[i]))  
chol(gr_list[[i]])
}
 # Make sure all calls to chol succeed

########################################### Loop over layer sizes, and make reciprocity parameters
for(layers in 1:N_time_top){
 N_layers = Time_steps[layers]

# Set effect sizes
sr_mu = rep(0, N_layers*2)  
sr_sigma = rgamma(N_layers*2, 2.5, 2)+1.1
sr_Rho = gr_list[[layers]]
dr_mu = rep(0, N_layers) 
dr_sigma = rgamma(N_layers, 2.5, 2)+1.1
dr_Rho = dr_list[[layers]]

############################### Now loop over different sample sizes, and simulate network datasets
for(samp in 1:N_nodes_top){
  N_id = Samp[samp]

########### Simple intercepts
groups_temp = data.frame(group = factor(rep(1,N_id)))

groups = B = list()

for(i in 1:N_layers){
  B[[i]] = list(matrix(runif(1,-4.5,-2),nrow=1,ncol=1))
  groups[[i]] = groups_temp
}

# Save generative parameters
G_list[[(layers-1)*N_nodes_top + samp]] = sr_Rho
D_list[[(layers-1)*N_nodes_top + samp]] = dr_Rho
N_list[[(layers-1)*N_nodes_top + samp]] = Samp[samp]
L_list[[(layers-1)*N_nodes_top + samp]] = Time_steps[layers]

########## Simulate network
G = simulate_longitudinal_network(
  N_id = N_id,            
  N_timesteps = N_layers,                   
  B = B,                       
  V = 1,       
  groups = groups,                     
  sr_mu = sr_mu,            
  sr_sigma = sr_sigma,                        
  sr_Rho = sr_Rho,                     
  dr_mu = dr_mu,                            
  dr_sigma = dr_sigma,                         
  dr_Rho = dr_Rho,                          
  outcome_mode="bernoulli",
  link_mode="logit",                  
  individual_predictors = NULL,    
  dyadic_predictors = NULL,        
  individual_effects = NULL,        
  dyadic_effects = NULL           
 )


# Create the STRAND data object
outcome = NULL
exposure = NULL

for(i in 1:N_layers){
  outcome[[i]] = G$network[i,,]
  exposure[[i]] = G$exposure[i,,]
  colnames(outcome[[i]]) = rownames(outcome[[i]]) = paste0("ID",c(1:N_id))
  colnames(exposure[[i]]) = rownames(exposure[[i]]) = paste0("ID",c(1:N_id))
}

names(outcome) = layer_labels[1:N_layers]
names(exposure) = layer_labels[1:N_layers]

long_dat = NULL
for(i in 1:N_layers){
long_dat[[i]] = make_strand_data(outcome = list(Playing=outcome[[i]]),
                       block_covariates = NULL, 
                       individual_covariates = NULL, 
                       dyadic_covariates = NULL,
                       exposure = list(Playing=exposure[[i]]),
                       outcome_mode="bernoulli",
                       link_mode="logit",
                       longitudinal = TRUE)
 }

 names(long_dat) = paste0("Time ",1:N_layers)

# Store the STRAND data
dat_list[[(layers-1)*N_nodes_top + samp]] = long_dat
dat_list_Rho[[(layers-1)*N_nodes_top + samp]]$sr_Rho = sr_Rho
dat_list_Rho[[(layers-1)*N_nodes_top + samp]]$dr_Rho = dr_Rho

##################################
# Store the STRAND data for cholsky mode
dat_list_cholesky[[(layers-1)*N_nodes_top + samp]] = long_dat
dat_list_cholesky[[(layers-1)*N_nodes_top + samp]]$sr_Rho = sr_Rho
dat_list_cholesky[[(layers-1)*N_nodes_top + samp]]$dr_Rho = dr_Rho
dat_list_cholesky[[(layers-1)*N_nodes_top + samp]]$model_to_deploy = 1
M_list[[(layers-1)*N_nodes_top + samp]]$model = "cholesky"

# Store the STRAND data for l2norm mode
dat_list_l2norm[[(layers-1)*N_nodes_top + samp]] = long_dat
dat_list_l2norm[[(layers-1)*N_nodes_top + samp]]$sr_Rho = sr_Rho
dat_list_l2norm[[(layers-1)*N_nodes_top + samp]]$dr_Rho = dr_Rho
dat_list_l2norm[[(layers-1)*N_nodes_top + samp]]$model_to_deploy = 2
M_list[[(layers-1)*N_nodes_top + samp]]$model = "l2norm"
}

print("Finished")
print(layers)
}

dat_list = c(dat_list_cholesky, rev(dat_list_l2norm))

####################################################
# Run models using a server

######## Test fit using Frobenius_norm
Frobenius_norm = function(x,y){
 z = sqrt(tr(t(x-y) %*% (x-y)))
 return(z)
}

# Get random Frobenius_norm fit from prior
test_F_norm_prior = function(data){
  scrap = rep(NA, 1000)
  lay = dim(data)[1]

  for(i in 1:1000){
  scrap[i] = Frobenius_norm(data, rlkjcorr( 1 , lay , eta=1.5 )) 
  }
  
  results = rep(NA,3)
  results[1] = mean(scrap) 
  results[2:3] = HPDI(scrap, 0.9) 

  return(results)
}

# Get posterior Frobenius_norm fit 
test_F_norm_post = function(data, post){
  scrap = rep(NA, 1000)
  lay = dim(data)[2]

  for(i in 1:1000){
  scrap[i] = Frobenius_norm(data, post[i,,]) 
  }
  
  results = rep(NA,3)
  results[1] = mean(scrap) 
  results[2:3] = HPDI(scrap, 0.9) 

  return(results)
}


# Model function
run_and_parse = function(XXX){
  if(XXX$model_to_deploy == 1){
  fit = fit_longitudinal_model(long_data=head(XXX,-3),
                          block_regression = ~ 1,
                          focal_regression = ~ 1,
                          target_regression = ~ 1,
                          dyad_regression = ~ 1,
                          coefficient_mode="varying",
                          random_effects_mode="fixed",
                          mode="mcmc",
                          bandage_penalty = -1,
                          stan_mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 1,
                                                        iter_warmup = 1000, iter_sampling = 1000,
                                                        max_treedepth = 12, adapt_delta = 0.98, init=0)
 )

  res_temp = as.data.frame(fit$fit$summary(variables ="D_corr"))
  res_temp$fit_time = fit$fit$time()$total
  res_temp$node_count = XXX[[1]]$N_id
  res_temp$method = "cholesky"
}

  if(XXX$model_to_deploy == 2){
  fit = fit_longitudinal_model(long_data=head(XXX,-3),
                          block_regression = ~ 1,
                          focal_regression = ~ 1,
                          target_regression = ~ 1,
                          dyad_regression = ~ 1,
                          coefficient_mode="varying",
                          random_effects_mode="fixed",
                          mode="mcmc",
                          bandage_penalty = 0.01,
                          stan_mcmc_parameters = list(chains = 1, parallel_chains = 1, refresh = 1,
                                                        iter_warmup = 1000, iter_sampling = 1000,
                                                        max_treedepth = 12, adapt_delta = 0.98, init=0)
)

  res_temp = as.data.frame(fit$fit$summary(variables ="D_corr"))
  res_temp$fit_time = fit$fit$time()$total
  res_temp$node_count = XXX[[1]]$N_id
  res_temp$method = "l2norm"
}

res = summarize_strand_results(fit)

G_corr = apply(res$samples$srm_model_samples$G_corr, 2:3, median)
D_corr = apply(res$samples$srm_model_samples$D_corr, 2:3, median)

G_frob = test_F_norm_post(XXX$sr_Rho, res$samples$srm_model_samples$G_corr)
D_frob = test_F_norm_post(XXX$dr_Rho, res$samples$srm_model_samples$D_corr)

return(list(G_corr=G_corr, D_corr=D_corr, G_frob=G_frob, D_frob=D_frob, res_temp=res_temp))
}

# Deploy models on a server in parallel
fit = mclapply(1:(2*N_datasets), function(z){
                run_and_parse(dat_list[[z]])
                },
  mc.cores = N_datasets)

######################################################### Test visuals
#################################  Runtime
df_viz_1 = data.frame(factor(unlist(L_list)), factor(unlist(N_list)), NA)
df_viz_2 = data.frame(factor(unlist(L_list)), factor(unlist(N_list)), NA)

for(i in 1:70){
 df_viz_1[i,3] = log10(fit[[i]]$res_temp$fit_time[1])
 df_viz_2[i,3] = log10(fit[[i+70]]$res_temp$fit_time[1])
}

df_viz_3 = data.frame(factor(unlist(L_list)), factor(unlist(N_list)), df_viz_1[,3]-df_viz_2[,3])

df_viz_1$mode = "cholesky"
df_viz_2$mode = "l2norm"

df_viz = rbind(df_viz_1, df_viz_2)

colnames(df_viz) = c("x", "y", "Runtime","mode")
p1 = ggplot(df_viz, aes(x = x, y = y, fill = Runtime)) + ylab("Network size (nodes)") + xlab("Network layers") + facet_wrap(vars(mode)) +
  geom_tile() + coord_fixed() +   scale_fill_gradientn(colors = c(colors2[c(1,5,6,8,7,3)])) + theme(legend.position="bottom") + 
  labs(fill="Runtime (log10 seconds)\n ") +    guides(fill = guide_colorbar(title.vjust = 0.4, barwidth = 12, nbin = 1000))
p1

colnames(df_viz_3) = c("x", "y", "Difference")
p2 = ggplot(df_viz_3, aes(x = x, y = y, fill = Difference)) + ylab("Network size (nodes)") + xlab("Network layers") + 
  geom_tile() +
  coord_fixed() +   scale_fill_gradientn(colors = c(colors2[c(1,5,6,8,7,3)])) + theme(legend.position="bottom")
p2

ggsave("RunTimeAbs_DR_Long.pdf", p1, height=5, width=6)

df_viz_time = df_viz


#################################  Effective samples bulk
df_viz_1 = data.frame(factor(unlist(L_list)), factor(unlist(N_list)), NA)
df_viz_2 = data.frame(factor(unlist(L_list)), factor(unlist(N_list)), NA)

for(i in 1:70){
 df_viz_1[i,3] = (mean(fit[[i]]$res_temp$ess_bulk,na.rm=TRUE))
 df_viz_2[i,3] = (mean(fit[[i+70]]$res_temp$ess_bulk,na.rm=TRUE))
}

df_viz_3 = data.frame(factor(unlist(L_list)), factor(unlist(N_list)), df_viz_1[,3]-df_viz_2[,3])

df_viz_1$mode = "cholesky"
df_viz_2$mode = "l2norm"

df_viz = rbind(df_viz_1, df_viz_2)

colnames(df_viz) = c("x", "y", "ESSbulk","mode")
p1 = ggplot(df_viz, aes(x = x, y = y, fill = ESSbulk)) + ylab("Network size (nodes)") + xlab("Network layers") + facet_wrap(vars(mode)) +
  geom_tile() + coord_fixed() +   scale_fill_gradientn(colors = c(colors2[c(1,5,6,8,7,3)])) + theme(legend.position="bottom") +
  guides(fill = guide_colorbar(title.vjust = 0.4, barwidth = 12, nbin = 1000)) + labs(fill="Effective sample size (Bulk)\n ")
p1

colnames(df_viz_3) = c("x", "y", "Difference")
p2 = ggplot(df_viz_3, aes(x = x, y = y, fill = Difference)) + ylab("Network size (nodes)") + xlab("Network layers") + 
  geom_tile() +
  coord_fixed() +   scale_fill_gradientn(colors = c(colors2[c(1,5,6,8,7,3)])) + theme(legend.position="bottom")
p2

ggsave("ESS_Bulk_DR_Long.pdf", p1, height=5, width=6)

df_viz_bulk = df_viz


#################################  Effective samples tail
df_viz_1 = data.frame(factor(unlist(L_list)), factor(unlist(N_list)), NA)
df_viz_2 = data.frame(factor(unlist(L_list)), factor(unlist(N_list)), NA)

for(i in 1:70){
 df_viz_1[i,3] = (mean(fit[[i]]$res_temp$ess_tail,na.rm=TRUE))
 df_viz_2[i,3] = (mean(fit[[i+70]]$res_temp$ess_tail,na.rm=TRUE))
}

df_viz_3 = data.frame(factor(unlist(L_list)), factor(unlist(N_list)), df_viz_1[,3]-df_viz_2[,3])

df_viz_1$mode = "cholesky"
df_viz_2$mode = "l2norm"

df_viz = rbind(df_viz_1, df_viz_2)

colnames(df_viz) = c("x", "y", "ESStail","mode")
p1 = ggplot(df_viz, aes(x = x, y = y, fill = ESStail)) + ylab("Network size (nodes)") + xlab("Network layers") + facet_wrap(vars(mode)) +
  geom_tile() + coord_fixed() +   scale_fill_gradientn(colors = c(colors2[c(1,5,6,8,7,3)])) + theme(legend.position="bottom") +
  guides(fill = guide_colorbar(title.vjust = 0.4, barwidth = 12, nbin = 1000)) + labs(fill="Effective sample size (Tail)\n ")
p1

colnames(df_viz_3) = c("x", "y", "Difference")
p2 = ggplot(df_viz_3, aes(x = x, y = y, fill = Difference)) + ylab("Network size (nodes)") + xlab("Network layers") + 
  geom_tile() +
  coord_fixed() +   scale_fill_gradientn(colors = c(colors2[c(1,5,6,8,7,3)])) + theme(legend.position="bottom")
p2


ggsave("ESS_Tail_DR_Long.pdf", p1, height=5, width=6)

df_viz_tail = df_viz


#################################  Frob norm
df_viz_1 = data.frame(factor(unlist(L_list)), factor(unlist(N_list)), NA, NA, NA)
df_viz_2 = data.frame(factor(unlist(L_list)), factor(unlist(N_list)), NA, NA, NA)

for(i in 1:70){
 df_viz_1[i,3] = fit[[i]]$D_frob[1]
 df_viz_2[i,3] = fit[[i+70]]$D_frob[1]

 df_viz_1[i,4] = test_F_norm_prior(D_list[[i]])[1]
 df_viz_2[i,4] = test_F_norm_prior(D_list[[i]])[1]

 df_viz_1[i,5] = df_viz_1[i,3]/df_viz_1[i,4]
 df_viz_2[i,5] = df_viz_2[i,3]/df_viz_2[i,4]
}

df_viz_1$mode = "cholesky"
df_viz_2$mode = "l2norm"

df_viz = rbind(df_viz_1, df_viz_2)

colnames(df_viz) = c("x", "y", "PostFrob","PriorFrob","Frobenius","mode")
p1 = ggplot(df_viz, aes(x = x, y = y, fill = Frobenius)) + ylab("Network size (nodes)") + xlab("Network layers") + facet_wrap(vars(mode)) +
  geom_tile() + coord_fixed() +   scale_fill_gradientn(colors = c(colors2[c(1,5,6,8,7,3)])) + theme(legend.position="bottom") + 
  guides(fill = guide_colorbar(title.vjust = 0.4, barwidth = 12, nbin = 1000)) + labs(fill="Frobenius norm ratio\n ")
p1

ggsave("Frob_DR_Long.pdf", p1, height=5, width=6)

################################### ESS per Sec
df_viz = cbind(df_viz_time, df_viz_bulk, df_viz_tail)
df_viz$EffTail = log10(df_viz$ESStail) - df_viz$Runtime
df_viz$EffBulk = log10(df_viz$ESSbulk) - df_viz$Runtime

df_viz_short = data.frame(Outcome = c(c(df_viz$EffBulk[1:70] - df_viz$EffBulk[71:140]), c(df_viz$EffTail[1:70] - df_viz$EffTail[71:140])),
                          Type = rep(c("Bulk","Tail"),each=70),
                          x = c(df_viz$x, df_viz$x),
                          y = c(df_viz$y, df_viz$y)
                          )

col_limit = max(
  abs(max(df_viz_short$Outcome)),
  abs(min(df_viz_short$Outcome))
 )

p1 = ggplot(df_viz_short, aes(x = x, y = y, fill = Outcome)) + ylab("Network size (nodes)") + xlab("Time-steps in data") + facet_wrap(vars(Type)) +
  geom_tile() + coord_fixed() + 
  scale_fill_gradient2(low = colors2[5], mid = "white", high = colors2[7], limits = c(-col_limit,col_limit)) +
   theme(legend.position="bottom") + 
  labs(fill="Efficiency") + guides(fill = guide_colorbar(title.vjust = 0.4, barwidth = 12, nbin = 1000, draw.ulim = TRUE, draw.llim = TRUE))
p1

ggsave("Efficiency_Long.pdf", p1, height=5, width=6)

#################################### Loop over data sets and make calculations
G_layers = G_SampSize = G_res_f_test_m = G_res_f_test_l = G_res_f_test_h = matrix(NA, nrow=N_time_top, ncol=8)
D_layers = D_SampSize = D_res_f_test_m = D_res_f_test_l = D_res_f_test_h = matrix(NA, nrow=N_time_top, ncol=8)

for(i in 1:N_datasets){
  n = which(Samp == N_list[[i]])
  l = which(Time_steps == L_list[[i]])

  G_res_f_test_m[l,n] = fit[[i]]$G_frob[1]
  D_res_f_test_m[l,n] = fit[[i]]$D_frob[1]

  G_res_f_test_l[l,n] = fit[[i]]$G_frob[2]
  D_res_f_test_l[l,n] = fit[[i]]$D_frob[2]

  G_res_f_test_h[l,n] = fit[[i]]$G_frob[3]
  D_res_f_test_h[l,n] = fit[[i]]$D_frob[3]

  G_layers[l,n] = Time_steps[l]
  D_layers[l,n] = Time_steps[l]

  G_SampSize[l,n] = Samp[n]
  D_SampSize[l,n] = Samp[n]
}

for(i in 1:N_time_top){
 scrap_g = test_F_norm_prior(G_list[[i*N_nodes_top]])
 scrap_d = test_F_norm_prior(D_list[[i*N_nodes_top]])

 G_res_f_test_m[i,8] = scrap_g[1]
 D_res_f_test_m[i,8] = scrap_d[1]

 G_res_f_test_l[i,8] = scrap_g[2]
 D_res_f_test_l[i,8] = scrap_d[2]

 G_res_f_test_h[i,8] = scrap_g[3]
 D_res_f_test_h[i,8] = scrap_d[3]

 G_layers[i,8] = Time_steps[i]
 D_layers[i,8] = Time_steps[i]

 G_SampSize[i,8] = "Random"
 D_SampSize[i,8] = "Random"
}

######## Parse data and plot
df_F_G = data.frame(SampleSize = c(G_SampSize), Layers = c(G_layers), Median=c(G_res_f_test_m), L = c(G_res_f_test_l), H=c(G_res_f_test_h))
df_F_D = data.frame(SampleSize = c(D_SampSize), Layers = c(D_layers), Median=c(D_res_f_test_m), L = c(D_res_f_test_l), H=c(D_res_f_test_h))

df_F_G$Outcome = "Generalized"
df_F_D$Outcome = "Dyadic"

df_F = rbind(df_F_G,df_F_D)

df_F$Layers = factor(df_F$Layers)

df_F$SampleSize = factor(df_F$SampleSize)
df_F$SampleSize = factor(df_F$SampleSize, levels = c(Samp,"Random"))

nbcols = 7
mycolors = colorRampPalette(plvs_vltra("dust_storm", rev=FALSE, elements=c(5,2), show=FALSE))(nbcols)
mycolors = c(mycolors,"black")

p = ggplot(df_F, aes(x = Layers, y = Median, group = SampleSize, color=SampleSize, ymin = L, ymax = H)) + 
        geom_linerange(size = 1,position = position_dodge(width = 0.3)) + 
        geom_point(size = 2,position = position_dodge(width = 0.3)) + 
        facet_grid(. ~ Outcome, scales = "free", space = "free") + 
        geom_hline(aes(yintercept = 0), color = "black", linetype = "dashed") + 
        labs(y = "Frobenius norm", x = "Number of network layers") + 
        theme(strip.text.x = element_text(size = 12, face = "bold"), 
              strip.text.y = element_text(size = 12, face = "bold"), 
              axis.text = element_text(size = 12), 
              axis.title = element_text(size = 14, face = "bold")) + 
        theme(panel.spacing = grid::unit(1, "lines")) + 
        scale_color_manual(values = mycolors, name = "Sample size:")+
        theme(legend.position="bottom") + guides(color=guide_legend(nrow=1))



ggsave("Frobenius_Norms_Bernoulli_cholesky_Long.pdf",p, width=9, height=5.5)



############################# Loop over data sets and make calculations
G_layers = G_SampSize = G_res_f_test_m = G_res_f_test_l = G_res_f_test_h = matrix(NA, nrow=N_time_top, ncol=8)
D_layers = D_SampSize = D_res_f_test_m = D_res_f_test_l = D_res_f_test_h = matrix(NA, nrow=N_time_top, ncol=8)

for(i in 1:N_datasets){
  n = which(Samp == N_list[[i]])
  l = which(Time_steps == L_list[[i]])

  G_res_f_test_m[l,n] = fit[[i+N_datasets]]$G_frob[1]
  D_res_f_test_m[l,n] = fit[[i+N_datasets]]$D_frob[1]

  G_res_f_test_l[l,n] = fit[[i+N_datasets]]$G_frob[2]
  D_res_f_test_l[l,n] = fit[[i+N_datasets]]$D_frob[2]

  G_res_f_test_h[l,n] = fit[[i+N_datasets]]$G_frob[3]
  D_res_f_test_h[l,n] = fit[[i+N_datasets]]$D_frob[3]

  G_layers[l,n] = Time_steps[l]
  D_layers[l,n] = Time_steps[l]

  G_SampSize[l,n] = Samp[n]
  D_SampSize[l,n] = Samp[n]
}

for(i in 1:N_time_top){
 scrap_g = test_F_norm_prior(G_list[[i*N_nodes_top]])
 scrap_d = test_F_norm_prior(D_list[[i*N_nodes_top]])

 G_res_f_test_m[i,8] = scrap_g[1]
 D_res_f_test_m[i,8] = scrap_d[1]

 G_res_f_test_l[i,8] = scrap_g[2]
 D_res_f_test_l[i,8] = scrap_d[2]

 G_res_f_test_h[i,8] = scrap_g[3]
 D_res_f_test_h[i,8] = scrap_d[3]

 G_layers[i,8] = Time_steps[i]
 D_layers[i,8] = Time_steps[i]

 G_SampSize[i,8] = "Random"
 D_SampSize[i,8] = "Random"
}

######## Parse data and plot
df_F_G = data.frame(SampleSize = c(G_SampSize), Layers = c(G_layers), Median=c(G_res_f_test_m), L = c(G_res_f_test_l), H=c(G_res_f_test_h))
df_F_D = data.frame(SampleSize = c(D_SampSize), Layers = c(D_layers), Median=c(D_res_f_test_m), L = c(D_res_f_test_l), H=c(D_res_f_test_h))

df_F_G$Outcome = "Generalized"
df_F_D$Outcome = "Dyadic"

df_F = rbind(df_F_G,df_F_D)

df_F$Layers = factor(df_F$Layers)

df_F$SampleSize = factor(df_F$SampleSize)
df_F$SampleSize = factor(df_F$SampleSize, levels = c(Samp,"Random"))

nbcols = 7
mycolors = colorRampPalette(plvs_vltra("dust_storm", rev=FALSE, elements=c(5,2), show=FALSE))(nbcols)
mycolors = c(mycolors,"black")

px = ggplot(df_F, aes(x = Layers, y = Median, group = SampleSize, color=SampleSize, ymin = L, ymax = H)) + 
        geom_linerange(size = 1,position = position_dodge(width = 0.3)) + 
        geom_point(size = 2,position = position_dodge(width = 0.3)) + 
        facet_grid(. ~ Outcome, scales = "free", space = "free") + 
        geom_hline(aes(yintercept = 0), color = "black", linetype = "dashed") + 
        labs(y = "Frobenius norm", x = "Number of network layers") + 
        theme(strip.text.x = element_text(size = 12, face = "bold"), 
              strip.text.y = element_text(size = 12, face = "bold"), 
              axis.text = element_text(size = 12), 
              axis.title = element_text(size = 14, face = "bold")) + 
        theme(panel.spacing = grid::unit(1, "lines")) + 
        scale_color_manual(values = mycolors, name = "Sample size:")+
        theme(legend.position="bottom") + guides(color=guide_legend(nrow=1))



ggsave("Frobenius_Norms_Bernoulli_l2norm_Long.pdf",px, width=9, height=5.5)


