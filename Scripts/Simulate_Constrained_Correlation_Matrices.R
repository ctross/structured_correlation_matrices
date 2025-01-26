################################################################################################## Simulate multiplex networks
########################################################################################### Setup
N_layers = 13

# mode options: "cholesky", "l2norm"
# setting options: "multiplex_dyadic_reciprocity", "longitudinal_dyadic_reciprocity", "longitudinal_generalized_reciprocity"

data = list(N_responses = N_layers, bandage_penalty=0.01)  

################################################################################################################################################# Multiplex_dyadic_reciprocity
set.seed(3)
####################################################################################################### Cholesky method
dr_cholesky_style = generate_multiplex_correlation_matrix(data = data, mode="cholesky", eta=2.0, setting="multiplex_dyadic_reciprocity", 
                                                                   stan_mcmc_parameters = list(
                                                                   chains = 1, parallel_chains = 1, refresh = 1,
                                                                   iter_warmup = 1000, iter_sampling = 35, init = 0,
                                                                   max_treedepth = 12, adapt_delta = 0.98))

dr_Rho_1a = matrix(matrix(dr_cholesky_style$fit$draws("D_corr", format = "matrix")[1,]), nrow=(2*N_layers), ncol=(2*N_layers))
image(dr_Rho_1a[,(2*N_layers):1])
dr_Rho_1a

m1 = matrix(matrix(dr_cholesky_style$fit$draws("D_corr", format = "matrix")[34,]), nrow=(2*N_layers), ncol=(2*N_layers))[,(2*N_layers):1]
colnames(m1) = c(1:(2*N_layers))
rownames(m1) = c(1:(2*N_layers))

####################################################################################################### Cholesky method
dr_l2norm_style = generate_multiplex_correlation_matrix(data = data, mode="l2norm", eta=2.0, setting="multiplex_dyadic_reciprocity", 
                                                                   stan_mcmc_parameters = list(
                                                                   chains = 1, parallel_chains = 1, refresh = 1,
                                                                   iter_warmup = 1000, iter_sampling = 10, init = 0,
                                                                   max_treedepth = 12, adapt_delta = 0.98))

dr_Rho_2a = matrix(matrix(dr_l2norm_style$fit$draws("D_corr", format = "matrix")[1,]), nrow=(2*N_layers), ncol=(2*N_layers))
image(dr_Rho_2a[,(2*N_layers):1])
dr_Rho_2a




################################################################################################################################################# Longitudinal_dyadic_reciprocity
set.seed(3)
####################################################################################################### Cholesky method
 dr_cholesky_style = generate_multiplex_correlation_matrix(data = data, mode="cholesky", eta=2.0, setting="longitudinal_dyadic_reciprocity", 
                                                                    stan_mcmc_parameters = list(
                                                                    chains = 1, parallel_chains = 1, refresh = 1,
                                                                    iter_warmup = 1000, iter_sampling = 35, init = 0,
                                                                    max_treedepth = 12, adapt_delta = 0.98))

 dr_Rho_1b = matrix(matrix(dr_cholesky_style$fit$draws("D_corr", format = "matrix")[1,]), nrow=(2*N_layers), ncol=(2*N_layers))
 image(dr_Rho_1b[,(2*N_layers):1])
 dr_Rho_1b

 m2 = matrix(matrix(dr_cholesky_style$fit$draws("D_corr", format = "matrix")[19,]), nrow=(2*N_layers), ncol=(2*N_layers))[,(2*N_layers):1]
 colnames(m2) = c(1:(2*N_layers))
 rownames(m2) = c(1:(2*N_layers))

####################################################################################################### Cholesky method
dr_l2norm_style = generate_multiplex_correlation_matrix(data = data, mode="l2norm", eta=2.0, setting="longitudinal_dyadic_reciprocity", 
                                                                   stan_mcmc_parameters = list(
                                                                   chains = 1, parallel_chains = 1, refresh = 1,
                                                                   iter_warmup = 1000, iter_sampling = 10, init = 0,
                                                                   max_treedepth = 12, adapt_delta = 0.98))

dr_Rho_2b = matrix(matrix(dr_l2norm_style$fit$draws("D_corr", format = "matrix")[1,]), nrow=(2*N_layers), ncol=(2*N_layers))
image(dr_Rho_2b[,(2*N_layers):1])
dr_Rho_2b



################################################################################################################################################# Longitudinal_generalized_reciprocity
set.seed(1)
####################################################################################################### Cholesky method
 dr_cholesky_style = generate_multiplex_correlation_matrix(data = data, mode="cholesky", eta=2.0, setting="longitudinal_generalized_reciprocity", 
                                                                    stan_mcmc_parameters = list(
                                                                    chains = 1, parallel_chains = 1, refresh = 1,
                                                                    iter_warmup = 1000, iter_sampling = 25, init = 0,
                                                                    max_treedepth = 12, adapt_delta = 0.98))

 dr_Rho_1c = matrix(matrix(dr_cholesky_style$fit$draws("D_corr", format = "matrix")[1,]), nrow=(2*N_layers), ncol=(2*N_layers))
 image(dr_Rho_1c[,(2*N_layers):1])
 dr_Rho_1c

 m3 = matrix(matrix(dr_cholesky_style$fit$draws("D_corr", format = "matrix")[24,]), nrow=(2*N_layers), ncol=(2*N_layers))[,(2*N_layers):1]
 colnames(m3) = c(1:(2*N_layers))
 rownames(m3) = c(1:(2*N_layers))

####################################################################################################### Cholesky method
dr_l2norm_style = generate_multiplex_correlation_matrix(data = data, mode="l2norm", eta=2.0, setting="longitudinal_generalized_reciprocity", 
                                                                   stan_mcmc_parameters = list(
                                                                   chains = 1, parallel_chains = 1, refresh = 1,
                                                                   iter_warmup = 1000, iter_sampling = 10, init = 0,
                                                                   max_treedepth = 12, adapt_delta = 0.98))

dr_Rho_2c = matrix(matrix(dr_l2norm_style$fit$draws("D_corr", format = "matrix")[3,]), nrow=(2*N_layers), ncol=(2*N_layers))
image(dr_Rho_2c[,(2*N_layers):1])
dr_Rho_2c


################################################################### Paper plots

# Transform the matrix in long format - Multiplex DR
df = reshape2::melt(m1, as.is = TRUE)
colnames(df) = c("x", "y", "value")
p1 = ggplot(df, aes(x = as.numeric(x), y = as.numeric(y), fill = value)) +
  geom_tile() +
  coord_fixed() +   scale_fill_gradientn(colors = rev(colors2[1:8])) + xlab("") + ylab("") + theme_void() + theme(legend.position="none") + 
theme(axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()) + theme(panel.background = element_blank()) 
ggsave("Multiplex_DR.pdf", p1, height=4, width=4)

# Transform the matrix in long format - Long DR
df = reshape2::melt(m2, as.is = TRUE)
colnames(df) = c("x", "y", "value")
p2 = ggplot(df, aes(x = as.numeric(x), y = as.numeric(y), fill = value)) +
  geom_tile() +
  coord_fixed() +   scale_fill_gradientn(colors = rev(colors2[1:8])) + xlab("") + ylab("") + theme_void() + theme(legend.position="none") + 
theme(axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()) + theme(panel.background = element_blank()) 
ggsave("Long_DR.pdf", p2, height=4, width=4)


# Transform the matrix in long format - Long GR
df = reshape2::melt(m3, as.is = TRUE)
colnames(df) = c("x", "y", "value")
p3 = ggplot(df, aes(x = as.numeric(x), y = as.numeric(y), fill = value)) +
  geom_tile() +
  coord_fixed() +   scale_fill_gradientn(colors = rev(colors2[1:8])) + xlab("") + ylab("") + theme_void() + theme(legend.position="none") + 
theme(axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()) + theme(panel.background = element_blank()) 
ggsave("Long_GR.pdf", p3, height=4, width=4)



