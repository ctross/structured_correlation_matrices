####################### Simulate multiplex networks
# install_github('ctross/PlvsVltra') # Takes way to long to install for a color palatte, sorry
 library(PlvsVltra) # For colors
 colors = plvs_vltra("dust_storm", rev=FALSE, elements=NULL, show=FALSE)
 colors2 = plvs_vltra("mystic_mausoleum", rev=FALSE, elements=NULL, show=FALSE)

################## Packages
library(STRAND)
library(rethinking)
library(cmdstanr)
library(reshape2)
library(stringr)
library(ggplot2)
library(psych)
library(parallel)
library(purrr)
library(RColorBrewer)
library(Matrix)

################## Setup
setwd("C:/Users/pound/Dropbox/Open Papers/Block Correlation Matrices")

source("Code/Functions/build_multiplex_bindings.R")
source("Code/Functions/generate_multiplex_correlation_matrix.R")

################## Scripts
source("Code/Scripts/Simulate_Constrained_Correlation_Matrices.R") # This is all you really need to look at, tracing back to Stan code

# Below assumes that you have a sever with 100+ cores for running the test examples, and you need the private version of STRAND
source("Code/Scripts/Multiplex_Empirical_Example.R")
source("Code/Scripts/Multiplex_Empirical_Example_Sweep.R")
source("Code/Scripts/Multiplex_Simulation_Sweep.R")
source("Code/Scripts/Longitudinal_Empirical_Example.R")
source("Code/Scripts/Longitudinal_Simulation_Sweep.R")
