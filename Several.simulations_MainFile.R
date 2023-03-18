
cluster = FALSE # if computer clusters are used or not

############################################################
# Load required packages
############################################################

packages <- c('rgdal','rgeos','spdep','ggplot2','colorspace','gridExtra','raster','psych',
              'disaggregation','parallel','INLABMA','patchwork','plyr', 'maptools',
              'viridis','ggthemes','broom','dplyr','reshape2','parallel')
lapply(packages, require, character.only = TRUE)
library(INLA)
inla.setOption(inla.mode='experimental')

############################################################
# Set working director, load functions / data / scenario settings
############################################################

if(!cluster) {
  setwd("/Users/stephenjunvillejo/Documents/Phd Studies - University of Glasgow/PhD Dissert/Model runs")
  source(here::here(getwd(), "/R/FinalFunctions_Paper1/DataSimulation.R"))
  source(here::here(getwd(), "/R/FinalFunctions_Paper1/Fitting.Stage1.Stage2.R"))
  source(here::here(getwd(), "/R/FinalFunctions_Paper1/FullSimulation.OneReplicate.R"))
  source(here::here(getwd(), "/R/FinalFunctions_Paper1/OtherHelperFunctions.R"))
  path = paste0(getwd(),"/FinalFigures/") # this is were figures are saved
  scenario.settings = read.csv(paste0(getwd(),"/Data/Scenarios.Settingsv2.csv"))
}else{
  source(here::here("/home/pgrad2/2510006v/Utils/DataSimulation.R"))
  source(here::here("/home/pgrad2/2510006v/Utils/Fitting.Stage1.Stage2.R"))
  source(here::here("/home/pgrad2/2510006v/Utils/FullSimulation.OneReplicate.R"))
  source(here::here("/home/pgrad2/2510006v/Utils/OtherHelperFunctions.R"))
  path = "/home/pgrad2/2510006v/FinalFigures/"
  scenario.settings = read.csv("/home/pgrad2/2510006v/Scenarios.Settingsv2.csv")
}


############################################################
# Set the parameters which do not change across simulations or scenarios
############################################################

kappa = 1.5; sigma2xi = 1.5; rho = sqrt(8)/kappa # Matern model parameters
b0 = 0; b1 = 2; a = 0.7 # exposure model parameters
RR = 1.2; gamma1_Poisson = log(RR); gamma0_Poisson = -3; sigma2iid_Poisson = 0.02; sigma2iid_time = 0.02 # poisson model parameters
highres.error.sigma2 = 1; highres.b0 = -1; highres.b1 = 1.5 # proxy data parameters
sigma2_e = .1; # classical error variance
res_pred_grid = 100 # resolution of the prediction grid 
highres_simulation_reggrid =100 # resolution of the simulation grid
unif.pop.a = 20; unif.pop.b = 50 # parameter of the uniform distribution for the population-at-risk per grid

prop.points.for.area = 0.05 # proportion of monitors in each area for the case of non-sparse data

n_random_fromposterior = 200 # how many values to sample from the parameter posterior distr of the exposure model
n_random_fromPPD = 10 # how many values to sample from the exposure posterior distribution (for each area)
n_simulations = 5 # number of simulations/replications for each scenario
save.plot = TRUE

subgroup = "A" # independent batches (just to make computing faster)
#SET HERE the scenarios you want to run!
run_scenarios = 1

n.cores = detectCores()

#################################
# Run the simulation! 
#################################

for (i in run_scenarios){
  
  Simulation.Run.Melding(x = i,
                         scenario.settings = scenario.settings,
                         subgroup = subgroup,
                         save.plot = save.plot)
  
}


