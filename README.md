# DataFusionINLASPDE

These are the R codes to implement the method discussed in the paper "Data Fusion in a Two-stage Spatio-Temporal Model using the INLA-SPDE Approach" by Stephen Jun Villejo, Janine Illian, and Ben Swallow.

The following provides a description of the different files:

1. DataSimulation.R

  This file contains five important functions:
    a) "Data_simulation()" simulates a data from the proposed spatio-temporal model. It gives you the Poisson counts at the area level, the observed values in a network of monitoring stations, a high-resolution data (can be an output of a numerical model or satellite data), and the covariate data.
    b) The "rspde()" function simulates a Matern field. This function is borrowed from [tutorials/spde/spde-tutorial-functions.R](https://github.com/grantbrown/inla/blob/master/tutorials/spde/spde-tutorial-functions.R).
    c) The "get.monitoring.stations()" function gets the point locations (coordinates) of the monitors. We use the same point locations for all the scenarios and replications. Hence, this function will either simulate the locations first if they have not yet been defined then compute the observed values at those location; or it will use previously defined point locations and simply update the observed values at those locations. This function is adapted from https://github.com/michelacameletti/INLA_COSP. 
    d) The "compute.intersection()" function computes the intersection between each area in the shapefile and the prediction grid and/or simulation grid. The output from this function will be used to compute spatial averages, i.e., to aggregate values from point-level to area-level. 
    e) The "compute.predpoints()" function identifies the prediction points in each area in the shapefile. The output from this function will also be used to compute spatial averages, i.e., to aggregate values from point-level to area-level. This function and the previous one are also adapted from https://github.com/michelacameletti/INLA_COSP. 

  
2. Fitting.Stage1.Stage2.R 

  This file contains two important functions: 
    a) The function "stage1.fit.ST.final.Melding()" performs data fusion using the Bayesian melding model. The default approach is to perform a two-stage process  which accounts for the uncertainty in the first stage model estimation. The actual simulation from the posterior predictive distributions of the fitted model is performed in another function. This function returns the model estimates both for the case of using informative and non-informative priors; it also returns some useful plots if one wishes. 
    b) The function "stage2.fit.ST()" fits the generalized linear mixed model to link the Poisson counts and the estimated exposures at the area level. This functions simulates several times from the posterior predictive distribution of the stage one model, and then fits the GLMM each time. This function returns the model estimates both for the case of using informative and non-informative priors; it also returns some useful plots if one wishes. 
  
  3. OtherHelperFunctions.R
  
  This file contains helper functions - mostly functions which produce useful plots and functions which computes the bias, RMSE, and coverage probabilities to assess the performance of the proposed method under different scenarios.
 
  4. FullSimulation.OneReplicate.R

  This file contains a single function called "Simulation.Run.Melding()". This function performs one full simulation, .i.e., from simulating a spatio-temporal data, to model fitting, saving important results, and saving some useful plots if one wishes. The user needs to specify the scenario of interest. Also, this function can perform several simulations (several independent replications) via the argument "n_simulations", which is a global parameter.
  
  5. Several.simulations_MainFile.R

  This functions loads the packages and sets the working directory and the filepath of the folder for the figures and results. This R file also sets the true values of the parameters for the simulation study. Other important aspects in the simulation set-up and model fitting are also in this file; this includes the number of samples that we want from the marginal posteriors, the number of times we sample from the posterior predictive distribution of the stage one model, the number of independent replicates for each scenario setting, and the subgroup identifier since the scenarios are performed in different batches. This function runs the function "Simulation.Run.Melding()".
  
  
  
