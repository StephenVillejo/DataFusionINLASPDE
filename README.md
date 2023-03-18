# DataFusionINLASPDE

These are the R codes to implement the method discussed in the paper "Data Fusion in a Two-stage Spatio-Temporal Model using the INLA-SPDE Approach" by Stephen jun Villejo, Janine Illian, and Ben Swallow.

The following provides a description of the different files:

1. DataSimulation.R

  This files contains five important functions:
    a) "Data_simulation()" simulates a data from the proposed spatio-temporal model. It gives you the Poisson counts at the area level, the observed values in a network of monitoring stations, a high-resolution data (can be an output of a numerical model or satellite data), and the covariate data.
    b) The "rspde()" function simulates a Matern field. This function is borrowed from [tutorials/spde/spde-tutorial-functions.R](https://github.com/grantbrown/inla/blob/master/tutorials/spde/spde-tutorial-functions.R).
    c) The get.monitoring.stations() function gets the locations of the monitors. We use the same point locations for all the scenarios and replication. Hence, this function will either simulate the locations first if they have not yet been defined then compute the observed values at those location; or it will use previously defined locations and simply update the observed values at those locations. This function is adapted from https://github.com/michelacameletti/INLA_COSP. 
    d) The compute.intersection() function computes the intersection between between each area in the shapefile and the prediction grid and/or simulation grid. The output from this function will be used to compute spatial averages, i.e., too aggregate values from point-level to area-level. 
    e) The compute.predpoints() function identifies the prediction points in each area in the shapefile. The output from this function will be used to compute spatial averages, i.e., too aggregate values from point-level to area-level. This function and the previous one are also adapted from https://github.com/michelacameletti/INLA_COSP. 

  
2. Fitting.Stage1.Stage2.R 

  This file contains two important functions: a) The function "stage1.fit.ST.final.Melding()" performs the data fusion using the Bayesian melding model. The default approach is to perform a two-stage process while accounting for the uncertainty in the first stage model estimation. The actual simulation from the posterior predictive distributions of the stage one model is performed in another function. This function returns the model estimates both for the case of having informative and non-informative priors and also some plots. b) The function "stage2.fit.ST()" fits the generalized linear mixed model to link the Poisson counts and the estimated exposures at the area level. This functions simulates several times from the psoterior predictive distribution of the stage one model, and then fits the GLMM each time. This function returns the model estimates both for the case of having informative and non-informative priors and also some plots.
  
  3. OtherHelperFunctions.R
  
  This file contains helper functions, mostly functions used to produce plots and functions which computes the bias, RMSE, and coverage probabiltiies to assess the performance of the proposed method under different scenarios.
  
  
