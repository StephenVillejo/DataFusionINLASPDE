# DataFusionINLASPDE

These are the R codes to implement the method discussed in the paper "Data Fusion in a Two-stage Spatio-Temporal Model using the INLA-SPDE Approach" by Stephen jun Villejo, Janine Illian, and Ben Swallow.

The following provides a description of the different files:

1. DataSimulation.R
  This files contains five important functions:
    i. "Data_simulation()" simulates a data from the proposed spatio-temporal model. It gives you the Poisson counts at the area level, the observed values in a network of monitoring stations, a high-resolution data (can be an output of a numerical model or satellite data), and the covariate data.
    ii. The "rspde()" function simulates a Matern field. This function is borrowed from [tutorials/spde/spde-tutorial-functions.R](https://github.com/grantbrown/inla/blob/master/tutorials/spde/spde-tutorial-functions.R).
   iii. The get.monitoring.stations() function gets the locations of the monitors. We use the same point locations for all the scenarios and replication. Hence, this function will either simulate the locations first if they have not yet been defined then compute the observed values at those location; or it will use a the previously defined locations and simply update the observed values at those locations. This function is adapted from https://github.com/michelacameletti/INLA_COSP. 
    iv. The compute.intersection() function computes the intersection between between each area in the shapefile and the prediction grid and/or simulation grid. The output from this function will be used to compute spatial averages, i.e., too aggregate values from point-level to area-level. This function is adapted from https://github.com/michelacameletti/INLA_COSP. 
    v. The compute.predpoints() function identifies the prediction points in each area in the shapefile. The output from this function will be used to compute spatial averages, i.e., too aggregate values from point-level to area-level. This function is adapted from https://github.com/michelacameletti/INLA_COSP. 

  
2. 
