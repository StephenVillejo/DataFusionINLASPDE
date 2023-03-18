# DataFusionINLASPDE

These are the R codes to implement the method discussed in the paper "Data Fusion in a Two-stage Spatio-Temporal Model using the INLA-SPDE Approach" by Stephen jun Villejo, Janine Illian, and Ben Swallow.

The following provides a description of the different files:

1. DataSimulation.R

  This files contains five important functions:
  
    i. "Data_simulation()" simulates a data from the proposed spatio-temporal model. It gives you the Poisson counts at the area level, the observed values in a network of monitoring stations, a high-resolution data (can be an output of a numerical model or satellite data), and the covariate data.
    ii. The "rspde()" function simulates a Matern field. This function is borrowed from [tutorials/spde/spde-tutorial-functions.R](https://github.com/grantbrown/inla/blob/master/tutorials/spde/spde-tutorial-functions.R).

  
2. 
