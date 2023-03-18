


#' Perform full simulation for a single scenario (FINAL version)
#' 
#' @param x - scenario number
#' @param scenario.settings -dataframe of all scenarios considered in the simulation
#' @param subgroup - run each scenario in different groups simulateneously 
#' @save.plot - TRUE/FALSE, whether to save plots or not
#' 

Simulation.Run.Melding = function(x,
                                  scenario.settings,
                                  subgroup,
                                  save.plot){
  
  scenario_n = x
  
  sparse = scenario.settings[scenario_n,"sparse"] 
  n.T = scenario.settings[scenario_n,"n.T"] 
  
  posterior_b0_compiled.inf = c()
  posterior_b1_compiled.inf = c()
  posterior_a0_compiled.inf = c()
  posterior_a1_compiled.inf = c()
  posterior_sigma2delta_compiled.inf = c()
  posterior_sigma2e_compiled.inf = c()
  posterior_sigma2xi_compiled.inf = c()
  posterior_range_compiled.inf =  c()
  posterior_AR_compiled.inf = c()
  cor_pred_obs_predgrid.inf = c()
  cor_pred_obs_area.inf =
  
  posterior_gamma0_compiled.m1.inf = c()
  posterior_gamma0_compiled.m2.inf = c()
  posterior_gamma1_compiled.m1.inf = c()
  posterior_gamma1_compiled.m2.inf = c()
  posterior_sigma2iid_compiled.m1.inf = c()
  posterior_sigma2iid_compiled.m2.inf = c()
  posterior_sigma2time_compiled.m1.inf = c()
  posterior_sigma2time_compiled.m2.inf = c()
  
  posterior_lincombs_method1.inf = vector("list",length=n_simulations)
  posterior_lincombs_method2.inf = vector("list",length=n_simulations)
  simulated_data_shapefile.inf = vector("list",length=n_simulations)
  
  time.stage1.inf = c()
  time.stage2.simfromppd.inf = c()
  time.stage2.M1.inf = c()
  time.stage2.M2.inf = c()
  
  posterior_b0_compiled.noninf = c()
  posterior_b1_compiled.noninf = c()
  posterior_a0_compiled.noninf = c()
  posterior_a1_compiled.noninf = c()
  posterior_sigma2delta_compiled.noninf = c()
  posterior_sigma2e_compiled.noninf = c()
  posterior_sigma2xi_compiled.noninf = c()
  posterior_range_compiled.noninf =  c()
  posterior_AR_compiled.noninf = c()
  cor_pred_obs_predgrid.noninf = c()
  cor_pred_obs_area.noninf = c()
  
  posterior_gamma0_compiled.m1.noninf = c()
  posterior_gamma0_compiled.m2.noninf = c()
  posterior_gamma1_compiled.m1.noninf = c()
  posterior_gamma1_compiled.m2.noninf = c()
  posterior_sigma2iid_compiled.m1.noninf = c()
  posterior_sigma2iid_compiled.m2.noninf = c()
  posterior_sigma2time_compiled.m1.noninf = c()
  posterior_sigma2time_compiled.m2.noninf = c()
  
  posterior_lincombs_method1.noninf = vector("list",length=n_simulations)
  posterior_lincombs_method2.noninf = vector("list",length=n_simulations)
  simulated_data_shapefile.noninf = vector("list",length=n_simulations)
  
  time.stage1.noninf = c()
  time.stage2.simfromppd.noninf = c()
  time.stage2.M1.noninf = c()
  time.stage2.M2.noninf = c()
  
  failure.informative = rep(0,n_simulations)
  failure.noninformative = rep(0,n_simulations)
  
  #sink(file = paste0(getwd(),'/Update.of.Sim.Scenario_n=',scenario_n,'.subgroup=',subgroup,'.txt'))
  #paste("Let's begin!")
  
  for (sim_i in 1:n_simulations){
    
    cat("Scenario_n = ",scenario_n,", Sim_i = ", sim_i," out of ",n_simulations, sep="","\n" )
    cat("Scenario_n = ",scenario_n,", Sim_i = ", sim_i," out of ",n_simulations," ----", " doing data simulation",sep="","\n" )
    #sink(file = paste0(getwd(),'/Update.of.Sim.Scenario_n=',scenario_n,'.subgroup=',subgroup,'.txt'), append = TRUE)
    #paste0("Scenario_n = ",scenario_n,", Sim_i = ", sim_i," out of ",n_simulations," ----", " doing data simulation")
    
    simulated.data <- Data_simulation(scenario_n = scenario_n,
                                      save.plot = save.plot,
                                      res_pred_grid = res_pred_grid,
                                      highres_simulation_reggrid = highres_simulation_reggrid,
                                      prop.points.for.area = prop.points.for.area,
                                      sigma2xi = sigma2xi,
                                      kappa = kappa,
                                      b0 = b0,
                                      b1 = b1,
                                      sigma2_e = sigma2_e,
                                      gamma0_Poisson = gamma0_Poisson,
                                      gamma1_Poisson = gamma1_Poisson,
                                      sigma2iid_Poisson = sigma2iid_Poisson,
                                      highres.error.sigma2 = highres.error.sigma2,
                                      highres.b0 = highres.b0,
                                      highres.b1 = highres.b1,
                                      n.T = n.T,
                                      a = a,
                                      unif.pop.a = unif.pop.a,
                                      unif.pop.b = unif.pop.b,
                                      sigma2iid_time = sigma2iid_time,
                                      sparse = sparse,
                                      sim_i = sim_i)
    
    #sink(file = paste0(getwd(),'/Update.of.Sim.Scenario_n=',scenario_n,'.subgroup=',subgroup,'.txt'), append = TRUE)
    #paste0("Scenario_n = ",scenario_n,", Sim_i = ", sim_i," out of ",n_simulations," ----", " doing stage 1 model fitting")
    cat("Scenario_n = ",scenario_n,", Sim_i = ", sim_i," out of ",n_simulations," ----", " doing stage 1 model fitting, priors = informative",sep="","\n" )
    
    stage1.res.informative <- tryCatch({
      stage1.fit.ST.final.Melding(inputdata = simulated.data,
                                  n.T = n.T,
                                  save.plot = save.plot,
                                  scenario_n = scenario_n,
                                  informative = TRUE,
                                  two.stage.process = TRUE)
    }, warning = function(w) {
      failure.informative[sim_i] <- 1 
      failure.noninformative[sim_i] <- 1 
      message("Encountered a warning in the stage 1 model (informative case) fitting!")
      next
    }, error = function(e) {
      failure.informative[sim_i] <- 1 
      failure.noninformative[sim_i] <- 1 
      message("Encountered an error in the stage 1 model (informative case) fitting!")
      next
    })
    
    if(any(is.na(stage1.res.informative$res$summary.hyperpar$mean)) == T |
       any(is.na(stage1.res.informative$res$summary.hyperpar$sd)) == T |
       any(is.na(stage1.res.informative$res$summary.hyperpar$mode)) == T |
       class(try(inla.rmarginal(2, stage1.res.informative$res$marginals.fixed$Intercept))) == "try-error" | 
       class(try(inla.rmarginal(2, stage1.res.informative$res$marginals.fixed$x))) == "try-error" |
       class(try(inla.rmarginal(2, stage1.res.informative$res$marginals.fixed$alpha))) == "try-error" |
       class(try(inla.rmarginal(2, stage1.res.informative$res$marginals.hyperpar$`Beta for z.highres`))) == "try-error" |
       class(try(inla.rmarginal(2,inla.tmarginal(function(x) 1/x, stage1.res.informative$res$marginals.hyperpar$`Precision for the Gaussian observations[2]`)))) == "try-error" |
       class(try(inla.rmarginal(2,inla.tmarginal(function(x) 1/x, stage1.res.informative$res$marginals.hyperpar$`Precision for the Gaussian observations`)))) == "try-error" |
       class(try(inla.rmarginal(2,inla.tmarginal(function(x) x^2, stage1.res.informative$res$marginals.hyperpar$`Stdev for spatial.field`)))) == "try-error" |
       class(try(inla.rmarginal(2, stage1.res.informative$res$marginals.hyperpar$`Range for spatial.field`))) == "try-error" |
       class(try(inla.rmarginal(2, stage1.res.informative$res$marginals.hyperpar$`GroupRho for spatial.field`))) == "try-error"){
      
      failure.informative[sim_i] <- 1
      failure.noninformative[sim_i] <- 1
      
      next
      
    }
    
    
    cat("Scenario_n = ",scenario_n,", Sim_i = ", sim_i," out of ",n_simulations," ----", " doing stage 1 model fitting, priors = non-informative",sep="","\n" )
    stage1.res.noninformative <- tryCatch({
      stage1.fit.ST.final.Melding(inputdata = simulated.data,
                                  n.T = n.T,
                                  save.plot = save.plot,
                                  scenario_n = scenario_n,
                                  informative = FALSE,
                                  two.stage.process = TRUE)
    }, warning = function(w) {
      failure.informative[sim_i] <- 1 
      failure.noninformative[sim_i] <- 1 
      message("Encountered a warning in the stage 1 model (non-informative case) fitting!")
      next
    }, error = function(e) {
      failure.informative[sim_i] <- 1 
      failure.noninformative[sim_i] <- 1 
      message("Encountered an error in the stage 1 model (non-informative case) fitting!")
      next
    })
    
    
    if(any(is.na(stage1.res.noninformative$res$summary.hyperpar$mean)) == T |
       any(is.na(stage1.res.noninformative$res$summary.hyperpar$sd)) == T |
       any(is.na(stage1.res.noninformative$res$summary.hyperpar$mode)) == T |
       class(try(inla.rmarginal(2, stage1.res.noninformative$res$marginals.fixed$Intercept))) == "try-error" | 
       class(try(inla.rmarginal(2, stage1.res.noninformative$res$marginals.fixed$x))) == "try-error" |
       class(try(inla.rmarginal(2, stage1.res.noninformative$res$marginals.fixed$alpha))) == "try-error" |
       class(try(inla.rmarginal(2, stage1.res.noninformative$res$marginals.hyperpar$`Beta for z.highres`))) == "try-error" |
       class(try(inla.rmarginal(2,inla.tmarginal(function(x) 1/x, stage1.res.noninformative$res$marginals.hyperpar$`Precision for the Gaussian observations[2]`)))) == "try-error" |
       class(try(inla.rmarginal(2,inla.tmarginal(function(x) 1/x, stage1.res.noninformative$res$marginals.hyperpar$`Precision for the Gaussian observations`)))) == "try-error" |
       class(try(inla.rmarginal(2,stage1.res.noninformative$transformed.params.res$marginals.variance.nominal$variance.nominal.1))) == "try-error" |
       class(try(inla.rmarginal(2,stage1.res.noninformative$transformed.params.res$marginals.range.nominal$range.nominal.1))) == "try-error"  |
       class(try(inla.rmarginal(2, stage1.res.noninformative$res$marginals.hyperpar$`GroupRho for spatial.field`))) == "try-error"){
      
      failure.informative[sim_i] <- 1
      failure.noninformative[sim_i] <- 1
      
      next
      
    }

    cat("Encountered no problems with stage 1 model fitting!",sep="","\n" )
    
    
    #sink(file = paste0(getwd(),'/Update.of.Sim.Scenario_n=',scenario_n,'.subgroup=',subgroup,'.txt'), append = TRUE)
    #paste0("Scenario_n = ",scenario_n,", Sim_i = ", sim_i," out of ",n_simulations," ----", " doing stage 2 model fitting")
    cat("Scenario_n = ",scenario_n,", Sim_i = ", sim_i," out of ",n_simulations," ----", " doing stage 2 model fitting, priors = informative",sep="","\n" )
    stage2.res.informative <- stage2.fit.ST(scenario_n = scenario_n,
                                            inputdata = stage1.res.informative,
                                            simulated.data = simulated.data,
                                            n_random_fromPPD = n_random_fromPPD,
                                            n_random_fromposterior = n_random_fromposterior,
                                            save.plot = save.plot,
                                            n.T = n.T,
                                            informative = TRUE)
    cat("Scenario_n = ",scenario_n,", Sim_i = ", sim_i," out of ",n_simulations," ----", " doing stage 2 model fitting, priors = non-informative",sep="","\n" )
    stage2.res.noninformative <-  stage2.fit.ST(scenario_n = scenario_n,
                                                inputdata = stage1.res.noninformative,
                                                simulated.data = simulated.data,
                                                n_random_fromPPD = n_random_fromPPD,
                                                n_random_fromposterior = n_random_fromposterior,
                                                save.plot = save.plot,
                                                n.T = n.T,
                                                informative = FALSE)
    
    
    
    if(stage2.res.informative$failure.m1 == TRUE | stage2.res.informative$failure.m2 == TRUE |
       stage2.res.noninformative$failure.m1 == TRUE | stage2.res.noninformative$failure.m2 == TRUE){
      
      failure.informative[sim_i] <- 1 
      failure.noninformative[sim_i] <- 1
      
      next
      
    }
    
    
    cat("Encountered no problems with stage 2 model fitting!",sep="","\n" )
    
    
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
    ##############################################################################################
    # SUMMARISE RESULTS FOR THE CASE OF INFORMATIVE PRIORS
    ##############################################################################################
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
    
    ##############################################################################################
    # Compile simulated block-level exposures from the PPD of stage 1 model (informative priors)
    ##############################################################################################
    
    temp_method1 = do.call(cbind,lapply(stage2.res.informative$list_samples_joint_posterior_method1,
                                        function(x) x$sample_joint_post_exposure_method1))
    colnames(temp_method1) = paste("Sim",seq(1,n_random_fromPPD))
    posterior_lincombs_method1.inf[[sim_i]] <- temp_method1
    
    temp_method2 = do.call(cbind,lapply(stage2.res.informative$list_samples_joint_posterior_method2,
                                        function(x) x$sample_joint_post_exposure_method2))
    colnames(temp_method2) = paste("Sim",seq(1,n_random_fromPPD))
    posterior_lincombs_method2.inf[[sim_i]] <- temp_method2
    
    simulated_data_shapefile.inf[[sim_i]] = simulated.data$shapefile_sim_complete@data
    
    ########################################################
    # Compile the stage 1 model parameter estimates
    ########################################################
    
    posterior_b0_compiled.inf <- c(posterior_b0_compiled.inf,
                                   inla.rmarginal(n_random_fromposterior, stage1.res.informative$res$marginals.fixed$Intercept))
    posterior_b1_compiled.inf <- c(posterior_b1_compiled.inf,
                                   inla.rmarginal(n_random_fromposterior, stage1.res.informative$res$marginals.fixed$x))
    posterior_a0_compiled.inf <- c(posterior_a0_compiled.inf,
                                   inla.rmarginal(n_random_fromposterior, stage1.res.informative$res$marginals.fixed$alpha))
    posterior_a1_compiled.inf <- c(posterior_a1_compiled.inf,
                                   inla.rmarginal(n_random_fromposterior, stage1.res.informative$res$marginals.hyperpar$`Beta for z.highres`))
    trans <- inla.tmarginal(function(x) 1/x, stage1.res.informative$res$marginals.hyperpar$`Precision for the Gaussian observations[2]`)
    posterior_sigma2delta_compiled.inf <- c(posterior_sigma2delta_compiled.inf,
                                            inla.rmarginal(n_random_fromposterior, trans))
    trans <- inla.tmarginal(function(x) 1/x, stage1.res.informative$res$marginals.hyperpar$`Precision for the Gaussian observations`)
    posterior_sigma2e_compiled.inf <- c(posterior_sigma2e_compiled.inf,
                                        inla.rmarginal(n_random_fromposterior, trans))
    posterior_AR_compiled.inf <- c(posterior_AR_compiled.inf,
                                   inla.rmarginal(n_random_fromposterior, stage1.res.informative$res$marginals.hyperpar$`GroupRho for spatial.field`))
    trans <- inla.tmarginal(function(x) x^2, stage1.res.informative$res$marginals.hyperpar$`Stdev for spatial.field`)
    posterior_sigma2xi_compiled.inf <- c(posterior_sigma2xi_compiled.inf,
                                         inla.rmarginal(n_random_fromposterior, trans))
    posterior_range_compiled.inf <- c(posterior_range_compiled.inf,
                                      inla.rmarginal(n_random_fromposterior, stage1.res.informative$res$marginals.hyperpar$`Range for spatial.field`))
    
    
    ########################################################
    # Compile the stage 2 model parameter estimates
    ########################################################
    
    posterior_gamma0_compiled.m1.inf = c(posterior_gamma0_compiled.m1.inf, stage2.res.informative$posterior_gamma0_FF_method1)
    posterior_gamma0_compiled.m2.inf = c(posterior_gamma0_compiled.m2.inf, stage2.res.informative$posterior_gamma0_FF_method2)
    posterior_gamma1_compiled.m1.inf = c(posterior_gamma1_compiled.m1.inf, stage2.res.informative$posterior_gamma1_FF_method1)
    posterior_gamma1_compiled.m2.inf = c(posterior_gamma1_compiled.m2.inf, stage2.res.informative$posterior_gamma1_FF_method2)
    posterior_sigma2iid_compiled.m1.inf = c(posterior_sigma2iid_compiled.m1.inf, stage2.res.informative$posterior_sigma2iid_FF_method1)
    posterior_sigma2iid_compiled.m2.inf = c(posterior_sigma2iid_compiled.m2.inf, stage2.res.informative$posterior_sigma2iid_FF_method2)
    posterior_sigma2time_compiled.m1.inf = c(posterior_sigma2time_compiled.m1.inf, stage2.res.informative$posterior_sigma2time_FF_method1)
    posterior_sigma2time_compiled.m2.inf = c(posterior_sigma2time_compiled.m2.inf, stage2.res.informative$posterior_sigma2time_FF_method2)
    
    #######################################################################################################
    # Compile the correlations of predicted exposures and true exposures (both at the grid and block level)
    #######################################################################################################
    
    index_pred <- c()
    for(i in 1:n.T){
      temp <- inla.stack.index(stage1.res.informative$stk.full, paste0("pred.",i,sep=""))$data
      index_pred <- c(index_pred, temp)
    }
    cor_pred_obs_predgrid.inf[sim_i] <- cor(simulated.data$compile.sim_grid_df$exposure,
                                            stage1.res.informative$res$summary.linear.predictor$`0.5quant`[index_pred])
    cor_pred_obs_area.inf[sim_i] <- cor(simulated.data$shapefile_sim_complete@data$true_area_exposure,
                                        stage1.res.informative$res$summary.lincomb.derived$mean)

    time.stage1.inf[sim_i] <- stage1.res.informative$time
    time.stage2.simfromppd.inf[sim_i] <- stage2.res.informative$time.simfromppd
    time.stage2.M1.inf[sim_i] <- stage2.res.informative$time.M1
    time.stage2.M2.inf[sim_i] <- stage2.res.informative$time.M2
    
    
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
    ##############################################################################################
    # SUMMARISE RESULTS FOR THE CASE OF NON-INFORMATIVE PRIORS
    ##############################################################################################
    #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^#
    
    ################################################################################################
    # Compile simulated block-level exposures from the PPD of stage 1 model (non-informative priors)
    ################################################################################################
    
    temp_method1 = do.call(cbind,lapply(stage2.res.noninformative$list_samples_joint_posterior_method1,
                                        function(x) x$sample_joint_post_exposure_method1))
    colnames(temp_method1) = paste("Sim",seq(1,n_random_fromPPD))
    posterior_lincombs_method1.noninf[[sim_i]] <- temp_method1
    
    temp_method2 = do.call(cbind,lapply(stage2.res.noninformative$list_samples_joint_posterior_method2,
                                        function(x) x$sample_joint_post_exposure_method2))
    colnames(temp_method2) = paste("Sim",seq(1,n_random_fromPPD))
    posterior_lincombs_method2.noninf[[sim_i]] <- temp_method2
    
    simulated_data_shapefile.noninf[[sim_i]] = simulated.data$shapefile_sim_complete@data
    
    ########################################################
    # Compile the stage 1 model parameter estimates
    ########################################################
    
    posterior_b0_compiled.noninf <- c(posterior_b0_compiled.noninf,
                                      inla.rmarginal(n_random_fromposterior, stage1.res.noninformative$res$marginals.fixed$Intercept))
    posterior_b1_compiled.noninf <- c(posterior_b1_compiled.noninf,
                                      inla.rmarginal(n_random_fromposterior, stage1.res.noninformative$res$marginals.fixed$x))
    posterior_a0_compiled.noninf <- c(posterior_a0_compiled.noninf,
                                      inla.rmarginal(n_random_fromposterior, stage1.res.noninformative$res$marginals.fixed$alpha))
    posterior_a1_compiled.noninf <- c(posterior_a1_compiled.noninf,
                                      inla.rmarginal(n_random_fromposterior, stage1.res.noninformative$res$marginals.hyperpar$`Beta for z.highres`))
    trans <- inla.tmarginal(function(x) 1/x, stage1.res.noninformative$res$marginals.hyperpar$`Precision for the Gaussian observations[2]`)
    posterior_sigma2delta_compiled.noninf <- c(posterior_sigma2delta_compiled.noninf,
                                               inla.rmarginal(n_random_fromposterior, trans))
    trans <- inla.tmarginal(function(x) 1/x, stage1.res.noninformative$res$marginals.hyperpar$`Precision for the Gaussian observations`)
    posterior_sigma2e_compiled.noninf <- c(posterior_sigma2e_compiled.noninf,
                                           inla.rmarginal(n_random_fromposterior, trans))
    posterior_AR_compiled.noninf <- c(posterior_AR_compiled.noninf,
                                      inla.rmarginal(n_random_fromposterior, stage1.res.noninformative$res$marginals.hyperpar$`GroupRho for spatial.field`))
    posterior_sigma2xi_compiled.noninf <- c(posterior_sigma2xi_compiled.noninf,
                                            inla.rmarginal(n_random_fromposterior, stage1.res.noninformative$transformed.params.res$marginals.variance.nominal$variance.nominal.1))
    posterior_range_compiled.noninf <- c(posterior_range_compiled.noninf,
                                         inla.rmarginal(n_random_fromposterior, stage1.res.noninformative$transformed.params.res$marginals.range.nominal$range.nominal.1))
    
    
    ########################################################
    # Compile the stage 2 model parameter estimates
    ########################################################
    
    posterior_gamma0_compiled.m1.noninf = c(posterior_gamma0_compiled.m1.noninf, stage2.res.noninformative$posterior_gamma0_FF_method1)
    posterior_gamma0_compiled.m2.noninf = c(posterior_gamma0_compiled.m2.noninf, stage2.res.noninformative$posterior_gamma0_FF_method2)
    posterior_gamma1_compiled.m1.noninf = c(posterior_gamma1_compiled.m1.noninf, stage2.res.noninformative$posterior_gamma1_FF_method1)
    posterior_gamma1_compiled.m2.noninf = c(posterior_gamma1_compiled.m2.noninf, stage2.res.noninformative$posterior_gamma1_FF_method2)
    posterior_sigma2iid_compiled.m1.noninf = c(posterior_sigma2iid_compiled.m1.noninf, stage2.res.noninformative$posterior_sigma2iid_FF_method1)
    posterior_sigma2iid_compiled.m2.noninf = c(posterior_sigma2iid_compiled.m2.noninf, stage2.res.noninformative$posterior_sigma2iid_FF_method2)
    posterior_sigma2time_compiled.m1.noninf = c(posterior_sigma2time_compiled.m1.noninf, stage2.res.noninformative$posterior_sigma2time_FF_method1)
    posterior_sigma2time_compiled.m2.noninf = c(posterior_sigma2time_compiled.m2.noninf, stage2.res.noninformative$posterior_sigma2time_FF_method2)
    
    
    #######################################################################################################
    # Compile the correlations of predicted exposures and true exposures (both at the grid and block level)
    #######################################################################################################
    
    index_pred <- c()
    for(i in 1:n.T){
      temp <- inla.stack.index(stage1.res.noninformative$stk.full, paste0("pred.",i,sep=""))$data
      index_pred <- c(index_pred, temp)
    }
    
    cor_pred_obs_predgrid.noninf[sim_i] <- cor(simulated.data$compile.sim_grid_df$exposure,
                                               stage1.res.noninformative$res$summary.linear.predictor$`0.5quant`[index_pred])
    cor_pred_obs_area.noninf[sim_i] <- cor(simulated.data$shapefile_sim_complete@data$true_area_exposure,
                                           stage1.res.noninformative$res$summary.lincomb.derived$mean)

    
    time.stage1.noninf[sim_i] <- stage1.res.noninformative$time
    time.stage2.simfromppd.noninf[sim_i] <- stage2.res.noninformative$time.simfromppd
    time.stage2.M1.noninf[sim_i] <- stage2.res.noninformative$time.M1
    time.stage2.M2.noninf[sim_i] <- stage2.res.noninformative$time.M2
    
    
  }
  
  
  cat("DONE with scenario_n = ",scenario_n," ----"," completed ",n_simulations," simulations", " ---- now saving the samples from the marginal posterior distributions",sep="","\n" )
  
  #sink(file =  paste0(getwd(),'/Update.of.Sim.Scenario_n=',scenario_n,'.subgroup=',subgroup,'.txt'), append = TRUE)
  #paste0("DONE with scenario_n = ",scenario_n," ----"," completed ",n_simulations," simulations", " ---- now compiling results")
  
  posterior_b0_compiled = posterior_b0_compiled.inf  
  posterior_b1_compiled = posterior_b1_compiled.inf  
  posterior_a0_compiled = posterior_a0_compiled.inf 
  posterior_a1_compiled = posterior_a1_compiled.inf 
  posterior_sigma2delta_compiled = posterior_sigma2delta_compiled.inf 
  posterior_sigma2e_compiled = posterior_sigma2e_compiled.inf 
  posterior_sigma2xi_compiled = posterior_sigma2xi_compiled.inf 
  posterior_range_compiled = posterior_range_compiled.inf 
  posterior_AR_compiled = posterior_AR_compiled.inf 
  cor_pred_obs_predgrid = cor_pred_obs_predgrid.inf
  cor_pred_obs_area = cor_pred_obs_area.inf 
  
  posterior_gamma0_compiled.m1 = posterior_gamma0_compiled.m1.inf 
  posterior_gamma0_compiled.m2 = posterior_gamma0_compiled.m2.inf
  posterior_gamma1_compiled.m1 = posterior_gamma1_compiled.m1.inf
  posterior_gamma1_compiled.m2 = posterior_gamma1_compiled.m2.inf 
  posterior_sigma2iid_compiled.m1 = posterior_sigma2iid_compiled.m1.inf 
  posterior_sigma2iid_compiled.m2 = posterior_sigma2iid_compiled.m2.inf 
  posterior_sigma2time_compiled.m1 = posterior_sigma2time_compiled.m1.inf 
  posterior_sigma2time_compiled.m2 = posterior_sigma2time_compiled.m2.inf 
  
  posterior_lincombs_method1 = posterior_lincombs_method1.inf 
  posterior_lincombs_method2 = posterior_lincombs_method2.inf
  simulated_data_shapefile = simulated_data_shapefile.inf 
  
  time.stage1 = time.stage1.inf
  time.stage2.simfromppd = time.stage2.simfromppd.inf
  time.stage2.M1 = time.stage2.M1.inf
  time.stage2.M2 = time.stage2.M2.inf
  
  save(posterior_gamma0_compiled.m1, posterior_gamma0_compiled.m2,
       posterior_gamma1_compiled.m1, posterior_gamma1_compiled.m2,
       posterior_sigma2iid_compiled.m1, posterior_sigma2iid_compiled.m2,
       posterior_sigma2time_compiled.m1, posterior_sigma2time_compiled.m2,
       posterior_b0_compiled, posterior_b1_compiled,
       posterior_a0_compiled, posterior_a1_compiled,
       posterior_sigma2delta_compiled, posterior_sigma2e_compiled,
       posterior_sigma2xi_compiled, posterior_range_compiled, posterior_AR_compiled,
       posterior_lincombs_method1,
       posterior_lincombs_method2,
       simulated_data_shapefile,
       cor_pred_obs_predgrid, 
       cor_pred_obs_area, 
       failure.informative,
       failure.noninformative,
       time.stage1,
       time.stage2.simfromppd,
       time.stage2.M1,
       time.stage2.M2,
       file=paste0(getwd(),"/FinalOutput/","AllposteriorsANDsimulatedData_scenarioN",scenario_n,".subgroup=",subgroup,"(InformativePriors).RData"))
  
  
  posterior_b0_compiled = posterior_b0_compiled.noninf  
  posterior_b1_compiled = posterior_b1_compiled.noninf  
  posterior_a0_compiled = posterior_a0_compiled.noninf 
  posterior_a1_compiled = posterior_a1_compiled.noninf 
  posterior_sigma2delta_compiled = posterior_sigma2delta_compiled.noninf 
  posterior_sigma2e_compiled = posterior_sigma2e_compiled.noninf 
  posterior_sigma2xi_compiled = posterior_sigma2xi_compiled.noninf 
  posterior_range_compiled = posterior_range_compiled.noninf 
  posterior_AR_compiled = posterior_AR_compiled.noninf 
  cor_pred_obs_predgrid = cor_pred_obs_predgrid.noninf
  cor_pred_obs_area = cor_pred_obs_area.noninf 
  
  posterior_gamma0_compiled.m1 = posterior_gamma0_compiled.m1.noninf 
  posterior_gamma0_compiled.m2 = posterior_gamma0_compiled.m2.noninf
  posterior_gamma1_compiled.m1 = posterior_gamma1_compiled.m1.noninf
  posterior_gamma1_compiled.m2 = posterior_gamma1_compiled.m2.noninf 
  posterior_sigma2iid_compiled.m1 = posterior_sigma2iid_compiled.m1.noninf 
  posterior_sigma2iid_compiled.m2 = posterior_sigma2iid_compiled.m2.noninf 
  posterior_sigma2time_compiled.m1 = posterior_sigma2time_compiled.m1.noninf 
  posterior_sigma2time_compiled.m2 = posterior_sigma2time_compiled.m2.noninf 
  
  posterior_lincombs_method1 = posterior_lincombs_method1.noninf 
  posterior_lincombs_method2 = posterior_lincombs_method2.noninf
  simulated_data_shapefile = simulated_data_shapefile.noninf 
  
  time.stage1 = time.stage1.noninf
  time.stage2.simfromppd = time.stage2.simfromppd.noninf
  time.stage2.M1 = time.stage2.M1.noninf
  time.stage2.M2 = time.stage2.M2.noninf
  
  
  save(posterior_gamma0_compiled.m1, posterior_gamma0_compiled.m2,
       posterior_gamma1_compiled.m1, posterior_gamma1_compiled.m2,
       posterior_sigma2iid_compiled.m1, posterior_sigma2iid_compiled.m2,
       posterior_sigma2time_compiled.m1, posterior_sigma2time_compiled.m2,
       posterior_b0_compiled, posterior_b1_compiled,
       posterior_a0_compiled, posterior_a1_compiled,
       posterior_sigma2delta_compiled, posterior_sigma2e_compiled,
       posterior_sigma2xi_compiled, posterior_range_compiled, posterior_AR_compiled,
       posterior_lincombs_method1,
       posterior_lincombs_method2,
       simulated_data_shapefile,
       cor_pred_obs_predgrid,
       cor_pred_obs_area,
       failure.informative,
       failure.noninformative,
       time.stage1,
       time.stage2.simfromppd,
       time.stage2.M1,
       time.stage2.M2,
       file=paste0(getwd(),"/FinalOutput/","AllposteriorsANDsimulatedData_scenarioN",scenario_n,".subgroup=",subgroup,"(Non-InformativePriors).RData"))
  
  
  cat("DONE saving the samples from the marginal posterior distributions!! :')",sep="","\n" )
  
  
  
  cat("This time, I will try to output the plots of the approximate posterior marginals.....",sep="","\n" )
  cat("However, this might fail if there are problems with the previous INLA outputs.",sep="","\n" )
  cat("Not a major problem though, since you already have the samples from the posterior marginals!",sep="","\n" )
  
  ######################################################################################
  # Plots of posterior distributions, combined for all simulations (informative priors)
  ######################################################################################
  
  posterior_b0_compiled = posterior_b0_compiled.inf  
  posterior_b1_compiled = posterior_b1_compiled.inf  
  posterior_a0_compiled = posterior_a0_compiled.inf 
  posterior_a1_compiled = posterior_a1_compiled.inf 
  posterior_sigma2delta_compiled = posterior_sigma2delta_compiled.inf 
  posterior_sigma2e_compiled = posterior_sigma2e_compiled.inf 
  posterior_sigma2xi_compiled = posterior_sigma2xi_compiled.inf 
  posterior_range_compiled = posterior_range_compiled.inf 
  posterior_AR_compiled = posterior_AR_compiled.inf 
 
  posterior_gamma0_compiled.m1 = posterior_gamma0_compiled.m1.inf 
  posterior_gamma0_compiled.m2 = posterior_gamma0_compiled.m2.inf
  posterior_gamma1_compiled.m1 = posterior_gamma1_compiled.m1.inf
  posterior_gamma1_compiled.m2 = posterior_gamma1_compiled.m2.inf 
  posterior_sigma2iid_compiled.m1 = posterior_sigma2iid_compiled.m1.inf 
  posterior_sigma2iid_compiled.m2 = posterior_sigma2iid_compiled.m2.inf 
  posterior_sigma2time_compiled.m1 = posterior_sigma2time_compiled.m1.inf 
  posterior_sigma2time_compiled.m2 = posterior_sigma2time_compiled.m2.inf 
  
  
  
  b0.plot <- post.dens.plot.stage1.ALLSIM(samples = posterior_b0_compiled,
                                          true.value = b0,
                                          title = expression(paste(hat(pi),"(",beta[0],"|y)",)),
                                          xlab = expression(beta[0]))
  b1.plot <- post.dens.plot.stage1.ALLSIM(samples = posterior_b1_compiled,
                                          true.value = b1,
                                          title = expression(paste(hat(pi),"(",beta[1],"|y)",)),
                                          xlab = expression(beta[1])) 
  a0.plot <- post.dens.plot.stage1.ALLSIM(samples = posterior_a0_compiled,
                                          true.value = highres.b0,
                                          title = expression(paste(hat(pi),"(",alpha[0],"|y)",)),
                                          xlab = expression(alpha[0])) 
  a1.plot <- post.dens.plot.stage1.ALLSIM(samples = posterior_a1_compiled,
                                          true.value = highres.b1,
                                          title = expression(paste(hat(pi),"(",alpha[1],"|y)",)),
                                          xlab = expression(alpha[1]))
  sigma2delta.plot <- post.dens.plot.stage1.ALLSIM(samples = posterior_sigma2delta_compiled,
                                                   true.value = highres.error.sigma2,
                                                   title = expression(paste(hat(pi),"(",sigma[delta]^2,"|y)",)),
                                                   xlab = expression(sigma[delta]^2))
  sigma2e.plot <- post.dens.plot.stage1.ALLSIM(samples = posterior_sigma2e_compiled,
                                               true.value = sigma2_e,
                                               title = expression(paste(hat(pi),"(",sigma[e]^2,"|y)",)),
                                               xlab = expression(sigma[e]^2))
  sigma2xi.plot <- post.dens.plot.stage1.ALLSIM(samples = posterior_sigma2xi_compiled,
                                                true.value = sigma2xi,
                                                title = expression(paste(hat(pi),"(",sigma[omega]^2,"|y)",)),
                                                xlab = expression(sigma[omega]^2))
  range.plot <- post.dens.plot.stage1.ALLSIM(samples = posterior_range_compiled,
                                             true.value = rho,
                                             title = expression(paste(hat(pi),"(",sigma[rho]^2,"|y)",)),
                                             xlab = expression(sigma[rho]^2))
  AR.plot <- post.dens.plot.stage1.ALLSIM(samples = posterior_AR_compiled,
                                          true.value = a,
                                          title = expression(paste(hat(pi),"(",varsigma,"|y)",)),
                                          xlab = expression(varsigma))
  
  png(paste0(path,'Scenario=',scenario_n,'.subgroup=',subgroup,'__ALLSIM(InformativePriors).Stage1.Posterior.Densities.png',sep=""), 
      width=40, height=20, units = 'cm', res = 300)
  print((b0.plot + b1.plot + a0.plot)/
          (a1.plot + sigma2xi.plot + range.plot)/
          (sigma2e.plot + sigma2delta.plot + AR.plot))  
  dev.off()
  
  
  
  gamma0.plot <- combined.post.dens.plot.FF(sample1 = posterior_gamma0_compiled.m1,
                                            sample2 = posterior_gamma0_compiled.m2,
                                            true.value = gamma0_Poisson,
                                            title = expression(paste(hat(pi),"(",gamma[0],"|y)",)),
                                            xlab = expression(gamma[0]))
  gamma1.plot <- combined.post.dens.plot.FF(sample1 = posterior_gamma1_compiled.m1,
                                            sample2 = posterior_gamma1_compiled.m2,
                                            true.value = gamma1_Poisson,
                                            title = expression(paste(hat(pi),"(",gamma[1],"|y)",)),
                                            xlab = expression(gamma[1]))
  sigma2iid.plot <- combined.post.dens.plot.FF(sample1 = posterior_sigma2iid_compiled.m1,
                                               sample2 = posterior_sigma2iid_compiled.m2,
                                               true.value = sigma2iid_Poisson,
                                               title = expression(paste(hat(pi),"(",sigma[phi]^2,"|y)",)),
                                               xlab = expression(sigma[phi]^2))
  sigma2time.plot <- combined.post.dens.plot.FF(sample1 = posterior_sigma2time_compiled.m1,
                                                sample2 = posterior_sigma2time_compiled.m2,
                                                true.value = sigma2iid_time,
                                                title = expression(paste(hat(pi),"(",sigma[nu]^2,"|y)",)),
                                                xlab = expression(sigma[nu]^2))
  
  png(paste0(path,'Scenario=',scenario_n,'.subgroup=',subgroup,'__ALLSIM(InformativePriors).Stage2.Posterior.Densities.png',sep=""), 
      width=40, height=20, units = 'cm', res = 300)
  print((gamma0.plot + gamma1.plot)/
          (sigma2iid.plot + sigma2time.plot))
  dev.off()
  

  
  
  posterior_b0_compiled = posterior_b0_compiled.noninf  
  posterior_b1_compiled = posterior_b1_compiled.noninf  
  posterior_a0_compiled = posterior_a0_compiled.noninf 
  posterior_a1_compiled = posterior_a1_compiled.noninf 
  posterior_sigma2delta_compiled = posterior_sigma2delta_compiled.noninf 
  posterior_sigma2e_compiled = posterior_sigma2e_compiled.noninf 
  posterior_sigma2xi_compiled = posterior_sigma2xi_compiled.noninf 
  posterior_range_compiled = posterior_range_compiled.noninf 
  posterior_AR_compiled = posterior_AR_compiled.noninf 
  
  posterior_gamma0_compiled.m1 = posterior_gamma0_compiled.m1.noninf 
  posterior_gamma0_compiled.m2 = posterior_gamma0_compiled.m2.noninf
  posterior_gamma1_compiled.m1 = posterior_gamma1_compiled.m1.noninf
  posterior_gamma1_compiled.m2 = posterior_gamma1_compiled.m2.noninf 
  posterior_sigma2iid_compiled.m1 = posterior_sigma2iid_compiled.m1.noninf 
  posterior_sigma2iid_compiled.m2 = posterior_sigma2iid_compiled.m2.noninf 
  posterior_sigma2time_compiled.m1 = posterior_sigma2time_compiled.m1.noninf 
  posterior_sigma2time_compiled.m2 = posterior_sigma2time_compiled.m2.noninf 
  
  posterior_lincombs_method1 = posterior_lincombs_method1.noninf 
  posterior_lincombs_method2 = posterior_lincombs_method2.noninf
  simulated_data_shapefile = simulated_data_shapefile.noninf 
  
  
  ######################################################################################
  # Plots of posterior distributions, combined for all simulations (informative priors)
  ######################################################################################
  
  b0.plot <- post.dens.plot.stage1.ALLSIM(samples = posterior_b0_compiled,
                                          true.value = b0,
                                          title = expression(paste(hat(pi),"(",beta[0],"|y)",)),
                                          xlab = expression(beta[0]))
  b1.plot <- post.dens.plot.stage1.ALLSIM(samples = posterior_b1_compiled,
                                          true.value = b1,
                                          title = expression(paste(hat(pi),"(",beta[1],"|y)",)),
                                          xlab = expression(beta[1])) 
  a0.plot <- post.dens.plot.stage1.ALLSIM(samples = posterior_a0_compiled,
                                          true.value = highres.b0,
                                          title = expression(paste(hat(pi),"(",alpha[0],"|y)",)),
                                          xlab = expression(alpha[0])) 
  a1.plot <- post.dens.plot.stage1.ALLSIM(samples = posterior_a1_compiled,
                                          true.value = highres.b1,
                                          title = expression(paste(hat(pi),"(",alpha[1],"|y)",)),
                                          xlab = expression(alpha[1]))
  sigma2delta.plot <- post.dens.plot.stage1.ALLSIM(samples = posterior_sigma2delta_compiled,
                                                   true.value = highres.error.sigma2,
                                                   title = expression(paste(hat(pi),"(",sigma[delta]^2,"|y)",)),
                                                   xlab = expression(sigma[delta]^2))
  sigma2e.plot <- post.dens.plot.stage1.ALLSIM(samples = posterior_sigma2e_compiled,
                                               true.value = sigma2_e,
                                               title = expression(paste(hat(pi),"(",sigma[e]^2,"|y)",)),
                                               xlab = expression(sigma[e]^2))
  sigma2xi.plot <- post.dens.plot.stage1.ALLSIM(samples = posterior_sigma2xi_compiled,
                                                true.value = sigma2xi,
                                                title = expression(paste(hat(pi),"(",sigma[omega]^2,"|y)",)),
                                                xlab = expression(sigma[omega]^2))
  range.plot <- post.dens.plot.stage1.ALLSIM(samples = posterior_range_compiled,
                                             true.value = rho,
                                             title = expression(paste(hat(pi),"(",sigma[rho]^2,"|y)",)),
                                             xlab = expression(sigma[rho]^2))
  AR.plot <- post.dens.plot.stage1.ALLSIM(samples = posterior_AR_compiled,
                                          true.value = a,
                                          title = expression(paste(hat(pi),"(",varsigma,"|y)",)),
                                          xlab = expression(varsigma))
  
  png(paste0(path,'Scenario=',scenario_n,'.subgroup=',subgroup,'__ALLSIM(Non-InformativePriors).Stage1.Posterior.Densities.png',sep=""), 
      width=40, height=20, units = 'cm', res = 300)
  print((b0.plot + b1.plot + a0.plot)/
          (a1.plot + sigma2xi.plot + range.plot)/
          (sigma2e.plot + sigma2delta.plot + AR.plot))  
  dev.off()
  
  gamma0.plot <- combined.post.dens.plot.FF(sample1 = posterior_gamma0_compiled.m1,
                                            sample2 = posterior_gamma0_compiled.m2,
                                            true.value = gamma0_Poisson,
                                            title = expression(paste(hat(pi),"(",gamma[0],"|y)",)),
                                            xlab = expression(gamma[0]))
  gamma1.plot <- combined.post.dens.plot.FF(sample1 = posterior_gamma1_compiled.m1,
                                            sample2 = posterior_gamma1_compiled.m2,
                                            true.value = gamma1_Poisson,
                                            title = expression(paste(hat(pi),"(",gamma[1],"|y)",)),
                                            xlab = expression(gamma[1]))
  sigma2iid.plot <- combined.post.dens.plot.FF(sample1 = posterior_sigma2iid_compiled.m1,
                                               sample2 = posterior_sigma2iid_compiled.m2,
                                               true.value = sigma2iid_Poisson,
                                               title = expression(paste(hat(pi),"(",sigma[phi]^2,"|y)",)),
                                               xlab = expression(sigma[phi]^2))
  sigma2time.plot <- combined.post.dens.plot.FF(sample1 = posterior_sigma2time_compiled.m1,
                                                sample2 = posterior_sigma2time_compiled.m2,
                                                true.value = sigma2iid_time,
                                                title = expression(paste(hat(pi),"(",sigma[nu]^2,"|y)",)),
                                                xlab = expression(sigma[nu]^2))
  
  png(paste0(path,'Scenario=',scenario_n,'.subgroup=',subgroup,'__ALLSIM(Non-InformativePriors).Stage2.Posterior.Densities.png',sep=""), 
      width=40, height=20, units = 'cm', res = 300)
  print((gamma0.plot + gamma1.plot)/
          (sigma2iid.plot + sigma2time.plot))
  dev.off()
  
  
  cat("The creation of the plots is successful!",sep="","\n" )
  
  
}


