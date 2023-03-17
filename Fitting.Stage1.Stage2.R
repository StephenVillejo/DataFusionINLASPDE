

#' Stage 1 Fitting of a spatio-temporal model - Bayesian Melding 
#' 
#' @param inputdata - simulated data
#' @param n.T - number of time points
#' @param save.plot - TRUE/FALSE, if plots are saved or not 
#' @scenario_n - scenario number
#' @param informative - TRUE/FALSE, if informative priors are used 
#' @param two.stage.process - TRUE/FALSE, if a two-stage process which accounts for uncertainty from first stage model will be done
#' 
#' @return - posterior estimates from stage 1 model 

stage1.fit.ST.final.Melding <- function(inputdata = simulated.data,
                                        n.T = 4,
                                        save.plot = TRUE,
                                        scenario_n = scenario_n,
                                        informative = TRUE,
                                        two.stage.process = TRUE){
  
  for (i in 1:n.T){
    inputdata$compile.exposures.highres.list[[i]]$Time = i
  }
  temp <- bind_rows(inputdata$compile.exposures.highres.list, .id = "column_label")
  highres.data <- merge(x = temp, y = inputdata$compile.sim_grid_inside_df[,c('x1','x2','z1','Time')],
                        by = c('x1','x2','Time'))
  
  n_stations <- nrow(inputdata$monitoringstations.coords) 
  n_data <- nrow(inputdata$compile_monitors) 
  n_time <- n_data/n_stations
  
  bnd = unionSpatialPolygons(inputdata$shapefile_sim_complete, rep(1,98))
  boundary <- list(as.inla.mesh.segment(bnd),NULL)
  mesh <- inla.mesh.2d(loc=cbind(inputdata$monitoringstations.coords$x1,inputdata$monitoringstations.coords$x2),
                       boundary = boundary,
                       offset=c(.17, .685),
                       max.edge=c(.15, .5),
                       cutoff= 0.05,
                       min.angle=c(30,21))
  
  if (informative == TRUE){
    spde <- inla.spde2.pcmatern(mesh, prior.range = c(rho,.05),
                                prior.sigma = c(sqrt(sigma2xi),.05))
    s_index <- inla.spde.make.index(name="spatial.field",
                                    n.spde=spde$n.spde,
                                    n.group=n_time)
  }else{
    spde <- inla.spde2.matern(mesh=mesh, alpha=2)
    s_index <- inla.spde.make.index(name="spatial.field",
                                    n.spde=spde$n.spde,
                                    n.group=n_time)
  }
  
  
  
  A.true.exposures <- inla.spde.make.A(mesh=mesh, 
                                       loc=rbind(cbind(inputdata$compile_monitors$x1,inputdata$compile_monitors$x2),cbind(highres.data$x1, highres.data$x2)),
                                       group=c(inputdata$compile_monitors$Time, highres.data$Time),
                                       n.group=n_time)
  
  #--- stack for the measurement error of the monitors
  stk_monitors <- inla.stack(data=list(y=cbind(as.vector(inputdata$compile_monitors$exposure_with_meas_err),NA,NA)),
                             A=list(1),
                             effects=list(
                               data.frame(z.monitors=1:nrow(inputdata$compile_monitors),
                                          z.monitors.weight = rep(1,nrow(inputdata$compile_monitors)))),
                             tag='est.monitors')
  
  #--- stack for the latent exposure field 
  stk_true.exposures <- inla.stack(data=list(y = cbind(NA,NA,rep(0, nrow(highres.data)+length(inputdata$compile_monitors$exposure_with_meas_err)))),
                                   A=list(A.true.exposures, 1),
                                   effects=list(c(s_index,list(Intercept=1)),
                                                data.frame(x = c(inputdata$compile_monitors$z1,highres.data$z1),
                                                           z.true = 1:(nrow(highres.data)+nrow(inputdata$compile_monitors)),
                                                           z.true.weight = rep(-1,nrow(highres.data)+nrow(inputdata$compile_monitors)))),
                                   tag="est.true.exposures")
  
  #--- stack for the high resolution data
  stk_highres <- inla.stack(data = list(y=cbind(NA, highres.data$exposure.highres, NA)),  
                            A = list(1),
                            effects = list(data.frame(alpha = 1,
                                                      z.highres = (nrow(inputdata$compile_monitors)+1):(nrow(highres.data)+nrow(inputdata$compile_monitors)))),
                            tag = "est.highres")                       
  
  #--- stacks for the prediction grid
  for(i in 1:n.T){
    A.pred <- inla.spde.make.A(mesh=mesh, loc=inputdata$pred_grid@coords,
                               group=i,
                               n.group=n_time)
    stk_pred <- inla.stack(data=list(y = cbind(NA,NA,rep(NA,length(inputdata$pred_grid)))),
                           A=list(A.pred,1),
                           effects=list(c(s_index,list(Intercept=1)),
                                        data.frame(x=inputdata$compile.sim_grid_df[which(inputdata$compile.sim_grid_df$Time==i),]$z1)),
                           tag=paste0("pred.",i,sep=""))
    assign(paste0("stk_pred.",i,sep=""),stk_pred)
  }
  
  temp.stk.full <- inla.stack(stk_monitors, stk_highres, stk_true.exposures)
  for(i in 1:n.T){
    temp <- inla.stack(temp.stk.full, get(paste0("stk_pred.",i,sep="")))  
    temp.stk.full <- temp
  }
  stk.full <- temp.stk.full
  
  if (informative == TRUE){
    
    control.family <- list(list(hyper = list(prec = list(initial = 4, fixed = FALSE,
                                                         prior = "loggamma",
                                                         param = c(100,10)))),
                           list(hyper = list(prec = list(initial = 4, fixed = FALSE,
                                                         prior = "loggamma",
                                                         param = c(10,10)))),
                           list(hyper = list(prec = list(initial = 15, fixed = TRUE))))
    
    form <- y ~ -1 + Intercept + x + alpha +
      f(spatial.field, model=spde,
        group=spatial.field.group, 
        control.group=list(model="ar1", hyper = list(theta = list(param=c(log((1+a)/(1-a)),0.05), initial=log((1+a)/(1-a)),fixed=FALSE)))) +
      f(z.monitors, z.monitors.weight, copy='z.true',
        hyper = list(theta = list(initial = 1, fixed = TRUE))) +
      f(z.true, z.true.weight, model='iid', values = 1:(nrow(inputdata$compile_monitors)+nrow(highres.data)),
        hyper = list(theta = list(initial = -5, fixed = TRUE))) +
      f(z.highres, copy='z.true',
        hyper = list(theta = list(initial = highres.b1, fixed = FALSE,
                                  prior = "normal",
                                  param = c(highres.b1, .5))))
  } else{
    
    control.family <- list(list(),
                           list(),
                           list(hyper = list(prec = list(initial = 15, fixed = TRUE))))
    
    form <- y ~ -1 + Intercept + x + alpha +
      f(spatial.field, model=spde,
        group=spatial.field.group, 
        control.group=list(model="ar1")) +
      f(z.monitors, z.monitors.weight, copy='z.true',
        hyper = list(theta = list(initial = 1, fixed = TRUE))) +
      f(z.true, z.true.weight, model='iid', values = 1:(nrow(inputdata$compile_monitors)+nrow(highres.data)),
        hyper = list(theta = list(initial = -5, fixed = TRUE))) +
      f(z.highres, copy='z.true',
        hyper = list(theta = list(initial = highres.b1, fixed = FALSE)))
    
  } 
  
  
  load(paste0(getwd(), 
              "/Data","/weights_points_areas_inside_closest_",res_pred_grid,".Rdata")) 
  load(paste0(getwd(), 
              "/Data","/weights_grid_areas_intersection_res_pred_grid_",res_pred_grid,".Rdata")) 
  
  ##############################################################
  # Define the linear combinations in the INLA style for method 1: intersections
  ##############################################################
  
  index_monitors <- inla.stack.index(stk.full,"est.monitors")$data
  index_highres <- inla.stack.index(stk.full,"est.highres")$data
  index_latent <- inla.stack.index(stk.full,"est.true.exposures")$data 
  for(i in 1:n.T){
    index_pred <- inla.stack.index(stk.full, paste0("pred.",i,sep=""))$data
    assign(paste0("index_pred.",i,sep=""),index_pred)
  }
  
  monitoringstations.coords = inputdata$monitoringstations.coords
  n_areas <- nrow(inputdata$shapefile_sim_complete)/n.T
  dim_lp = length(index_monitors) + length(index_highres) + length(index_latent) + (length(index_pred.1)*n.T)
  #each row is an area, each column is a component of the linear combination
  lincombs.weights_method1 = matrix(0, nrow = n_areas, ncol = length(inputdata$pred_grid) * n.T)  
  
  lc_all_method1 = c() 
  
  for (time.ind in 1:n.T){
    
    for(i in seq(1, length(unique(weights_grid_areas_intersection$AREAKEY)))) { #this corresponds to n_area
      
      index_pred <- get(paste0("index_pred.",time.ind,sep=""))
      #cat("------------- Area n.",i,"----", "\n" )
      pixel_area = weights_grid_areas_intersection[which(weights_grid_areas_intersection$AREAKEY==unique(weights_grid_areas_intersection$AREAKEY)[i]),]
      
      #Extract from the SP_ID string only the number related to the pixel (remove g and other numbers at the end)  
      pixel_area$index = as.numeric(substring(gsub( " .*$", "", pixel_area$pixel ),2))
      pixel_area = pixel_area[order(pixel_area$index),] 
      
      #Fill in the matrix with the weights (sum by row = 1)
      lincombs.weights_method1[i,(index_pred[pixel_area$index] - (length(index_monitors) + length(index_highres) + length(index_latent) ))]  = pixel_area$weights  
      
      lc_vec <- rep(NA,times=dim_lp)
      lc_vec[index_pred][pixel_area$index] <- pixel_area$weights
      lc_predictor <- inla.make.lincomb(APredictor = lc_vec)
      names(lc_predictor) <- paste0("lc_area_ID=",pixel_area$AREAKEY[1],",Time=",time.ind)
      assign(paste0("lc_area_ID=",pixel_area$AREAKEY[1],",Time=",time.ind), lc_predictor)
      
      lc_all_method1 = c(lc_all_method1, get(paste0("lc_area_ID=",pixel_area$AREAKEY[1],",Time=",time.ind)))   
    }
    
  }
  
  
  ##############################################################
  # Define the linear combinations in the INLA style for method 2: points inside/closest
  ##############################################################
  #each row is an area, each column is a component of the linear combination
  lincombs.weights_method2 = matrix(0, nrow = n_areas, ncol = length(inputdata$pred_grid) * n.T)  
  lc_all_method2 = c() 
  
  for (time.ind in 1:n.T){
    
    for(i in seq(1, length(unique(weights_points_areas_inside_closest$AREAKEY)))) { #this corresponds to n_area
      
      index_pred <- get(paste0("index_pred.",time.ind,sep=""))
      
      #cat("------------- Area n.",i,"----", "\n" )
      pixel_area = weights_points_areas_inside_closest[which(weights_points_areas_inside_closest$AREAKEY==unique(weights_points_areas_inside_closest$AREAKEY)[i]),]
      
      #Fill in the matrix with the weights (sum by row = 1)
      pixel_area$index = pixel_area$pixel
      lincombs.weights_method2[i,(index_pred[pixel_area$index] - (length(index_monitors) + length(index_highres) + length(index_latent) ))]  = pixel_area$weights  
      
      lc_vec <- rep(NA,times=dim_lp)
      lc_vec[index_pred][pixel_area$index] <- pixel_area$weights
      lc_predictor <- inla.make.lincomb(APredictor = lc_vec)
      names(lc_predictor) <- paste0("lc_area_ID=",pixel_area$AREAKEY[1],",Time=",time.ind)
      assign(paste0("lc_area_ID=",pixel_area$AREAKEY[1],",Time=",time.ind), lc_predictor)
      
      lc_all_method2 = c(lc_all_method2, get(paste0("lc_area_ID=",pixel_area$AREAKEY[1],",Time=",time.ind)))   
    }
    
  }
  
  
  if(two.stage.process == TRUE){
    
    cat("-----------------------------------------", " now fitting the INLA on stage 1 model",sep="","\n" )
    
    start <- Sys.time()
    
    res <- inla(formula = form,
                family=c("gaussian","gaussian","gaussian"),
                data=inla.stack.data(stk.full),
                control.predictor=list(compute=TRUE,link=1, A=inla.stack.A(stk.full)),
                control.compute = list(dic = TRUE,cpo=TRUE, config=TRUE, waic=TRUE),
                control.inla = list(strategy="adaptive",int.strategy='eb'),
                control.family = control.family,
                lincomb = lc_all_method2,
                verbose = FALSE)
    
    end <- Sys.time()
    
    cat("-----------------------------------------", " done fitting the INLA on stage 1 model",sep="","\n" )
    
  }else{
    
    cat("-----------------------------------------", " now fitting the INLA on stage 1 model using method 1",sep="","\n" )
    
    start <- Sys.time()
    
    res1 <- inla(formula = form,
                 family=c("gaussian","gaussian","gaussian"),
                 data=inla.stack.data(stk.full),
                 control.predictor=list(compute=TRUE,link=1, A=inla.stack.A(stk.full)),
                 control.compute = list(dic = TRUE,cpo=TRUE, config=TRUE, waic=TRUE),
                 control.inla = list(strategy="adaptive",int.strategy='eb'),
                 control.family = control.family,
                 lincomb = lc_all_method1,
                 verbose = FALSE)
    cat("-----------------------------------------", " now fitting the INLA on stage 1 model using method 2",sep="","\n" )
    res2 <- inla(formula = form,
                 family=c("gaussian","gaussian","gaussian"),
                 data=inla.stack.data(stk.full),
                 control.predictor=list(compute=TRUE,link=1, A=inla.stack.A(stk.full)),
                 control.compute = list(dic = TRUE,cpo=TRUE, config=TRUE, waic=TRUE),
                 control.inla = list(strategy="adaptive",int.strategy='eb'),
                 control.family = control.family,
                 lincomb = lc_all_method2,
                 verbose = FALSE)
    
    end <- Sys.time()
    
    cat("-----------------------------------------", " done fitting the INLA on stage 1 model",sep="","\n" )
    
  }
  
  
  
  
  
  if (save.plot == TRUE & two.stage.process == FALSE){
    
    b0.plot <- post.dens.plot.stage1(res1.param = res1$marginals.fixed$Intercept,
                                     res2.param = res2$marginals.fixed$Intercept,
                                     true.value = b0,
                                     title = expression(paste(hat(pi),"(",beta[0],"|y)",)),
                                     xlab = expression(beta[0]),
                                     trans = 'reciprocal',
                                     prec.param = FALSE,
                                     n_random_fromposterior = 200)
    
    b1.plot <- post.dens.plot.stage1(res1.param = res1$marginals.fixed$x,
                                     res2.param = res2$marginals.fixed$x,
                                     true.value = b1,
                                     title = expression(paste(hat(pi),"(",beta[1],"|y)",)),
                                     xlab = expression(beta[1]),
                                     trans = 'reciprocal',
                                     prec.param = FALSE,
                                     n_random_fromposterior = 200)
    
    a0.plot <- post.dens.plot.stage1(res1.param = res1$marginals.fixed$alpha,
                                     res2.param = res2$marginals.fixed$alpha,
                                     true.value = highres.b0,
                                     title = expression(paste(hat(pi),"(",alpha[0],"|y)",)),
                                     xlab = expression(alpha[0]),
                                     trans = 'reciprocal',
                                     prec.param = FALSE,
                                     n_random_fromposterior = 200)
    
    if (informative == TRUE){
      sigma2xi.plot <- post.dens.plot.stage1(res1.param = res1$marginals.hyperpar$`Stdev for spatial.field`,
                                             res2.param = res2$marginals.hyperpar$`Stdev for spatial.field`,
                                             true.value = sigma2xi,
                                             title = expression(paste(hat(pi),"(",sigma[omega]^2,"|y)",)),
                                             xlab = expression(sigma[omega]^2),
                                             trans = 'squared',
                                             prec.param = TRUE,
                                             n_random_fromposterior = 200)
      
      rho.plot <- post.dens.plot.stage1(res1.param = res1$marginals.hyperpar$`Range for spatial.field`,
                                        res2.param = res2$marginals.hyperpar$`Range for spatial.field`,
                                        true.value = rho,
                                        title = expression(paste(hat(pi),"(",rho,"|y)",)),
                                        xlab = expression(rho),
                                        trans = 'squared',
                                        prec.param = FALSE,
                                        n_random_fromposterior = 200)
    }else{
      transformed.params.res1 <- inla.spde2.result(inla=res1,
                                                   name="spatial.field", spde=spde, do.transf=TRUE)
      transformed.params.res2 <- inla.spde2.result(inla=res2,
                                                   name="spatial.field", spde=spde, do.transf=TRUE)
      sigma2xi.plot <- post.dens.plot.stage1(res1.param = transformed.params.res1$marginals.variance.nominal$variance.nominal.1,
                                             res2.param =transformed.params.res2$marginals.variance.nominal$variance.nominal.1,
                                             true.value = sigma2xi,
                                             title = expression(paste(hat(pi),"(",sigma[omega]^2,"|y)",)),
                                             xlab = expression(sigma[omega]^2),
                                             trans = 'reciprocal',
                                             prec.param = FALSE,
                                             n_random_fromposterior = 200)
      rho.plot <- post.dens.plot.stage1(res1.param = transformed.params.res1$marginals.range.nominal$range.nominal.1,
                                        res2.param = transformed.params.res2$marginals.range.nominal$range.nominal.1,
                                        true.value = rho,
                                        title = expression(paste(hat(pi),"(",rho,"|y)",)),
                                        xlab = expression(rho),
                                        trans = 'squared',
                                        prec.param = FALSE,
                                        n_random_fromposterior = 200)
    }
    
    
    sigma2_e.plot <- post.dens.plot.stage1(res1.param = res1$marginals.hyperpar$`Precision for the Gaussian observations`,
                                           res2.param = res2$marginals.hyperpar$`Precision for the Gaussian observations`,
                                           true.value = sigma2_e,
                                           title = expression(paste(hat(pi),"(",sigma[phi]^2,"|y)",)),
                                           xlab = expression(sigma[e]^2),
                                           trans = 'reciprocal',
                                           prec.param = TRUE,
                                           n_random_fromposterior = 200)
    
    sigma2_highres.plot <- post.dens.plot.stage1(res1.param = res1$marginals.hyperpar$`Precision for the Gaussian observations[2]`,
                                                 res2.param = res2$marginals.hyperpar$`Precision for the Gaussian observations[2]`,
                                                 true.value = highres.error.sigma2,
                                                 title = expression(paste(hat(pi),"(",sigma[delta]^2,"|y)",)),
                                                 xlab = expression(sigma[delta]^2),
                                                 trans = 'reciprocal',
                                                 prec.param = TRUE,
                                                 n_random_fromposterior = 200)
    
    a1.plot <- post.dens.plot.stage1(res1.param = res1$marginals.hyperpar$`Beta for z.highres`,
                                     res2.param = res2$marginals.hyperpar$`Beta for z.highres`,
                                     true.value = highres.b1,
                                     title = expression(paste(hat(pi),"(",alpha[1],"|y)",)),
                                     xlab = expression(alpha[1]),
                                     trans = 'nothing',
                                     prec.param = FALSE,
                                     n_random_fromposterior = 200)
    
    AR.param <- post.dens.plot.stage1(res1.param = res1$marginals.hyperpar$`GroupRho for spatial.field`,
                                      res2.param = res2$marginals.hyperpar$`GroupRho for spatial.field`,
                                      true.value = a,
                                      title = expression(paste(hat(pi),"(",varsigma,"|y)",)),
                                      xlab = expression(varsigma),
                                      trans = 'nothing',
                                      prec.param = FALSE,
                                      n_random_fromposterior = 200)
    
    
    png(paste0(path,'Scenario=',scenario_n,'.S1.Posterior.Densities.ST.png',sep=""), width=25, height=25, units = 'cm', res = 300)
    print((b0.plot + b1.plot + a0.plot)/
            (a1.plot + sigma2xi.plot + rho.plot)/
            (sigma2_e.plot + sigma2_highres.plot + AR.param))  
    dev.off()
    
    
    
    if (save.plot == TRUE){
      
      scatterplots <- scatterplot(inla.m1 = res1, inla.m2 = res2, ST = TRUE, n.T = n.T, inputdata = inputdata)
      
      if (n.T == 3){
        ncols = n.T
      }else{
        ncols = n.T/2
      }
      
      png(paste(path,'Scenario=',scenario_n,'.S1.ScatterPlot.M1.ST.png',sep=""), width=50, height=30, units = 'cm', res = 300)
      print(do.call("grid.arrange", c(scatterplots$compile.m1, ncol=ncols)))
      dev.off()
      
      png(paste(path,'Scenario=',scenario_n,'.S1.ScatterPlot.M2.ST.png',sep=""), width=50, height=30, units = 'cm', res = 300)
      print(do.call("grid.arrange", c(scatterplots$compile.m2, ncol=ncols)))
      dev.off()
      
      exp.surface.plots <- plot.exposures.surface.ST(inlamodel1 = res1,
                                                     inlamodel2 = res2,
                                                     stk = stk.full,
                                                     n.T = n.T,
                                                     inputdata = inputdata)
      
      for(time.ind in 1:n.T){
        i <- (time.ind - 1) * 3 + 1
        png(paste0(path,'Scenario=',scenario_n,'.S1.Exposures.Surface.Time=',time.ind,'.png',sep=""), width=30, height=20, units = 'cm', res = 300)
        #print(exp.surface.plots$compile.plots[[i]] + exp.surface.plots$compile.plots[[i+1]]  + exp.surface.plots$compile.plots[[i+2]] )
        print(exp.surface.plots$compile.plots[[i]] + exp.surface.plots$compile.plots[[i+1]] )
        dev.off()
      }
      
      
      #--- Summaries from the posterior distributions of the lincombs:
      lincomb_marginals_method1 = res1$marginals.lincomb.derived
      lincomb_marginals_method2 = res2$marginals.lincomb.derived 
      
      #--- Define prior mean and precision for exposure for each area (taken from the lincombs marginal posterior distributions)
      all_post_mean_method1 = unlist(lapply(lincomb_marginals_method1, function(x) inla.zmarginal(x,silent=T)$mean))
      all_post_sd_method1 = unlist(lapply(lincomb_marginals_method1, function(x) inla.zmarginal(x,silent=T)$sd))
      all_post_mean_method2 = unlist(lapply(lincomb_marginals_method2, function(x) inla.zmarginal(x,silent=T)$mean))
      all_post_sd_method2 = unlist(lapply(lincomb_marginals_method2, function(x) inla.zmarginal(x,silent=T)$sd))
      
      all_vals = c(all_post_mean_method1, all_post_mean_method1, inputdata$shapefile_sim_complete$true_area_exposure)
      
      for (i in 1:n.T){
        
        sub.all_post_mean_method1 <- all_post_mean_method1[(98*(i-1)+1):(98*i)]
        sub.all_post_mean_method2 <- all_post_mean_method2[(98*(i-1)+1):(98*i)]
        
        #--- Prepare the data
        shapefile = inputdata$shapefile_sim_complete
        temp = inputdata$shapefile_sim_complete@data[(98*(i-1)+1):(98*i),]
        shapefile@data <- temp
        forplot <- shapefile
        forplot@data$value1 = sub.all_post_mean_method1
        forplot@data$value2 = sub.all_post_mean_method2
        
        forplot_tidy <- tidy(forplot)
        forplot$id <- row.names(forplot)
        colnames(forplot@data)[which(colnames(forplot@data) == 'true_area_exposure')] <- 'true'
        forplot <- left_join(forplot_tidy, forplot@data)
        
        estimatedcounts.m1 <- ggplot(forplot, aes(x = long, y = lat, group = group, fill = value1)) +
          geom_polygon(color = "black", size = 0.1) +
          coord_equal(ratio=1) +
          labs(x="", y="", title=paste0("Estimated block-level exposures (method 1), Time = ",i,sep="")) +
          theme(legend.key.width=unit(1, "cm")) + 
          theme_map() +
          scale_fill_viridis(option = 'D', limits = c(min(all_vals), max(all_vals))) +
          theme(legend.position="bottom")
        estimatedcounts.m2 <- ggplot(forplot, aes(x = long, y = lat, group = group, fill = value2)) +
          geom_polygon(color = "black", size = 0.1) +
          coord_equal(ratio=1) +
          labs(x="", y="", title=paste0("Estimated block-level exposures (method 2), Time = ",i,sep="")) +
          theme(legend.key.width=unit(1, "cm")) + 
          theme_map() +
          scale_fill_viridis(option = 'D', limits = c(min(all_vals), max(all_vals))) +
          theme(legend.position="bottom")
        true <- ggplot(forplot, aes(x = long, y = lat, group = group, fill = true)) +
          geom_polygon(color = "black", size = 0.1) +
          coord_equal(ratio=1) +
          labs(x="", y="", title=paste0("True block-level exposures, Time = ",i,sep="")) +
          theme(legend.key.width=unit(1, "cm")) + 
          theme_map() +
          scale_fill_viridis(option = 'D', limits = c(min(all_vals), max(all_vals))) +
          theme(legend.position="bottom")
        
        png(paste0(path,'Scenario=',scenario_n,'.S1.Exposures.Blocks.Time=',i,'.png',sep=""), width=30, height=20, units = 'cm', res = 300)
        print(true + estimatedcounts.m1 + estimatedcounts.m2)
        dev.off()
        
      }
      
      index_monitors <- inla.stack.index(stk.full,"est.monitors")$data
      temp1 <- res1$summary.linear.predictor[index_monitors,'mean']
      temp2 <- res2$summary.linear.predictor[index_monitors,'mean']
      
      for.plot <- inputdata$compile_monitors
      for.plot$method.1 <- temp1
      for.plot$method.2 <- temp2
      
      for.plot <- for.plot[,c('exposure_with_meas_err','Time','monitor','exposure','method.1','method.2')]
      names(for.plot)[1] <- "observed"
      names(for.plot)[4] <- "true value"
      
      
      for (i in 1:4){
        
        sub.for.plot <- subset(for.plot, subset=(monitor==paste0("M",i,sep="")))
        meltdf <- melt(sub.for.plot,id.vars=c("Time","monitor"),
                       measure.vars=c('observed','true value','method.1','method.2'))
        png(paste0(path,'Scenario=',scenario_n,".S1.Preds.Monitor=M",i,".png",sep=""), width=20, height=20, units = 'cm', res = 300)
        print(ggplot(meltdf,aes(x=Time,y=value,colour=variable,group=variable)) + geom_line() + 
                ggtitle(paste0("Monitor = M",i,sep=""))+geom_point())
        dev.off()
        
      }
      
    }
  }
  
  if(two.stage.process == FALSE){
    
    if (informative == TRUE){
      transformed.params.res1 = c()
      transformed.params.res2 = c()
    }else{
      transformed.params.res1 <- inla.spde2.result(inla=res1,
                                                   name="spatial.field", spde=spde, do.transf=TRUE)
      transformed.params.res2 <- inla.spde2.result(inla=res2,
                                                   name="spatial.field", spde=spde, do.transf=TRUE)
    }
    
  }else{
    
    if (informative == TRUE){
      transformed.params.res = c()
    }else{
      transformed.params.res = inla.spde2.result(inla=res,
                                                 name="spatial.field", spde=spde, do.transf=TRUE)
    }
    
  }
  
  
  time <- end - start
  
  if(two.stage.process == TRUE){
    
    out <- list(res = res,
                stk.full = stk.full,
                lc_all_method1 = lc_all_method1,
                lc_all_method2 = lc_all_method2,
                weights_grid_areas_intersection = weights_grid_areas_intersection,
                weights_points_areas_inside_closest = weights_points_areas_inside_closest,
                transformed.params.res = transformed.params.res,
                time = time)
    
  }else{
    
    out <- list(res1 = res1,
                res2 = res2,
                stk.full = stk.full,
                lc_all_method1 = lc_all_method1,
                lc_all_method2 = lc_all_method2,
                weights_grid_areas_intersection = weights_grid_areas_intersection,
                weights_points_areas_inside_closest = weights_points_areas_inside_closest,
                transformed.params.res1 = transformed.params.res1,
                transformed.params.res2 = transformed.params.res2,
                time = time)
    
  }
  
  
  return(out)
  
}




#' Stage 2 Fitting of a spatio-temporal model
#' 
#' @param inputdata - stage 1 model results
#' @param simulated.data - simulated data
#' @param n_random_fromPPD - number of simulations from the posterior predictive distribution of stage 1
#' @param n_random_fromposterior - number of simulations from the posterior marginals of stage 2 parameters
#' @param save.plot - TRUE/FALSE, if plots are saved or not
#' @param n.T - number of time points
#' @scenario_n - scenario number
#' @param informative - T/F, if priors are informative or not
#' 
#' @return - posterior estimates from stage 1 model 

stage2.fit.ST <- function(inputdata = stage1.res,
                          simulated.data = simulated.data,
                          n_random_fromPPD = n_random_fromPPD,
                          n_random_fromposterior = n_random_fromposterior,
                          save.plot = TRUE,
                          n.T = n.T,
                          scenario_n = scenario_n,
                          informative = informative){
  
  posterior_gamma0_FF_method1 = posterior_gamma0_FF_method2 = c() 
  posterior_gamma1_FF_method1 = posterior_gamma1_FF_method2 = c() 
  posterior_sigma2iid_FF_method1 = posterior_sigma2iid_FF_method2 = c() 
  posterior_sigma2time_FF_method1 = posterior_sigma2time_FF_method2 = c()
  counts_compile_method1 = counts_compile_method2 = c()
  
  #--- Prepare the data
  data_Poisson =simulated.data$shapefile_sim_complete@data
  data_Poisson$y_Poisson = as.integer(data_Poisson$y_Poisson)
  
  list_samples_joint_posterior_method1 = list()
  list_samples_joint_posterior_method2 = list()
  
  index_monitors <- inla.stack.index(inputdata$stk.full,"est.monitors")$data
  index_highres <- inla.stack.index(inputdata$stk.full,"est.highres")$data
  index_latent <- inla.stack.index(inputdata$stk.full,"est.true.exposures")$data 
  for(i in 1:n.T){
    index_pred <- inla.stack.index(inputdata$stk.full, paste0("pred.",i,sep=""))$data
    assign(paste0("index_pred.",i,sep=""),index_pred)
  }
  
  monitoringstations.coords = inputdata$monitoringstations.coords
  n_areas <- nrow(simulated.data$shapefile_sim_complete)/n.T
  dim_lp = length(index_monitors) + length(index_highres) + length(index_latent) + (length(index_pred.1)*n.T)
  #each row is an area, each column is a component of the linear combination
  
  cat("-----------------------------------------", " simulating from the PPD of stage 1 model ", n_random_fromPPD, " times. This might take a while!", sep="","\n" )
  
  start = Sys.time()
  
  for(jj in 1:n_random_fromPPD){
    
    sample_joint_posterior_method = inla.posterior.sample(1, inputdata$res)

    sample_joint_posterior_area_method1 = c()
    sample_joint_posterior_area_method2 = c()
    
    for (ii in 1:n.T){
      
      index_pred <- get(paste0("index_pred.",ii,sep=""))
      
      # Extract the part about the predictions (at the regular grid resolution)
      sample_joint_posterior_grid_method1 = sample_joint_posterior_method[[1]]$latent[1:dim_lp][index_pred]
      sample_joint_posterior_grid_method2 = sample_joint_posterior_method[[1]]$latent[1:dim_lp][index_pred]
      
      temp = c()
      for(j in 1:nrow(data_Poisson[which(data_Poisson$Time==ii),])){
        index_area = inputdata$lc_all_method1[[j+(98*(ii-1))]][[1]]$APredictor$idx - (length(index_monitors) + length(index_highres) + length(index_latent) + 
                                                                                        (length(index_pred) * (ii-1)))
        weight_area = inputdata$lc_all_method1[[j+(98*(ii-1))]][[1]]$APredictor$weight
        temp[j] = sum(sample_joint_posterior_grid_method1[index_area] * weight_area) #posterior area exposure mean
      }
      
      temp2 <-  data.frame(sample_joint_post_exposure_method1 = temp,
                           "AREAKEY"=unique(inputdata$weights_grid_areas_intersection$AREAKEY),
                           Time = ii)
      sample_joint_posterior_area_method1 <- rbind(sample_joint_posterior_area_method1, temp2)
      
      temp = c()
      for(j in 1:nrow(data_Poisson[which(data_Poisson$Time==ii),])){
        index_area = inputdata$lc_all_method2[[j+(98*(ii-1))]][[1]]$APredictor$idx - (length(index_monitors) + length(index_highres) + length(index_latent) + 
                                                                                        (length(index_pred) * (ii-1)))
        weight_area = inputdata$lc_all_method2[[j+(98*(ii-1))]][[1]]$APredictor$weight
        temp[j] = sum(sample_joint_posterior_grid_method2[index_area] * weight_area) #posterior area exposure mean
      }
      
      temp2 <-  data.frame(sample_joint_post_exposure_method2 = temp,
                           "AREAKEY"=unique(inputdata$weights_grid_areas_intersection$AREAKEY),
                           Time = ii)
      sample_joint_posterior_area_method2 <- rbind(sample_joint_posterior_area_method2, temp2)
      
    }
    
    list_samples_joint_posterior_method1[[jj]] <- sample_joint_posterior_area_method1
    list_samples_joint_posterior_method2[[jj]] <- sample_joint_posterior_area_method2
    
  }
  
  end = Sys.time()
  
  time.simfromppd = end - start
  
  ##############################################################
  # SECOND STAGE: Poisson model with FF approach (method1 and method2)
  ##############################################################
  #--- Method1
  #sample_joint_post_exposure_method1
  if (informative == TRUE){
    formula.Poisson.FF.method1 = y_Poisson ~ 1 + sample_joint_post_exposure_method1 +
      f(AREAKEY, model="iid", hyper=list(prec = list(prior = "pc.prec", param = c(sqrt(sigma2iid_Poisson), .05)))) +
      f(Time, model="iid", hyper=list(prec = list(prior = "pc.prec", param = c(sqrt(sigma2iid_time), .05))))
  }else{
    formula.Poisson.FF.method1 = y_Poisson ~ 1 + sample_joint_post_exposure_method1 + 
      f(AREAKEY, model="iid") +
      f(Time, model="iid")
  }
  
  
  compile.temporal.RW.mean.m1 <- c()
  compile.temporal.RW.p.25.m1 <- c()
  compile.temporal.RW.p.975.m1 <- c()
  
  cat("-----------------------------------------", " I am fitting the Poisson model - method 1 -- ","\n" )
  
  start = Sys.time()
  
  stop = FALSE
  for(j_FF in 1:n_random_fromPPD){
    
    # cat("-----------------------------------------", " I am fitting the Poisson model - method 1 -- ", j_FF, " out of ", n_random_fromPPD, sep="","\n" )
    
    output.Poisson.FF.method1 <- try(inla(formula.Poisson.FF.method1,family="poisson",
                                          data=data.frame(data_Poisson,
                                                          list_samples_joint_posterior_method1[[j_FF]]),
                                          offset=log(data_Poisson$Pop),
                                          control.predictor=list(compute=TRUE),
                                          control.compute=list(config = TRUE),
                                          verbose = FALSE))
    if (class(output.Poisson.FF.method1) == "try-error"){
      stop = TRUE
      break}

    if(any(is.na(output.Poisson.FF.method1$summary.hyperpar$mean)) == T |
       any(is.na(output.Poisson.FF.method1$summary.hyperpar$sd)) == T |
       any(is.na(output.Poisson.FF.method1$summary.hyperpar$mode)) == T |
       class(try(inla.rmarginal(2, output.Poisson.FF.method1$marginals.fixed$`(Intercept)`))) == "try-error" |
       class(try(inla.rmarginal(2, output.Poisson.FF.method1$marginals.fixed$sample_joint_post_exposure_method1))) == "try-error" |
       class(try(inla.rmarginal(2,inla.tmarginal(function(x) 1/x, output.Poisson.FF.method1$marginals.hyperpar$`Precision for AREAKEY`)))) == "try-error" |
       class(try(inla.rmarginal(2,inla.tmarginal(function(x) 1/x, output.Poisson.FF.method1$marginals.hyperpar$`Precision for Time`)))) == "try-error"){
     
      stop = TRUE
      break
       
    }
          
    posterior_gamma0_FF_method1 = c(posterior_gamma0_FF_method1,
                                    inla.rmarginal(n_random_fromposterior, output.Poisson.FF.method1$marginals.fixed$`(Intercept)`))
    posterior_gamma1_FF_method1 = c(posterior_gamma1_FF_method1,
                                    inla.rmarginal(n_random_fromposterior, output.Poisson.FF.method1$marginals.fixed$sample_joint_post_exposure_method1))
    sigma2iid_marg_post_FF.method1 = inla.tmarginal(function(x) 1/x,
                                                    output.Poisson.FF.method1$marginals.hyperpar$`Precision for AREAKEY`)
    posterior_sigma2iid_FF_method1 = c(posterior_sigma2iid_FF_method1,
                                       inla.rmarginal(n_random_fromposterior, sigma2iid_marg_post_FF.method1))
    sigma2time_marg_post_FF.method1 = inla.tmarginal(function(x) 1/x,
                                                     output.Poisson.FF.method1$marginals.hyperpar$`Precision for Time`)
    posterior_sigma2time_FF_method1 = c(posterior_sigma2time_FF_method1,
                                        inla.rmarginal(n_random_fromposterior, sigma2time_marg_post_FF.method1))
    
    counts_compile_method1 = cbind(counts_compile_method1, output.Poisson.FF.method1$summary.fitted.values$mean)      
          
  }
  
  end = Sys.time()
  
  time.M1 = end - start
  
  if(stop){failure.m1 = TRUE}else{failure.m1 = FALSE}
  
  #--- Method2
  #sample_joint_post_exposure_method2
  if (informative == TRUE){
    formula.Poisson.FF.method2 = y_Poisson ~ 1 + sample_joint_post_exposure_method2 +
      f(AREAKEY,model="iid", hyper=list(prec = list(prior = "pc.prec", param = c(sqrt(sigma2iid_Poisson), .05)))) +
      f(Time,model="iid", hyper=list(prec = list(prior = "pc.prec", param = c(sqrt(sigma2iid_time), .05))))
  }else{
    formula.Poisson.FF.method2 = y_Poisson ~ 1 + sample_joint_post_exposure_method2 +
      f(AREAKEY,model="iid") +
      f(Time,model="iid")
  }
  
  
  compile.temporal.RW.mean.m2 <- c()
  compile.temporal.RW.p.25.m2 <- c()
  compile.temporal.RW.p.975.m2 <- c()
  
  cat("-----------------------------------------", " I am fitting the Poisson model - method 2 -- ","\n" )
  
  start = Sys.time()
  
  stop = FALSE
  for(j_FF in 1:n_random_fromPPD){
    
    #cat("-----------------------------------------", " I am fitting the Poisson model - method 2 -- ", j_FF, " out of ", n_random_fromPPD, sep="","\n" )
    
    output.Poisson.FF.method2 <- try(inla(formula.Poisson.FF.method2,family="poisson",
                                          data=data.frame(data_Poisson,
                                                          list_samples_joint_posterior_method2[[j_FF]]),
                                          offset=log(data_Poisson$Pop),
                                          control.predictor=list(compute=TRUE),
                                          control.compute=list(config = TRUE),
                                          verbose = FALSE))
    
    if (class(output.Poisson.FF.method2) == "try-error"){
      stop = TRUE
      break}
    
    if(any(is.na(output.Poisson.FF.method2$summary.hyperpar$mean)) == T |
       any(is.na(output.Poisson.FF.method2$summary.hyperpar$sd)) == T |
       any(is.na(output.Poisson.FF.method2$summary.hyperpar$mode)) == T |
       class(try(inla.rmarginal(2, output.Poisson.FF.method2$marginals.fixed$`(Intercept)`))) == "try-error" |
       class(try(inla.rmarginal(2, output.Poisson.FF.method2$marginals.fixed$sample_joint_post_exposure_method2))) == "try-error" |
       class(try(inla.rmarginal(2,inla.tmarginal(function(x) 1/x, output.Poisson.FF.method2$marginals.hyperpar$`Precision for AREAKEY`)))) == "try-error" |
       class(try(inla.rmarginal(2,inla.tmarginal(function(x) 1/x, output.Poisson.FF.method2$marginals.hyperpar$`Precision for Time`)))) == "try-error"){
      
      stop = TRUE
      break
      
    }
    
    posterior_gamma0_FF_method2 = c(posterior_gamma0_FF_method2,
                                    inla.rmarginal(n_random_fromposterior, output.Poisson.FF.method2$marginals.fixed$`(Intercept)`))
    posterior_gamma1_FF_method2 = c(posterior_gamma1_FF_method2,
                                    inla.rmarginal(n_random_fromposterior, output.Poisson.FF.method2$marginals.fixed$sample_joint_post_exposure_method2))
    sigma2iid_marg_post_FF.method2 = inla.tmarginal(function(x) 1/x,
                                                    output.Poisson.FF.method2$marginals.hyperpar$`Precision for AREAKEY`)
    posterior_sigma2iid_FF_method2 = c(posterior_sigma2iid_FF_method2,
                                       inla.rmarginal(n_random_fromposterior, sigma2iid_marg_post_FF.method2))
    sigma2time_marg_post_FF.method2 = inla.tmarginal(function(x) 1/x,
                                                     output.Poisson.FF.method2$marginals.hyperpar$`Precision for Time`)
    posterior_sigma2time_FF_method2 = c(posterior_sigma2time_FF_method2,
                                        inla.rmarginal(n_random_fromposterior, sigma2time_marg_post_FF.method2))
    
    counts_compile_method2 = cbind(counts_compile_method2, output.Poisson.FF.method2$summary.fitted.values$mean)  
      
  
  }
  
  end = Sys.time()
  
  time.M2 = end - start
  
  if(stop){failure.m2 = TRUE}else{failure.m2 = FALSE}
  
  
  cat("-----------------------------------------", " DONE with stage 2 model fitting --- now saving plots and outputs ", sep="","\n" )
  
  
  if (save.plot == TRUE){
    
    gamma0.plot <- combined.post.dens.plot.FF(sample1 = posterior_gamma0_FF_method1,
                                              sample2 = posterior_gamma0_FF_method2,
                                              true.value = gamma0_Poisson,
                                              title = expression(paste(hat(pi),"(",gamma[0],"|y)",)),
                                              xlab = expression(gamma[0]))
    gamma1.plot <- combined.post.dens.plot.FF(sample1 = posterior_gamma1_FF_method1,
                                              sample2 = posterior_gamma1_FF_method2,
                                              true.value = gamma1_Poisson,
                                              title = expression(paste(hat(pi),"(",gamma[1],"|y)",)),
                                              xlab = expression(gamma[1]))
    sigma2iid.plot <- combined.post.dens.plot.FF(sample1 = posterior_sigma2iid_FF_method1,
                                                 sample2 = posterior_sigma2iid_FF_method2,
                                                 true.value = sigma2iid_Poisson,
                                                 title = expression(paste(hat(pi),"(",sigma[phi]^2,"|y)",)),
                                                 xlab = expression(sigma^2))
    sigma2time.plot <- combined.post.dens.plot.FF(sample1 = posterior_sigma2time_FF_method1,
                                                  sample2 = posterior_sigma2time_FF_method2,
                                                  true.value = sigma2iid_time,
                                                  title = expression(paste(hat(pi),"(",sigma[nu]^2,"|y)",)),
                                                  xlab = expression(sigma[nu]^2))
    
    
    png(paste0(path,'Scenario=',scenario_n,'.S2.Results.posteriors.ST.png',sep=""), width=30, height=20, units = 'cm', res = 300)
    print((gamma0.plot + gamma1.plot)/
            (sigma2iid.plot + sigma2time.plot))  
    dev.off()
    
    counts.m1 <- rowMeans(counts_compile_method1)
    counts.m2 <- rowMeans(counts_compile_method2)
    
    forplot <- simulated.data$shapefile_sim_complete
    colnames(forplot@data)[which(colnames(forplot@data) == 'y_Poisson')] <- 'true'
    forplot$M1 <- counts.m1
    forplot$M2 <- counts.m2
    
    all_vals = c(forplot$true, forplot$M1, forplot$M2)
    
    for (i in 1:n.T){
      
      sub.forplot <- forplot@data[which(forplot$Time == i),]
      
      plotdata <- simulated.data$shapefile_sim_complete
      plotdata@data <- sub.forplot
      plotdata_tidy <- tidy(plotdata)
      plotdata$id <- row.names(plotdata)
      plotdata <- left_join(plotdata_tidy, plotdata@data)
      
      counts.true.map <- ggplot(plotdata, aes(x = long, y = lat, group = group, fill = true)) +
        geom_polygon(color = "black", size = 0.1) +
        coord_equal(ratio=1) +
        labs(x="", y="", title=paste0("True Poisson counts, Time = ",i,sep="")) +
        theme(legend.key.width=unit(1, "cm")) + 
        theme_map() +
        scale_fill_viridis(option = 'C', limits = c(min(all_vals), max(all_vals))) +
        theme(legend.position="bottom")
      counts.M1.map <- ggplot(plotdata, aes(x = long, y = lat, group = group, fill = M1)) +
        geom_polygon(color = "black", size = 0.1) +
        coord_equal(ratio=1) +
        labs(x="", y="", title=paste0("Estimated counts (Method 1), Time = ",i,sep="")) +
        theme(legend.key.width=unit(1, "cm")) + 
        theme_map() +
        scale_fill_viridis(option = 'C', limits = c(min(all_vals), max(all_vals))) +
        theme(legend.position="bottom")
      counts.M2.map <- ggplot(plotdata, aes(x = long, y = lat, group = group, fill = M2)) +
        geom_polygon(color = "black", size = 0.1) +
        coord_equal(ratio=1) +
        labs(x="", y="", title=paste0("Estimated counts (Method 2), Time = ",i,sep="")) +
        theme(legend.key.width=unit(1, "cm")) + 
        theme_map() +
        scale_fill_viridis(option = 'C', limits = c(min(all_vals), max(all_vals))) +
        theme(legend.position="bottom")
      
      png(paste0(path,'Scenario=',scenario_n,'.S2.ST.Results.maps.Time=',i,'.png',sep=""), width=30, height=20, units = 'cm', res = 300)
      print(counts.true.map + counts.M1.map + counts.M2.map)
      dev.off()
      
    }
    
    temp <- simulated.data$shapefile_sim_complete@data[,c('y_Poisson','Time','AREAKEY')]
    subset.areas <- unique(simulated.data$shapefile_sim_complete$AREAKEY)[1:5]
    for.plot <- cbind(temp, rowMeans(counts_compile_method1), rowMeans(counts_compile_method2))
    names(for.plot)[c(1,4,5)] <- c('true','method 1', 'method 2')
    
    for (i in 1:length(subset.areas)){
      
      sub.for.plot <- subset(for.plot, subset = (AREAKEY == subset.areas[i]))
      meltdf <- melt(sub.for.plot,id.vars=c("Time","AREAKEY"),measure.vars=c('true','method 1','method 2'))
      
      png(paste0(path,'Scenario=',scenario_n,'S2.Obs.vs.Preds.Count.',subset.areas[i],'.png',sep=""), width=20, height=20, units = 'cm', res = 300)
      print(ggplot(meltdf,aes(x=Time,y=value,colour=variable,group=variable)) + geom_line() + 
              ggtitle(paste0("Observed and predicted counts at AREA = ",subset.areas[i],sep=""))+geom_point())
      dev.off()
      
    }
    
  }
  
  
  out <- list(posterior_gamma0_FF_method1 = posterior_gamma0_FF_method1,
              posterior_gamma1_FF_method1 = posterior_gamma1_FF_method1,
              posterior_sigma2iid_FF_method1 = posterior_sigma2iid_FF_method1,
              posterior_sigma2time_FF_method1 = posterior_sigma2time_FF_method1,
              posterior_gamma0_FF_method2 = posterior_gamma0_FF_method2,
              posterior_gamma1_FF_method2 = posterior_gamma1_FF_method2,
              posterior_sigma2iid_FF_method2 = posterior_sigma2iid_FF_method2,
              posterior_sigma2time_FF_method2 = posterior_sigma2time_FF_method2,
              list_samples_joint_posterior_method1 = list_samples_joint_posterior_method1,
              list_samples_joint_posterior_method2 = list_samples_joint_posterior_method2,
              time.simfromppd = time.simfromppd,
              time.M1 = time.M1,
              time.M2 = time.M2,
              failure.m1 = failure.m1,
              failure.m2 = failure.m2)
  
  return(out)
  
}







