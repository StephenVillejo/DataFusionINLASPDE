

#' DATA GENERATION FUNCTION
#' 
#' @param scenario_n - scenario number 
#' @param save.plot - T/F, if plots are produced or not
#' @param res_pred_grid - resolution of prediction grid
#' @param highres_simulation_reggrid - resolution of simulation grid
#' @param prop.points.for.area - proportion of points to simulated in each area
#' @param sigma2xi - spatial variance of Matern function
#' @param kappa - kappa parameter of Matern function 
#' @param b0 - intercept in the exposure model
#' @param b1 - coefficient of the covariate in the exposures model
#' @param sigma2_e - measurement error variance
#' @param gamma0_Poisson - intercept of the poisson model
#' @param gamma1_Poisson - coefficient of the true exposures in the poisson model
#' @param sigma2iid_Poisson - area-specific iid effect
#' @param highres.error.sigma2 - error variance for the highres data
#' @param highres.b0 - intercept of the high res data melding model
#' @param highres.b1 - slope of the high res data melding model
#' @param n.T - number of time points
#' @param a - autocorrelation parameter
#' @param unif.pop.a - uniform dist paramater (lowerbound) for the grid population size
#' @param unif.pop.b - uniform dist paramater (upperbound) for the grid population sizea
#' @param sigma2iid_time - time-specific iid effect
#' @param sparse - TRUE/FALSE, if true, then monitoring stations are sparsely located in the map
#' @param sim_i - simulation/replication number
#' 
#' @return A matrix of size nrow(locs) x n, each column corresponds to a realisation of the random field at the locations specified in locs.

Data_simulation = function(scenario_n,
                           save.plot,
                           res_pred_grid,
                           highres_simulation_reggrid,
                           prop.points.for.area,
                           sigma2xi,
                           kappa,
                           b0,
                           b1,
                           sigma2_e,
                           gamma0_Poisson,
                           gamma1_Poisson,
                           sigma2iid_Poisson,
                           highres.error.sigma2,
                           highres.b0,
                           highres.b1,
                           n.T,
                           a,
                           unif.pop.a,
                           unif.pop.b,
                           sigma2iid_time,
                           sparse,
                           sim_i) {
  
  #--- Import the shapefile
  shapefile_sim.orig = readOGR(system.file("etc/shapes/bhicv.shp",package="spdep")[1]) #coordinates in km 
  colnames(shapefile_sim.orig@data)[2] = "AREAKEY" 
  
  #--- Get the area of each polygon
  areas <- sapply(slot(shapefile_sim.orig, "polygons"), slot, "area")
  areas <- as.data.frame(cbind(areas, shapefile_sim.orig@data$AREAKEY))
  areas$areas <- as.numeric(areas$areas)
  colnames(areas)[2] <- c("AREAKEY")
  
  #--- Get the centroids of each polygon
  n_areas = nrow(shapefile_sim.orig@data)
  shapefile_sim_centroids = gCentroid(shapefile_sim.orig, byid = T)
  centroids_df = data.frame(AREAKEY = shapefile_sim.orig@data$AREAKEY, centroids = shapefile_sim_centroids)
  colnames(centroids_df)[2:3] = c("x1", "x2")
  
  ############################################################
  # Create a high resolution grid of points  for exposure simulation
  # The resolution is highres_simulation_reggrid
  ############################################################
  
  sim_grid_x1 = seq(shapefile_sim.orig@bbox[1, 1] - 0.1, shapefile_sim.orig@bbox[1, 2] + 0.1, l = highres_simulation_reggrid) 
  sim_grid_x2 = seq(shapefile_sim.orig@bbox[2, 1] - 0.1, shapefile_sim.orig@bbox[2, 2] + 0.1, l = highres_simulation_reggrid) 
  
  sim_grid = expand.grid(x = sim_grid_x1, y = sim_grid_x2)
  coordinates(sim_grid) <- ~ x + y #this converts the dataframe to a spatialpoints object
  colnames(sim_grid@coords) = c("x1", "x2")
  proj4string(sim_grid) = proj4string(shapefile_sim.orig) #this line is not necessary
  
  ############################################################
  # Exposure simulation grid and parameters
  ############################################################
  
  bnd = unionSpatialPolygons(shapefile_sim.orig, rep(1, nrow(shapefile_sim.orig@data)))
  boundary <- list(
    as.inla.mesh.segment(bnd),
    NULL)
  mesh <- inla.mesh.2d(boundary=boundary,
                       max.edge=c(0.075, 0.5),
                       min.angle=c(30, 21),
                       max.n=c(48000, 16000), ## Safeguard against large meshes.
                       max.n.strict=c(128000, 128000), ## Don't build a huge mesh!
                       cutoff=0.022, ## Filter away adjacent points.
                       offset=c(0.17, 0.685)) ## Offset for extra boundaries, if needed.
  
  z1.matrix <- matrix(NA,nrow=highres_simulation_reggrid^2, ncol=n.T)
  cov2.matrix <- matrix(NA,nrow=highres_simulation_reggrid^2, ncol=n.T)
  
  for (i in 1:n.T){
    z1.matrix[,i] <- rnorm(highres_simulation_reggrid^2,0,1)
  }
  
  ### Simulate the exposures field
  
  omega <- rspde(coords = as.matrix(sim_grid@coords),
                 kappa = kappa,
                 variance = sigma2xi,
                 n = n.T,
                 mesh = mesh)
  xi <- omega
  for (j in 2:n.T)
    xi[,j] <- a*xi[,j-1] + sqrt(1-a^2)*omega[,j]
  
  exposure = b0 + b1 * z1.matrix + xi 
  
  # To save all rasters (proxy data for exposure, covariate raster, population raster)
  raster.exposures.highres.list <- vector("list", n.T)
  raster.covariate.2.list <- vector("list", n.T)
  raster.population.list <- vector("list", n.T)
  compile.exposures.highres.list <- vector("list", n.T)
  compile.population.list <- vector("list", n.T)
  compile.risk.grid <- vector("list",n.T)
  
  # To compile all datasets (area level and monitors)
  data_Poisson <- c()
  compile_monitors <- c()
  compile.sim_grid_df <- c()
  compile.sim_grid_inside_df <- c()
  
  # Simulate first the nus (time-specific effects)
  #temp <- rnorm(n = n.T, 0, sd = sqrt(sigma2iid_time))
  #nu <- cumsum(temp)
  nu <- rnorm(n = n.T, mean = 0, sd = sqrt(sigma2iid_time))
  
  # Simulate first the phis for non-aggregated poisson count (area-specific iid effect)
  phi = rnorm(n_areas, mean = 0, sd = sqrt(sigma2iid_Poisson))
  
  for (i in 1:n.T){
    
    shapefile_sim <- shapefile_sim.orig
    
    #--- Associate each point to an area of the shapefile and select points inside the shapefile
    sim_grid_df = over(sim_grid, shapefile_sim)
    sim_grid_df$z1 = z1.matrix[,i]
    sim_grid_df = cbind(sim_grid_df,sim_grid@coords)
    # Note that the arrangement of the data is arranged by row first - first row values, then second row values, so on and so forth
    sim_grid_df$exposure = exposure[,i]
    sim_grid_df$Time = i
    
    #--- Just select the points inside the shapefile
    which_inside = which(!is.na(sim_grid_df$AREAKEY))
    sim_grid_inside_df = sim_grid_df[which_inside, ] #the data/information on the spatial points inside the polygon
    sim_grid_inside = sim_grid[which_inside] #the spatial points inside the polygon
    
    ############################################################
    # Compute the TRUE *AREA* exposure as mean of the sim_grid points inside each area
    ############################################################
    true_area_exposure = aggregate(exposure ~ AREAKEY,
                                   data=sim_grid_inside_df,
                                   FUN=mean)
    colnames(true_area_exposure)[2] = "true_area_exposure"
    
    #--- Include the "true" area (mean) value in the shapefile
    shapefile_sim = merge(shapefile_sim,
                          true_area_exposure,
                          by = "AREAKEY",
                          all.x = T)
    
    #--- Compute the variability of the true area points inside each area
    true_area_exposure_var = aggregate(exposure ~ AREAKEY,
                                       data=sim_grid_inside_df,
                                       FUN=var)
    colnames(true_area_exposure_var)[2] = "var_true_area_exposure"
    shapefile_sim = merge(shapefile_sim,
                          true_area_exposure_var,
                          by = "AREAKEY",
                          all.x = T)
    
    #--- Compute the average of z1 inside each area
    area_z1 = aggregate(z1 ~ AREAKEY,
                        data=sim_grid_inside_df,
                        FUN=mean)
    shapefile_sim = merge(shapefile_sim,
                          area_z1,
                          by = "AREAKEY",
                          all.x = T)
    
    shapefile_sim@data$Time = i
    sim_grid_inside_df$Time = i
    sim_grid_df$Time = i
    
    
    compile.sim_grid_df = rbind(compile.sim_grid_df,sim_grid_df)
    compile.sim_grid_inside_df = rbind(compile.sim_grid_inside_df, sim_grid_inside_df)
    
    compute_Y_Poisson = merge(sim_grid_inside_df, true_area_exposure, by="AREAKEY")
    compute_Y_Poisson$pop.raster = round(runif(n = nrow(compute_Y_Poisson), min = unif.pop.a, max = unif.pop.b))
    
    area_pop = as.data.frame(describeBy(compute_Y_Poisson$pop.raster, group=compute_Y_Poisson$AREAKEY, mat=TRUE))
    area_pop = area_pop[,c('group1','n','mean')]
    area_pop$Pop = area_pop$n * area_pop$mean
    area_pop = area_pop[,c('group1','Pop')]
    colnames(area_pop)[1] = "AREAKEY"
    
    temp.data_Poisson = merge(shapefile_sim, area_pop, by = "AREAKEY")
    
    log_rate_mean_exposure = gamma0_Poisson + (gamma1_Poisson * shapefile_sim@data$true_area_exposure) + phi +
      rep(nu[i],nrow(shapefile_sim@data))
    temp.data_Poisson$y_Poisson = rpois(n_areas, lambda = exp(log_rate_mean_exposure) * temp.data_Poisson@data$Pop)
    
    data_Poisson = rbind(data_Poisson,temp.data_Poisson@data)
    
    compile.risk.grid[[i]] <- compute_Y_Poisson
    
    
    #############################################################
    # Sample the monitoring station sites 
    #############################################################
    # The proportion of monitoring stations is a subset of the high resolution simulation grid which is fixed
    # We use the same location of monitoring stations for all scenarios with the same prop.points.for.area in the simulation settings
    
    monitoringstations.data <- get.monitoring.stations(prop.points.for.area = prop.points.for.area, compute_Y_Poisson = compute_Y_Poisson,
                                                       shapefile_sim = shapefile_sim.orig,
                                                       sim_grid_df = sim_grid_df,
                                                       sparse = sparse)
    monitoringstations.coords = monitoringstations.data[,c("x1","x2")]
    
    ############################################################
    # Simulate the exposure for the monitoring stations (just add the measurement error)
    ############################################################
    meas_error = rnorm(length(monitoringstations.data$exposure),mean = 0,sd = sqrt(sigma2_e))
    monitoringstations.data$exposure_with_meas_err = monitoringstations.data$exposure + meas_error
    
    final.data = monitoringstations.data[,c('x1','x2','AREAKEY','z1','exposure','pop.raster','exposure_with_meas_err')]
    final.data$Time = i
    
    compile_monitors <- rbind(compile_monitors, final.data)
    
    #############################################################
    # Rasters (covariate rasters and aggregation/population raster)
    #############################################################
    
    exposures.highres = sim_grid_inside_df[,c('x1','x2','exposure')]
    exposures.highres$exposure.highres <- highres.b0 + (exposures.highres$exposure * highres.b1) +  rnorm(nrow(exposures.highres),mean = 0,sd = sqrt(highres.error.sigma2))
    exposures.highres.data = exposures.highres
    exposures.highres = exposures.highres[,c('x1','x2','exposure.highres')]
    coordinates(exposures.highres) <- ~ x1 + x2
    gridded(exposures.highres) <- TRUE
    raster.exposures.highres <- raster(exposures.highres)
    
    population.raster = compute_Y_Poisson[,c('x1','x2','pop.raster')]
    coordinates(population.raster) <- ~ x1 + x2
    gridded(population.raster) <- TRUE
    raster.population <- raster(population.raster)
    
    raster.exposures.highres.list[[i]] <- raster.exposures.highres
    raster.population.list[[i]] <- raster.population
    compile.exposures.highres.list[[i]] <- exposures.highres.data
    
  }
  
  shapefile_sim_complete <- shapefile_sim.orig
  shapefile_sim_complete@data <- data_Poisson
  
  tag <- compile_monitors[which(compile_monitors$Time==1),c('x1','x2')]
  tag$monitor <- paste0("M",c(1:nrow(tag)),sep="")
  temp <- merge(x = compile_monitors, y = tag, by = c('x1','x2'))
  compile_monitors <- temp
  
  ############################################################
  # Create a regular grid for the prediction
  # grid resolution = highres_simulation_reggrid/res_pred_grid
  ############################################################
  # Note that if highres_simulation_reggrid = res_pred_grid, then the prediction grid is the same as the high res simulation grid
  pred_grid_x1 = unique(sim_grid@coords[, 1])[seq(1, highres_simulation_reggrid,
                                                  by = highres_simulation_reggrid / res_pred_grid)]
  pred_grid_x2 = unique(sim_grid@coords[, 2])[seq(1, highres_simulation_reggrid,
                                                  by = highres_simulation_reggrid / res_pred_grid)]
  
  pred_grid = expand.grid(x = pred_grid_x1, y = pred_grid_x2)
  coordinates(pred_grid) <- ~ x + y
  proj4string(pred_grid) = proj4string(sim_grid)
  # Transform the regular grid into a SpatialPolygons object
  pred_grid_poly = as.SpatialPolygons.GridTopology(points2grid(pred_grid))
  proj4string(pred_grid_poly) = proj4string(shapefile_sim)
  #plot(shapefile_sim)
  #plot(pred_grid_poly,add=T,border="grey")
  
  # Change the name of the pixels of the prediction grid (from left to right from bottom to top)
  ncol_pred_grid = res_pred_grid
  nrow_pred_grid = res_pred_grid
  mydesiredseq = c()
  for (i in (nrow_pred_grid - 1):0) {
    temp = seq((ncol_pred_grid * i) + 1, l = ncol_pred_grid)
    mydesiredseq = c(mydesiredseq, temp)
  }
  for (i in 1:length(pred_grid_poly@polygons)) {
    pred_grid_poly@polygons[[i]]@ID = paste("g", mydesiredseq[i], sep = "")
  }
  
  weights_grid_areas_intersection <- compute.intersections(res_pred_grid = 100, compute_Y_Poisson=compute_Y_Poisson)
  weights_points_areas_inside_closest <- compute.predpoints(res_pred_grid = 100, compute_Y_Poisson=compute_Y_Poisson)
  
  
  ########################################################################################################################
  #--- Prediction dataset
  ########################################################################################################################
  which_in_pred_grid = over(pred_grid, sim_grid)
  compile <- c()
  for (i in 1:n.T){
    sim_grid_inside_df_pred = compile.sim_grid_df[which(compile.sim_grid_df$Time==i),][which_in_pred_grid,]
    compile <- rbind(compile, sim_grid_inside_df_pred)
  }
  sim_grid_inside_df_pred <- compile
  
  
  if(save.plot == TRUE){
  
    all.exposures <- bind_rows(compile.exposures.highres.list, .id = "column_label")
    all_vals = c(all.exposures$exposure, all.exposures$exposure.highres)
    scl_1 = scale_colour_viridis_c(limits = c(min(all_vals), max(all_vals)))
    scl_2 = scale_fill_viridis_c(limits = c(min(all_vals), max(all_vals)))
    
    plots.proxy.data <- vector("list", n.T)      
    plots.exposuresfield <- vector("list", n.T)
    plots.monitors <- vector("list", n.T)
    plots.healthdata <- vector("list", n.T)
    plots.risk.surface <- vector("list", n.T)
    
    
    for (i in 1:n.T) {
      
      exposures.highres.data <- compile.exposures.highres.list[[i]]
      colnames(exposures.highres.data)[4] <- 'value'
      proxydata <- ggplot() +  
        geom_tile(data=exposures.highres.data, aes(x=x1, y=x2, fill=value),
                  alpha=0.8) + 
        coord_equal() +
        theme_map() +
        theme(legend.position="bottom") +
        theme(legend.key.width=unit(1, "cm")) + 
        ggtitle(paste("Proxy Data, time =",i,sep=" ")) +
        scl_2
      plots.proxy.data[[i]] <- proxydata
      
      exposurefield <- ggplot() +  
        geom_tile(data=exposures.highres.data, aes(x=x1, y=x2, fill=exposure),
                  alpha=0.8) + 
        coord_equal() +
        theme_map() +
        theme(legend.position="bottom") +
        theme(legend.key.width=unit(1, "cm")) + 
        ggtitle(paste("Exposures Field, time =",i,sep=" ")) +
        scl_2
      plots.exposuresfield[[i]] <- exposurefield
      
      
      forplot2 <- compile_monitors[which(compile_monitors$Time==i),]
      colnames(forplot2)[which(colnames(forplot2) == 'exposure_with_meas_err')] <- 'value'
      monitors <- ggplot() + geom_polygon(data=shapefile_sim, aes(x=long, y=lat, group=group),fill='grey99',color='gray',alpha=1) +
        geom_point(aes(x=x1, y=x2, color=value), data=forplot2, alpha=1, size=1) +
        coord_equal(ratio=1) +
        theme_map() +
        ggtitle(paste("Monitoring stations, time =",i,sep=" ")) +
        scale_colour_viridis_c(limits = c(min(compile_monitors$exposure_with_meas_err), max(compile_monitors$exposure_with_meas_err))) + 
        theme(legend.position="bottom") +
        theme(legend.key.width=unit(1, "cm"))
      plots.monitors[[i]] <- monitors
      
      forplot3 <- shapefile_sim_complete@data[which(shapefile_sim_complete@data$Time==i),]
      temp.shapefile <- shapefile_sim_complete
      temp.shapefile@data <- forplot3
      forplot3 <- temp.shapefile
      colnames(forplot3@data)[which(colnames(forplot3@data) == 'y_Poisson')] <- 'value'
      forplot3_tidy <- tidy(forplot3)
      forplot3$id <- row.names(forplot3)
      forplot3 <- left_join(forplot3_tidy, forplot3@data)
      
      healthdata <- ggplot(forplot3, aes(x = long, y = lat, group = group, fill = value)) +
        geom_polygon(color = "black", size = 0.1) +
        coord_equal(ratio=1) +
        ggtitle(paste("Health Data, time =",i,sep=" ")) +
        theme(legend.key.width=unit(1, "cm")) + 
        theme_map() +
        scale_fill_viridis(option = 'C', limits = c(min(shapefile_sim_complete@data$y_Poisson), max(shapefile_sim_complete@data$y_Poisson))) +
        theme(legend.position="bottom")
      plots.healthdata[[i]] <- healthdata
      
    }
    
    
    if(n.T<6){
      png(paste0(path,'Scenario=',scenario_n,'.DATA=Proxy.Data.All.Times.png',sep=""), width=30, height=30, units = 'cm', res = 300)
      print(do.call("grid.arrange", c(plots.proxy.data, ncol=n.T)))
      dev.off()
      
      png(paste0(path,'Scenario=',scenario_n,'.DATA=Exposures.Field.All.Times.png',sep=""), width=30, height=30, units = 'cm', res = 300)
      print(do.call("grid.arrange", c(plots.exposuresfield, ncol=n.T)))
      dev.off()
      
      png(paste0(path,'Scenario=',scenario_n,'.DATA=Monitors.All.Times.png',sep=""), width=30, height=30, units = 'cm', res = 300)
      print(do.call("grid.arrange", c(plots.monitors, ncol=n.T)))
      dev.off()
      
      png(paste0(path,'Scenario=',scenario_n,'.DATA=Health.Data.All.Times.png',sep=""), width=30, height=30, units = 'cm', res = 300)
      print(do.call("grid.arrange", c(plots.healthdata, ncol=n.T)))
      dev.off()
      
    }else{
      
      png(paste0(path,'Scenario=',scenario_n,'.DATA=Proxy.Data.All.Times.png',sep=""), width=60, height=30, units = 'cm', res = 300)
      print(do.call("grid.arrange", c(plots.proxy.data, ncol=n.T/2)))
      dev.off()
      
      png(paste0(path,'Scenario=',scenario_n,'.DATA=Exposures.Field.All.Times.png',sep=""), width=60, height=30, units = 'cm', res = 300)
      print(do.call("grid.arrange", c(plots.exposuresfield, ncol=n.T/2)))
      dev.off()
      
      png(paste0(path,'Scenario=',scenario_n,'.DATA=Monitors.All.Times.png',sep=""), width=60, height=30, units = 'cm', res = 300)
      print(do.call("grid.arrange", c(plots.monitors, ncol=n.T/2)))
      dev.off()
      
      png(paste0(path,'Scenario=',scenario_n,'.DATA=Health.Data.All.Times.png',sep=""), width=60, height=30, units = 'cm', res = 300)
      print(do.call("grid.arrange", c(plots.healthdata, ncol=n.T/2)))
      dev.off()
      
    }
    
    for.plot <- subset(compile_monitors, subset=(monitor=='M1'|monitor=='M2'|
                                                   monitor=='M3'|monitor=='M4'|
                                                   monitor=='M5'|monitor=='M6'|
                                                   monitor=='M7'|monitor=='M8'|
                                                   monitor=='M9'|monitor=='M10'))
    for.plot <- for.plot[,c('exposure_with_meas_err','Time','monitor')]
    meltdf <- melt(for.plot,id.vars=c("Time","monitor"),measure.vars='exposure_with_meas_err')
    meltdf <- meltdf[order(for.plot$monitor,for.plot$Time),]
    png(paste0(path,'Scenario=',scenario_n,'.DATA=Time.Plot.Monitors.png',sep=""), width=30, height=20, units = 'cm', res = 300)
    print(ggplot(meltdf,aes(x=Time,y=value,colour=monitor,group=monitor)) + geom_line() + 
            ggtitle("Time series plot of observed values at selected monitors")+geom_point())
    dev.off()
    
    for.plot <- shapefile_sim_complete@data[,c('y_Poisson','Time','AREAKEY')]
    subset.areas <- shapefile_sim.orig$AREAKEY[1:10]
    for.plot <- subset(for.plot, subset=(AREAKEY==subset.areas[1]|AREAKEY==subset.areas[2]|
                                           AREAKEY==subset.areas[3]|AREAKEY==subset.areas[4]|
                                           AREAKEY==subset.areas[5]|AREAKEY==subset.areas[6]|
                                           AREAKEY==subset.areas[7]|AREAKEY==subset.areas[8]|
                                           AREAKEY==subset.areas[9]|AREAKEY==subset.areas[10]))
    meltdf <- melt(for.plot,id.vars=c("Time","AREAKEY"),measure.vars='y_Poisson')
    meltdf <- meltdf[order(for.plot$AREAKEY,for.plot$Time),]
    
    png(paste0(path,'Scenario=',scenario_n,'.DATA=Time.Plot.Health.Data.png',sep=""), width=30, height=20, units = 'cm', res = 300)
    print(ggplot(meltdf,aes(x=Time,y=value,colour=AREAKEY,group=AREAKEY)) + geom_line() + 
            ggtitle("Time series plot of observed counts at selected areas")+geom_point())
    dev.off()
    
    
  }
  
  
  ############################################################
  # Save all the relevant data for model fitting
  ############################################################
  
  return(list(shapefile_sim_complete=shapefile_sim_complete,
              monitoringstations.data=monitoringstations.data,
              monitoringstations.coords=monitoringstations.coords,
              pred_grid_x1 = pred_grid_x1,
              pred_grid_x2 = pred_grid_x2,
              pred_grid=pred_grid,
              n.monitoringstats = nrow(monitoringstations.data),                
              final.data = compile_monitors,
              raster.exposures.highres.list = raster.exposures.highres.list,
              raster.population.list = raster.population.list,
              centroids_df = centroids_df,
              compile_monitors = compile_monitors,
              compile.exposures.highres.list = compile.exposures.highres.list,
              compile.sim_grid_inside_df = compile.sim_grid_inside_df,
              compile.sim_grid_df = compile.sim_grid_df,
              compile.risk.grid = compile.risk.grid,
              phi = phi,
              nu = nu))
  
}


#' Simulate Matern Field 
#' @param coords - spatial coordinate
#' @param kappa - kappa parameter of the Matern function
#' @param variance - spatial variance term
#' @param alppha - alpha parameter (function of nu) of the Matern, this is typically fixed
#' @param n - this refers to the number of samples
#' @param mesh - mesh for the discretization
#' @param verbose - 
#' @param seed - seed number
#' @param return.attributes - if attributes of the objects are returned 
#' 

rspde <- function (coords, kappa, variance = 1, alpha = 2, n = 1, mesh,
                   verbose = FALSE, seed, return.attributes = FALSE){
  t0 <- Sys.time()
  theta <- c(-0.5 * log(4 * pi * variance * kappa^2), log(kappa))
  if (verbose)
    cat("theta =", theta, "\n")
  if (missing(mesh)) {
    mesh.pars <- c(0.5, 1, 0.1, 0.5, 1) * sqrt(alpha - ncol(coords)/2)/kappa
    if (verbose)
      cat("mesh.pars =", mesh.pars, "\n")
    attributes <- list(mesh = inla.mesh.2d(, coords[chull(coords),
    ], max.edge = mesh.pars[1:2], cutoff = mesh.pars[3],
    offset = mesh.pars[4:5]))
    if (verbose)
      cat("n.mesh =", attributes$mesh$n, "\n")
  }
  else attributes <- list(mesh = mesh)
  attributes$spde <- inla.spde2.matern(attributes$mesh, alpha = alpha)
  attributes$Q <- inla.spde2.precision(attributes$spde, theta = theta)
  attributes$A <- inla.mesh.project(mesh = attributes$mesh,
                                    loc = coords)$A
  if (n == 1)
    result <- drop(attributes$A %*% inla.qsample(Q = attributes$Q,
                                                 constr = attributes$spde$f$extraconstr))
  t1 <- Sys.time()
  result <- inla.qsample(n, attributes$Q, seed = ifelse(missing(seed),
                                                        0, seed), constr = attributes$spde$f$extraconstr)
  if (nrow(result) < nrow(attributes$A)) {
    result <- rbind(result, matrix(NA, nrow(attributes$A) -
                                     nrow(result), ncol(result)))
    dimnames(result)[[1]] <- paste("x", 1:nrow(result), sep = "")
    for (j in 1:ncol(result)) result[, j] <- drop(attributes$A %*%
                                                    result[1:ncol(attributes$A), j])
  }
  else {
    for (j in 1:ncol(result)) result[1:nrow(attributes$A),
                                     j] <- drop(attributes$A %*% result[, j])
    result <- result[1:nrow(attributes$A), ]
  }
  t2 <- Sys.time()
  attributes$cpu <- c(prep = t1 - t0, sample = t2 - t1, total = t2 -
                        t0)
  if (return.attributes)
    attributes(result) <- c(attributes(result), attributes)
  return(drop(result))
}


#' MONITORING STATIONS DATA - Sample from the random field locations to get monitoring stations 
#' 
#' @param prop.points.for.area - proportion of points to get per area
#' @param compute_Y_Poisson - simulated Poisson data
#' @param shapefile_sim - the shapefile
#' @param sim_grid_df - simulation grid
#' @param sparse - T/F, if monitoring stations data is sparse or not
#'
#' @return locations (coordinates) of the monitoring stations

get.monitoring.stations <- function(prop.points.for.area,
                                    compute_Y_Poisson,
                                    shapefile_sim,
                                    sim_grid_df,
                                    sparse){
  
  if (sparse == FALSE){
    
    if (!file.exists(paste0(getwd(),
                            "/Data",
                            "/monitoringstationsdata.sparse=FALSE.Rdata"))) {
      
      
      #--- Select the areas according to the sampling method ("all areas")
      #--- Each area has at least one monitoring station
      selected_area = shapefile_sim@data$AREAKEY 
      
      #Compute how many monitoring stations we have for each area (according to prop.points.for.area)
      npoints.for.area = round(table(compute_Y_Poisson$AREAKEY) * prop.points.for.area,0)
      
      #--- Select randomly a given number of locations within each selected area
      compute_Y_Poisson_reduced = compute_Y_Poisson[which(compute_Y_Poisson$AREAKEY %in% selected_area), ]
      compute_Y_Poisson_splitted = split(compute_Y_Poisson_reduced,
                                         list(compute_Y_Poisson_reduced$AREAKEY))
      #just keep elements with nrow>0 (remove areas with no points inside)
      area.nopoints = which(table(sim_grid_df$AREAKEY) == 0) #area without grid points inside
      compute_Y_Poisson_splitted = compute_Y_Poisson_splitted[sapply(compute_Y_Poisson_splitted, nrow) >0]
      
      #all the areas, sample a given % of points
      if (length(area.nopoints) > 0) { #if there are some areas without points inside
        samples <- lapply(1:length(compute_Y_Poisson_splitted),
                          function(x) {
                            compute_Y_Poisson_splitted[[x]][sample(1:nrow(compute_Y_Poisson_splitted[[x]]),
                                                                   npoints.for.area[-area.nopoints][x],
                                                                   replace = FALSE), ]})
      } else {
        samples <- lapply(1:length(compute_Y_Poisson_splitted),
                          function(x) {
                            compute_Y_Poisson_splitted[[x]][sample(1:nrow(compute_Y_Poisson_splitted[[x]]),
                                                                   npoints.for.area[x],
                                                                   replace = FALSE), ]})
      }
      
      #--- Define the df with the monitoring stations
      selected_loc_df <- do.call(rbind, samples)
      monitoringstations.coords = selected_loc_df[, c("x1", "x2")]
      monitoringstations.data = selected_loc_df
      
      save(file = paste0(getwd(),
                         "/Data",
                         "/monitoringstationsdata.sparse=FALSE.Rdata"),
           monitoringstations.coords, monitoringstations.data)
      
    } else { #retrieve just the stations coordinates and update the exposure values
      
      load(paste0(getwd(),
                  "/Data",
                  "/monitoringstationsdata.sparse=FALSE.Rdata"))
      
      #-- Merge the dataframe with all the exposure simulation points (sim_grid_inside_df2)
      #-- with the selected monitoring stats dataframe containing just the x1, x2 coordinates
      #-- This is required to select again the same monitoring stations (obviously with different exposure values)
      monitoringstations.data = merge(compute_Y_Poisson,
                                      monitoringstations.data[, c("x1","x2")], #just x1, x2
                                      by=c("x1","x2"), sort=F)
      
      monitoringstations.coords = monitoringstations.data[,c("x1","x2")]
      
    }
    
  }else{
    
    if (!file.exists(paste0(getwd(),
                            "/Data",
                            "/monitoringstationsdata.sparse=TRUE.Rdata"))) {
      
      
      #--- Select the areas according to the sampling method ("all areas")
      #--- Each area has at least one monitoring station
      selected_area <- sample(shapefile_sim@data$AREAKEY , size = 28)
      new_prop.points.for.area <- c(rep(prop.points.for.area,14),
                                    rep(prop.points.for.area*.5, 10),
                                    rep(prop.points.for.area*.25, 4))
      
      #--- Select randomly a given number of locations within each selected area
      compute_Y_Poisson_reduced = compute_Y_Poisson[which(compute_Y_Poisson$AREAKEY %in% selected_area), ]
      compute_Y_Poisson_splitted = split(compute_Y_Poisson_reduced,
                                         list(compute_Y_Poisson_reduced$AREAKEY))
      
      #Compute how many monitoring stations we have for each area (according to prop.points.for.area)
      npoints.for.area = round(table(compute_Y_Poisson_reduced$AREAKEY) * new_prop.points.for.area,0)
      
      samples <- lapply(1:length(compute_Y_Poisson_splitted),
                        function(x) {
                          compute_Y_Poisson_splitted[[x]][sample(1:nrow(compute_Y_Poisson_splitted[[x]]),
                                                                 npoints.for.area[x],
                                                                 replace = FALSE), ]})
      
      #--- Define the df with the monitoring stations
      selected_loc_df <- do.call(rbind, samples)
      monitoringstations.coords = selected_loc_df[, c("x1", "x2")]
      monitoringstations.data = selected_loc_df
      
      save(file = paste0(getwd(),
                         "/Data",
                         "/monitoringstationsdata.sparse=TRUE.Rdata"),
           monitoringstations.coords, monitoringstations.data)
      
    }else{
      
      
      load(paste0(getwd(),
                  "/Data",
                  "/monitoringstationsdata.sparse=TRUE.Rdata"))
      
      #-- Merge the dataframe with all the exposure simulation points (sim_grid_inside_df2)
      #-- with the selected monitoring stats dataframe containing just the x1, x2 coordinates
      #-- This is required to select again the same monitoring stations (obviously with different exposure values)
      monitoringstations.data = merge(compute_Y_Poisson,
                                      monitoringstations.data[, c("x1","x2")], #just x1, x2
                                      by=c("x1","x2"), sort=F)
      
      monitoringstations.coords = monitoringstations.data[,c("x1","x2")]
    
    }
  }
  
  return(monitoringstations.data)
}



#' COMPUTE INTERSECTIONS between the shapefile and the grid 
#' 
#' @param res_pred_grid - resolution of the prediction grid
#' @param compute_Y_Poisson - simulated Poisson data
#'
#' @return intersections of the grids to the areas

compute.intersections <- function(res_pred_grid, 
                                  compute_Y_Poisson){
  
  if (!file.exists(
    paste0(getwd(),"/Data","/weights_grid_areas_intersection_res_pred_grid_",res_pred_grid,".Rdata"))) {
    
    weights_grid_areas_intersection = data.frame()
    cat("--- I'm computing the intersection weight matrix ... it takes a while" , "\n")
    
    for (i in 1:length(shapefile_sim@data$AREAKEY)) {
      cat("--- Area n.",i,"---- of", length(shapefile_sim@data$AREAKEY), "\n" )
      # Select an area
      sel_area = shapefile_sim[i, ]
      if (!gIsValid(sel_area)) {
        sel_area = gBuffer(sel_area, width = 0, byid = TRUE)
      }
      
      # Compute the intersection between the selected area and the grid
      inter = gIntersection(pred_grid_poly, sel_area, byid = TRUE)
      if (is.null(inter))
        next #Skip if there is no intersection between the shapefile and the grid 
      
      # Define a data frame with the name+area intersection pixels and the weights given as area/total_area
      inter_w = data.frame(
        pixel = sapply(inter@polygons, function(x) x@ID),
        area = sapply(inter@polygons, function(x) x@area))
      inter_w$weights = inter_w$area / sum(inter_w$area)
      inter_w$AREAKEY = sel_area@data$AREAKEY
      # Check the areas (total area of the polygon - area given by the inserctions --> should be close to 0)
      abs_diff_area = abs(sum(sapply(inter@polygons, function(x) x@area)) - sel_area@polygons[[1]]@area)
      
      if (abs_diff_area  > 0.0005) {
        warning(
          paste("----------- ATTENTION: there's something wrong with the areas! for area i=",i,"The abs_diff_area is ", round(abs_diff_area, 6)))
      }
      weights_grid_areas_intersection = rbind(weights_grid_areas_intersection, inter_w)
    }
    save(file = paste0(getwd(),
                       "/Data",
                       "/weights_grid_areas_intersection_res_pred_grid_",res_pred_grid,".Rdata"),
         weights_grid_areas_intersection)
    
  } else {
    
    load(
      paste0(getwd(),"/Data",
             "/weights_grid_areas_intersection_res_pred_grid_",
             res_pred_grid,".Rdata"))
  }
  
  return(weights_grid_areas_intersection)
  
}



#' COMPUTE THE PRED POINTS INSIDE EACH AREA OR THE CLOSEST ONE 
#' 
#' @param res_pred_grid - resolution of the prediction grid
#' @param compute_Y_Poisson - simulated Poisson data
#'
#' @return values of pred points for each area

compute.predpoints <- function(res_pred_grid,
                               compute_Y_Poisson){
  
  if (!file.exists(
    paste0(getwd(),"/Data","/weights_points_areas_inside_closest_",res_pred_grid,".Rdata"))) {
    
    weights_points_areas_inside_closest = data.frame()
    cat("--- I'm computing the inside points weight matrix ... it takes a while" , "\n")
    
    for (i in 1:length(shapefile_sim@data$AREAKEY)) {
      
      sel_area = shapefile_sim[i, ]
      #Pred grid points inside the area
      pointsinside = which(!is.na(over(pred_grid,sel_area)$AREAKEY))
      
      if(length(pointsinside)>0){
        points_df = data.frame(pixel = pointsinside,
                               weights = rep(1/length(pointsinside),length(pointsinside)),
                               AREAKEY = shapefile_sim@data$AREAKEY[i])
      } else {
        nearestpoints = apply(gDistance(pred_grid, gCentroid(sel_area), byid=TRUE), 1, which.min)
        points_df = data.frame(pixel = nearestpoints,
                               weights = rep(1/length(nearestpoints),length(nearestpoints)),
                               AREAKEY = shapefile_sim@data$AREAKEY[i])
        
      }
      
      weights_points_areas_inside_closest = rbind(weights_points_areas_inside_closest, points_df)
    }
    
    save(file = paste0(getwd(),
                       "/Data",
                       "/weights_points_areas_inside_closest_",res_pred_grid,".Rdata"),
         weights_points_areas_inside_closest)
    
  } else {
    
    load(
      paste0(getwd(),"/Data",
             "/weights_points_areas_inside_closest_",
             res_pred_grid,".Rdata"))
  }
  
  return(weights_points_areas_inside_closest)
  
}



