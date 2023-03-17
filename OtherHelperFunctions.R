
#' PLOTS - produce plot of predicted exposures surface from the first stage spatio-temporal model fit 
#' 
#' @param inlamodel1 - output from inla method 1
#' @param inlamodel2 - output from inla method 2
#' @param stk - inla stack
#' @param n.T - number of timepoints
#' @param inputdata - input data
#' 
#' @return - plot of exposures surface

plot.exposures.surface.ST <- function(inlamodel1,
                                      inlamodel2,
                                      stk = stk,
                                      n.T = 1,
                                      inputdata = inputdata){
  
  pred_grid = inputdata$pred_grid
  
  for(i in 1:n.T){
    index_pred <- inla.stack.index(stk, paste0("pred.",i,sep=""))$data
    assign(paste0("index_pred.",i,sep=""),index_pred)
  }
  
  all <- bind_rows(inputdata$compile.risk.grid, .id = "column_label")
  all.preds <- c()
  for (i in 1:n.T){
    index.pred <- get(paste0("index_pred.",i,sep=""))
    temp1 <- inlamodel1$summary.linear.predictor[index.pred,'mean']
    temp2 <- inlamodel2$summary.linear.predictor[index.pred,'mean']
    all.preds <- c(all.preds,temp1,temp2)
  }
  
  all_vals = c(all.preds, all$exposure)
  scl = scale_fill_viridis_c(limits = c(min(all_vals), max(all_vals)))
  
  compile.plots <- vector("list",n.T*3)
  
  for (time.ind in 1:n.T){
    
    index.pred <- get(paste0("index_pred.",time.ind,sep=""))
    
    predicted.exposure.m1 <- cbind(pred_grid@coords,inlamodel1$summary.linear.predictor[index.pred,'mean'])
    colnames(predicted.exposure.m1) <- c("x1","x2","exposure")
    predicted.exposure_inside.m1 <- merge(y = predicted.exposure.m1,
                                          x = inputdata$compile.sim_grid_inside_df[which(inputdata$compile.sim_grid_inside_df$Time==time.ind),c('x1','x2')], by = c('x1','x2'), all.x = TRUE)
    
    predicted.exposure.m2 <- cbind(pred_grid@coords,inlamodel2$summary.linear.predictor[index.pred,'mean'])
    colnames(predicted.exposure.m2) <- c("x1","x2","exposure")
    predicted.exposure_inside.m2 <- merge(y = predicted.exposure.m2,
                                          x = inputdata$compile.sim_grid_inside_df[which(inputdata$compile.sim_grid_inside_df$Time==time.ind),c('x1','x2')], by = c('x1','x2'), all.x = TRUE)
    
    true <- ggplot() +
      geom_tile(aes(x = x1, y = x2, fill = exposure), data = inputdata$compile.risk.grid[[time.ind]]) +
      theme_bw() + 
      geom_path(aes(x = long, y = lat, group = group), data = inputdata$shapefile_sim, alpha=2) +
      coord_equal() + 
      theme_map() +
      labs(x="", y="", title=paste0("True exposures surface, Time = ",time.ind,sep="")) +
      theme(legend.key.width=unit(1, "cm")) + 
      theme(legend.position="bottom") +
      scl
    
    predicted.m1 <- ggplot() +
      geom_tile(aes(x = x1, y = x2, fill = exposure), data = predicted.exposure_inside.m1) +
      theme_map() +
      geom_path(aes(x = long, y = lat, group = group), data = inputdata$shapefile_sim, alpha=2) +
      coord_equal() + 
      labs(x="", y="", title=paste0("Predicted exposures surface (method 1), Time = ",time.ind)) +
      theme(legend.key.width=unit(1, "cm")) + 
      theme(legend.position="bottom") +
      scl
    
    predicted.m2 <- ggplot() +
      geom_tile(aes(x = x1, y = x2, fill = exposure), data = predicted.exposure_inside.m2) +
      theme_map() +
      geom_path(aes(x = long, y = lat, group = group), data = inputdata$shapefile_sim, alpha=2) +
      coord_equal() + 
      labs(x="", y="", title=paste0("Predicted exposures surface (method 2), Time = ",time.ind)) +
      theme(legend.key.width=unit(1, "cm")) + 
      theme(legend.position="bottom") +
      scl
    
    i <- (time.ind - 1) * 3 + 1
    compile.plots[[i]] <- true
    compile.plots[[i+1]] <- predicted.m1
    compile.plots[[i+2]] <- predicted.m2
    
  }
  
  out <- list(compile.plots = compile.plots)
  return(out)
}



#' PLOTS - produce scatterplot of true block exposures versus predicted block exposures
#' 
#' @param inla.m1 - output from inla using method 1
#' @param inla.m2 - output from inla using method 2
#' @param ST - TRUE/FALSE, if spatiotemporal model or not
#' @param n.T - number of timepoints
#' @param inputdata - input data
#' 
#' @return - scatterplots

scatterplot <- function(inla.m1, 
                        inla.m2, 
                        ST = FALSE, 
                        n.T = 1, 
                        inputdata){
  
  if (ST == FALSE){
    forscatter <- as.data.frame(cbind(inputdata$shapefile_sim_complete@data$true_area_exposure,
                                      inla.m1$summary.lincomb.derived$mean,
                                      inla.m2$summary.lincomb.derived$mean))
    names(forscatter) <- c("truepres", "predicted.m1", "predicted.m2")
    
    scatterplot.m1 <- ggplot(forscatter, aes(x=truepres, y=predicted.m1)) + 
      geom_point() +
      geom_abline(slope=1, intercept=0, col='red') +
      labs(x ="true area exposure", y = "predicted exposure") +
      ggtitle("Plot of true versus predicted area exposures \n (Method 1)") +
      theme(plot.title = element_text(hjust=0.5))
    scatterplot.m2 <- ggplot(forscatter, aes(x=truepres, y=predicted.m2)) + 
      geom_point() +
      geom_abline(slope=1, intercept=0, col='red') +
      labs(x ="true area exposure", y = "predicted exposure") +
      ggtitle("Plot of true versus predicted area exposures \n (Method 2)") +
      theme(plot.title = element_text(hjust=0.5))
    
    out <- list(scatterplot.m1 = scatterplot.m1,
                scatterplot.m2 = scatterplot.m2)
    
    return(out)
  }
  
  if (ST == TRUE){
    
    compile.m1 <- vector("list",n.T)
    compile.m2 <- vector("list",n.T)
    
    for (i in 1:n.T){
      
      forscatter <- as.data.frame(cbind(inputdata$shapefile_sim_complete@data[which(inputdata$shapefile_sim_complete$Time==i),]$true_area_exposure,
                                        inla.m1$summary.lincomb.derived$mean[((i-1)*98+1):(98*i)],
                                        inla.m2$summary.lincomb.derived$mean[((i-1)*98+1):(98*i)]))
      names(forscatter) <- c("truepres", "predicted.m1", "predicted.m2")
      
      scatterplot.m1 <- ggplot(forscatter, aes(x=truepres, y=predicted.m1)) + 
        geom_point() +
        geom_abline(slope=1, intercept=0, col='red') +
        labs(x = paste0("true area exposure, Time = ",i,sep=""), y = paste0("predicted exposure, Time = ",i,sep="")) +
        ggtitle(paste0("Plot of true versus predicted area exposures, Time = ",i, "\n (Method 1)",sep="")) +
        theme(plot.title = element_text(hjust=0.5))
      scatterplot.m2 <- ggplot(forscatter, aes(x=truepres, y=predicted.m2)) + 
        geom_point() +
        geom_abline(slope=1, intercept=0, col='red') +
        labs(x =paste0("true area exposure, Time = ",i,sep=""), y = paste0("predicted exposure, Time = ",i,sep="")) +
        ggtitle(paste0("Plot of true versus predicted area exposures, Time = ",i, "\n (Method 2)",sep="")) +
        theme(plot.title = element_text(hjust=0.5))
      
      compile.m1[[i]] <- scatterplot.m1
      compile.m2[[i]] <- scatterplot.m2
      
    }
    
    out <- list(compile.m1 = compile.m1,
                compile.m2 = compile.m2)
    return(out)
  }
  
  
}



#' PLOTS - produce plots of estimated marginal posterior distributions for stage 2 model parameters - method 1 vs method 2  
#' 
#' @param sample1 - posterior samples from method 1
#' @param sample2 - posterior samples from method 2
#' @param true.value - true value of the parameter 
#' @param title - plot title
#' @param xlab - x-axis label 
#' 
#' @return - posterior density plots 

combined.post.dens.plot.FF <- function(sample1 = posterior_gamma0_FF_method1,
                                       sample2 = posterior_gamma0_FF_method2,
                                       true.value = gamma0_Poisson,
                                       title = expression(paste(hat(pi),"(",gamma[0],"|y)",)),
                                       xlab = expression(gamma[0])){
  
  sample1 <- as.data.frame(cbind(as.numeric(sample1),method = 1))
  names(sample1)[1] <- 'x'
  sample1['method'][sample1['method'] == 1] <- "Method 1"
  
  sample2 <- as.data.frame(cbind(as.numeric(sample2),method = 2))
  names(sample2)[1] <- 'x'
  sample2['method'][sample2['method'] == 2] <- "Method 2"
  
  combined <- rbind(sample1, sample2)
  
  median <- ddply(combined, "method", summarise, value=quantile(x,probs=0.5))
  q2.5 <- ddply(combined, "method", summarise, value=quantile(x,probs=0.025))
  q97.5 <- ddply(combined, "method", summarise, value=quantile(x,probs=0.975))
  
  ggplot(combined, aes(x=x, color=method)) +
    geom_density() + 
    theme_bw() +
    geom_vline(aes(xintercept = true.value),color='black') +
    geom_vline(data=median, aes(xintercept=value, color=method)) + 
    geom_vline(data=q2.5, aes(xintercept=value, color=method),
               linetype="dashed") + 
    geom_vline(data=q97.5, aes(xintercept=value, color=method),
               linetype="dashed") +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab(xlab) +
    ylab('')
  
}



#' PLOTS - produce plots of estimated marginal posterior distributions for stage 1 model parameters
#' 
#' @param res1.param - paramater result in res1 object - method 1
#' @param res2.param - paramater result in res1 object - method 2
#' @param true.value - true value of the parameter 
#' @param title - plot title
#' @param xlab - x-axis label 
#' @param trans - transformation to be applied on the parameter
#' @param prec.param - TRUE/FALSE, whether the parameter is a precision parameter or not
#' @param n_random_fromposterior - number of samples to simulate from the posterior distribution
#' 
#' @return - posterior density plots 

post.dens.plot.stage1 <- function(res1.param = res1$marginals.hyperpar$`Stdev for s`,
                                  res2.param = res2$marginals.hyperpar$`Stdev for s`,
                                  true.value = sigma2xi,
                                  title = expression(paste(hat(pi),"(",sigma[xi]^2,"|y)",)),
                                  xlab = expression(sigma[xi]^2),
                                  trans = 'reciprocal',
                                  prec.param = FALSE,
                                  n_random_fromposterior = 200){
  
  if (prec.param == FALSE){
    res1.samples <- inla.rmarginal(n_random_fromposterior, res1.param)
    res2.samples <- inla.rmarginal(n_random_fromposterior, res2.param)
    combined <- c(res1.samples, res2.samples)
  }else{
    if(trans == 'reciprocal'){
      res1.trans <- inla.tmarginal(function(x) 1/x, res1.param)
      res2.trans <- inla.tmarginal(function(x) 1/x, res2.param)
      res1.samples <- inla.rmarginal(n_random_fromposterior, res1.trans)
      res2.samples <- inla.rmarginal(n_random_fromposterior, res2.trans)
      combined <- c(res1.samples, res2.samples)
    }else{
      res1.trans <- inla.tmarginal(function(x) x^2, res1.param)
      res2.trans <- inla.tmarginal(function(x) x^2, res2.param)
      res1.samples <- inla.rmarginal(n_random_fromposterior, res1.trans)
      res2.samples <- inla.rmarginal(n_random_fromposterior, res2.trans)
      combined <- c(res1.samples, res2.samples)
    }
  }
  
  median <- quantile(combined,probs=0.5)
  q2.5 <- quantile(combined,probs=0.025)
  q97.5 <- quantile(combined,probs=0.975)
  
  
  ggplot(as.data.frame(combined), aes(x=combined)) +
    geom_density() + 
    theme_bw() +
    geom_vline(aes(xintercept = true.value),color='blue') +
    geom_vline(aes(xintercept=median),color='red') + 
    geom_vline(aes(xintercept=q2.5),color='red',
               linetype="dashed") + 
    geom_vline(aes(xintercept=q97.5),color='red',
               linetype="dashed") +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab(xlab) +
    ylab('')
  
}



#This function computes the mean bias, RMSE and coverage given some values generated
#from the posterior distribution of a PARAMETER and the true values
#' 
#' @param vec_simulatedvalues - compiled posterior estimates
#' @param n_simulations - number of simulations 
#' @param pseudotrue_vec - vector of true values
#' 
#' @return - performance indices

my_bias_RMSE_coverage = function(vec_simulatedvalues, n_simulations, pseudotrue_vec){
  
  matrix_simulatedvalues = matrix(vec_simulatedvalues, ncol=n_simulations)
  bias = rmse = cov = c()
  
  for(i in 1:ncol(matrix_simulatedvalues)){
    #Relative indexes:
    #bias[i] = mean( (matrix_simulatedvalues[,i] - pseudotrue_vec[i]) /pseudotrue_vec[i] )
    #rmse[i] = sqrt(mean(((matrix_simulatedvalues[,i] - pseudotrue_vec[i])/pseudotrue_vec[i])^2))
    
    #Absolute indexes:
    bias[i] = mean( matrix_simulatedvalues[,i] - pseudotrue_vec[i] )
    rmse[i] = sqrt(mean( (matrix_simulatedvalues[,i] - pseudotrue_vec[i])^2 ))
    
    q0025 = quantile(matrix_simulatedvalues[,i],0.025)
    q0975 = quantile(matrix_simulatedvalues[,i],0.975)
    cov[i] = (q0025 < pseudotrue_vec[i]  & pseudotrue_vec[i] < q0975)
  }
  return(list(bias=mean(bias),
              rmse=mean(rmse),
              cov=sum(cov)/n_simulations*100))
}



# This function computes the mean bias and rmse in a vector form (length given by the number of simulations)
#'
#' @param vec_simulatedvalues - compiled posterior estimates
#' @param n_simulations - number of simulations 
#' @param pseudotrue_vec - vector of true values
#' 
#' @return - performance indices in vector

my_bias_RMSE_vec = function(vec_simulatedvalues, n_simulations, pseudotrue_vec){
  matrix_simulatedvalues = matrix(vec_simulatedvalues, ncol=n_simulations)
  bias = rmse = c()
  
  for(i in 1:ncol(matrix_simulatedvalues)){
    bias[i] = mean( matrix_simulatedvalues[,i] - pseudotrue_vec[i] )
    rmse[i] = sqrt(mean( (matrix_simulatedvalues[,i] - pseudotrue_vec[i])^2 ))
  }
  return(data.frame(bias=bias,
                    rmse=rmse))
}



#' PLOTS - this is just another plot for posterior distributions of stage 1 model parameters
#' 
#' @param samples - compiled samples (could be from several replications)
#' @param true.value - true value of the parameter 
#' @param title - plot title
#' @param xlab - x-axis label 
#' 
#' @return - posterior density plots 

post.dens.plot.stage1.ALLSIM <- function(samples = posterior_sigma2xi_compiled,
                                         true.value = sigma2xi,
                                         title = expression(paste(hat(pi),"(",sigma[xi]^2,"|y)",)),
                                         xlab = expression(sigma[xi]^2)){
  
  median <- quantile(samples,probs=0.5)
  q2.5 <- quantile(samples,probs=0.025)
  q97.5 <- quantile(samples,probs=0.975)
  
  ggplot(as.data.frame(samples), aes(x=samples)) +
    geom_density() + 
    theme_bw() +
    geom_vline(aes(xintercept = true.value),color='blue') +
    geom_vline(aes(xintercept=median),color='red') + 
    geom_vline(aes(xintercept=q2.5),color='red',
               linetype="dashed") + 
    geom_vline(aes(xintercept=q97.5),color='red',
               linetype="dashed") +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab(xlab) +
    ylab('')
  
}



#' PLOTS - plot bias + rmse + coverage for exposure area prediction
#' 
#' @param values.method1 - performance values from using method 1
#' @param values.method1 - performance values from using method 2
#' @param simulated.data - simulated data
#' @param title - plot title
#' @param style - color style of plot
#' 
#' @return - plots
#' 

plot.performance.area.preds <- function(values.method1,
                                        values.method2,
                                        simulated.data = simulated.data,
                                        title,
                                        style,
                                        n.T,
                                        scenario_n){
  
  all_vals = c(values.method1,values.method2)
  
  plots.method1 <- vector("list", n.T)      
  plots.method2 <- vector("list", n.T)
  
  for (i in 1:n.T){
    
    sub.values.method1 <- values.method1[(98*(i-1)+1):(98*i)]
    sub.values.method2 <- values.method2[(98*(i-1)+1):(98*i)]
    
    #--- Prepare the data
    shapefile = simulated.data$shapefile_sim_complete
    temp = simulated.data$shapefile_sim_complete@data[(98*(i-1)+1):(98*i),]
    shapefile@data <- temp
    forplot <- shapefile
    forplot@data$value1 = sub.values.method1
    forplot@data$value2 = sub.values.method2
    
    forplot_tidy <- tidy(forplot)
    forplot$id <- row.names(forplot)
    forplot <- left_join(forplot_tidy, forplot@data)
    
    method1 <- ggplot(forplot, aes(x = long, y = lat, group = group, fill = value1)) +
      geom_polygon(color = "black", size = 0.1) +
      coord_equal(ratio=1) +
      labs(x="", y="", title=paste0(title,", Scenario = ",scenario_n,", (method 1), Time = ",i,sep="")) +
      theme(legend.key.width=unit(1, "cm")) + 
      theme_map() +
      scale_fill_viridis(option = style, limits = c(min(all_vals), max(all_vals))) +
      theme(legend.position="bottom")
    
    method2 <- ggplot(forplot, aes(x = long, y = lat, group = group, fill = value2)) +
      geom_polygon(color = "black", size = 0.1) +
      coord_equal(ratio=1) +
      labs(x="", y="", title=paste0(title,", Scenario = ",scenario_n,", (method 2), Time = ",i,sep="")) +
      theme(legend.key.width=unit(1, "cm")) + 
      theme_map() +
      scale_fill_viridis(option = style, limits = c(min(all_vals), max(all_vals))) +
      theme(legend.position="bottom")
    
    plots.method1[[i]] <- method1
    plots.method2[[i]] <- method2 
    
  }
  
  return(list(plots.method1 = plots.method1,
              plots.method2= plots.method2))
}

