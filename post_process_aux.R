###
### This file contains most of the BEESTS post processing routines.
###
# Original Authors: Matzke et al
# Substantial Update: Daniel Scott
#
library(grid)
library(gtable)

# ------------------------------------------------------------------------------ #
# A function for checking silly input values
# ------------------------------------------------------------------------------ #
check_silly_things <- function(samples, burn.in, thinning, number.of.chains, posterior.predictor.samples, int_lower, int_upper, priors){
   if (samples < 1){ stop("The total number of MCMC samples must be greater than zero.") }
   if (samples <= burn.in){ stop("The total number of MCMC samples must be greater than the number of burn-in samples.")}
   if (thinning < 1){ stop("The thinning factor must be greater than 0.")}
   if (((samples-burn.in)/thinning) < 1){ stop("No MCMC samples will be retained. Increse the number of retained samples or decrese the thinning factor.")}
   if (post.preds == T){
      if ((((samples-burn.in)/thinning)*number.of.chains) < post.preds){
           stop("The number of posterior predictive samples cannot be greater than the number of retained MCMC samples.")
      }
   }
   if (int_lower>=int_upper){stop("The lower limit of integration must be lower than the upper limit.")}

   #for (par in priors){
   #   if  (as.numeric(as.vector(vars[vars[,1]==paste(par, " upper",sep=""),2])) <= as.numeric(as.vector(vars[vars[,1]==paste(par, " lower",sep=""),2]))) {
   #     if (par=="stop pf"){
   #              stop("Silly prior: Upper bound for P(TF) must be higher than lower bound")
   #      } else {
   #              stop(paste(paste("Silly prior: Upper bound for ", par,sep=""), " must be higher than lower bound.", sep=""))
   #      }
   #   }

   #   if  ((as.numeric(as.vector(vars[vars[,1]==paste(par, " start",sep=""),2])) < as.numeric(as.vector(vars[vars[,1]==paste(par, " lower",sep=""),2]))) | (as.numeric(as.vector(vars[vars[,1]==paste(par, " start",sep=""),2])) > as.numeric(as.vector(vars[vars[,1]==paste(par, " upper",sep=""),2])))){
   #   if (par=="stop pf"){
   #              stop("Silly start value: Start value for P(TF) must be within prior range")
   #      } else {
   #              stop(paste(paste("Start value for ", par,sep=""), " must be within prior range.", sep=""))
   #    }
   #   }
   #}

   #if (ncol==4){
   #  if (pf_lower < 0 | pf_upper > 1){
   #        stop("The range of the uniform prior for P(TF) must be between 0 and 1")
   #    }
   #}
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
### Load python traces in R and make array for further use
# ------------------------------------------------------------------------------ #
read_prep = function(dir, n_chains) {

  for(i in 1:n_chains){

    # Read data, convert to matrix:
    data      <- read.csv(file = paste(dir, "/parameters", i, ".csv", sep=""), head=TRUE, sep=";")
    col_names <- names(data)

    # The data contains an indexing column, assigned colname 'X'...
    #data <- data[setdiff(col_names, 'X')]
    #col_names <- names(data)

    # Remove 'pt' from 'pf_stop_subjpt.X' column names
    # It's a hack, sorry...
    chain <- as.matrix(data[,order(col_names)])
    col_names <- colnames(chain)
    colnames(chain) <- gsub('pt', '', col_names)

    regex <- paste('*_subj.*', sep='')
    subj_params <- grep(glob2rx(regex), col_names, value=TRUE)
    subjects <- sort(strtoi(unique( gsub('.+_subj.', '', subj_params) )))
    n_subject <- length(subjects)

    cat(sprintf('read_prep() has determined there are %s subjects in parameter file %d', toString(n_subject), i), '\n')

    first_subj <- min(subjects)
    last_subj  <- max(subjects)
    cat(sprintf('... with subject id %s:%s less %s', first_subj, last_subj, toString(setdiff(first_subj:last_subj,subjects))), '\n')

    if (i == 1) {
      col_names <- dimnames(chain)[[2]]
      n_sample  <- nrow(data)
      traces    <- array(NA, dim=c(dim(data),n_chains), dimnames=list(1:n_sample, col_names,paste("chain", 1:n_chains)))
    }
    traces[,,i] = chain

  }
  return(list(traces = traces, n_subject = n_subject, my_names = col_names, n_samples = n_sample, n_chains = n_chains, subjects = subjects))
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# Computes summary stats for each parameter (by participant, collapsed over chains) and save output to csv.
# ------------------------------------------------------------------------------ #
summary_stats <- function(traces, params, n_subj, subj_idx = NULL){
  library(abind)

  # In both subject and group cases, we need a row for each param
  summary_rows <- params

  # In subject level data we also want some exta summary stuff
  if (!is.null(subj_idx)) {summary_rows <- c(params, "mean go", "sd go", "mean SSRT", "sd SSRT")}

  # Create the matrix to hold summary statistics:
  prctiles     <- c("2.5%", "25%", "50%", "75%", "97.5%")
  summary_cols <- c("Mean", "Sd", prctiles)

  summary = matrix(NA,length(summary_rows), length(summary_cols))
  colnames(summary) = summary_cols
  rownames(summary) = summary_rows

  if (!is.null(subj_idx)) {
    # Extract lists of param names by type:
    mu_go_subj      = grep(glob2rx("mu_go*"   ), params, value=TRUE)
    sigma_go_subj   = grep(glob2rx("sigma_go*"), params, value=TRUE)
    tau_go_subj     = grep(glob2rx("tau_go*"  ), params, value=TRUE)

    mu_stop_subj    = grep(glob2rx("mu_stop*"   ), params, value=TRUE)
    sigma_stop_subj = grep(glob2rx("sigma_stop*"), params, value=TRUE)
    tau_stop_subj   = grep(glob2rx("tau_stop*"  ), params, value=TRUE)

    go     <-      traces[,     mu_go_subj,]   + traces[,  tau_go_subj,]
    goSd   <- sqrt(traces[,  sigma_go_subj,]^2 + traces[,  tau_go_subj,]^2)
    ssrt   <-      traces[,   mu_stop_subj,]   + traces[,tau_stop_subj,]
    ssrtSd <- sqrt(traces[,sigma_stop_subj,]^2 + traces[,tau_stop_subj,]^2)

    traces <- abind( traces , go, ssrt, goSd, ssrtSd, along = 2 )
    new_names <- dimnames(traces)
    addnl_names <- paste(c('go_subj', 'ssrt_subj', 'goSd_subj', 'ssrtSd_subj'), subj_idx, sep='.')

    ncols <- dim(traces)[2]
    new_names[[2]][(ncols-3):ncols] <- addnl_names
    dimnames(traces) <- new_names

    meanGo   =      as.vector(traces[,mu_go_subj,])        + as.vector(traces[,tau_go_subj,])
    sdGo     = sqrt(as.vector(traces[,sigma_go_subj,])^2   + as.vector(traces[,tau_go_subj,])^2)
    meanSRRT =      as.vector(traces[,mu_stop_subj ,])     + as.vector(traces[,tau_stop_subj ,])
    sdSRRT   = sqrt(as.vector(traces[,sigma_stop_subj,])^2 + as.vector(traces[,tau_stop_subj ,])^2)

    summary["mean go",1]   = round(mean(meanGo),4)
    summary["mean go",2]   = round(sd(meanGo),4)
    summary["mean go",3:7] = round(quantile(meanGo),4)

    summary["sd go",1] = round(mean(sdGo),4)
    summary["sd go",2] = round(sd(sdGo),4)
    summary["sd go",3:7] = round(quantile(sdGo),4)

    summary["mean SSRT",1] = round(mean(meanSRRT),4)
    summary["mean SSRT",2] = round(sd(meanSRRT),4)
    summary["mean SSRT",3:7] = round(quantile(meanSRRT),4)

    summary["sd SSRT",1] = round(mean(sdSRRT),4)
    summary["sd SSRT",2] = round(sd(sdSRRT),4)
    summary["sd SSRT",3:7] = round(quantile(sdSRRT),4)
  }

  pf_stop_subj = grep(glob2rx("pf_stop*"), params, value=TRUE)

  row = 0
  for (param in params){
    row <- row + 1
    if (!any(param == colnames(traces[,,1]))) {next}

    #if (col == pf_stop_subj & !is.null(subj_idx)) {val <- pnorm(traces[,col,])} else {val <- traces[,col,]}
    #if (param == pf_stop_subj) {next}
    val <- traces[,param,]
    summary[row,1]   <- round(    mean(as.vector(val)),4)
    summary[row,2]   <- round(      sd(as.vector(val)),4)
    summary[row,3:7] <- round(quantile(as.vector(val), prob = c(0.025,0.25,0.5,0.75,0.975)),4)
  }

  # Write table
  table <- tableGrob(summary)

  #title <- textGrob("Summary statistics group level parameters", y=unit(0.5,"npc") + 0.5*h, vjust=0, gp=gpar(fontsize=20))
  if (is.null(subj_idx)) {
    level <- 'group'
  } else {
    level <- paste('subject ', toString(subj_idx), sep='')
  }

  title    <- paste("Summary statistics for", level, "parameters", sep = ' ')
  title    <- textGrob(title, gp=gpar(fontsize=20))
  footnote <- textGrob("footnote", x=0, hjust=0, gp=gpar( fontface="italic"))

  padding <- unit(1,"line")

  table <- gtable_add_rows(table, heights = grobHeight(title) + padding, pos = 0)
  table <- gtable_add_rows(table, heights = grobHeight(footnote)+ padding)

  table <- gtable_add_grob(table, list(title, footnote), t=c(1, nrow(table)), l=c(1,2), r=ncol(table))

  #gt <- gTree(children=gList(table, title))
  grid.newpage()
  grid.draw(table)

  return(traces)
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
#
# ------------------------------------------------------------------------------ #
plot_posteriors_wrapper = function(mcmc_samples, priors){

  regex <- paste('*_subj.', '*', sep='')
  params_subj  <- grep(glob2rx(regex), colnames(mcmc_samples$traces), value=TRUE)
  params_group <- setdiff(colnames(mcmc_samples$traces), c(params_subj,'X'))

  if (mcmc_samples$n_subject==1) {
    # Individual posteriors
    # 'params_group' is used because there are no *_subj.[num] suffixes in non
    # hierarchical model.
    cat(sprintf('Plotting posteriors'), '\n')
    plot_posteriors(mcmc_samples, priors, params_group, plot_priors = TRUE)

    cat(sprintf('Plotting chains'), '\n')
    chains <- mcmc_samples$traces[,params_group,]
    plot_chains(chains, priors, params_group, probit = FALSE)

  }
  else {
    # Group posteriors
    cat(sprintf('Plotting group posteriors'), '\n')
    plot_posteriors(mcmc_samples, priors, params_group, plot_priors = TRUE, title_addendum = 'for Group')

    cat(sprintf('Plotting chains for group params'), '\n')
    chains <- mcmc_samples$traces[,params_group,]
    plot_chains(chains, priors, params_group, probit = FALSE)

    cat(sprintf('Plotting diagnostics for group params'), '\n')
    plot_diagnostics(mcmc_samples, params_group, nrows = 6, ncols = 3)

    #cat(sprintf('Plotting histogram of individual posterior means'), '\n')
    #plot_individual_posterior_means(chains, params_group)

    for (subj_num in mcmc_samples$subjects){
      cat(sprintf('Plotting chains, posteriors, and diagnostics for subject %s', subj_num), '\n')

      # Individual posteriors without prior plotting
      regex <- paste('*_subj.', toString(subj_num), sep='')
      params_subj_n  <- grep(glob2rx(regex), params_subj, value=TRUE)
      #cat(sprintf('Parameters found for subject %s: %s', subj_num, toString(params_subj_n)), '\n')

      plot_posteriors(mcmc_samples, priors, params_subj_n, plot_priors = FALSE, title_addendum = c('for Subject ', toString(subj_num)))

      chains <- mcmc_samples$traces[,params_subj_n,]
      plot_chains(chains, priors, params_subj_n, subject_idx = subj_num, probit = FALSE)

      plot_diagnostics(mcmc_samples, params_subj_n, nrows = 4, ncols = 4)
    }
  }
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# Geweke, Gelman-Rubin Statistics, and ...
# ------------------------------------------------------------------------------ #
plot_diagnostics = function(mcmc_samples, params, nrows, ncols){
  library(coda)
  mcmc_object <- as.mcmc.list(lapply(as.data.frame(mcmc_samples$traces[,params,]),mcmc))

  # Geweke
  layout(matrix(1:(nrows*ncols), nrows, ncols, byrow = T))
  par(cex.main=1.4)
  #browser()
  #geweke.plot(mcmc_object, auto.layout = FALSE)

  # Gelman-rubin statistics
  layout(matrix(1:(nrows*ncols), nrows, ncols, byrow = T))
  par(cex.main=1.4)
  for (param in params) {
    mcmc_object <- as.mcmc.list(lapply(as.data.frame(mcmc_samples$traces[,param,]),mcmc))
    gelman.plot(mcmc_object, auto.layout = FALSE)
    title(param)
  }

  # Auto-Correlations
  layout(matrix(1:(nrows*ncols), nrows, ncols, byrow = T))
  par(cex.main=1.4)
  for (param in params) {
    mcmc_object <- as.mcmc.list(lapply(as.data.frame(mcmc_samples$traces[,param,]),mcmc))
    maxlag = 500 #dim(mcmc_samples$traces)[1]
    acf(mcmc_samples$traces[,param,1], lag.max = maxlag, main = '')
    title(param)
  }
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# Posteriors
# ------------------------------------------------------------------------------ #
plot_posteriors = function(mcmc_samples, priors, params, plot_priors, title_addendum = ''){
  library(msm)

  if (plot_priors) {
    # Set prior densities
    prior_densities <- list()
    for (name in params) {

      # Model contains SD terms, but output on parameter distribution uncertainties is variances.
      name <- gsub("_var", "_sd", name)
      lower_lim <- priors[[name]][1]
      upper_lim <- priors[[name]][2]

      name_var  <- gsub("_sd", "_var", name)

      cat(sprintf('Setting %s prior density using start:%s, lower:%s, upper:%s', name_var, toString(lower_lim), toString(lower_lim), toString(upper_lim) ), '\n')
      prior_densities[[name_var]] <- dunif(lower_lim, lower_lim, upper_lim)
    }

    # An exception to the above uniform dists
    p_lims_prior = seq(priors$pf_stop[1],priors$pf_stop[2],0.1)
    prior_densities$pf_stop = dtnorm(p_lims_prior, priors$pf_stop[3], priors$pf_stop[4], lower=priors$pf_stop[1], upper=priors$pf_stop[2])
  }

  n_params <- length(params)
  n_rows   <- ceiling(n_params/4)
  # Setup plotting
  #par(mar = c(6, 5, 6, 3) + 0.1) # group
  #par(mar = c(8, 5, 6, 3) + 0.1) # individual
  layout(matrix(1:(n_rows*4), n_rows, 4, byrow=T))
  par(cex.main=1.2)

  # Plot posteriors
  for(par in params){

    if (par == 'pf_stop_var' | par == 'pf_stop') {next}

    #title   <- paste("Posterior", par, title_addendum, sep=" ")
    x_lower <- min(mcmc_samples$traces[,par,])-30
    x_upper <- max(mcmc_samples$traces[,par,])+30
    xlimits <- c(x_lower, x_upper)

    posterior <- density(mcmc_samples$traces[,par,])
    plot(posterior, xlim = xlimits, main = par, xlab=paste(par,"(ms)",sep=" "), ylab = "Density")

    # We plot param priors for (individual & not hierarchical) or group get plotted
    if (plot_priors) {
      # These are mostly just uniform lines
      curves_to_plot <- c(prior_densities[params == par], prior_densities[params == par])
      lines(c(x_lower, x_upper), curves_to_plot , lty=2)
    }
  }

  #plot(density(mcmc_samples$traces[,"pf_stop",]),xlim = c(priors$pf_stop[1], priors$pf_stop[2]), main = "Posterior pf_stop",xlab="pf_stop on probit scale",ylab = "Density")
  #plot(density(mcmc_samples$traces[,"pf_stop_sd",]),xlim = c(priors$pf_stop_sd[1],priors$pf_stop_sd[2]), main = "Posterior pf_stop_sd",xlab="pf_stop_sd on probit scale",ylab = "Density")

  if (plot_priors) {
    # These are mostly just uniform lines
    #lines(p_lims_prior,as.numeric(unlist(prior_den_group[names(prior_den_group)=="pf_stop"])),lty=2)
    #lines(c(priors$pf_stop_sd[1],priors$pf_stop_sd[2]),c(prior_den_group[names(prior_den_group)=="pf_stop_sd"],prior_den_group[names(prior_den_group)=="pf_stop_sd"]),lty=2)
  }
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# Posterior means...
# ------------------------------------------------------------------------------ #
plot_individual_posterior_means = function(chains, n_subj){

  a_big_enough_num <- 20
  post_means_ind   <- array(NA,dim=c(n_subj, a_big_enough_num))

  # Get posterior means
  for (subj_num in 1:n_subj) {
    regex <- paste('*_subj.', toString(subj_num), sep='')
    params = grep(glob2rx(regex), colnames(chains), value=TRUE)

    n_params <- length(params)
    for(index in 1:n_params){
        post_means_ind[subj_num, index] = mean(chains[,params[index],])
    }
  }
  post_means_ind <- post_means_ind[, 1:n_params]
  colnames(post_means_ind) = gsub('_subj', params)

  #par(mar = c(8, 5, 6, 3) + 0.1)
  #layout(matrix(1:n_params, 4, byrow=T))
  par(cex.main=1.4)

  for (index in 1:n_params){
    data  <- post_means_ind[,index]
    param <- params[index]
    hist(data, main=param, xlab=param)
  }
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# MCMC Chains
# ------------------------------------------------------------------------------ #
plot_chains = function(chains, priors, params, subject_idx = NULL, probit){

  n_samples <- dim(chains)[1]
  n_chains  <- dim(chains)[3]
  n_params  <- length(params)

  # Plot setup for group data
  #if (is.null(subject_idx)) {
    #par(mar = c(6, 5, 6, 3) + 0.1)
  #} else {
   # par(mar = c(8, 5, 6, 3) + 0.1)
  #}

  n_rows <- ceiling(n_params/4)
  layout(matrix(1:(n_rows*4), n_rows, 4, byrow=T))
  par(cex.main=1.4)

  pf_stop = c(priors$pf_stop[1],priors$pf_stop[2])
  #if (probit) {
  #  pf_stop_subj = grep(glob2rx("pf_stop*"), params, value=TRUE)
  #}

  for(par in params){
    if (par == 'pf_stop' | par == 'pf_stop_sd') {next}
    lim = c(min(chains[,par,])-10, max(chains[,par,])+10)


    plot(1:n_samples, chains[,par,1], ylim = lim,
      xlim = c(1,n_samples), main = par, xlab="Iteration",
      ylab = par,type="l")

    if(n_chains > 1){
       #if multiple chains, draw lines
      for(j in 2:n_chains){
        lines(1:n_samples,chains[,par,j],col=j)
      }
    }
  }

  #if (!probit) {
  #  plot(1:n_samples, chains[,"pf_stop",1], ylim = c(pf_stop[1],pf_stop[2]), xlim = c(1,n_samples), main = "MCMC chains pf_stop",xlab="Iteration",ylab = "pf_stop",type="l")
  #} else {
  #  plot(1:n_samples, pnorm(chains[,pf_stop_subj,1]),ylim = c(0,1), xlim = c(1,n_samples), main = paste("MCMC chains",pf_stop_subj,sep=" "), xlab="Iteration",ylab = pf_stop_subj,type="l")
  #}

  #if(n_chains>1){  #if multiple chains, draw lines
  #    for(j in 2:n_chains){
  #      if (!probit) {
  #        lines(1:n_samples,chains[,"pf_stop",j],col=j)
  #      } else {
  #        lines(1:n_samples,pnorm(chains[,pf_stop_subj,j]),col=j)
  #      }
  #    }
  #}
}
# ------------------------------------------------------------------------------ #
