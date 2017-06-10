###
### This file contains most of the BEESTS post processing routines.
###
# Original Authors: Matzke et al
# Substantial Update: Daniel Scott
#

# ------------------------------------------------------------------------------ #
# A function for checking silly input values
# ------------------------------------------------------------------------------ #
check_silly_things <- function(samples, burn.in, thinning, number.of.chains, posterior.predictor.samples, int_lower, int_upper, priors){
   if (samples < 1){ stop("The total number of MCMC samples must be greater than zero.") }
   if (samples <= burn.in){ stop("The total number of MCMC samples must be greater than the number of burn-in samples.")}
   if (thinning < 1){ stop("The thinning factor must be greater than 0.")}
   if (((samples-burn.in)/thinning) < 1){ stop("No MCMC samples will be retained. Increse the number of retained samples or decrese the thinning factor.")}
   if (posterior.predictors == T){
      if ((((samples-burn.in)/thinning)*number.of.chains) < posterior.predictors.samples){
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
read_prep = function(dir, n_chains, n_params) {

  for(i in 1:n_chains){

      # Read data, convert to matrix:
    data  <- read.csv(file = paste(dir, "/parameters", i, ".csv", sep=""),head=TRUE,sep=";")
    chain <- as.matrix(data[,order(names(data))])

      # Get number of subjects
      # Magic numbers:
      # 1: The number of columns which are not parameter values
      # 2: The number of columns used to account for the group mean and std of each param
    n_subject <-(ncol(chain) - 1) / n_params - 2

    if (i == 1) {
      col_names <- dimnames(chain)[[2]]
         n_sample  <- nrow(data)
      traces    <- array(NA, dim=c(dim(data),n_chains), dimnames=list(1:n_sample, col_names,paste("chain", 1:n_chains)))
    }
    traces[,,i] = chain

  }
  return(list(traces = traces, n_subject = n_subject, my_names = col_names, n_samples = n_sample, n_chains = n_chains))
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# Computes summary stats for each parameter (by participant, collapsed over chains) and save output to csv.
# ------------------------------------------------------------------------------ #
summary_stats <- function(pars, params, n_subj, subj_idx = NULL){

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

    meanGo   =      as.vector(pars$traces[,mu_go_subj,])        + as.vector(pars$traces[,tau_go_subj,])
    sdGo     = sqrt(as.vector(pars$traces[,sigma_go_subj,])^2   + as.vector(pars$traces[,tau_go_subj,])^2)
    meanSRRT =      as.vector(pars$traces[,mu_stop_subj ,])     + as.vector(pars$traces[,tau_stop_subj ,])
    sdSRRT   = sqrt(as.vector(pars$traces[,sigma_stop_subj,])^2 + as.vector(pars$traces[,tau_stop_subj ,])^2)

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
  for(row in 1:length(params)){
    col <- params[[row]]

    if (col == pf_stop_subj & !is.null(subj_idx)) {val <- pnorm(pars$traces[,col,])} else {val <- pars$traces[,col,]}
    summary[row,1]   <- round(    mean(as.vector(val)),4)
    summary[row,2]   <- round(      sd(as.vector(val)),4)
    summary[row,3:7] <- round(quantile(as.vector(val), prob = c(0.025,0.25,0.5,0.75,0.975)),4)
  }
  # Write table
  table = tableGrob(summary)
  h <- grobHeight(table)
  w <- grobWidth(table)
  title <- textGrob("Summary statistics group level parameters", y=unit(0.5,"npc") + 0.5*h, vjust=0, gp=gpar(fontsize=20))
  gt <- gTree(children=gList(table, title))
  grid.newpage()
  grid.draw(gt)
}
# ------------------------------------------------------------------------------ #



# ------------------------------------------------------------------------------ #
# Posterior means...
# ------------------------------------------------------------------------------ #
plot_individual_posterior_means = function(pars, priors){

  n_subj <- pars$n_subject
  a_big_enough_num <- 20
  post_means_ind <- array(NA,dim=c(n_subj, a_big_enough_num))

  # Get posterior means
  for (subj_num in 1:n_subj) {
    regex <- paste('*_subj.', toString(subj_num), '*', sep='')
    params = grep(glob2rx(regex), colnames(pars$traces), value=TRUE)

    n_params <- length(params)
    for(index in 1:n_params){
        post_means_ind[subj_num, index] = mean(pars$traces[,params[index],])
    }
  }
  post_means_ind <- post_means_ind[, 1:n_params]
  colnames(post_means_ind) = params

  par(mar = c(8, 5, 6, 3) + 0.1)
  layout(matrix(1:n_params, 4, byrow=T))
  par(cex.main=1.4)

  for (index in 1:n_params){
    data  <- post_means_ind[,index]
    param <- params[index]

    hist(data, main=param, xlab=param)
  }
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# Posteriors
# ------------------------------------------------------------------------------ #
plot_posteriors = function(pars, priors, plot_priors, title_addendum = ''){
  library(msm)

  # Model contains SD terms, but output on parameter distribution uncertainties is variances.
  prior_names     <- names(priors)
  prior_names_var <- gsub("_sd", "_var", prior_names)

  if (plot_priors) {
    # Set prior densities
    prior_densities <- list()
    for (name in prior_names) {
      lower_lim <- priors[[name]][1]
      upper_lim <- priors[[name]][2]

      name_var  <- gsub("_sd", "_var", name)
      prior_densities[[name_var]] <- dunif(lower_lim, lower_lim, upper_lim)
    }

    # An exception to the above uniform dists
    p_lims_prior = seq(priors$pf_stop[1],priors$pf_stop[2],0.1)
    prior_densities$pf_stop = dtnorm(p_lims_prior, priors$pf_stop[3], priors$pf_stop[4], lower=priors$pf_stop[1], upper=priors$pf_stop[2])
  }

  # Setup plotting
  par(mar = c(4.5, 5, 2, 3) + 0.1) # group
  #par(mar = c(8, 5, 6, 3) + 0.1) # individual
  layout(matrix(1:14,,2,byrow=T))
  par(cex.main=1.2)

  # Plot posteriors
  for(par in prior_names_var){

    if (par == 'pf_stop_var' | par == 'pf_stop') {next}

    title   <- paste("Posterior", par, title_addendum, sep=" ")
    xlimits <- c(min(pars$traces[,par,])-30, max(pars$traces[,par,])+30)

    posterior <- density(pars$traces[,par,])
    plot(posterior, xlim = xlimits, main = title, xlab=paste(par,"(ms)",sep=" "), ylab = "Density")

    # We plot param priors for (individual & not hierarchical) or group get plotted
    if (plot_priors) {
      # These are mostly just uniform lines
      curves_to_plot <- c(prior_densities[prior_names_var == par], prior_densities[prior_names_var == par])
      lines(c(0, xlimits[2]), curves_to_plot , lty=2)
    }
  }

  #plot(density(pars$traces[,"pf_stop",]),xlim = c(priors$pf_stop[1], priors$pf_stop[2]), main = "Posterior pf_stop",xlab="pf_stop on probit scale",ylab = "Density")
  #plot(density(pars$traces[,"pf_stop_sd",]),xlim = c(priors$pf_stop_sd[1],priors$pf_stop_sd[2]), main = "Posterior pf_stop_sd",xlab="pf_stop_sd on probit scale",ylab = "Density")

  if (plot_priors) {
    # These are mostly just uniform lines
    #lines(p_lims_prior,as.numeric(unlist(prior_den_group[names(prior_den_group)=="pf_stop"])),lty=2)
    #lines(c(priors$pf_stop_sd[1],priors$pf_stop_sd[2]),c(prior_den_group[names(prior_den_group)=="pf_stop_sd"],prior_den_group[names(prior_den_group)=="pf_stop_sd"]),lty=2)
  }
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# Posteriors
# ------------------------------------------------------------------------------ #
plot_posteriors_wrapper = function(pars, priors){

  if (pars$n_subject==1) {
    # Individual posteriors
    plot_posteriors(pars, priors, plot_priors = TRUE)
    #plot_chains(pars, priors, probit = FALSE)
  }
  else {
    # Group posteriors
    flog.info('Plotting group posteriors')
    plot_posteriors(pars, priors, plot_priors = TRUE, title_addendum = 'for Group')
    plot_individual_posterior_means(pars,priors)

    #plot_chains(pars, priors, probit = FALSE)

    for (n in 1:pars$n_subject){
      # Individual posteriors without prior plotting

      flog.info('Plotting group posteriors for subject %s', n)
      param_names = paste(names(priors),n,sep="")
      plot_posteriors(pars, priors, plot_priors = FALSE, title_addendum = c('for Subject ', toString(n)))

      flog.info('Plotting chains for subject %s', n)
      #plot_chains(pars, priors, subject_idx = n, probit = TRUE)
    }
  }
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# MCMC Chains
# ------------------------------------------------------------------------------ #
plot_chains = function(pars, priors, subject_idx = NULL, probit){

  # Extract parameter names
  params   <- colnames(pars$traces)
  n_params <- length(params)

  # Plot setup for group data
  if (is.null(subject_idx)) {
    par(mar = c(4.5, 5, 2, 3) + 0.1)
  } else {
    par(mar = c(8, 5, 6, 3) + 0.1)
  }

  layout(matrix(1:n_params, 4, byrow=T))
  par(cex.main=1.4)

  pf_stop = c(priors$pf_stop[1],priors$pf_stop[2])
  if (probit) {
    pf_stop_subj = grep(glob2rx("pf_stop*"), params, value=TRUE)
  }

  for(par in params){
    if (par == 'pf_stop' | par == 'pf_stop_sd') {next}
    lim = c(min(pars$traces[,par,])-10, max(pars$traces[,par,])+10)

    plot(1:pars$n_samples, pars$traces[,par,1],ylim = lim,
      xlim = c(1,pars$n_samples), main = paste("MCMC chains",par,sep=" "), xlab="Iteration",
      ylab = par,type="l")

    if(pars$n_chains>1){
       #if multiple chains, draw lines
      for(j in 2:pars$n_chains){
        lines(1:pars$n_samples,pars$traces[,par,j],col=j)
      }
    }
  }

  if (!probit) {
    plot(1:pars$n_samples,pars$traces[,"pf_stop",1],ylim = c(pf_stop[1],pf_stop[2]), xlim = c(1,pars$n_samples), main = "MCMC chains pf_stop",xlab="Iteration",ylab = "pf_stop",type="l")
  } else {
    plot(1:pars$n_samples,pnorm(pars$traces[,pf_stop_subj,1]),ylim = c(0,1), xlim = c(1,pars$n_samples), main = paste("MCMC chains",pf_stop_subj,sep=" "), xlab="Iteration",ylab = pf_stop_subj,type="l")
  }

  if(pars$n_chains>1){  #if multiple chains, draw lines
      for(j in 2:pars$n_chains){
        if (!probit) {
          lines(1:pars$n_samples,pars$traces[,"pf_stop",j],col=j)
        } else {
          lines(1:pars$n_samples,pnorm(pars$traces[,pf_stop_subj,j]),col=j)
        }
      }
  }
}
# ------------------------------------------------------------------------------ #
