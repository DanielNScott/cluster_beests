###
### This file contains most of the BEESTS post processing routines.
###
# Original Authors: Matzke et al
# Substantial Update: Daniel Scott
#

# A function for checking silly input values
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


### Load python traces in R and make array for further use
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

# Computes summary stats for each parameter (by participant, collapsed over chains) and save output to csv.
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

   	pf_stop_subj       = grep(glob2rx("pf_stop*"), params, value=TRUE)

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

      # Get across subject parameter means, standard deviations, and quantiles:
      for(par in params){
         browser()
         if (par == pf_stop_subj) {val <- pnorm(pars$traces[,pf_stop_subj,])} else {val <- pars$traces[,par,]}

         # Summary statistics
         summary[par,1  ] = round(    mean(as.vector(val)), 4)
         summary[par,2  ] = round(      sd(as.vector(val)), 4)
         summary[par,3:7] = round(quantile(as.vector(val), prob = c(0.025,0.25,0.5,0.75,0.975)), 4)
      }

     # Write table
      table = tableGrob(summary)
      h <- grobHeight(table)
      w <- grobWidth(table)
      title <- textGrob(paste("Summary statistics Subject",subject_idx,sep=" "), y=unit(0.5,"npc") + 0.5*h,
                     vjust=0, gp=gpar(fontsize=20))
      gt <- gTree(children=gList(table, title))
      grid.newpage()
      grid.draw(gt)
   }
   else {
      # Group data:
      for(row in 1:length(params)){
         col <- params[[row]]
         summary[row,1]   <- round(mean(as.vector(pars$traces[,col,])),4)
         summary[row,2]   <- round(sd(as.vector(pars$traces[,col,])),4)
         summary[row,3:7] <- round(quantile(as.vector(pars$traces[,col,]), prob = c(0.025,0.25,0.5,0.75,0.975)),4)
      }
      browser()
      # Write table
      table = tableGrob(summary)
      h <- grobHeight(table)
      w <- grobWidth(table)
      title <- textGrob("Summary statistics group level parameters", y=unit(0.5,"npc") + 0.5*h,
                     vjust=0, gp=gpar(fontsize=20))
      gt <- gTree(children=gList(table, title))
      grid.newpage()
      grid.draw(gt)
   }
}


#--------------------------------------------------------------------------------------------------------------------------------------------------Plots posteriors for individual participants for the individual model
plot_individual_posteriors = function(pars,params,priors){

  prior_den = list(mu_go = dunif(priors$mu_go[1],priors$mu_go[1],priors$mu_go[2]),
                    mu_stop = dunif(priors$mu_stop[1],priors$mu_stop[1],priors$mu_stop[2]),
                    sigma_go = dunif(priors$sigma_go[1],priors$sigma_go[1],priors$sigma_go[2]),
                    sigma_stop = dunif(priors$sigma_stop[1],priors$sigma_stop[1],priors$sigma_stop[2]),
                    tau_go = dunif(priors$tau_go[1],priors$tau_go[1],priors$tau_go[2]),
                    tau_stop = dunif(priors$tau_stop[1],priors$tau_stop[1],priors$tau_stop[2]),
                    pf_stop = dunif(priors$pf_stop[1],priors$pf_stop[1],priors$pf_stop[2]))

	par(mar = c(8, 5, 6, 3) + 0.1)
	layout(matrix(1:8,4,byrow=T))
	par(cex.main=1.4)

	for(par in params[1:6]){
    lim = c(min(pars$traces[,par,])-30,max(pars$traces[,par,])+30)
    plot(density(pars$traces[,par,]),xlim = lim, main = paste("Posterior",par,sep=" "),xlab=paste(par,"(ms)",sep=" "),ylab = "Density")
    lines(c(0,lim[2]),c(prior_den[names(prior_den)==par],prior_den[names(prior_den)==par]),lty=2)
  }

  plot(density(pars$traces[,"pf_stop",]),xlim = c(priors$pf_stop[1],priors$pf_stop[2]), main = "Posterior pf_stop",xlab="Probability",ylab = "Density")
  lines(c(priors$pf_stop[1],priors$pf_stop[2]),c(prior_den[names(prior_den)=="pf_stop"],prior_den[names(prior_den)=="pf_stop"]),lty=2)
 }

#--------------------------------------------------------------------------------------------------------------------------------------------------Plots posteriors for individual participants for the hierarchical model
plot_individual_posteriors_for_hier = function(pars,params,subject_idx){

	par(mar = c(8, 5, 6, 3) + 0.1)
	layout(matrix(1:8,4,byrow=T))
	par(cex.main=1.4)

  pf_stop_subj = grep(glob2rx("pf_stop*"), params, value=TRUE)

	for(par in params[1:6]){
    lim = c(min(pars$traces[,par,])-30,max(pars$traces[,par,])+30)
    plot(density(pars$traces[,par,]),xlim = lim, main = paste("Posterior",par,sep=" "),xlab=paste(par,"(ms)",sep=" "),ylab = "Density")
  }
  plot(density(pnorm(pars$traces[,pf_stop_subj,])),xlim = c(0,1), main = paste("Posterior", pf_stop_subj,sep=" "),xlab="Probability",ylab = "Density")
}

#--------------------------------------------------------------------------------------------------------------------------------------------------Hist of posterior means for individual participants for the hierarchical model
plot_individual_posterior_means = function(pars,priors){
  nsub = pars$n_subject
  post_means_ind = array(NA,dim=c(nsub, 7))
  nam = c("mu_go","mu_stop","sigma_go","sigma_stop","tau_go","tau_stop","pf_stop")
  colnames(post_means_ind) = nam

  # Get posterior means
  for (n in 1:nsub){

    params = c(paste("mu_go_subj.",n,sep=""),paste("mu_stop_subj.",n,sep=""),
                        paste("sigma_go_subj.",n,sep=""),paste("sigma_stop_subj.",n,sep=""),
                        paste("tau_go_subj.",n,sep=""),paste("tau_stop_subj.",n,sep=""),
                        paste("pf_stop_subjpt.",n,sep=""))

    for(p in 1:7){
        post_means_ind[n,p] = mean(pars$traces[,params[p],])
    }
  }

  par(mar = c(8, 5, 6, 3) + 0.1)
  layout(matrix(1:8,4,byrow=T))
	par(cex.main=1.4)
  #p=1
	#print(priors[names(priors)==nam[p]])

	for (p in 1:7){
    #print(as.numeric(unlist(priors[names(priors)==nam[p]])))
    hist(post_means_ind[,p],main=nam[p],xlab=nam[p],xlim=c(as.numeric(unlist(priors[names(priors)==nam[p]]))[1],as.numeric(unlist(priors[names(priors)==nam[p]]))[2]))
  }
}
#---------------------------------------------------------------------------------------------------------------------------------------------------Plots posterior of group level parameters for hierarchical model
plot_group_posteriors = function(pars,priors){

library(msm)
p_lims_prior = seq(priors$pf_stop[1],priors$pf_stop[2],0.1)

  prior_den_group = list(mu_go = dunif(priors$mu_go[1],priors$mu_go[1],priors$mu_go[2]),
                    mu_stop = dunif(priors$mu_stop[1],priors$mu_stop[1],priors$mu_stop[2]),
                    sigma_go = dunif(priors$sigma_go[1],priors$sigma_go[1],priors$sigma_go[2]),
                    sigma_stop = dunif(priors$sigma_stop[1],priors$sigma_stop[1],priors$sigma_stop[2]),
                    tau_go = dunif(priors$tau_go[1],priors$tau_go[1],priors$tau_go[2]),
                    tau_stop = dunif(priors$tau_stop[1],priors$tau_stop[1],priors$tau_stop[2]),
                    pf_stop = dtnorm(p_lims_prior,priors$pf_stop[3],priors$pf_stop[4],lower=priors$pf_stop[1],upper=priors$pf_stop[2]),

                    mu_go_var = dunif(priors$mu_go_var[1],priors$mu_go_var[1],priors$mu_go_var[2]),
                    mu_stop_var = dunif(priors$mu_stop_var[1],priors$mu_stop_var[1],priors$mu_stop_var[2]),
                    sigma_go_var = dunif(priors$sigma_go_var[1],priors$sigma_go_var[1],priors$sigma_go_var[2]),
                    sigma_stop_var = dunif(priors$sigma_stop_var[1],priors$sigma_stop_var[1],priors$sigma_stop_var[2]),
                    tau_go_var = dunif(priors$tau_go_var[1],priors$tau_go_var[1],priors$tau_go_var[2]),
                    tau_stop_var = dunif(priors$tau_stop_var[1],priors$tau_stop_var[1],priors$tau_stop_var[2]),
                    pf_stop_var = dunif(priors$pf_stop_var[1],priors$pf_stop_var[1],priors$pf_stop_var[2]))

	par(mar = c(4.5, 5, 2, 3) + 0.1)
	layout(matrix(1:14,,2,byrow=T))
	par(cex.main=1.2)

	group_pars = c("mu_go","mu_go_var","mu_stop","mu_stop_var","sigma_go","sigma_go_var","sigma_stop", "sigma_stop_var","tau_go","tau_go_var","tau_stop","tau_stop_var")

	for(par in group_pars){
    		lim = c(min(pars$traces[,par,])-30,max(pars$traces[,par,])+30)
    		plot(density(pars$traces[,par,]),xlim = lim, main = paste("Posterior",par,sep=" "),xlab=paste(par,"(ms)",sep=" "),ylab = "Density")
    		lines(c(0,lim[2]),c(prior_den_group[names(prior_den_group)==par],prior_den_group[names(prior_den_group)==par]),lty=2)
  	}
  plot(density(pars$traces[,"pf_stop",]),xlim = c(priors$pf_stop[1],priors$pf_stop[2]), main = "Posterior pf_stop",xlab="pf_stop on probit scale",ylab = "Density")
  lines(p_lims_prior,as.numeric(unlist(prior_den_group[names(prior_den_group)=="pf_stop"])),lty=2)

  plot(density(pars$traces[,"pf_stop_var",]),xlim = c(priors$pf_stop_var[1],priors$pf_stop_var[2]), main = "Posterior pf_stop_var",xlab="pf_stop_var on probit scale",ylab = "Density")
  lines(c(priors$pf_stop_var[1],priors$pf_stop_var[2]),c(prior_den_group[names(prior_den_group)=="pf_stop_var"],prior_den_group[names(prior_den_group)=="pf_stop_var"]),lty=2)
}

#---------------------------------------------------------------------------------------------------------------------------------------------------Plots posteriors using plot_individual_posteriors() and plot_group_posteriors() functions
plot_posteriors = function(pars,all_pars=T,priors){
	if (all_pars){
		if (pars$n_subject==1) ###Individual model
			plot_individual_posteriors(pars, params = c("mu_go","mu_stop","sigma_go","sigma_stop","tau_go","tau_stop", "pf_stop"), priors)

		else {
			plot_group_posteriors(pars,priors)
      plot_individual_posterior_means(pars,priors)

			for (n in 1:pars$n_subject){ #Seperate figure for each participant
				param_name = c(paste("mu_go_subj.",n,sep=""),paste("mu_stop_subj.",n,sep=""),
                        paste("sigma_go_subj.",n,sep=""),paste("sigma_stop_subj.",n,sep=""),
                        paste("tau_go_subj.",n,sep=""),paste("tau_stop_subj.",n,sep=""),
                        paste("pf_stop_subjpt.",n,sep=""))

				plot_individual_posteriors_for_hier(pars, params = param_name, subject_idx = n)
			}
		}
	} else {
         plot_group_posteriors(pars,priors)
         plot_individual_posterior_means(pars,priors)
    }
}

#---------------------------------------------------------------------------------------------------------------------------------------------------Plots chains for individual participants for the individual model as well as the hierarchical model
plot_individual_chains = function(pars,params,priors,subject_idx = NULL){

  pf_stop = c(priors$pf_stop[1],priors$pf_stop[2])

  par(mar = c(8, 5, 6, 3) + 0.1)
	layout(matrix(1:8,4,byrow=T))
	par(cex.main=1.4)

	for(par in params[1:6]){
    lim = c(min(pars$traces[,par,])-10,max(pars$traces[,par,])+10)
    plot(1:pars$n_samples,pars$traces[,par,1],ylim = lim, xlim = c(1,pars$n_samples), main = paste("MCMC chains",par,sep=" "),xlab="Iteration",ylab = par,type="l")
		if(pars$n_chains>1){  #if multiple chains, draw lines
			for(j in 2:pars$n_chains){
				lines(1:pars$n_samples,pars$traces[,par,j],col=j)
			}
		}
	}

	plot(1:pars$n_samples,pars$traces[,"pf_stop",1],ylim = c(pf_stop[1],pf_stop[2]), xlim = c(1,pars$n_samples), main = "MCMC chains pf_stop",xlab="Iteration",ylab = "pf_stop",type="l")

	if(pars$n_chains>1){  #if multiple chains, draw lines
			for(j in 2:pars$n_chains){
				lines(1:pars$n_samples,pars$traces[,"pf_stop",j],col=j)
	  	}
	}
}

#---------------------------------------------------------------------------------------------------------------------------------------------------Plots chains for individual participants for the individual model as well as the hierarchical model
plot_individual_chains_probit = function(pars,params,subject_idx = NULL){

	par(mar = c(8, 5, 6, 3) + 0.1)
	layout(matrix(1:8,4,byrow=T))
	par(cex.main=1.4)

  pf_stop_subj = grep(glob2rx("pf_stop*"), params, value=TRUE)

	for(par in params[1:6]){
    lim = c(min(pars$traces[,par,])-10,max(pars$traces[,par,])+10)
    plot(1:pars$n_samples,pars$traces[,par,1],ylim = lim, xlim = c(1,pars$n_samples), main = paste("MCMC chains",par),xlab="Iteration",ylab = par,type="l")

		if(pars$n_chains>1){  #if multiple chains, draw lines
			for(j in 2:pars$n_chains){
				lines(1:pars$n_samples,pars$traces[,par,j],col=j)
			}
		}
	}

	plot(1:pars$n_samples,pnorm(pars$traces[,pf_stop_subj,1]),ylim = c(0,1), xlim = c(1,pars$n_samples), main = paste("MCMC chains",pf_stop_subj,sep=" "), xlab="Iteration",ylab = pf_stop_subj,type="l")

	if(pars$n_chains>1){  #if multiple chains, draw lines
			for(j in 2:pars$n_chains){
				lines(1:pars$n_samples,pnorm(pars$traces[,pf_stop_subj,j]),col=j)
			}
	}
}

#-------------------------------------------------------------------------------------------------------------------------------------------------Plots chains of group level parameters for hierarchical model
plot_group_chains = function(pars,priors){

	par(mar = c(4.5, 5, 2, 3) + 0.1)
	layout(matrix(1:14,,2,byrow=T))
	par(cex.main=1.4)

	group_pars = c("mu_go","mu_go_var","mu_stop","mu_stop_var",
                "sigma_go","sigma_go_var", "sigma_stop","sigma_stop_var",
                 "tau_go","tau_go_var","tau_stop","tau_stop_var")

	for(par in group_pars){
    		lim = c(min(pars$traces[,par,])-10,max(pars$traces[,par,])+10)
    		plot(1:pars$n_samples,pars$traces[,par,1],ylim = lim, xlim = c(1,pars$n_samples), main = paste("MCMC chains",par,sep=" "),xlab="Iteration",ylab = par,type="l")

		if(pars$n_chains>1){ #if multiple chains, draw lines
			for(j in 2:pars$n_chains){
				lines(1:pars$n_samples,pars$traces[,par,j],col=j)
			}
		}
  }

  plot(1:pars$n_samples,pars$traces[,"pf_stop",1],ylim = c(priors$pf_stop[1],priors$pf_stop[2]), xlim = c(1,pars$n_samples), main = "MCMC chains pf_stop",xlab="Iteration",ylab = "pf_stop",type="l")

	if(pars$n_chains>1){ #if multiple chains, draw lines
			for(j in 2:pars$n_chains){
				lines(1:pars$n_samples,pars$traces[,"pf_stop",j],col=j)
			}
	}

  plot(1:pars$n_samples,pars$traces[,"pf_stop_var",1],ylim = c(priors$pf_stop_var[1],priors$pf_stop_var[2]), xlim = c(1,pars$n_samples), main = "MCMC chains pf_stop_var",xlab="Iteration",ylab = "pf_stop_var",type="l")

	if(pars$n_chains>1){ #if multiple chains, draw lines
			for(j in 2:pars$n_chains){
				lines(1:pars$n_samples,pars$traces[,"pf_stop_var",j],col=j)
			}
	}
}

#------------------------------------------------------------------------------------------------------------------------------------------------------Plots chains using plot_individual_chains() and plot_group_chains() functions
plot_chains = function(pars,all_pars=T,priors){
	if (all_pars){
		if (pars$n_subject==1) #Individual model
			plot_individual_chains(pars, params = c("mu_go","mu_stop","sigma_go","sigma_stop","tau_go","tau_stop","pf_stop"),priors)
		else {                 #Group model
			plot_group_chains(pars,priors)

			for (n in 1:pars$n_subject){ #Seperate figure for each participant
				param_name = c(paste("mu_go_subj.",n,sep=""),paste("mu_stop_subj.",n,sep=""),
                        paste("sigma_go_subj.",n,sep=""),paste("sigma_stop_subj.",n,sep=""),
                        paste("tau_go_subj.",n,sep=""),paste("tau_stop_subj.",n,sep=""),
                        paste("pf_stop_subjpt.",n,sep=""))

				plot_individual_chains_probit(pars, params = param_name, subject_idx = n)
			}
		}
	} else plot_group_chains(pars,priors)
}

#-------------------------------------------------------------------------------------------------------------------------------------------------------Posterior predictions
# Sample from the joint distribution of the parameters
sample_joint_posterior = function(pars,n_post_samples,subject_idx = NULL){

	#n_post_samples must to be smaller than pars$n_samples*pars$n_chains!!!
	if(pars$n_subject==1){
		params = c("mu_go","mu_stop","sigma_go","sigma_stop","tau_go","tau_stop","pf_stop")
	} else {
		params = c(paste("mu_go_subj.",subject_idx,sep=""),paste("mu_stop_subj.",subject_idx,sep=""),
				paste("sigma_go_subj.",subject_idx,sep=""),paste("sigma_stop_subj.",subject_idx,sep=""),
				paste("tau_go_subj.",subject_idx,sep=""),paste("tau_stop_subj.",subject_idx,sep=""),paste("pf_stop_subjpt.",subject_idx,sep=""))
	}

	par_vectors = array(NA,dim = c(n_post_samples,7),list(NULL,params))
	it_ind = sample(1:(pars$n_samples*pars$n_chains),n_post_samples,replace=F)

	for(p in params){
		par_vectors[,p] = as.vector(pars$traces[,p,])[it_ind]
	}

	return(list(par_vectors = par_vectors, params = params))
}

#------------------------------------------------------------------------------------------------------Go RTs
load_prep_observed_go= function(pars,subject_idx,data){

  n_subj = pars$n_subject
	if(n_subj == 1){
		data = cbind(subj_idx = rep(1,nrow(data)),data)
	}

  ###select go RTs
  go_rt = data$rt[data$subj_idx==subject_idx&data$ss_presented==0&data$inhibited==-999&data$ssd==-999&data$rt!=-999]
  n_go = length(go_rt)

  return(list(go_rt=go_rt,n_go=n_go,n_subj=n_subj))
}

#-------------------------------------------------------------------------------------------------------Median SRRTs
load_prep_observed_srrt= function(pars,subject_idx,data){

  n_subj = pars$n_subject
	if(n_subj == 1){
		data = cbind(subj_idx = rep(1,nrow(data)),data)
	}

  # determine number of srrts and inhibitions for EACH delay
  delays_all = sort(unique(data$ssd[data$subj_idx==subject_idx&data$ss_presented==1]))
  n_delays_all = length(delays_all)
  n_observed_srrt_all = rep(0,n_delays_all)
  n_observed_inhibit_all = rep(0,n_delays_all)

  for(d2 in 1:n_delays_all){
    n_observed_srrt_all[d2] = nrow(data[data$subj_idx==subject_idx&data$ss_presented==1&data$inhibited==0&data$ssd==delays_all[d2],])
    n_observed_inhibit_all[d2] = nrow(data[data$subj_idx==subject_idx&data$ss_presented==1&data$inhibited==1&data$ssd==delays_all[d2],])
  }
  n_observed_all = n_observed_srrt_all + n_observed_inhibit_all

	# select delays with at least one srrt for posterior predictive checks for median srrt and srrt distributions
	# and get observed SRRTs per delay; and compute median
	delays = sort(unique(data$ssd[data$subj_idx==subject_idx&data$ss_presented==1&data$inhibited==0&data$rt!=-999]))
	n_delays = length(delays)

	median_observed_srrt = rep(NA,length(delays))
	n_observed_srrt = rep(NA,length(delays))
	n_observed_inhibit = rep(NA,length(delays))
  srrt_obs = list()

	for (d in 1:n_delays){
		srrt_temp = data$rt[data$subj_idx==subject_idx&data$ss_presented==1&data$inhibited==0&data$ssd==delays[d]]
		srrt_obs[[d]] = srrt_temp
		median_observed_srrt[d] = median(srrt_temp)
		n_observed_srrt[d] = length(srrt_temp)
		n_observed_inhibit[d] = nrow(data[data$subj_idx==subject_idx&data$ss_presented==1&data$inhibited==1&data$rt==-999&data$ssd==delays[d],])
	}
	n_observed = n_observed_srrt + n_observed_inhibit

  # response rate per delay
  RR_obs = n_observed_srrt_all/n_observed_all

	return(list(delays = delays, n_delays = n_delays, delays_all = delays_all, n_delays_all = n_delays_all,median_observed_srrt = median_observed_srrt, n_observed_srrt = n_observed_srrt, n_observed = n_observed,
      n_observed_srrt_all = n_observed_srrt_all, n_observed_all = n_observed_all,n_subj = n_subj, RR_obs=RR_obs,srrt_obs=srrt_obs))
}

#---------------------------Generate predicted go RTs
generate_posterior_predictions_go = function(n_go,n_post_samples,par_vectors,params){

	go_pred = matrix(NA,n_go,n_post_samples)
	for(j in 1:n_post_samples){
    go_pred[,j] = rexGAUS(n_go,mu = par_vectors[j,params[1]],sigma = par_vectors[j,params[3]],nu = par_vectors[j,params[5]])
	}

	return(list(go_pred=go_pred))
}

#--------------------------Plot go RTs and posterior predictions
plot_posterior_predictions_go = function(go_rt,n_go,go_pred,n_subj,subject_idx=NULL,n_post_samples){

	par(mar = c(13, 11, 11, 10) + 0.1,cex.main=1.5)
  x = max(go_rt)+600
  y = max(hist(go_rt,plot=F)$density)
	y = y+0.6*y
  if(n_subj==1){
    my_main = "Posterior predictive model check for go RT distribution"
  } else {
    my_main = paste("Posterior predictive model check for go RT distribution \nSubject",subject_idx,sep=" ")
  }

  layout(matrix(1))
  hist(go_rt,freq=F,xlim=c(0,x),ylim=c(0,y),xlab="go RT (ms)",las=1, main = my_main,ylab="Density")

  for(j in 1:n_post_samples){
		lines(density(go_pred[,j]),col="gray")
	}
  hist(go_rt,freq=F,add=T)
  #grid.newpage()
}

#--------------------------Generates posterior predictions for median SRRTs, SRRT distributions, and inhibition functions
generate_posterior_predictions_srrt = function(delays_all,n_delays_all,n_observed_all,n_observed_srrt_all,par_vectors,n_post_samples,params,n_subj,subject_idx=NULL){

	median_post_pred_srrt = matrix(NA,n_delays_all,n_post_samples)
	RR_pred = matrix(NA,n_delays_all,n_post_samples)
	#post_pred_srrt = sapply(as.character(delays_all),function(x) NULL)

	if(n_subj==1){
		prob = par_vectors[,"pf_stop"]
	} else {
    prob = pnorm(par_vectors[,paste("pf_stop_subjpt.",subject_idx,sep="")])
  }
  #print(prob)
  #print(n_delays_all)
  #print(n_observed_all)

	for(d in 1:n_delays_all){
    n_tf = n_observed_all[d]*prob
    n_tf = round(n_tf,0)
    n_notf = n_observed_all[d]-n_tf
    #print(n_tf)
    #print(n_notf)
    #print(n_observed_all)

    for(j in 1:n_post_samples){
      if(n_tf[j]>0){
        ###not triggered trials
        srrt_tf = rexGAUS(n_tf[j],mu = par_vectors[j,params[1]],sigma = par_vectors[j,params[3]],nu = par_vectors[j,params[5]])
      } else srrt_tf = NULL

      #print(n_notf[j]==0)
      if(n_notf[j]>0){
        ###triggered trials
        go_rts = rexGAUS(n_notf[j],mu = par_vectors[j,params[1]],sigma = par_vectors[j,params[3]],nu = par_vectors[j,params[5]])
        stop_rt = delays_all[d] + rexGAUS(n_notf[j],mu=par_vectors[j,params[2]],sigma = par_vectors[j,params[4]],nu = par_vectors[j,params[6]])
        srrt = go_rts[go_rts<=stop_rt]
       } else srrt = NULL

      all_pred_srrt = c(srrt,srrt_tf)
      post_pred_srrt[[d]][[j]] = list(all_pred_srrt)
      median_post_pred_srrt[d,j] = median(all_pred_srrt)
      RR_pred[d,j] = (length(all_pred_srrt))/n_observed_all[d]
    }
   }
   median_post_pred_srrt = median_post_pred_srrt[n_observed_srrt_all>0,]
   post_pred_srrt = post_pred_srrt[n_observed_srrt_all>0]

   return(list(median_post_pred_srrt = median_post_pred_srrt, RR_pred=RR_pred, post_pred_srrt = post_pred_srrt))
}

####--------------------------Generates output for posterior predictions for SRRTs and corresponding p values
post_pred_output_median_srrt = function(n_delays,delays,n_subj,median_post_pred_srrt,median_observed_srrt,n_observed_srrt,subject_idx=NULL){

  pvalue_one = rep(NA,n_delays)
	pvalue_two = rep(NA,n_delays)

  # Compute posterior predictive p values/delay
	for(d in 1:n_delays){
		pvalue_one[d] = mean(median_post_pred_srrt[d,]> median_observed_srrt[d],na.rm=T)
	}
	ipvalue_one = 1-pvalue_one

	for(d in 1:n_delays){
			pvalue_two[d] = 2*min(pvalue_one[d],ipvalue_one[d])
	}
	average_post_pred = apply(median_post_pred_srrt,1,mean,na.rm=T)

	# Make output table
	post_pred_summary = rbind(round(rbind(n_observed_srrt,median_observed_srrt,average_post_pred),2),
                            round(rbind(pvalue_one,pvalue_two),3))
	colnames(post_pred_summary) = paste(rep("SSD=",each=n_delays),delays,sep="")
	rownames(post_pred_summary) = c("Number of observed SRRT", "Observed median SRRT","Average predicted SRRT","One-sided p value","Two-sided p value")

  # Write table
  nc = ncol(post_pred_summary)
  fsize=19-nc
  if (fsize<1)  {
      fsize=1
  }
  table = tableGrob(post_pred_summary)#, gpar.coretext = gpar(fontsize=fsize),equal.width=F,gpar.coltext = gpar(fontsize=fsize),gpar.rowtext = gpar(fontsize=fsize))
	h <- grobHeight(table)
  w <- grobWidth(table)

  if(n_subj==1){
    title <- textGrob("Posterior predictive p values for median SRRT", y=unit(0.5,"npc") + 0.5*h,
                  vjust=0, gp=gpar(fontsize=20))
    my_main = "Posterior predictive model check for median SRRT"
  } else {
    title <- textGrob(paste("Posterior predictive p values for median SRRT \nSubject",subject_idx,sep=" "), y=unit(0.5,"npc") + 0.5*h,
                  vjust=0, gp=gpar(fontsize=20))
    my_main = paste("Posterior predictive model check for median SRRT \nSubject",subject_idx,sep=" ")
  }

  gt <- gTree(children=gList(table, title))
  grid.newpage()
  grid.draw(gt)

  # Make violin plot
  my_min = min(c(min(median_post_pred_srrt,na.rm=T),min(median_observed_srrt)))-80
	my_max = max(c(max(median_post_pred_srrt,na.rm=T),max(median_observed_srrt)))+80
	ran = round(c(my_min,my_max))
	par(cex=1.2,cex.lab = 1.2,mar=c(11, 10, 10, 8),cex.main=1.4)
	plot(1:length(delays),median_observed_srrt,pch=17,xlab="SSD (ms)",axes=F,
		ylim=ran,xlim=c(0,n_delays+1),ylab ="",main=my_main,cex=1.3)
	axis(1,at=0:(n_delays+1),labels=c(NA,delays,NA))
	axis(2,at=ran,las=2)
	mtext("Median SRRT (ms)",2, cex=1.4,line=1)
	for(i in 1:n_delays){
    preds_med = median_post_pred_srrt[i,][!is.na(median_post_pred_srrt[i,])]
    if (length(preds_med) >0) {
      vioplot(preds_med,add=T,at=i,col="gray",cex=2)
    }
	}
	points(1:n_delays,median_observed_srrt,cex=1.5,pch=17)
	lines(1:n_delays,median_observed_srrt,lwd=2,lty=2)
	lines(1:n_delays,average_post_pred,lwd=2,lty=1)
}

#--------------------------Generates output for posterior predictions for inhibition function and corresponding p values
post_pred_output_inhibition_function = function(params,pars,n_subj,n_delays_all,delays_all,RR_pred,RR_obs,n_observed_all,subject_idx=NULL){

  pvalue_one_if = rep(NA,n_delays_all)
	pvalue_two_if = rep(NA,n_delays_all)

	if(n_subj==1){
		prob = pars$traces[,params[7],]
	} else {
    prob = pnorm(pars$traces[,grep(glob2rx("pf_stop*"), params, value=TRUE),])
  }
  q_tf = quantile(prob,prob=c(0.025,0.975))
  mean_prob = mean(prob)

  # Compute posterior predictive p values/delay
	for(d in 1:n_delays_all){
    if (RR_obs[d]==1)
      pvalue_one_if[d] = mean(RR_pred[d,]>=RR_obs[d],na.rm=T)
    else
      pvalue_one_if[d] = mean(RR_pred[d,]> RR_obs[d],na.rm=T)
	}
	ipvalue_one_if = 1-pvalue_one_if

	for(d in 1:n_delays_all){
			pvalue_two_if[d] = 2*min(pvalue_one_if[d],ipvalue_one_if[d])
	}
	average_pred_if = apply(RR_pred,1,mean,na.rm=T)

	# Make output table
	post_pred_if_summary = rbind(round(rbind(n_observed_all,RR_obs,average_pred_if),2),
                            round(rbind(pvalue_one_if,pvalue_two_if),3))
	colnames(post_pred_if_summary) = paste(rep("SSD=",each=n_delays_all),delays_all,sep="")
	rownames(post_pred_if_summary) = c("Number of stop-signal trials", "Observed response rate","Average predicted response rate","One-sided p value","Two-sided p value")

  # Write table
  nc = ncol(post_pred_if_summary)
  fsize=19-nc
  if (fsize<1)  {
      fsize=1
  }
  table = tableGrob(post_pred_if_summary) #,name="IF", gpar.coretext = gpar(fontsize=fsize),equal.width=F,gpar.coltext = gpar(fontsize=fsize),gpar.rowtext = gpar(fontsize=fsize))
 	h <- grobHeight(table)
  w <- grobWidth(table)

  if(n_subj==1){
		title <- textGrob("Posterior predictive p values for inhibition function", y=unit(0.5,"npc") + 0.5*h,
                  vjust=0, gp=gpar(fontsize=16))
    my_main = "Posterior predictive model check for inhibition function"
	} else {
    title <- textGrob(paste("Posterior predictive p values for inhibition function \nSubject",subject_idx,sep=" "), y=unit(0.5,"npc") + 0.5*h,
                  vjust=0, gp=gpar(fontsize=16))
    my_main = paste("Posterior predictive model check for inhibition function \nSubject",subject_idx,sep=" ")
  }

  gt <- gTree(children=gList(table, title))
  grid.newpage()
  grid.draw(gt)

  # Make violin plot
	par(cex=1.2,cex.lab = 1.2,mar=c(11, 10, 10, 8),cex.main=1.4)
	plot(1:length(delays_all),RR_obs,pch=17,xlab="SSD (ms)",axes=F,
		ylim=c(0,1),xlim=c(0,n_delays_all+1),ylab ="",main=my_main,cex=1.3)
	axis(1,at=0:(n_delays_all+1),labels=c(NA,delays_all,NA))
	axis(2,at=c(0,1),las=2)
	mtext("Response Rate",2, cex=1.4,line=1)

	for(i in 1:n_delays_all){
    if(length(unique(RR_pred[i,]))!=1){
      vioplot(RR_pred[i,],add=T,at=i,col="gray",cex=2)
    }
  }

	for(i in 1:n_delays_all){
    vioplot(RR_pred[i,],add=T,at=i,col="gray",cex=2)
	}
	points(1:n_delays_all,RR_obs,cex=1.5,pch=17)
	lines(1:n_delays_all,RR_obs,lwd=2,lty=2)
	lines(1:n_delays_all,average_pred_if,lwd=2,lty=1)
	# add TF parameter for reference
  lines(c(1,n_delays_all),rep(q_tf[1],2),lwd=2,lty=4)
  lines(c(1,n_delays_all),rep(q_tf[2],2),lwd=2,lty=4)
  lines(c(1,n_delays_all),rep(mean_prob,2),lwd=2,lty=3)
}

#--------------------------Generates output for posterior predictions for srrt distributions
# only for SSDs with at least 10 observed SSRTs
post_pred_output_srrt_distribution = function(delays,n_observed_srrt,n_subj,post_pred_srrt,srrt_obs,n_post_samples,subject_idx=NULL){

  ind = n_observed_srrt>=10
  dels = delays[ind]
  n_dels = length(dels)
  srrt_obs_use = srrt_obs[ind]
  post_pred_srrt_use =  post_pred_srrt[ind]

  if(n_dels>0){
    for(d in 1:n_dels){
      SRRT_obs = unlist(srrt_obs_use[d])
      ssd_obs = dels[d]
      if(n_subj == 1){
        my_main = paste("Posterior predictive model check for SRRT distribution \nat SSD =",ssd_obs,sep=" ")
      } else {
        my_main = paste(paste(paste("Posterior predictive model check for SRRT distribution \nSubject",subject_idx, sep=" "), "\nat SSD =",sep=" "),ssd_obs,sep=" ")
      }

      SRRT_pred = post_pred_srrt_use[d]

      x = max(SRRT_obs)+600
      y = max(hist(SRRT_obs,plot=F)$density)
      y = y+0.6*y
      par(mar = c(12, 10, 10, 8) + 0.1,cex.main=1.3)
      hist(SRRT_obs,freq=F,xlim=c(0,x),ylim=c(0,y),xlab="SRRT (ms)",las=1,main = my_main,ylab="Density")
      for(j in 1:n_post_samples){
        pdat = unlist(SRRT_pred[[1]][j])
          if(length(pdat)>1){
            lines(density(pdat),col="gray")
          }
      }
    hist(SRRT_obs,freq=F,add=T)
    }
  }
}
# produce posterior predictive model checks
posterior_predictions_computation = function(pars,n_post_samples,data){
   suppressPackageStartupMessages(library("gamlss.dist"))
   suppressPackageStartupMessages(library("vioplot"))

	if(pars$n_subject==1){
    post_pars = sample_joint_posterior(pars = pars,n_post_samples = n_post_samples)

    obs_go = load_prep_observed_go(pars = pars,subject_idx = 1,data=data)
		obs_srrt = load_prep_observed_srrt(pars = pars,subject_idx = 1,data = data)

		post_predictions_go = generate_posterior_predictions_go(n_go=obs_go$n_go,n_post_samples=n_post_samples,
                                                            par_vectors = post_pars$par_vectors,params = post_pars$params)

		post_predictions_srrt = generate_posterior_predictions_srrt(delays_all = obs_srrt$delays_all,n_delays_all = obs_srrt$n_delays_all,
                                                                n_observed_all = obs_srrt$n_observed_all,
                                                                n_observed_srrt_all = obs_srrt$n_observed_srrt_all,
                                                                par_vectors = post_pars$par_vectors,n_post_samples = n_post_samples,
                                                                params = post_pars$params,n_subj=obs_srrt$n_subj)

    plot_posterior_predictions_go(go_rt=obs_go$go_rt,n_go=obs_go$n_go,n_subj=obs_go$n_subj,go_pred=post_predictions_go$go_pred,
                                  subject_idx=n,n_post_samples=n_post_samples)

    post_pred_output_median_srrt(n_delays = obs_srrt$n_delays,delays = obs_srrt$delays,n_subj=obs_srrt$n_subj,
                                 median_post_pred_srrt = post_predictions_srrt$median_post_pred_srrt,
                                 median_observed_srrt = obs_srrt$median_observed_srrt,
                                 n_observed_srrt = obs_srrt$n_observed_srrt)

    post_pred_output_inhibition_function(params = post_pars$params,pars = mcmc.samples,n_subj = obs_srrt$n_subj,
                                         n_delays_all = obs_srrt$n_delays_all,delays_all = obs_srrt$delays_all,
                                         RR_pred = post_predictions_srrt$RR_pred,RR_obs = obs_srrt$RR_obs,
                                         n_observed_all = obs_srrt$n_observed_all,subject_idx=n)

    post_pred_output_srrt_distribution(delays=obs_srrt$delays,n_observed_srrt=obs_srrt$n_observed_srrt,n_subj = obs_srrt$n_subj,
                                       post_pred_srrt=post_predictions_srrt$post_pred_srrt,
                                       srrt_obs=obs_srrt$srrt_obs,n_post_samples=n_post_samples)

    } else{

      	for(n in 1:pars$n_subject){
          post_pars = sample_joint_posterior(pars = pars,n_post_samples = n_post_samples,subject_idx=n)

          obs_go = load_prep_observed_go(pars = pars,subject_idx = n,data=data)
		      obs_srrt = load_prep_observed_srrt(pars = pars,subject_idx = n,data = data)

          post_predictions_go = generate_posterior_predictions_go(n_go=obs_go$n_go,n_post_samples=n_post_samples,
                                                                  par_vectors = post_pars$par_vectors,params = post_pars$params)

          post_predictions_srrt = generate_posterior_predictions_srrt(delays_all = obs_srrt$delays_all,n_delays_all = obs_srrt$n_delays_all,
                                                                      n_observed_all = obs_srrt$n_observed_all,
                                                                      n_observed_srrt_all = obs_srrt$n_observed_srrt_all,
                                                                      par_vectors = post_pars$par_vectors,n_post_samples = n_post_samples,
                                                                      params = post_pars$params,n_subj=obs_srrt$n_subj,subject_idx=n)

          plot_posterior_predictions_go(go_rt=obs_go$go_rt,n_go=obs_go$n_go,n_subj = obs_go$n_subj,go_pred=post_predictions_go$go_pred,
                                        subject_idx=n,n_post_samples=n_post_samples)


          post_pred_output_median_srrt(n_delays = obs_srrt$n_delays,delays = obs_srrt$delays,n_subj = obs_srrt$n_subj,
                                       median_post_pred_srrt = post_predictions_srrt$median_post_pred_srrt,
                                       median_observed_srrt = obs_srrt$median_observed_srrt,
                                       n_observed_srrt = obs_srrt$n_observed_srrt,subject_idx=n)

          post_pred_output_inhibition_function(params = post_pars$params,pars = mcmc.samples,n_subj = obs_srrt$n_subj,
                                               n_delays_all = obs_srrt$n_delays_all,delays_all = obs_srrt$delays_all,
                                               RR_pred = post_predictions_srrt$RR_pred,RR_obs = obs_srrt$RR_obs,
                                               n_observed_all = obs_srrt$n_observed_all,subject_idx=n)

          post_pred_output_srrt_distribution(delays=obs_srrt$delays,n_observed_srrt=obs_srrt$n_observed_srrt,n_subj = obs_srrt$n_subj,
                                             post_pred_srrt=post_predictions_srrt$post_pred_srrt,
                                             srrt_obs=obs_srrt$srrt_obs,n_post_samples=n_post_samples,subject_idx=n)

      	}
    	}
}
