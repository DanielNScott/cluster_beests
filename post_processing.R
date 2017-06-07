##
## This file coordinates the R post processing routines.
##

# Read Arguments:
args             <- commandArgs(trailingOnly=TRUE)
path_analysisDir <- args[1]
data.file        <- paste(args[1], 'sst_data.csv', sep='')
path__analysisDescription <- paste(args[1], 'analysis.txt', sep='')

dat <- read.csv(data.file,head=TRUE,sep=",")
ncol <- ncol(dat)
vars <- read.table(path__analysisDescription, sep=",")

samples	<- as.numeric(as.vector(vars[vars[,1]=="samples",2]))
burn.in	<- as.numeric(as.vector(vars[vars[,1]=="burn-in",2]))
number.of.chains	<- as.numeric(as.vector(vars[vars[,1]=="number of chains",2]))
thinning	<- as.numeric(as.vector(vars[vars[,1]=="thinning",2]))
estimates.for.all	<- (vars[vars[,1]=="estimates for subjects or groups",2]=="All")
summary.statistics	<- (vars[vars[,1]=="summary statistics",2]=="1")
posterior.distributions	<- (vars[vars[,1]=="posterior distributions",2]=="1")
mcmc.chains	<- (vars[vars[,1]=="mcmc chains",2]=="1")
deviance	<- (vars[vars[,1]=="deviance",2]=="1")
posterior.predictors	<- (vars[vars[,1]=="posterior predictors",2]=="1")
posterior.predictors.samples	<- as.numeric(as.vector(vars[vars[,1]=="posterior predictor samples",2]))
num.cores	= as.numeric(as.vector(vars[vars[,1]=="cpu cores",2]))
int_lower	<- as.numeric(as.vector(vars[vars[,1]=="limits of integration lower",2]))
int_upper	<- as.numeric(as.vector(vars[vars[,1]=="limits of integration upper",2]))
version	<- (vars[vars[,1]=="model trigger failure",2]=="0")

#flog.info(samples)
#flog.info(burn.in)
#flog.info(number.of.chains)
#flog.info(thinning)
#flog.info(estimates.for.all)
#flog.info(summary.statistics)
#flog.info(posterior.distributions)
#flog.info(mcmc.chains)
#flog.info(deviance)
#flog.info(posterior.predictors)
#flog.info(posterior.predictors.samples)
#flog.info(num.cores)

#input vars
#print(vars)
go_mu_lower <- as.numeric(as.vector(vars[vars[,1]=="go mu lower",2]))
go_mu_upper <- as.numeric(as.vector(vars[vars[,1]=="go mu upper",2]))

go_sigma_lower <- as.numeric(as.vector(vars[vars[,1]=="go sigma lower",2]))
go_sigma_upper <- as.numeric(as.vector(vars[vars[,1]=="go sigma upper",2]))

go_tau_lower <- as.numeric(as.vector(vars[vars[,1]=="go tau lower",2]))
go_tau_upper <- as.numeric(as.vector(vars[vars[,1]=="go tau upper",2]))

stop_mu_lower <- as.numeric(as.vector(vars[vars[,1]=="stop mu lower",2]))
stop_mu_upper <- as.numeric(as.vector(vars[vars[,1]=="stop mu upper",2]))

stop_sigma_lower <- as.numeric(as.vector(vars[vars[,1]=="stop sigma lower",2]))
stop_sigma_upper <- as.numeric(as.vector(vars[vars[,1]=="stop sigma upper",2]))

stop_tau_lower <- as.numeric(as.vector(vars[vars[,1]=="stop tau lower",2]))
stop_tau_upper <- as.numeric(as.vector(vars[vars[,1]=="stop tau upper",2]))

pf_lower <- as.numeric(as.vector(vars[vars[,1]=="stop pf lower",2]))
pf_upper <- as.numeric(as.vector(vars[vars[,1]=="stop pf upper",2]))
pf_mean <- as.numeric(as.vector(vars[vars[,1]=="stop pf mean",2]))
pf_sd <- as.numeric(as.vector(vars[vars[,1]=="stop pf sd",2]))

go_mu_sd_lower <- as.numeric(as.vector(vars[vars[,1]=="go mu sd lower",2]))
go_mu_sd_upper <- as.numeric(as.vector(vars[vars[,1]=="go mu sd upper",2]))

go_sigma_sd_lower <- as.numeric(as.vector(vars[vars[,1]=="go sigma sd lower",2]))
go_sigma_sd_upper <- as.numeric(as.vector(vars[vars[,1]=="go sigma sd upper",2]))

go_tau_sd_lower <- as.numeric(as.vector(vars[vars[,1]=="go tau sd lower",2]))
go_tau_sd_upper <- as.numeric(as.vector(vars[vars[,1]=="go tau sd upper",2]))

stop_mu_sd_lower <- as.numeric(as.vector(vars[vars[,1]=="stop mu sd lower",2]))
stop_mu_sd_upper <- as.numeric(as.vector(vars[vars[,1]=="stop mu sd upper",2]))

stop_sigma_sd_lower <- as.numeric(as.vector(vars[vars[,1]=="stop sigma sd lower",2]))
stop_sigma_sd_upper <- as.numeric(as.vector(vars[vars[,1]=="stop sigma sd upper",2]))

stop_tau_sd_lower <- as.numeric(as.vector(vars[vars[,1]=="stop tau sd lower",2]))
stop_tau_sd_upper <- as.numeric(as.vector(vars[vars[,1]=="stop tau sd upper",2]))

pf_sd_lower <- as.numeric(as.vector(vars[vars[,1]=="stop pf sd lower",2]))
pf_sd_upper <- as.numeric(as.vector(vars[vars[,1]=="stop pf sd upper",2]))

if (version){
      prior_ran <- list(mu_go = c(go_mu_lower,go_mu_upper),sigma_go = c(go_sigma_lower,go_sigma_upper),
                 tau_go = c(go_tau_lower,go_tau_upper),mu_stop = c(stop_mu_lower,stop_mu_upper),
                 sigma_stop= c(stop_sigma_lower,stop_sigma_upper),tau_stop = c(stop_tau_lower,stop_tau_upper),
                 mu_go_var = c(go_mu_sd_lower,go_mu_sd_upper), sigma_go_var = c(go_sigma_sd_lower,go_sigma_sd_upper),
                 tau_go_var = c(go_tau_sd_lower,go_tau_sd_upper),mu_stop_var = c(stop_mu_sd_lower,stop_mu_sd_upper),
                 sigma_stop_var = c(stop_sigma_sd_lower,stop_sigma_sd_upper),tau_stop_var = c(stop_tau_sd_lower,stop_tau_sd_upper))
      } else {

      prior_ran <- list(mu_go = c(go_mu_lower,go_mu_upper),sigma_go = c(go_sigma_lower,go_sigma_upper),
                 tau_go = c(go_tau_lower,go_tau_upper),mu_stop = c(stop_mu_lower,stop_mu_upper),
                 sigma_stop= c(stop_sigma_lower,stop_sigma_upper),tau_stop = c(stop_tau_lower,stop_tau_upper),
                 p_tf = c(pf_lower,pf_upper,pf_mean,pf_sd),
                 mu_go_var = c(go_mu_sd_lower,go_mu_sd_upper), sigma_go_var = c(go_sigma_sd_lower,go_sigma_sd_upper),
                 tau_go_var = c(go_tau_sd_lower,go_tau_sd_upper),mu_stop_var = c(stop_mu_sd_lower,stop_mu_sd_upper),
                 sigma_stop_var = c(stop_sigma_sd_lower,stop_sigma_sd_upper),tau_stop_var = c(stop_tau_sd_lower,stop_tau_sd_upper),
                 p_tf_var = c(pf_sd_lower,pf_sd_upper))
      }

if (samples < 1){
        stop("The total number of MCMC samples must be greater than zero.")
}

if (samples <= burn.in){
        stop("The total number of MCMC samples must be greater than the number of burn-in samples.")
}

if (thinning < 1){
       stop("The thinning factor must be greater than 0.")
}

if  (((samples-burn.in)/thinning) < 1){
        stop("No MCMC samples will be retained. Increse the number of retained samples or decrese the thinning factor.")
}

if (posterior.predictors == T){
  if ((((samples-burn.in)/thinning)*number.of.chains) < posterior.predictors.samples){
        stop("The number of posterior predictive samples cannot be greater than the number of retained MCMC samples.")
  }
}

if (int_lower>=int_upper){
        stop("The lower limit of integration must be lower than the upper limit.")
}

check_priors_inits <- function(pars){
       for (par in pars){
            if  (as.numeric(as.vector(vars[vars[,1]==paste(par, " upper",sep=""),2])) <= as.numeric(as.vector(vars[vars[,1]==paste(par, " lower",sep=""),2]))) {
              if (par=="stop pf"){
                       stop("Silly prior: Upper bound for P(TF) must be higher than lower bound")
               } else {
                       stop(paste(paste("Silly prior: Upper bound for ", par,sep=""), " must be higher than lower bound.", sep=""))
               }
            }

           if  ((as.numeric(as.vector(vars[vars[,1]==paste(par, " start",sep=""),2])) < as.numeric(as.vector(vars[vars[,1]==paste(par, " lower",sep=""),2]))) | (as.numeric(as.vector(vars[vars[,1]==paste(par, " start",sep=""),2])) > as.numeric(as.vector(vars[vars[,1]==paste(par, " upper",sep=""),2])))){
            if (par=="stop pf"){
                       stop("Silly start value: Start value for P(TF) must be within prior range")
               } else {
                       stop(paste(paste("Start value for ", par,sep=""), " must be within prior range.", sep=""))
             }
            }
      }
}

if (version){
      check_priors_inits(pars=c("go mu", "go sigma", "go tau", "stop mu", "stop sigma", "stop tau","go mu sd", "go sigma sd", "go tau sd", "stop mu sd", "stop sigma sd", "stop tau sd"))
} else {
       if (ncol==4){
            if (pf_lower < 0 | pf_upper > 1){
                stop("The range of the uniform prior for P(TF) must be between 0 and 1")
            }
        }
        check_priors_inits(pars=c("go mu", "go sigma", "go tau", "stop mu", "stop sigma", "stop tau", "stop pf","go mu sd", "go sigma sd", "go tau sd", "stop mu sd", "stop sigma sd", "stop tau sd", "stop pf sd"))
}

if (version) source("Aux_BEEST.R")  else source("Aux_BEEST_wtf.R") # <- the wd for these will be the location of this run.R file

setwd(path_analysisDir)  # set the wd to the file location
suppressPackageStartupMessages(library("grid"))
suppressWarnings(library("gridExtra"))

if (summary.statistics == T|posterior.distributions==T|mcmc.chains==T|posterior.predictors==T){
  mcmc.samples <- read_prep(n_chains = number.of.chains)

  pdf("output.pdf", height=11, width=8.5)

    if (summary.statistics==T){
      summary_stats(pars = mcmc.samples,all_pars=estimates.for.all)
      print("Summary statistics are saved to file.")
    }

    if (posterior.distributions==T){
     #suppressWarnings(library(msm))
     plot_posteriors(pars = mcmc.samples,all_pars=estimates.for.all,priors = prior_ran)

     print("Posterior distributions are saved to file.")
    }

    if (mcmc.chains==T){
      plot_chains(pars = mcmc.samples,all_pars=estimates.for.all,priors = prior_ran)
      print("MCMC chains are saved to file.")
    }

    if (posterior.predictors==T){
      print("Running posterior predictive model checks. This might take a while...")
      suppressPackageStartupMessages(library("gamlss.dist"))
      suppressPackageStartupMessages(library("vioplot"))
      posterior_predictions_computation(pars = mcmc.samples,n_post_samples = posterior.predictors.samples,data=dat)
      print("Results of posterior predictive model checks are saved to file.")
    }
  cl <- dev.off()
  print("BEESTS is done!")
}















