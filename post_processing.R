##
## This file coordinates the R post processing routines.
##

# Load some required libraries
suppressPackageStartupMessages(library("grid"))
suppressWarnings(library("gridExtra"))
library(futile.logger)

# Load most of the helper functions
source("post_processing_wtf_aux.R")
source("posterior_predictive_checks.R")

# Read Arguments:
args       <- commandArgs()
analy_dir  <- args[1]
data.file  <- paste(args[1], 'sst_data.csv', sep='')
analy_file <- paste(args[1], 'analysis.txt', sep='')

flog.info('Analysis directory provided: %s', args[1])

# Read data file:
data <- read.csv(data.file,head=TRUE,sep=",")
ncol <- ncol(data)
vars <- read.table(analy_file, sep=",")

# Extract variables:
samples	        <- as.numeric(as.vector(vars[vars[,1]=="samples",2]))
burn.in          <- as.numeric(as.vector(vars[vars[,1]=="burn-in",2]))
number.of.chains <- as.numeric(as.vector(vars[vars[,1]=="number of chains",2]))
thinning	        <- as.numeric(as.vector(vars[vars[,1]=="thinning",2]))
estimates.for.all<- (vars[vars[,1]=="estimates for subjects or groups",2]=="All")
summary.statistics<- (vars[vars[,1]=="summary statistics",2]=="1")
posterior.distributions<- (vars[vars[,1]=="posterior distributions",2]=="1")
mcmc.chains	     <- (vars[vars[,1]=="mcmc chains",2]=="1")
deviance         <- (vars[vars[,1]=="deviance",2]=="1")
posterior.predictors	<- (vars[vars[,1]=="posterior predictors",2]=="1")
posterior.predictors.samples	<- as.numeric(as.vector(vars[vars[,1]=="posterior predictor samples",2]))
num.cores        <- as.numeric(as.vector(vars[vars[,1]=="cpu cores",2]))
int_lower        <- as.numeric(as.vector(vars[vars[,1]=="limits of integration lower",2]))
int_upper        <- as.numeric(as.vector(vars[vars[,1]=="limits of integration upper",2]))
version          <- (vars[vars[,1]=="model trigger failure",2]=="0")

# Log input fidelity:
flog.info('Vars read')
flog.info('Samples      : %s', samples)
flog.info('Burn in      : %s', burn.in)
flog.info('Num Chains   : %s', number.of.chains)
flog.info('Thinning     : %s', thinning)
flog.info('Estimates    : %s', estimates.for.all)
flog.info('Stats        : %s', summary.statistics)
flog.info('Posteriors   : %s', posterior.distributions)
flog.info('Chains       : %s', mcmc.chains)
flog.info('Deviance     : %s', deviance)
flog.info('Predictors   : %s', posterior.predictors)
flog.info('Pred. Samples: %s', posterior.predictors.samples)
flog.info('Num Cores    : %s', num.cores)


### ------- Set up params to look at -------- ###
# Specify the pieces defining a set of variables to look at:
param_bases <- list('mu_', 'sigma_', 'tau_', 'shift_')
param_dists <- list('go', 'stop')
param_var   <- list('', '_sd')
param_lims  <- list('_lower', '_upper')

param_names <- as.vector(outer(param_bases, param_dists, function(x,y) paste(x,y,sep='')))
param_wsds  <- as.vector(outer(param_names, param_var  , function(x,y) paste(x,y,sep='')))

# Exceptions to rules above
param_names <- c(param_names, 'pf_stop')
param_wsds  <- c(param_wsds , 'pf_stop', 'pf_stop_sd')

# Everything gets a lower lim and an upper lim
terms_wlims <- as.vector(outer(param_wsds , param_lims , function(x,y) paste(x,y,sep='')))

# Read upper and lower limits, otherwise tons of copy paste
for (term in terms_wlims) {
  term_string_spaced <- gsub("_", " ", term)
  assign(term, as.numeric(as.vector(vars[vars[,1] == term_string_spaced, 2])) )
}

# Pack prior bounds into a list of lists
priors <- list()
for (name in param_wsds) {
  lims <- as.vector(outer(name, param_lims , function(x,y) paste(x,y,sep='')))
  priors[[name]] <- c(get(lims[1]), get(lims[2]))
}

# More exceptions...
priors$pf_stop[3:4] <- c(as.numeric(as.vector(vars[vars[,1] == 'pf stop mean', 2])), as.numeric(as.vector(vars[vars[,1] == 'pf stop sd', 2])))


### ------------- Some Misc Stuff ------------- ###
check_silly_things(samples, burn.in, thinning, number.of.chains, posterior.predictor.samples, int_lower, int_upper, priors)

# Open pdf file to save output
pdf(paste(analy_dir, "/output.pdf", sep=''), height=11, width=8.5)

# Get MCMC output.
mcmc.samples <- read_prep(dir = analy_dir, n_chains = number.of.chains, n_params = length(param_names))


### --------- Get Summary Statistics ---------- ###
print('Calculating summary statistics...')

# Group summary statistics
flog.info('Calculating summary stats for group')
summary_stats(pars = mcmc.samples, params = param_names, n_subj = mcmc.samples$n_subject)

# Individual summary statistics
for (n in 1:mcmc.samples$n_subject){
  flog.info('Calculating summary stats for subject %s.', n)
  subj_param_names <- paste(param_names, '_subj.', n, sep='')
  summary_stats(pars = mcmc.samples, params = subj_param_names, n_subj = mcmc.samples$n_subject, subj_idx = n)
}


### ------------- Get posteriors --------------- ###
print('Calculating posterior distributions...')
plot_posteriors_wrapper(pars = mcmc.samples, priors = priors)


### -------- Posterior Pred. Checks ------------ ###
print("Running posterior predictive model checks. This might take a while...")
posterior_predictions(pars = mcmc.samples, n_post_samples = posterior.predictors.samples, data=data, n_subj = mcmc.samples$n_subject)


cl <- dev.off()
print("BEESTS is done!")















