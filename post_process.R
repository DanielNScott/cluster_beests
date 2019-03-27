##
## This file coordinates the R post processing routines.
##

# Load some required libraries
suppressPackageStartupMessages(library("grid"))
suppressWarnings(library("gridExtra"))
#library(futile.logger)
#flog.appender(appender.file('post_process.log'))

# Load most of the helper functions
source("post_process_aux.R")
source("posterior_predictive_checks.R")

# Check if Rscript invocation or called from R
# and determine analysis directory
args <- commandArgs()
if (any(args == '--args')) {
   analy_dir  <- args[which(args == '--args')+1]
}
cat('\npost_process.R believes the analysis dir is:\n')
cat(paste(analy_dir, '\n'))

data.file  <- paste(analy_dir, 'sst_data.csv', sep = '')
analy_file <- paste(analy_dir, 'analysis.txt', sep = '')

# Read data file:
data <- read.csv(data.file,head = TRUE,sep = ",")
ncol <- ncol(data)
vars <- read.table(analy_file, sep = ",")

# Extract variables:
extract  <- function(str) {vars[vars[,1] == str,2]}
convert  <- function(x) {as.numeric(as.vector(x))}

samples	          <- convert(extract("samples"))
burn.in            <- convert(extract("burn-in"))
number.of.chains   <- convert(extract("number of chains"))
thinning	          <- convert(extract("thinning"))
estimates.for.all  <- extract("estimates for subjects or groups") == "All"
summary.stats      <- extract("summary stats") == "1"
posterior.dists    <- extract("posterior distributions") == "1"
mcmc.chains	       <- extract("mcmc chains") == "1"
deviance           <- extract("deviance") == "1"
post.preds	       <- extract("posterior predictors") == "1"
post.preds.samples <- convert(extract("posterior predictor samples"))
num.cores          <- convert(extract("cpu cores"))
int_lower          <- convert(extract("limits of integration lower"))
int_upper          <- convert(extract("limits of integration upper"))
version            <- extract("model trigger failure") == "0"

# Log input fidelity:
cat(sprintf('Vars read'), '\n')
cat(sprintf('Samples      : %s', samples), '\n')
cat(sprintf('Burn in      : %s', burn.in), '\n')
cat(sprintf('Num Chains   : %s', number.of.chains), '\n')
cat(sprintf('Thinning     : %s', thinning), '\n')
cat(sprintf('Estimates    : %s', estimates.for.all), '\n')
cat(sprintf('Stats        : %s', summary.stats), '\n')
cat(sprintf('Posteriors   : %s', posterior.dists), '\n')
cat(sprintf('Chains       : %s', mcmc.chains), '\n')
cat(sprintf('Deviance     : %s', deviance), '\n')
cat(sprintf('Predictors   : %s', post.preds), '\n')
cat(sprintf('Pred. Samples: %s', post.preds.samples), '\n')
cat(sprintf('Num Cores    : %s', num.cores), '\n')


### ------- Set up params to look at -------- ###
# Specify the pieces defining a set of variables to look at:
param_bases <- list('mu_', 'sigma_', 'tau_', 'shift_')
param_dists <- list('go', 'stop')
param_var   <- list('', '_sd')
param_lims  <- list('_lower', '_upper')

param_names <- as.vector(outer(param_bases, param_dists, function(x,y) paste(x,y,sep = '')))
param_wsds  <- as.vector(outer(param_names, param_var  , function(x,y) paste(x,y,sep = '')))

# Exceptions to rules above
param_names <- c(param_names, 'pf_stop')
param_wsds  <- c(param_wsds , 'pf_stop', 'pf_stop_sd')

# Everything gets a lower lim and an upper lim
terms_wlims <- as.vector(outer(param_wsds , param_lims , function(x,y) paste(x,y,sep = '')))

# Read upper and lower limits, otherwise tons of copy paste
for (term in terms_wlims) {
  term_string_spaced <- gsub("_", " ", term)
  assign(term, convert(vars[vars[,1] == term_string_spaced, 2]) )
}

# Pack prior bounds into a list of lists
priors <- list()
for (name in param_wsds) {
  lims <- as.vector(outer(name, param_lims , function(x,y) paste(x,y,sep = '')))
  priors[[name]] <- c(get(lims[1]), get(lims[2]))
}

# More exceptions...
priors$pf_stop[3:4] <- c(convert(extract("pf stop mean")), convert(extract('pf stop sd')))
### ----------------------------------------- ###


### ------------- Some Misc Stuff ------------- ###
check_silly_things(samples, burn.in, thinning, number.of.chains, post.preds.samples, int_lower, int_upper, priors)

# Open pdf file to save output
pdf(paste(analy_dir, "/output.pdf", sep = ''), height = 11, width = 8.5)

# Get MCMC output.
cat('Retrieving MCMC samples...', '\n')
mcmc.samples <- read_prep(dir = analy_dir, n_chains = number.of.chains)
### ---


### --------- Get Summary Statistics ---------- ###
# Group summary statistics
cat('Calculating summary stats for group', '\n')
mcmc.samples$traces <- summary_stats(traces = mcmc.samples$traces, params = param_names, n_subj = mcmc.samples$n_subject)

# Individual summary statistics
for (subj_num in mcmc.samples$subjects){
  cat(sprintf('Calculating summary stats for subject %s.', toString(subj_num)), '\n')
  subj_param_names <- paste(param_names, '_subj.', subj_num, sep = '')
  mcmc.samples$traces <- summary_stats(traces = mcmc.samples$traces, params = subj_param_names, n_subj = mcmc.samples$n_subject, subj_idx = subj_num)
}
### ---

# Write modified traces back out (GoRT, SSRT, GoSd, SSRTSd have been appended)
for(i in 1:number.of.chains){
  write.table(x = mcmc.samples$traces[,,i], file = paste(analy_dir, "/new_parameters", i, ".csv", sep=""),
   col.names=TRUE, row.names=FALSE, sep=";", quote=FALSE)
}

### ------ Chains, posteriors, and PPCs ------ ###
cat('Calculating posterior distributions...', '\n')
plot_posteriors_wrapper(mcmc.samples, priors)

cat("Running posterior predictive model checks. This might take a while...", '\n')
posterior_predictions(mcmc.samples, post.preds.samples, data, mcmc.samples$n_subject)
### ---


dev.off()
cat("BEESTS is done!", '\n')
