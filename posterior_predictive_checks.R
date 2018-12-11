
# ------------------------------------------------------------------------------ #
# produce posterior predictive model checks
# ------------------------------------------------------------------------------ #
posterior_predictions = function(samples, n_post_samples, data, n_subject){
   suppressPackageStartupMessages(library("gamlss.dist"))
   suppressPackageStartupMessages(library("vioplot"))

  for (subj_num in unique(data$subj_idx)) {
    post_pars = sample_joint_posterior(samples, n_post_samples, subj_num)
    obs_go    = load_prep_observed_go(samples, subj_num, data)
    obs_srrt  = load_prep_observed_srrt(samples, subj_num, data)

    cat(sprintf('Processing subject %s', toString(subj_num)), '\n')
    cat('Generating predictions: Go...', '\n')
    post_preds_go = generate_post_preds_go(obs_go$n_go, n_post_samples, post_pars$par_vectors, subj_num)

    cat('Generating predictions: SRRT...', '\n')
    post_preds_srrt = generate_post_preds_srrt(obs_srrt, post_pars$par_vectors, n_post_samples, n_subject, subj_num)

    cat('Plotting predictions: Go...', '\n')
    plot_post_preds_go(go_rt=obs_go$go_rt, go_pred=post_preds_go$go_pred, subject_idx=subj_num)

    cat('Plotting predictions: SRRT...', '\n')
    plot_post_preds_srrt(delays=obs_srrt$delays, n_observed_srrt=obs_srrt$n_observed_srrt,
      post_pred_srrt=post_preds_srrt$post_pred_srrt,
      srrt_obs=obs_srrt$srrt_obs,subject_idx=subj_num)

    cat('Getting median SRRTs...', '\n')
    post_pred_output_median_srrt(n_delays = obs_srrt$n_delays,delays = obs_srrt$delays,
      n_subj = obs_srrt$n_subj,
      median_post_pred_srrt = post_preds_srrt$median_post_pred_srrt,
      median_observed_srrt = obs_srrt$median_observed_srrt,
      n_observed_srrt = obs_srrt$n_observed_srrt,subject_idx=subj_num)

    cat('Getting inhibition function...', '\n')
    post_pred_output_inhibition_function(params = post_pars$params,pars = mcmc.samples, obs_srrt = obs_srrt,
      RR_pred = post_preds_srrt$RR_pred,subject_idx=subj_num)

    cat('\n')
  }
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# Sample from the joint distribution of the parameters
# ------------------------------------------------------------------------------ #
sample_joint_posterior = function(pars, n_post_samples, subject_idx = NULL){

  regex  <- paste('*_subj.', subject_idx, '*', sep='')
  params <- grep(glob2rx(regex), colnames(pars$traces), value=TRUE)
  n_params <- length(params)

  par_vectors = array(NA, dim = c(n_post_samples, n_params), list(NULL, params))
  it_ind = sample(1:(pars$n_samples*pars$n_chains), n_post_samples, replace=F)


  for(p in params){
    par_vectors[,p] = as.vector(pars$traces[,p,])[it_ind]
  }

  return(list(par_vectors = par_vectors, params = params))
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# Go
# ------------------------------------------------------------------------------ #
load_prep_observed_go = function(pars, subject_idx, data){

  n_subj = pars$n_subject
  if(n_subj == 1){
    data = cbind(subj_idx = rep(1,nrow(data)),data)
  }

  ###select go RTs
  go_rt = data$rt[data$subj_idx==subject_idx&data$ss_presented==0&data$inhibited==-999&data$ssd==-999&data$rt!=-999]
  n_go = length(go_rt)

  return(list(go_rt=go_rt,n_go=n_go,n_subj=n_subj))
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# DFASDLKFJASDFKL
# ------------------------------------------------------------------------------ #
load_prep_observed_srrt = function(pars, subject_idx, data){

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
  # and get observed srrts per delay; and compute median
  delays = sort(unique(data$ssd[data$subj_idx==subject_idx&data$ss_presented==1&data$inhibited==0&data$rt!=-999]))
  n_delays = length(delays)

  median_observed_srrt = rep(NA,length(delays))
  n_observed_srrt = rep(NA,length(delays))
  n_observed_inhibit = rep(NA,length(delays))
  srrt_obs = list()

  for (d in 1:n_delays) {
    srrt_temp = data$rt[data$subj_idx == subject_idx & data$ss_presented == 1 & data$inhibited == 0 & data$ssd == delays[d]]
    srrt_obs[[d]] = srrt_temp
    median_observed_srrt[d] = median(srrt_temp)
    n_observed_srrt[d] = length(srrt_temp)
    n_observed_inhibit[d] = nrow(data[data$subj_idx == subject_idx & data$ss_presented == 1 & data$inhibited == 1 & data$rt == -999 & data$ssd == delays[d],])
  }
  n_observed = n_observed_srrt + n_observed_inhibit

  # response rate per delay
  RR_obs = n_observed_srrt_all/n_observed_all

  return(list(delays = delays, n_delays = n_delays, delays_all = delays_all, n_delays_all = n_delays_all,median_observed_srrt = median_observed_srrt, n_observed_srrt = n_observed_srrt, n_observed = n_observed,
      n_observed_srrt_all = n_observed_srrt_all, n_observed_all = n_observed_all,n_subj = n_subj, RR_obs = RR_obs,srrt_obs = srrt_obs))
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# Posterior predictions for GO RTs
# ------------------------------------------------------------------------------ #
generate_post_preds_go = function(n_go, n_post_samples, par_vectors, subj_num){

  go_pred <- matrix(NA, n_go, n_post_samples)
  params  <- paste(c('mu_go_subj', 'sigma_go_subj', 'tau_go_subj'), subj_num, sep = '.')
  for (j in 1:n_post_samples) {
    mu    <- par_vectors[j, params[1]]
    sigma <- par_vectors[j, params[2]]
    nu    <- par_vectors[j, params[3]]

    go_pred[,j] = rexGAUS(n_go, mu, sigma, nu)
  }

  return(list(go_pred = go_pred))
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# Plot posterior predictions for go dists and the GO RTs
# ------------------------------------------------------------------------------ #
plot_post_preds_go = function(go_rt, go_pred, subject_idx = NULL){

  # Layout
  par(mar = c(13, 11, 11, 10) + 0.1,cex.main=1.5)
  x = max(go_rt)+600
  y = max(hist(go_rt, plot=F)$density)
  y = y + 0.6*y

  # Plot the go RT density and histogram
  layout(matrix(1))
  subj_str <- ifelse(is.null(subject_idx), '', ', Subject')
  main <- paste("PPC for go RT dist", subject_idx, sep=" ")

  hist(go_rt, freq=F, xlim=c(0,x), ylim=c(0,y), xlab="go RT (ms)", las=1, main=main, ylab="Density")
  hist(go_rt, freq=F, add=T)

  # Add posterior density predictions
  n_samples <- dim(go_pred)[2]
  for(j in 1:n_samples) {
    lines(density(go_pred[,j]),col="gray")
  }
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# Generates posterior predictions for median srrts, srrt distributions, and inhibition functions
# ------------------------------------------------------------------------------ #
generate_post_preds_srrt = function(obs_srrt, par_vectors, n_post_samples, n_subj, subj_num = NULL){

  params_go    <- paste(c('mu_go_subj'  , 'sigma_go_subj'  , 'tau_go_subj'  ), subj_num, sep = '.')
  params_stop  <- paste(c('mu_stop_subj', 'sigma_stop_subj', 'tau_stop_subj'), subj_num, sep = '.')

  median_post_pred_srrt <- matrix(NA, obs_srrt$n_delays_all, n_post_samples)
  RR_pred               <- matrix(NA, obs_srrt$n_delays_all, n_post_samples)
  post_pred_srrt        <- matrix(NA, obs_srrt$n_delays_all, n_post_samples)

  if (n_subj == 1) {
    prob = par_vectors[, "pf_stop"]
  } else {
    prob = pnorm(par_vectors[, paste("pf_stop_subj.", subj_num, sep = "")])
  }

  for (d in 1:obs_srrt$n_delays_all) {

    n_tf = obs_srrt$n_observed_all[d]*prob
    n_tf = round(n_tf, 0)
    n_notf = obs_srrt$n_observed_all[d] - n_tf

    for (j in 1:n_post_samples) {

      go_mu    <- par_vectors[j, params_go[1]]
      go_sigma <- par_vectors[j, params_go[2]]
      go_nu    <- par_vectors[j, params_go[3]]

      stop_mu    <- par_vectors[j, params_stop[1]]
      stop_sigma <- par_vectors[j, params_stop[2]]
      stop_nu    <- par_vectors[j, params_stop[3]]

      if (n_tf[j] > 0) {
        ### not triggered trials
        srrt_tf = rexGAUS(n_tf[j], go_mu, go_sigma, go_nu)
      } else srrt_tf = NULL


      if (n_notf[j] > 0) {
        ###triggered trials
        go_rts  <- rexGAUS(n_notf[j], go_mu, go_sigma, go_nu)
        stop_rt <- obs_srrt$delays_all[d] + rexGAUS(n_notf[j], stop_mu, stop_sigma, stop_nu)

        srrt = go_rts[go_rts <= stop_rt]
       } else srrt = NULL

      all_pred_srrt = c(srrt, srrt_tf)
      post_pred_srrt[[d]][[j]] = list(all_pred_srrt)
      median_post_pred_srrt[d,j] = median(all_pred_srrt)
      RR_pred[d,j] = (length(all_pred_srrt))/obs_srrt$n_observed_all[d]
    }
   }
   median_post_pred_srrt = median_post_pred_srrt[obs_srrt$n_observed_srrt_all>0,]
   post_pred_srrt = post_pred_srrt[obs_srrt$n_observed_srrt_all>0]

   return(list(median_post_pred_srrt = median_post_pred_srrt, RR_pred=RR_pred, post_pred_srrt = post_pred_srrt))
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# Generates output for posterior predictions for srrts and corresponding p values
# ------------------------------------------------------------------------------ #
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
  rownames(post_pred_summary) = c("Number of observed srrt", "Observed median srrt","Average predicted srrt","One-sided p value","Two-sided p value")

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
    table_title <- "Post. pred. p vals for median srrt"
    plot_title  <- "Posterior predictive model check for median srrt"
  } else {
    table_title <- paste("Post. pred. p vals, median srrt, Subject",subject_idx,sep=" ")
    plot_title  <- paste("PPC, median srrt, Subject",subject_idx,sep=" ")
  }

  padding     <- unit(1,"line")
  table_title <- textGrob(table_title, gp=gpar(fontsize=20))
  footnote    <- textGrob("footnote", x=0, hjust=0, gp=gpar( fontface="italic"))

  table <- gtable_add_rows(table, heights = grobHeight(table_title) + padding, pos = 0)
  table <- gtable_add_rows(table, heights = grobHeight(footnote)+ padding)

  table <- gtable_add_grob(table, list(table_title, footnote), t=c(1, nrow(table)), l=c(1,2), r=ncol(table))

  #gt <- gTree(children=gList(table, table_title))
  grid.newpage()
  grid.draw(table)

  # Make violin plot
  my_min = min(c(min(median_post_pred_srrt,na.rm=T),min(median_observed_srrt)))-80
  my_max = max(c(max(median_post_pred_srrt,na.rm=T),max(median_observed_srrt)))+80
  ran = round(c(my_min,my_max))
  par(cex=1.2,cex.lab = 1.2,mar=c(11, 10, 10, 8),cex.main=1.4)

  # These are the circles I think...
  plot(1:length(delays),median_observed_srrt,pch=17,xlab="SSD (ms)",axes=F,
    ylim=ran,xlim=c(0,n_delays+1),ylab ="",main=plot_title,cex=1.3)
  axis(1,at=0:(n_delays+1),labels=c(NA,delays,NA))
  axis(2,at=ran,las=2)
  mtext("Median srrt (ms)",2, cex=1.4,line=1)
  for(i in 1:n_delays){
    preds_med = median_post_pred_srrt[i,][!is.na(median_post_pred_srrt[i,])]
    if (length(preds_med) >0) {
      vioplot(preds_med,add=T,at=i,col="gray",cex=2)
    }
  }

  # These are the observations, the dashed line and triangles
  points(1:n_delays,median_observed_srrt,cex=1.5,pch=17)
  lines(1:n_delays,median_observed_srrt,lwd=2,lty=2)

  # These are mean values, the black line
  lines(1:n_delays,average_post_pred,lwd=2,lty=1)
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# Generates output for posterior predictions for inhibition function and corresponding p values
# ------------------------------------------------------------------------------ #
post_pred_output_inhibition_function = function(params, pars, obs_srrt, RR_pred, subject_idx=NULL){

  n_subj         <- obs_srrt$n_subj
  n_delays_all   <- obs_srrt$n_delays_all
  delays_all     <- obs_srrt$delays_all
  RR_obs         <- obs_srrt$RR_obs
  n_observed_all <- obs_srrt$n_observed_all

  pvalue_one_if = rep(NA,n_delays_all)
  pvalue_two_if = rep(NA,n_delays_all)

  if(n_subj == 1) {
    prob = pars$traces[,params[7],]
  } else {
    prob = pnorm(pars$traces[,grep(glob2rx("pf_stop*"), params, value=TRUE),])
  }
  q_tf = quantile(prob, prob = c(0.025,0.975))
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

  colnames(post_pred_if_summary) = paste(rep("",each=n_delays_all),delays_all,sep="")
  rownames(post_pred_if_summary) = c("# Stop trials", "Observed RR", "Avg. predicted RR","One-sided p-val","Two-sided p-val")

  # Write table
  nc = ncol(post_pred_if_summary)
  fsize=19-nc
  if (fsize<1)  {
      fsize=1
  }
  table = tableGrob(post_pred_if_summary) #,name="IF", gpar.coretext = gpar(fontsize=fsize),equal.width=F,gpar.coltext = gpar(fontsize=fsize),gpar.rowtext = gpar(fontsize=fsize))
  h <- grobHeight(table)
  w <- grobWidth(table)

  # Set titles for table and plots
  subj_str    <- ifelse(is.null(subject_idx), '', ', Subject')
  table_title <- paste("Response rate data", subj_str, subject_idx, sep=" ")
  plot_title  <- paste("PPC for inhibition function", subj_str, subject_idx, sep=" ")

  padding     <- unit(1,"line")
  table_title <- textGrob(table_title, gp=gpar(fontsize=20))
  footnote    <- textGrob("footnote", x=0, hjust=0, gp=gpar( fontface="italic"))

  table <- gtable_add_rows(table, heights = grobHeight(table_title) + padding, pos = 0)
  table <- gtable_add_rows(table, heights = grobHeight(footnote)+ padding)

  table <- gtable_add_grob(table, list(table_title, footnote), t=c(1, nrow(table)), l=c(1,2), r=ncol(table))

  grid.newpage()
  grid.draw(table)

  # Make violin plot
  par(cex=1.2,cex.lab = 1.2,mar=c(11, 10, 10, 8),cex.main=1.4)
  plot(1:length(delays_all),RR_obs,pch=17,xlab="SSD (ms)",axes=F,
    ylim=c(0,1),xlim=c(0,n_delays_all+1),ylab ="",main=plot_title,cex=1.3)
  axis(1,at=0:(n_delays_all+1),labels=c(NA,delays_all,NA))
  axis(2,at=c(0,1),las=2)
  mtext("Response Rate",2, cex=1.4,line=1)

  for (i in 1:n_delays_all) {
    if (length(unique(RR_pred[i,])) != 1) {
      vioplot(RR_pred[i,],add=T,at=i,col="gray",cex=2)
    }
  }

  points(1:n_delays_all, RR_obs, cex=1.5, pch=17)
  lines( 1:n_delays_all, RR_obs, lwd=2, lty=2)
  lines( 1:n_delays_all, average_pred_if, lwd=2, lty=1)

  # add TF parameter for reference
  lines(c(1, n_delays_all), rep(q_tf[1],   2), lwd=2, lty=4)
  lines(c(1, n_delays_all), rep(q_tf[2],   2), lwd=2, lty=4)
  lines(c(1, n_delays_all), rep(mean_prob, 2), lwd=2, lty=3)
}
# ------------------------------------------------------------------------------ #


# ------------------------------------------------------------------------------ #
# Generates output for posterior predictions for srrt distributions
# only for SSDs with at least 10 observed SRRTs
# ------------------------------------------------------------------------------ #
plot_post_preds_srrt = function(delays,n_observed_srrt,post_pred_srrt,srrt_obs,subject_idx=NULL){

  ind = n_observed_srrt>=5
  dels = delays[ind]
  n_dels = length(dels)
  srrt_obs_use = srrt_obs[ind]
  post_pred_srrt_use =  post_pred_srrt[ind]

  if(n_dels>0){
    for(d in 1:n_dels){
      srrt_obs = unlist(srrt_obs_use[d])
      ssd_obs = dels[d]

      subj_str <- ifelse(is.null(subject_idx), '', ', Subject')
      main <- paste("PPC for SRRT dist", subj_str, subject_idx, "at SSD =", ssd_obs, sep=" ")

      srrt_pred = post_pred_srrt_use[d]

      x = max(srrt_obs)+600
      y = max(hist(srrt_obs,plot=F)$density)
      y = y+0.6*y
      par(mar = c(12, 10, 10, 8) + 0.1,cex.main=1.3)
      hist(srrt_obs, col='black', freq=F,xlim=c(0,x),ylim=c(0,y),xlab="srrt (ms)",las=1,main = main,ylab="Density")

      preds <- unlist(srrt_pred[[1]])
      n_post_samples <- length(preds)
      for(j in 1:n_post_samples){
        pdat = unlist(srrt_pred[[1]][j])
          if(length(pdat)>1){
            lines(density(pdat),col="gray")
          }
      }
    hist(srrt_obs,freq=F,add=T)
    }
  }
}
# ------------------------------------------------------------------------------ #

