#!/usr/bin/python
# -*- coding: utf-8 -*-
from __future__ import division
import os
import sys
import numpy as np
import numpy.lib.recfunctions as rec
import pymc as pm
from copy import copy
from scipy.stats import norm
from scipy.integrate import quad
try:
    from IPython.Debugger import Tracer; debug_here = Tracer()
except:
    def debug_here(): pass

#for import settings file
from kabuki import Hierarchical, Knode
from kabuki.utils import load_csv
import stop_likelihoods_wtf
from collections import OrderedDict

#settings file
path__analysisDescription = sys.argv[3]

#make dict for input arguments
vars = dict()
#load analysis.txt with input arguments and automatically set wd
with open(path__analysisDescription) as f:
    for line in f:
        eq_index = line.find(',')
        var_name = line[:eq_index].strip()
        var_name = var_name.replace('"', '').strip()
        number = line[eq_index + 1:].strip()
        number = number.replace('"', '').strip()
        vars[var_name] = number

#input vars
go_mu_lower = float(vars["go mu lower"])
go_mu_upper = float(vars["go mu upper"])
go_mu_start = float(vars["go mu start"]) 

go_sigma_lower = float(vars["go sigma lower"])
go_sigma_upper = float(vars["go sigma upper"])
go_sigma_start = float(vars["go sigma start"])

go_tau_lower = float(vars["go tau lower"])
go_tau_upper = float(vars["go tau upper"])
go_tau_start = float(vars["go tau start"])

stop_mu_lower = float(vars["stop mu lower"])
stop_mu_upper = float(vars["stop mu upper"])
stop_mu_start = float(vars["stop mu start"]) 

stop_sigma_lower = float(vars["stop sigma lower"])
stop_sigma_upper = float(vars["stop sigma upper"])
stop_sigma_start = float(vars["stop sigma start"])

stop_tau_lower = float(vars["stop tau lower"])
stop_tau_upper = float(vars["stop tau upper"])
stop_tau_start = float(vars["stop tau start"])

pf_lower = float(vars["stop pf lower"]) 
pf_upper = float(vars["stop pf upper"]) 
pf_mean = float(vars["stop pf mean"]) 
pf_sd = float(vars["stop pf sd"]) 
pf_start = float(vars["stop pf start"]) 

go_mu_sd_lower = float(vars["go mu sd lower"]) 
go_mu_sd_upper = float(vars["go mu sd upper"])
go_mu_sd_start = float(vars["go mu sd start"])

go_sigma_sd_lower = float(vars["go sigma sd lower"]) 
go_sigma_sd_upper = float(vars["go sigma sd upper"])
go_sigma_sd_start = float(vars["go sigma sd start"])

go_tau_sd_lower = float(vars["go tau sd lower"]) 
go_tau_sd_upper = float(vars["go tau sd upper"])
go_tau_sd_start = float(vars["go tau sd start"])

stop_mu_sd_lower = float(vars["stop mu sd lower"]) 
stop_mu_sd_upper = float(vars["stop mu sd upper"])
stop_mu_sd_start = float(vars["stop mu sd start"])

stop_sigma_sd_lower = float(vars["stop sigma sd lower"])
stop_sigma_sd_upper = float(vars["stop sigma sd upper"])
stop_sigma_sd_start = float(vars["stop sigma sd start"])

stop_tau_sd_lower = float(vars["stop tau sd lower"])
stop_tau_sd_upper = float(vars["stop tau sd upper"])
stop_tau_sd_start = float(vars["stop tau sd start"])

pf_sd_lower = float(vars["stop pf sd lower"])
pf_sd_upper = float(vars["stop pf sd upper"])
pf_sd_start = float(vars["stop pf sd start"])

iinteg_lower = int(vars["limits of integration lower"])
iinteg_upper = int(vars["limits of integration upper"])

def cython_Go(value, imu_go, isigma_go, itau_go):
    """Ex-Gaussian log-likelihood of GoRTs"""
    return stop_likelihoods_wtf.Go(value, imu_go, isigma_go, itau_go)

def cython_SRRT(value, issd, imu_go, isigma_go, itau_go, imu_stop, isigma_stop, itau_stop, ip_tf):
    """Censored ExGaussian log-likelihood of SRRTs"""
    return stop_likelihoods_wtf.SRRT(value, issd, imu_go, isigma_go, itau_go, imu_stop, isigma_stop, itau_stop, ip_tf)

def cython_Inhibitions(value, imu_go, isigma_go, itau_go, imu_stop, isigma_stop, itau_stop, ip_tf):
    """Censored ExGaussian log-likelihood of inhibitions"""
    return stop_likelihoods_wtf.Inhibitions(value, imu_go, isigma_go, itau_go, imu_stop, isigma_stop, itau_stop, ip_tf)

Go_like = pm.stochastic_from_dist(name="Ex-Gauss GoRT",
                                  logp=cython_Go,
                                  dtype=np.float,
                                  mv=False)

SRRT_like = pm.stochastic_from_dist(name="CensoredEx-Gauss SRRT",
                                    logp=cython_SRRT,
                                    dtype=np.float,
                                    mv=False)

Inhibitions_like = pm.stochastic_from_dist(name="CensoredEx-Gauss Inhibittions",
                                  logp=cython_Inhibitions,
                                  dtype=np.int32,
                                  mv=False)
class KnodeGo(Knode):
    def create_node(self, node_name, kwargs, data):
        new_data = data[data['ss_presented'] == 0]
        kwargs['value'] = new_data['rt']
        return self.pymc_node(name=node_name, **kwargs)

class KnodeSRRT(Knode):
    def create_node(self, node_name, kwargs, data):
        new_data = data[(data['ss_presented'] == 1) & (data['inhibited'] == 0)]
        kwargs['value'] = new_data['rt']
        kwargs['issd'] = np.array(new_data['ssd'], dtype=np.int32)
        return self.pymc_node(name=node_name, **kwargs)

class KnodeInhibitions(Knode):
    def create_node(self, node_name, kwargs, data):
        new_data = data[(data['ss_presented'] == 1) & (data['inhibited'] == 1)]
        uniq_ssds = np.unique(new_data['ssd'])
        ssd_inhib_trials = []
        for uniq_ssd in uniq_ssds:
            ssd_inhib_trials.append((uniq_ssd, len(new_data[new_data['ssd'] == uniq_ssd]), iinteg_lower, iinteg_upper))
        ssd_inhib_trials = np.array(ssd_inhib_trials, dtype=np.int32)
        #print(ssd_inhib_trials)        
        kwargs['value'] = ssd_inhib_trials
        return self.pymc_node(name=node_name, **kwargs)

class StopSignal(Hierarchical):
    def __init__(self, data, **kwargs):
        self.group_only_nodes = kwargs.pop('group_only_nodes', ())
        super(StopSignal, self).__init__(data, **kwargs)

    def create_ss_knode(self, knodes):
        ss_parents = OrderedDict()
        ss_parents['imu_go'] = knodes['mu_go_bottom']
        ss_parents['isigma_go'] = knodes['sigma_go_bottom']
        ss_parents['itau_go'] = knodes['tau_go_bottom']
        ss_parents['imu_stop'] = knodes['mu_stop_bottom']
        ss_parents['isigma_stop'] = knodes['sigma_stop_bottom']
        ss_parents['itau_stop'] = knodes['tau_stop_bottom']
        ss_parents['ip_tf'] = knodes['p_tf_bottom']

        go_like = KnodeGo(Go_like, 'go_like', col_name='rt', observed=True, imu_go=ss_parents['imu_go'], isigma_go=ss_parents['isigma_go'], itau_go=ss_parents['itau_go'])
        srrt_like = KnodeSRRT(SRRT_like, 'srrt_like', col_name='rt', observed=True, **ss_parents)
        inhibitions_like = KnodeInhibitions(Inhibitions_like, 'inhibitions_like', col_name='rt', observed=True, **ss_parents)

        return [go_like,srrt_like,inhibitions_like]

    def create_knodes(self):
        knodes = OrderedDict()
        knodes.update(self.create_family_trunc_normal('mu_go', lower=go_mu_lower, upper=go_mu_upper, value=go_mu_start,var_lower=go_mu_sd_lower, var_upper=go_mu_sd_upper, var_value=go_mu_sd_start))
        knodes.update(self.create_family_trunc_normal('sigma_go', lower=go_sigma_lower, upper=go_sigma_upper, value=go_sigma_start,var_lower=go_sigma_sd_lower, var_upper=go_sigma_sd_upper, var_value=go_sigma_sd_start))
        knodes.update(self.create_family_trunc_normal('tau_go', lower=go_tau_lower, upper=go_tau_upper, value=go_tau_start,var_lower=go_tau_sd_lower, var_upper=go_tau_sd_upper, var_value=go_tau_sd_start))
        knodes.update(self.create_family_trunc_normal('mu_stop', lower=stop_mu_lower, upper=stop_mu_upper, value=stop_mu_start,var_lower=stop_mu_sd_lower, var_upper=stop_mu_sd_upper, var_value=stop_mu_sd_start))
        knodes.update(self.create_family_trunc_normal('sigma_stop', lower=stop_sigma_lower, upper=stop_sigma_upper, value=stop_sigma_start,var_lower=stop_sigma_sd_lower, var_upper=stop_sigma_sd_upper, var_value=stop_sigma_sd_start))
        knodes.update(self.create_family_trunc_normal('tau_stop', lower=stop_tau_lower, upper=stop_tau_upper, value=stop_tau_start,var_lower=stop_tau_sd_lower, var_upper=stop_tau_sd_upper, var_value=stop_tau_sd_start))
        knodes.update(self.create_family_trunc_normal_probit('p_tf', lower_ind=pf_lower, upper_ind=pf_upper, mean_gr=pf_mean, var_gr=pf_sd, lower=pf_lower, upper=pf_upper, value=pf_start, var_lower=pf_sd_lower, var_upper=pf_sd_upper, var_value=pf_sd_start))
     
        likelihoods = self.create_ss_knode(knodes)
      
        return knodes.values() + likelihoods

    def create_family_trunc_normal_probit(self, name, value=None, lower_ind=None, upper_ind=None, mean_gr=None, var_gr=None, lower=None, upper=None, var_lower=1e-10,var_upper=1, var_value=.1):
        knodes = OrderedDict()
        probit = norm.cdf
        
        if self.is_group_model and name not in self.group_only_nodes:
            tau_gr = var_gr**-2
            g = Knode(pm.TruncatedNormal, '%s' % name, mu=mean_gr,tau=tau_gr, a=lower, b=upper, value=value, depends=self.depends[name])  
            var = Knode(pm.Uniform, '%s_var' % name, lower=var_lower,upper=var_upper, value=var_value)
            tau = Knode(pm.Deterministic, '%s_tau' % name,doc='%s_tau' % name, eval=lambda x: x**-2, x=var,plot=False, trace=False, hidden=True)
            subjpt = Knode(pm.TruncatedNormal, '%s_subjpt' % name, mu=g,tau=tau, a=lower, b=upper, value=value, depends=('subj_idx',), subj=True,plot=self.plot_subjs)         
            subj = Knode(pm.Deterministic, '%s_subj' % name, doc='%s_subj' % name, eval=lambda x: probit(x), x=subjpt,plot=False, trace=False, hidden=True)

            knodes['%s'%name] = g
            knodes['%s_var'%name] = var
            knodes['%s_tau'%name] = tau
            knodes['%s_subjpt'%name] = subjpt
            knodes['%s_bottom'%name] = subj

        else:
            subj = Knode(pm.Uniform, name, lower=lower_ind,upper=upper_ind, value=value,depends=self.depends[name])
            
            knodes['%s_bottom'%name] = subj

        return knodes