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
import pandas as pd
#import ipdb
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
path__analysisDescription = sys.argv[1] + 'analysis.txt'

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
mu_go_lower = float(vars["mu go lower"])
mu_go_upper = float(vars["mu go upper"])
mu_go_start = float(vars["mu go start"])

sigma_go_lower = float(vars["sigma go lower"])
sigma_go_upper = float(vars["sigma go upper"])
sigma_go_start = float(vars["sigma go start"])

tau_go_lower = float(vars["tau go lower"])
tau_go_upper = float(vars["tau go upper"])
tau_go_start = float(vars["tau go start"])

shift_go_lower = float(vars["shift go lower"])
shift_go_upper = float(vars["shift go upper"])
shift_go_start = float(vars["shift go start"])

mu_stop_lower = float(vars["mu stop lower"])
mu_stop_upper = float(vars["mu stop upper"])
mu_stop_start = float(vars["mu stop start"])

sigma_stop_lower = float(vars["sigma stop lower"])
sigma_stop_upper = float(vars["sigma stop upper"])
sigma_stop_start = float(vars["sigma stop start"])

tau_stop_lower = float(vars["tau stop lower"])
tau_stop_upper = float(vars["tau stop upper"])
tau_stop_start = float(vars["tau stop start"])

shift_stop_lower = float(vars["shift stop lower"])
shift_stop_upper = float(vars["shift stop upper"])
shift_stop_start = float(vars["shift stop start"])

pf_stop_lower = float(vars["pf stop lower"])
pf_stop_upper = float(vars["pf stop upper"])
pf_stop_mean  = float(vars["pf stop mean"])
pf_stop_sd    = float(vars["pf stop sd"])
pf_stop_start = float(vars["pf stop start"])

mu_go_sd_lower = float(vars["mu go sd lower"])
mu_go_sd_upper = float(vars["mu go sd upper"])
mu_go_sd_start = float(vars["mu go sd start"])

sigma_go_sd_lower = float(vars["sigma go sd lower"])
sigma_go_sd_upper = float(vars["sigma go sd upper"])
sigma_go_sd_start = float(vars["sigma go sd start"])

tau_go_sd_lower = float(vars["tau go sd lower"])
tau_go_sd_upper = float(vars["tau go sd upper"])
tau_go_sd_start = float(vars["tau go sd start"])

shift_go_sd_lower = float(vars["shift go sd lower"])
shift_go_sd_upper = float(vars["shift go sd upper"])
shift_go_sd_start = float(vars["shift go sd start"])

mu_stop_sd_lower = float(vars["mu stop sd lower"])
mu_stop_sd_upper = float(vars["mu stop sd upper"])
mu_stop_sd_start = float(vars["mu stop sd start"])

sigma_stop_sd_lower = float(vars["sigma stop sd lower"])
sigma_stop_sd_upper = float(vars["sigma stop sd upper"])
sigma_stop_sd_start = float(vars["sigma stop sd start"])

tau_stop_sd_lower = float(vars["tau stop sd lower"])
tau_stop_sd_upper = float(vars["tau stop sd upper"])
tau_stop_sd_start = float(vars["tau stop sd start"])

shift_stop_sd_lower = float(vars["shift stop sd lower"])
shift_stop_sd_upper = float(vars["shift stop sd upper"])
shift_stop_sd_start = float(vars["shift stop sd start"])

pf_stop_sd_lower = float(vars["pf stop sd lower"])
pf_stop_sd_upper = float(vars["pf stop sd upper"])
pf_stop_sd_start = float(vars["pf stop sd start"])

iinteg_lower = int(vars["limits of integration lower"])
iinteg_upper = int(vars["limits of integration upper"])

def cython_Go(value, ips, imu_go, isigma_go, itau_go, ishift_go):
    """Ex-Gaussian log-likelihood of GoRTs"""
    return stop_likelihoods_wtf.Go(value, ips, imu_go, isigma_go, itau_go, ishift_go)

def cython_SRRT(value, ips, issd, imu_go, isigma_go, itau_go, ishift_go, imu_stop, isigma_stop, itau_stop, ishift_stop, ipf_stop):
    """Censored ExGaussian log-likelihood of SRRTs"""
    return stop_likelihoods_wtf.SRRT(value, ips, issd, imu_go, isigma_go, itau_go, ishift_go, imu_stop, isigma_stop, itau_stop, ishift_stop, ipf_stop)

def cython_Inhibitions(value, ips, imu_go, isigma_go, itau_go, ishift_go, imu_stop, isigma_stop, itau_stop, ishift_stop, ipf_stop):
    """Censored ExGaussian log-likelihood of inhibitions"""
    return stop_likelihoods_wtf.Inhibitions(value, ips, imu_go, isigma_go, itau_go, ishift_go, imu_stop, isigma_stop, itau_stop, ishift_stop, ipf_stop)

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
        msk_ss          = data['ss_presented'] == 0
        kwargs['value'] = data[msk_ss]['rt']

        ps_bools      = abs(np.roll(data['ss_presented'],1) - 1)
        ps_bools[0]   = 0.
        ps_bools      = np.array(ps_bools, dtype=np.float)
        kwargs['ips'] = ps_bools[np.where(msk_ss)]

        #print('KnodeGo:')
        #print(kwargs['ips'].shape  , type(kwargs['ips']))
        #print(kwargs['value'].shape, type(kwargs['value']))
        return self.pymc_node(name=node_name, **kwargs)

class KnodeSRRT(Knode):
    def create_node(self, node_name, kwargs, data):
        msk_ss          = data['ss_presented'] == 1
        msk_respond     = data['inhibited'] == 0
        relevant_data   = data[(msk_ss) & (msk_respond)]

        kwargs['value'] = relevant_data['rt']
        kwargs['issd']  = np.array(relevant_data['ssd'], dtype=np.int32)

        ps_bools        = np.roll(data['ss_presented'],1)
        ps_bools[0]     = 0.
        kwargs['ips']   = np.array(ps_bools[np.where((msk_ss) & (msk_respond))], dtype=np.float)

        #print('KnodeSRRT:')
        #print(kwargs['ips'].shape  , type(kwargs['ips']))
        #print(kwargs['value'].shape, type(kwargs['value']))
        #print(kwargs['issd'].shape , type(['issd']))
        return self.pymc_node(name=node_name, **kwargs)

class KnodeInhibitions(Knode):
    def create_node(self, node_name, kwargs, data):
        msk_ss          = data['ss_presented'] == 1
        msk_respond     = data['inhibited'] == 1
        relevant_data   = data[(msk_ss) & (msk_respond)]

        ps_bools    = np.roll(data['ss_presented'],1)
        ps_bools[0] = 0.
        ps_bools    = np.array(ps_bools[np.where((msk_ss) & (msk_respond))], dtype=np.float)

        uniq_ssds = np.unique(relevant_data['ssd'])
        ssd_inhib_trials = []
        ips = []
        for uniq_ssd in uniq_ssds:
            this_ssd_msk = relevant_data['ssd'] == uniq_ssd

            n_not_ps = len(relevant_data[this_ssd_msk & (ps_bools == 0)])
            n_ps     = len(relevant_data[this_ssd_msk & (ps_bools == 1)])

            if n_not_ps > 0:
                ssd_inhib_trials.append((uniq_ssd, n_not_ps))
                ips.append(0)

            if n_ps > 0:
                ssd_inhib_trials.append((uniq_ssd, n_ps    ))
                ips.append(1)

        ssd_inhib_trials = np.array(ssd_inhib_trials, dtype=np.int32)
        ips = np.array(ips, dtype=np.float)

        kwargs['value'] = ssd_inhib_trials
        kwargs['ips']   = ips
        #kwargs['ips']   = np.array([ssd_inhib_trials.shape[0],1], dtype=np.float)
        #kwargs['ips'][1:len(uniq_ssds*2):2] = 1

        #kwargs['ips']   = ps_bools
        #kwargs['value'] = np.array(relevant_data['ssd'], dtype=np.int32)

        #print('KnodeInhibitions:')
        #print('ips:'  , kwargs['ips'].shape  )
        #print('value:', kwargs['value'].shape)
        #print(kwargs['value'])
        #print(kwargs['ips'])
        return self.pymc_node(name=node_name, **kwargs)

class StopSignal(Hierarchical):
    def __init__(self, data, **kwargs):
        self.group_only_nodes = kwargs.pop('group_only_nodes', ())
        super(StopSignal, self).__init__(data, **kwargs)

    def create_ss_knode(self, knodes):
        ss_parents = OrderedDict()
        ss_parents['imu_go']      = knodes['mu_go_bottom']
        ss_parents['isigma_go']   = knodes['sigma_go_bottom']
        ss_parents['itau_go']     = knodes['tau_go_bottom']
        ss_parents['ishift_go']   = knodes['shift_go_bottom']
        ss_parents['imu_stop']    = knodes['mu_stop_bottom']
        ss_parents['isigma_stop'] = knodes['sigma_stop_bottom']
        ss_parents['itau_stop']   = knodes['tau_stop_bottom']
        ss_parents['ishift_stop'] = knodes['shift_stop_bottom']
        ss_parents['ipf_stop']    = knodes['pf_stop_bottom']

        go_like = KnodeGo(Go_like, 'go_like', col_name='rt', observed=True, imu_go=ss_parents['imu_go'], isigma_go=ss_parents['isigma_go'], itau_go=ss_parents['itau_go'], ishift_go=ss_parents['ishift_go'])
        srrt_like = KnodeSRRT(SRRT_like, 'srrt_like', col_name='rt', observed=True, **ss_parents)
        inhibitions_like = KnodeInhibitions(Inhibitions_like, 'inhibitions_like', col_name='rt', observed=True, **ss_parents)

        return [go_like,srrt_like,inhibitions_like]

    def create_knodes(self):
        knodes = OrderedDict()
        knodes.update(self.create_family_trunc_normal('mu_go', lower=mu_go_lower, upper=mu_go_upper, value=mu_go_start,var_lower=mu_go_sd_lower, var_upper=mu_go_sd_upper, var_value=mu_go_sd_start))
        knodes.update(self.create_family_trunc_normal('sigma_go', lower=sigma_go_lower, upper=sigma_go_upper, value=sigma_go_start,var_lower=sigma_go_sd_lower, var_upper=sigma_go_sd_upper, var_value=sigma_go_sd_start))
        knodes.update(self.create_family_trunc_normal('tau_go', lower=tau_go_lower, upper=tau_go_upper, value=tau_go_start,var_lower=tau_go_sd_lower, var_upper=tau_go_sd_upper, var_value=tau_go_sd_start))
        knodes.update(self.create_family_trunc_normal('shift_go', lower=shift_go_lower, upper=shift_go_upper, value=shift_go_start,var_lower=shift_go_sd_lower, var_upper=shift_go_sd_upper, var_value=shift_go_sd_start))
        knodes.update(self.create_family_trunc_normal('mu_stop', lower=mu_stop_lower, upper=mu_stop_upper, value=mu_stop_start,var_lower=mu_stop_sd_lower, var_upper=mu_stop_sd_upper, var_value=mu_stop_sd_start))
        knodes.update(self.create_family_trunc_normal('sigma_stop', lower=sigma_stop_lower, upper=sigma_stop_upper, value=sigma_stop_start,var_lower=sigma_stop_sd_lower, var_upper=sigma_stop_sd_upper, var_value=sigma_stop_sd_start))
        knodes.update(self.create_family_trunc_normal('tau_stop', lower=tau_stop_lower, upper=tau_stop_upper, value=tau_stop_start,var_lower=tau_stop_sd_lower, var_upper=tau_stop_sd_upper, var_value=tau_stop_sd_start))
        knodes.update(self.create_family_trunc_normal('shift_stop', lower=shift_stop_lower, upper=shift_stop_upper, value=shift_stop_start,var_lower=shift_stop_sd_lower, var_upper=shift_stop_sd_upper, var_value=shift_stop_sd_start))
        knodes.update(self.create_family_trunc_normal_probit('pf_stop', lower_ind=pf_stop_lower, upper_ind=pf_stop_upper, mean_gr=pf_stop_mean, var_gr=pf_stop_sd, lower=pf_stop_lower, upper=pf_stop_upper, value=pf_stop_start, var_lower=pf_stop_sd_lower, var_upper=pf_stop_sd_upper, var_value=pf_stop_sd_start))

        likelihoods = self.create_ss_knode(knodes)

        return knodes.values() + likelihoods

    def create_family_trunc_normal_probit(self, name, value=None, lower_ind=None, upper_ind=None, mean_gr=None, var_gr=None, lower=None, upper=None, var_lower=1e-10,var_upper=1, var_value=.1):
        knodes = OrderedDict()
        probit = norm.cdf

        if self.is_group_model and name not in self.group_only_nodes:
            tau_gr = var_gr**-2
            g      = Knode(pm.TruncatedNormal, '%s'        % name, mu=mean_gr,tau=tau_gr, a=lower, b=upper, value=value, depends=self.depends[name])
            var    = Knode(pm.Uniform        , '%s_var'    % name, lower=var_lower,upper=var_upper, value=var_value)
            tau    = Knode(pm.Deterministic  , '%s_tau'    % name, doc='%s_tau' % name, eval=lambda x: x**-2, x=var,plot=False, trace=False, hidden=True)
            subjpt = Knode(pm.TruncatedNormal, '%s_subjpt' % name, mu=g,tau=tau, a=lower, b=upper, value=value, depends=('subj_idx',), subj=True,plot=self.plot_subjs)
            subj   = Knode(pm.Deterministic  , '%s_subj'   % name, doc='%s_subj' % name, eval=lambda x: probit(x), x=subjpt,plot=False, trace=False, hidden=True)

            knodes['%s'%name]        = g
            knodes['%s_var'%name]    = var
            knodes['%s_tau'%name]    = tau
            knodes['%s_subjpt'%name] = subjpt
            knodes['%s_bottom'%name] = subj

        else:
            subj = Knode(pm.Uniform, name, lower=lower_ind,upper=upper_ind, value=value,depends=self.depends[name])

            knodes['%s_bottom'%name] = subj

        return knodes
