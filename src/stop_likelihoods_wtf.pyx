#cython: embedsignature=True
#cython: cdivision=True
#cython: wraparound=False
#cython: boundscheck=False

from copy import copy
import numpy as np
cimport numpy as np

from scipy.stats import norm
from scipy.integrate import quad

cimport cython
from cython.parallel import *

from libc.math cimport log, exp, sqrt, pow

cdef extern from "math.h":
    double INFINITY, NAN

from cython_gsl cimport *

ctypedef double * double_ptr
ctypedef void * void_ptr

cdef extern from "stdlib.h":
    void free(void* ptr)
    void* malloc(size_t size)
    void* realloc(void* ptr, size_t size)


cdef inline double ExGauss_pdf(double value, double ips, double mu, double sigma, double tau, double shift) nogil:
    """ ExGaussian log pdf"""
    cdef double z
    
    # Standard deviation or mu shift on the basis of valence. Should not be both.
    sigma = sigma + shift*ips
    #mu = mu + shift*ips

    if tau > 0.05*sigma:
        z = value - mu - ((sigma**2)/tau)
        return -log(tau) - (z+(sigma**2/(2*tau)))/tau + log( gsl_cdf_gaussian_P(z/sigma, 1) )
    else:
        return log(gsl_ran_gaussian_pdf(value-mu, sigma))

cdef inline double ExGauss_cdf(double value, double ips, double mu, double sigma, double tau, double shift) nogil:
    """ ExGaussian log cdf upper tail"""
    cdef double z

    # Standard deviation or mu shift on the basis of valence. Should not be both.
    sigma = sigma + shift*ips
    #mu = mu + shift*ips

    if tau > 0.05*sigma:
        z = value - mu - ((sigma**2)/tau)
        return log(1-(gsl_cdf_gaussian_P((value-mu)/sigma,1)-gsl_cdf_gaussian_P(z/sigma,1)*exp((((mu+((sigma**2)/tau))**2)-
        (mu**2)-2*value *((sigma**2)/tau))/(2*(sigma**2)))))
    else:
        return log((1-(gsl_cdf_gaussian_P(((value-mu)/sigma), 1))))

def Go(np.ndarray[double, ndim=1] value, np.ndarray[double, ndim=1] ips, double imu_go, double isigma_go, double itau_go, double ishift_go):
    """Ex-Gaussian log-likelihood of GoRTs"""
    assert imu_go > 0
    assert isigma_go > 0
    assert itau_go > 0
    assert np.all(value != -999.), 'RT values on go trails cannot be -999.\n' + \
      'That is, you must remove "misses" from the data.\n GoRT values:\n' + \
      np.array_str(value) +'\n Valid:\n' + np.array_str(value != -999.)

    cdef Py_ssize_t size = value.shape[0]
    cdef Py_ssize_t i
    cdef double p
    cdef double sum_logp = 0

    for i in range(size):
        p = ExGauss_pdf(value[i], ips[i], imu_go, isigma_go, itau_go, ishift_go)

        if np.isinf(p) or np.isnan(p):
            return -np.inf

        sum_logp += p

    return sum_logp

def SRRT(np.ndarray[double, ndim=1] value, np.ndarray[double, ndim=1] ips, np.ndarray[int, ndim=1] issd, double imu_go, double isigma_go, double itau_go, double ishift_go, double imu_stop, double isigma_stop, double itau_stop, double ishift_stop, double ipf_stop):
    """Censored ExGaussian log-likelihood of SRRTs"""
    assert imu_go > 0
    assert isigma_go > 0
    assert itau_go > 0
    assert imu_stop > 0
    assert isigma_stop > 0
    assert itau_stop > 0
    assert ipf_stop >= 0
    assert ipf_stop <= 1

    assert np.all(value != -999.)

    cdef Py_ssize_t size = value.shape[0]
    cdef Py_ssize_t i
    cdef double p
    cdef double p1
    cdef double p2
    cdef double sum_logp = 0

    for i in range(size):
        p1= exp(ExGauss_pdf(value[i], ips[i], imu_go, isigma_go, itau_go, ishift_go))*ipf_stop
        p2= exp(ExGauss_pdf(value[i], ips[i], imu_go, isigma_go, itau_go, ishift_go))*exp(ExGauss_cdf(value[i]-issd[i], ips[i], imu_stop, isigma_stop, itau_stop, ishift_stop))*(1-ipf_stop)
        p=log(p1+p2)
        if np.isinf(p) or np.isnan(p):
            return -np.inf

        sum_logp += p

    return sum_logp

def Inhibitions(np.ndarray[int, ndim=2] value, np.ndarray[double, ndim=1] ips, double imu_go, double isigma_go, double itau_go, double ishift_go, double imu_stop, double isigma_stop, double itau_stop, double ishift_stop, double ipf_stop):
    """Censored ExGaussian log-likelihood of inhibitions"""
    assert imu_go > 0
    assert isigma_go > 0
    assert itau_go > 0
    assert imu_stop > 0
    assert isigma_stop > 0
    assert itau_stop > 0
    assert ipf_stop >= 0
    assert ipf_stop <= 1

    cdef Py_ssize_t size = value.shape[0]
    cdef Py_ssize_t i #, unique_ssd
    cdef double p

    cdef double sum_logp = 0
    cdef double p_ssd

    cdef int ssd, n_trials
    cdef int iinteg_lower, iinteg_upper

    #cdef np.ndarray[int, ndim=1] unique_ssds
    #cdef np.ndarray[int, ndim=1] this_ssd_msk
    #cdef np.ndarray[int, ndim=2] ssd_inhib_trials #= np.empty([,], dtype=np.int32)

    #uniq_ssds = np.unique(value)
    #ssd_inhib_trials = np.empty([0,3], dtype=np.int32)
    #for uniq_ssd in uniq_ssds:
    #    this_ssd_msk = (value == uniq_ssd)
    #    ssd_inhib_trials.append((uniq_ssd, len(value[this_ssd_msk & (ips == 0)]), 0))
    #    ssd_inhib_trials.append((uniq_ssd, len(value[this_ssd_msk & (ips == 1)]), 1))

    #ssd_inhib_trials = np.array(ssd_inhib_trials, dtype=np.int32)

    size = value.shape[0]
    for i in range(size):
        ssd          = value[i, 0]
        n_trials     = value[i, 1]
        #iinteg_lower = value[i, 2]
        #iinteg_upper = value[i, 3]
        #ps = ssd_inhib_trials[i,2]

        #ssd = value[i]
        ps  = ips[i]

        iinteg_lower = 0
        iinteg_upper = 10000

        assert (ssd != -999)
        #assert (iinteg_lower > 0)

        # Compute probability of single SSD
        p_ssd = log(integrate_cexgauss(iinteg_lower, iinteg_upper, imu_go, isigma_go, itau_go, ishift_go, imu_stop, isigma_stop, itau_stop, ishift_stop, ipf_stop, ps, ssd))
        if np.isinf(p_ssd) or np.isnan(p_ssd):
            return -np.inf
        # Multiply with number of trials and add to overall p
        sum_logp += n_trials * p_ssd

    return sum_logp

########################
# Integration routines #
########################
cdef double eval_cexgauss(double x, void * params) nogil:
    cdef double imu_go, isigma_go, itau_go, ishift_go, imu_stop, isigma_stop, itau_stop, ishift_stop, ipf_stop, ips, issd
    cdef double p
    imu_go = (<double_ptr> params)[0]
    isigma_go = (<double_ptr> params)[1]
    itau_go = (<double_ptr> params)[2]
    ishift_go = (<double_ptr> params)[3]
    imu_stop = (<double_ptr> params)[4]
    isigma_stop = (<double_ptr> params)[5]
    itau_stop = (<double_ptr> params)[6]
    ishift_stop = (<double_ptr> params)[7]
    ipf_stop = (<double_ptr> params)[8]
    ips = (<double_ptr> params)[9]
    issd = (<double_ptr> params)[10]

    p = exp(ExGauss_cdf(x, ips, imu_go, isigma_go, itau_go, ishift_go)) * exp(ExGauss_pdf(x-issd, ips, imu_stop, isigma_stop, itau_stop, ishift_stop)) * (1-ipf_stop)

    return p

cdef double integrate_cexgauss(double lower, double upper, double imu_go, double isigma_go, double itau_go, double ishift_go, double imu_stop, double isigma_stop, double itau_stop, double ishift_stop, double ipf_stop, double ips, int issd):

    cdef double alpha, result, error, expected
    cdef gsl_function F
    cdef double params[11]
    cdef size_t neval
    cdef gsl_integration_workspace * W
    W = gsl_integration_workspace_alloc(5000)

    params[0] = imu_go
    params[1] = isigma_go
    params[2] = itau_go
    params[3] = ishift_go
    params[4] = imu_stop
    params[5] = isigma_stop
    params[6] = itau_stop
    params[7] = ishift_stop
    params[8] = ipf_stop
    params[9] = ips
    params[10] = issd

    F.function = &eval_cexgauss
    F.params = params

    gsl_integration_qag(&F, lower, upper, 1e-4, 1e-4, 5000, GSL_INTEG_GAUSS41, W, &result, &error)
    gsl_integration_workspace_free(W)

    return result
