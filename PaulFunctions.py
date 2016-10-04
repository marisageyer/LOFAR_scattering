# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 17:20:15 2016

@author: marisa
"""
import numpy as np
import matplotlib.pyplot as plt
from math import exp,sqrt,pi,erf,log
from scipy.special import iv, ive
from lmfit import Model, conf_interval, printfuncs
from lmfit.models import GaussianModel



def read_Paul(filename):
    #function for reading /Paul/1937+21_mjd_57397.txt at the moment
    #Single column of ascii data
    pulsar_name = "B1937+21"
    datas = np.loadtxt(filename)
    nch = 1
    nbins = np.shape(datas)[0]
    mainpulse = datas[nbins/2:nbins]
    nbinsmain = nbins/2
    freq = 1.4 #(1.4GHz)
    return pulsar_name, nch, nbinsmain, mainpulse, freq
    
def makeprofile(nbins = 2**9, ncomps = 1, amps = 1, means = 100, sigmas = 10):
    if ncomps == 1:
        npamps = np.array([amps])
        npmeans = np.array([means])
        npsigmas = np.array([sigmas])
    else:    
        npamps = np.array(amps)
        npmeans = np.array(means)
        npsigmas = np.array(sigmas)
   
    profile = np.zeros(nbins)
    x = np.linspace(1,nbins,nbins)
    
    for i in range(ncomps):
#        print npmeans[i]
        profile = profile + \
        npamps[i]*np.exp(-pow(x-npmeans[i],2)/(2*pow(npsigmas[i],2)))
    return x, profile            

def pulsetrain_bins(npulses, numberofbins, profile):
    binsrange = np.linspace(1,numberofbins,numberofbins)    
    nbins = np.max(binsrange)
#    print nbins
    train = np.zeros(npulses*int(nbins))

    for i in range(npulses):
        startbin = i*nbins
        train[startbin:startbin + nbins] = profile
    return train
    
def broadfunc(x,tau):
    broadfunc = (1/tau)*np.exp(-x/tau)
    return broadfunc

def psrscatter(brfunc, profile):    
    scattered = np.convolve(profile,brfunc)
    profint = np.sum(profile)    
    scint = np.sum(scattered)
    scatterednorm = scattered / scint * profint
    bins = profile.shape[0]
    out = scatterednorm[0:bins]    
    return out

def extractpulse(train, pulsesfromend, binsperpulse):
    if pulsesfromend == 0:
        start = 0
        end = binsperpulse
        zerobpulse = train[start:end]-np.min(train[start:end])
        rectangle = np.min(train[start:end])*binsperpulse
        flux = np.sum(train[start:end]) - rectangle
        return train[start:end], zerobpulse, rectangle, flux
    
    else:     
        start = -pulsesfromend*binsperpulse
        end = start + binsperpulse 
        zerobpulse = train[start:end]-np.min(train[start:end])
        rectangle = np.min(train[start:end])*binsperpulse
        flux = np.sum(train[start:end]) - rectangle
        return train[start:end], zerobpulse, rectangle, flux







def GxETrainMulti(x, ncomps, mus1, mus2, sigmas1, sigmas2, As1, As2, tau, dc, nbins):
#    As = np.array(As)
#    mus = np.array(mus)
#    sigmas = np.array(sigmas)
#    As1,As2,mus1,mus2,sigmas1,sigmas2 = As[0],As[1],mus[0],mus[1],sigmas[0],sigmas[1]    
    bins, profile = makeprofile(nbins = nbins, ncomps = 2, amps = [As1,As2], means = [mus1,mus2], sigmas = [sigmas1,sigmas2])
    binstau = np.linspace(1,nbins,nbins)
#    binstau = x
    scat = psrscatter(broadfunc(binstau,tau),pulsetrain_bins(3, nbins, profile))   
    climb, observed_nonoise, rec, flux = extractpulse(scat, 2, nbins)
    return observed_nonoise + dc

   


def tau_fittermulti(data,nbins,guess_mean1, guess_mean2, guess_peak1, guess_peak2, maxmean):
#    binlen = len(data)
#    profile_peak1 = np.max(data[0:binlen/2])
#    profile_peak2 = np.max(data[binlen/2:binlen])
#    profile_peak1 = np.max(data[40:55])
#    profile_peak2 = np.max(data[55:60])
#    print profile_peak1
#    print profile_peak2
#    binpeak = np.argmax(data)
       
    modelname = GxETrainMulti
    model = Model(modelname)
            
    model.set_param_hint('nbins', value=nbins, vary=False)
    model.set_param_hint('ncomps', value=2, vary=False)       
    model.set_param_hint('sigmas1', value=15, vary=True, min =0, max = nbins/4)
    model.set_param_hint('sigmas2', value=15, vary=True, min =0, max = nbins/4)
    model.set_param_hint('mus1', value=guess_mean1, vary=True, min=660, max =720)
    model.set_param_hint('mus2', value=guess_mean2, vary=True, min=725, max = 750)
    model.set_param_hint('As1',value=guess_peak1, vary=True,min=0)
    model.set_param_hint('As2',value=guess_peak2, vary=True,min=0)
    model.set_param_hint('tau',value=10, vary=True, min=0)
    model.set_param_hint('dc',value = 0.1, vary = True)
#    val1, val2, val3, val4, val5, val6 = 20.25, 7.58, 689.5, 739.8, 1.05, 0.412  ##These are the mean values from the free fit
#    kt = 0.05
#    min1, min2, min3, min4, min5, min6 = (1-kt)*val1, (1-kt)*val2, (1-kt)*val3, (1-kt)*val4, (1-kt)*val5, (1-kt)*val6 
#    max1, max2, max3, max4, max5, max6 = (1+kt)*val1, (1+kt)*val2, (1+kt)*val3, (1+kt)*val4, (1+kt)*val5, (1+kt)*val6 
#    model.set_param_hint('nbins', value=nbins, vary=False)
#    model.set_param_hint('ncomps', value=2, vary=False)       
#    model.set_param_hint('sigmas1', value=val1*1.2, vary=True, min =min1, max = max1)
#    model.set_param_hint('sigmas2', value=val2*0.8, vary=True, min =min2, max = max2)
#    model.set_param_hint('mus1', value=val3*1.01, vary=True, min=min3, max =max3)
#    model.set_param_hint('mus2', value=val4*0.9, vary=True, min=min4, max = max4)
#    model.set_param_hint('As1',value=val5*0.6, vary=True,min=min5, max=max5)
#    model.set_param_hint('As2',value=val6*0.3, vary=True,min=min6, max=max6)
#    model.set_param_hint('tau',value=10, vary=True, min=0)
#    model.set_param_hint('dc',value = 0.1, vary = True)
    pars = model.make_params()
    xax=np.linspace(1,nbins,nbins)
    #"""Fit data"""
    result = model.fit(data,pars,x=xax)
    print(result.fit_report(show_correl = False))
    
    noiselessmodel = result.best_fit
    besttau = result.best_values['tau']
    taustd = result.params['tau'].stderr  ##estimated 1 sigma error

    bestsig1 = result.best_values['sigmas1']
    bestsig2 = result.best_values['sigmas2']
    bestmu1 = result.best_values['mus1']
    bestmu2 = result.best_values['mus2']
    bestA1 = result.best_values['As1']
    bestA2 = result.best_values['As2']
    
    bestdc = result.best_values['dc']    
    
    bestsig_std1 = result.params['sigmas1'].stderr
    bestsig_std2 = result.params['sigmas2'].stderr
    bestmu_std1 = result.params['mus1'].stderr
    bestmu_std2 = result.params['mus2'].stderr
    bestA_std1 = result.params['As1'].stderr
    bestA_std2 = result.params['As2'].stderr
    
    bestdc_std = result.params['dc'].stderr    
    
    bestparams = np.array([bestsig1,bestsig2,bestmu1,bestmu2, bestA1, bestA2,bestdc])
    bestparams_std = np.array([bestsig_std1,bestsig_std2, bestmu_std1,bestmu_std2, bestA_std1,bestA_std2,bestdc_std])
    
    rchi = result.redchi    
    
    return noiselessmodel, besttau, taustd, bestparams, bestparams_std, rchi
    
def TwoPeaksModel(x,nbins,A1,A2,mu1,mu2,sig1,sig2,dc):
    bins, twoprofile = makeprofile(nbins = nbins, ncomps = 2, amps = [A1,A2], means = [mu1,mu2], sigmas = [sig1,sig2])
    return twoprofile + dc
