# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 15:08:23 2016

@author: marisa
"""

import argparse
import os, sys
import pypsr_standalone as psr
import matplotlib.pyplot as plt
from lmfit import Model, conf_interval, printfuncs
from lmfit import minimize, Parameter, Parameters, fit_report
from lmfit.models import LinearModel, PowerLawModel, ExponentialModel, QuadraticModel
import numpy as np
from scipy import special
import DataReadIn as dri


"""This script is used to analyse the tau, alpha, DM, dist. etc of the pulsar set 
- and is therefore aimed at using data previous obtained from
1. Either using results from Standalone_Taufit_simuBF.py stored as txt files
2. Or having data as txt files from the literature """


"""Tau at 1 GHz vs DM plots"""


LTab_noIT = np.loadtxt('/Users/marisa/Dropbox/Aris/LOFARCommCensus/Lewandowski15/Lewand15Tau_noIT.txt',skiprows=1)
LTab = np.loadtxt('/Users/marisa/Dropbox/Aris/LOFARCommCensus/Lewandowski15/Lewand15Tau.txt')
LohTau01 = np.loadtxt('/Users/marisa/Dropbox/Aris/LOFARCommCensus/Lewandowski15/Lohmer01Tau.txt')
LohTau04 = np.loadtxt('/Users/marisa/Dropbox/Aris/LOFARCommCensus/Lewandowski15/Lohmer04Tau.txt')
GTau = np.loadtxt('/Users/marisa/Dropbox/Aris/LOFARCommCensus/Lewandowski15/GeyerTau.txt',skiprows=1)
L13_noIT = np.loadtxt(dri.skip_first_col('/Users/marisa/Dropbox/Aris/LOFARCommCensus/Lewandowski15/Lewan13Table_noIT.txt'),skiprows=1)
L13 = np.loadtxt(dri.skip_first_col('/Users/marisa/Dropbox/Aris/LOFARCommCensus/Lewandowski15/Lewan13Table.txt'),skiprows=1)

##Combine all the data sets to be able to do tau-fitting

## All the DMs:

#dm1 = LTab_noIT[:,0]
#Lewandowski has excluded more data points in their fits
dm1 = LTab_noIT[:,0]
dm2 = LohTau04[:,0]
dm3 = LohTau01[:,0]
dm4 = GTau[:,0]
dm5 = L13_noIT[:,0]
DMs = np.concatenate((dm1,dm2,dm3,dm4))


## All the Taus at 1 GHz (in ms)

#t1 = LTab_noIT[:,4]
#Lewandowski has excluded more data points in their fits
t1 = LTab_noIT[:,4]
t2 = LohTau04[:,4]
t3 = LohTau01[:,3]
t4 = GTau[:,4]*1000
t5 = L13_noIT[:,4]

Taus = np.concatenate((t1,t2,t3,t4))

DMTauGrid = np.column_stack((DMs,Taus))
G_dmtgrid = np.column_stack((dm4,t4))

# Sort the grid by DM:
DMTauGridorder = DMTauGrid[DMTauGrid[:,0].argsort()]
indDMlow = np.where(DMTauGridorder[:,0] < 300)
indtauhigh = np.where(DMTauGridorder[:,1] > 45)
indDMlow2 = np.where(DMTauGridorder[:,0] < 100)
indtauhigh2 = np.where(DMTauGridorder[:,1] > 8e-1)
indtoremove = np.intersect1d(indDMlow, indtauhigh)
indtoremove2 = np.intersect1d(indDMlow2, indtauhigh2)
rm_combo = np.concatenate((indtoremove,indtoremove2))
print rm_combo
DMTau = np.delete(DMTauGridorder,rm_combo,axis=0)
dm4order = G_dmtgrid[G_dmtgrid[:,0].argsort()]



## Fit with Ramachandran et al. 1997 model (power law)
Rfit, RBestP, Rchi = psr.Ramanch_fit(DMTau[:,0],DMTau[:,1])
G_Rfit, G_BestP, G_Rchi = psr.Ramanch_fit(dm4order[:,0],dm4order[:,1])
DMspace = np.linspace(10,1000,10)

## Fit with Bhat 2004 parabolic model
logDM = np.log10(DMTau[:,0])
logtau = np.log10(DMTau[:,1])
Bfit, BBestP, BRchi = psr.Bhat_fit(logDM,logtau)
LewParab_logtau = psr.Bhat_model(np.log10(DMspace),-6.344,1.467,0.509)
BhatParab_logtau = psr.Bhat_model(np.log10(DMspace),-6.46,0.154,1.07)

plt.figure()
plt.plot(LTab[:,0],LTab[:,4],'bo',alpha = 0.3,markersize=9,label='L15 Full set')
plt.plot(LTab_noIT[:,0],LTab_noIT[:,4],'r*',markersize=11,label='L15 Selected Set')
plt.plot(L13[:,0],L13[:,4],'co',alpha = 0.3,markersize=9,label='L13 Full set')
plt.plot(L13_noIT[:,0],L13_noIT[:,4],'r*',markersize=11,markerfacecolor='None',label='L13 Selected Set')
plt.plot(LohTau04[:,0],LohTau04[:,4],'k^',markersize=10,label='Lohmer 2004')
plt.plot(LohTau01[:,0],LohTau01[:,3],'k^',markersize=10,markerfacecolor='None',label='Lohmer 2001') 
plt.plot(GTau[:,0],GTau[:,4]*1000.,'go',markersize=10,label='LOFAR')
plt.plot(DMTau[:,0],Rfit,'m-',lw=1.5,label=r'$\gamma = %.3f$' % (RBestP[0]-2.2))
#plt.plot(dm4order[:,0],G_Rfit,'r--',label=r'$\gamma_{L} = %.3f$' % (G_BestP[0]-2.2))
plt.plot(DMTau[:,0],np.power(10.,Bfit),'b--',lw=1.5,label=r'$k1,k2,k3 = %.1f,%.1f,%.1f$' % (BBestP[0],BBestP[1],BBestP[2]))
plt.plot(DMspace,psr.Ramanch_model(DMspace,2.26e-7,2.26e-7*0.00205,1.74+2.2),'k-',alpha=0.5, label='L15 power law')
plt.plot(DMspace,10.**LewParab_logtau,'g--',alpha=0.7,label = 'L15 parabola')
plt.plot(DMspace,10.**BhatParab_logtau,'r:',alpha=0.7,label = 'Bhat04 parabola')
plt.plot(DMTau[:,0],DMTau[:,1],'yo',markersize=4.0)
plt.xscale('log')
plt.yscale('log')
plt.legend(loc = 'best',fontsize=10)

plt.rc('text', usetex = True)
plt.rc('font',family='serif')
plt.ylabel(r'$\tau$ (ms)',fontsize=16)
plt.xlabel(r'DM (pc/cm$^3$)',fontsize=16)
plt.xlim(5,2000)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)



plt.figure()
plt.plot(GTau[:,0],GTau[:,4]*1000.,'go',markersize=10,label='LOFAR')
plt.plot(dm4order[:,0],G_Rfit,'r--',label=r'$\gamma_{L} = %.3f$' % (G_BestP[0]-2.2))
plt.rc('text', usetex = True)
plt.rc('font',family='serif')
plt.xscale('log')
plt.yscale('log')
plt.legend(loc = 'best',fontsize=10)
plt.ylabel(r'$\tau$ (ms)',fontsize=16)
plt.xlabel(r'DM (pc/cm$^3$)',fontsize=16)
plt.xlim(5,2000)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)




"""Alpha vs DM plots"""

Lohm01 = np.loadtxt(dri.skip_first_col('/Users/marisa/Dropbox/Aris/LOFARCommCensus/Lewandowski15/Lohmer2001DMAlphaDist.txt'),skiprows=1)
Lohm04 = np.loadtxt('/Users/marisa/Dropbox/Aris/LOFARCommCensus/Lewandowski15/LohmerDMDistAlpha.txt',skiprows=1)
GTab = np.loadtxt('/Users/marisa/Dropbox/Aris/LOFARCommCensus/DMDistAlphaGeyer.txt', skiprows = 1)
John98 = np.loadtxt(dri.skip_first_col('/Users/marisa/Dropbox/Aris/LOFARCommCensus/Lewandowski15/Johnston98.txt'), skiprows = 1) 


gdm = GTab[:,0]
galf = GTab[:,2]
galferr = GTab[:,3]

plt.figure()

plt.errorbar(Lohm01[:,0],Lohm01[:,1], yerr=Lohm01[:,2], fmt = 'k^',markersize=10, markerfacecolor='None',alpha=0.8,label='Lohmer 2001')

plt.errorbar(Lohm04[:,0],Lohm04[:,2], yerr=Lohm04[:,3], fmt = 'k^',markersize=10, alpha=0.8,label='Lohmer 2004')

plt.errorbar(John98[:,0],John98[:,1], yerr=John98[:,2], fmt = 'k>',markersize=10, alpha=0.8,label='Jonhston 1998 et al.')

plt.errorbar(gdm,galf, yerr=galferr, fmt = 'go',markersize=10, alpha=0.8,label='LOFAR')

plt.errorbar(LTab[:,0], LTab[:,2], yerr=LTab[:,3],  fmt = 'bo',markersize =10,alpha=0.2,label = 'L15 Full set')

plt.errorbar(L13[:,0], L13[:,2], yerr=L13[:,3],  fmt = 'co',markersize =10,alpha=0.2,label = 'L13 Full set')

plt.errorbar(LTab_noIT[:,0], LTab_noIT[:,2], yerr=LTab_noIT[:,3], fmt = 'r*',markersize=14, alpha=0.8, label = 'L15 Selected set')

plt.errorbar(L13_noIT[:,0], L13_noIT[:,2], yerr=L13_noIT[:,3], fmt = 'k*',markersize=14,markerfacecolor='None', alpha=0.8, label = 'L13 Selected set')

plt.axhline(y=4.0,color='k',linestyle=':')
plt.axhline(y=4.4,color='k',linestyle='dashed')

plt.xscale('log')
plt.xlabel(r'DM (pc/cm$^3$)',fontsize=16)
plt.ylabel(r'$\alpha$',fontsize=16)


plt.text(240,2,"B1920+21")
plt.text(180,2.5,"B2255+58")

plt.plot(gdm[8],galf[8],'k+',markersize=16,lw=2.0)
plt.plot(gdm[12],galf[12],'k+',markersize=16, lw=2.0)
plt.plot(gdm[13],galf[13],'k+',markersize=16,lw=2.0)
plt.plot(gdm[14],galf[14],'k+',markersize=16,lw=2.0)
plt.plot(gdm[8],galf[8],'k+',markersize=16,lw=2.0)
plt.plot(gdm[12],galf[12],'k+',markersize=16, lw=2.0)
plt.plot(gdm[13],galf[13],'k+',markersize=16,lw=2.0)
plt.plot(gdm[14],galf[14],'k+',markersize=16,lw=2.0)
plt.plot(gdm[8],galf[8],'k+',markersize=16,lw=2.0)
plt.plot(gdm[12],galf[12],'k+',markersize=16, lw=2.0)
plt.plot(gdm[13],galf[13],'k+',markersize=16,lw=2.0)
plt.plot(gdm[14],galf[14],'k+',markersize=16,lw=2.0)

plt.legend(fontsize=10,loc='best')

"""Alpha vs Dist plots"""

gdist = GTab[:,1]

plt.figure()
plt.errorbar(gdist,galf, yerr=galferr, fmt = 'go',markersize=10, alpha=0.8,label='LOFAR')
plt.plot(gdist[-4:-1],galf[-4:-1],'y*',markersize=8.0)
plt.plot(gdist[-8],galf[-8],'y*',markersize=8.0)
#plt.xscale('log')
plt.xlabel(r'Dist (kpc)',fontsize=16)
plt.ylabel(r'$\alpha$',fontsize=16)


""" Alpha vs Tau (at 1Ghz)"""

plt.figure()
plt.errorbar(GTau[:,4],GTau[:,2], GTau[:,3],fmt = 'go',markersize=10, alpha=0.8,label='LOFAR')
plt.xlabel('tau')
plt.ylabel('alpha')

