# -*- coding: utf-8 -*-
"""
Created on Tue Sep  6 13:21:59 2016

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


pathd = "/home/geyer"

"""Path to the freq, tau, tauerr txtfiles"""
txtpath_iso = os.path.join(pathd, "Dropbox/Aris/LOFARCommCensus/Lewandowski15/FreqTau")
txtpath_1D = "/Users/marisa/python-workingdir/Calibrated_ArticleData/Lists/FreqTaulists/FreqTau_onedim"

freqtaus = []
pulsars = []

for filename in os.listdir(txtpath_iso):
    pulsar = filename[0:10]
    freqtaufile = np.loadtxt(os.path.join(txtpath_iso,filename))
    freqtaus.append(freqtaufile)
    pulsars.append(pulsar)

freqtaus_1D = []
pulsars_1D = [] 

for filename in os.listdir(txtpath_1D):
    pulsar_1D = filename[0:10]
    freqtaufile_1D = np.loadtxt(os.path.join(txtpath_1D,filename))
    freqtaus_1D.append(freqtaufile_1D)
    pulsars_1D.append(pulsar_1D)



"""Do tau fits with power law model"""

powmod = PowerLawModel()
alphas = []
alphaerrs = []
fits = []

alphas1D = []
alphaerrs1D = []
fits1D = []

for i in range(len(pulsars)): 
    powpars = powmod.guess(freqtaus[i][:,1], x=freqtaus[i][:,0])
    powout = powmod.fit(freqtaus[i][:,1], powpars, x=freqtaus[i][:,0], weights=1/((freqtaus[i][:,2])**2))
    fit = powout.best_fit    
    alpha = -powout.best_values['exponent']
    alphaerr = powout.params['exponent'].stderr
    alphas.append(alpha)
    alphaerrs.append(alphaerr)
    fits.append(fit)
    
    powpars1D = powmod.guess(freqtaus_1D[i][:,1], x=freqtaus_1D[i][:,0])
    powout1D = powmod.fit(freqtaus_1D[i][:,1], powpars1D, x=freqtaus_1D[i][:,0], weights=1/((freqtaus_1D[i][:,2])**2))
    fit1D = powout1D.best_fit    
    alpha1D = -powout1D.best_values['exponent']
    alphaerr1D = powout1D.params['exponent'].stderr
    alphas1D.append(alpha1D)
    alphaerrs1D.append(alphaerr1D)
    fits1D.append(fit1D)


for i in range(len(pulsars)):
    if pulsars[i] != pulsars_1D[i]:
        print "Iso and 1D pulsar arrays do not agree"
 

#fig = plt.figure(1,figsize=(10,8))
#for i in range(4):
#    plt.subplot(2,2,i+1)
#    fig.subplots_adjust(left = 0.08, right = 0.97, wspace=0.2, hspace = 0.3, bottom=0.1)
#    plt.rc('text', usetex=True)
#    plt.rc('font', family='serif')
#    plt.errorbar(1000*freqtaus[i][:,0],freqtaus[i][:,1],yerr=freqtaus[i][:,2],fmt='k*',markersize=9.0,capthick=2,linewidth=1.5,label=r'$\alpha = %.2f \pm %.2f$' %(alphas[i],alphaerrs[i]))
#    plt.errorbar(1000*freqtaus_1D[i][:,0],freqtaus_1D[i][:,1],yerr=freqtaus_1D[i][:,2],fmt='k^',markersize=9.0,alpha = 0.5, capthick=2,linewidth=1.5,label=r'$\alpha = %.2f \pm %.2f$' %(alphas1D[i],alphaerrs1D[i]))
#    plt.plot(1000*freqtaus[i][:,0],fits[i],'k--',linewidth=1.5)
#    plt.plot(1000*freqtaus_1D[i][:,0],fits1D[i],'k--',linewidth=1.5)
#    plt.title(pulsars[i])
#    plt.yticks(fontsize=12)
#    ticksMHz = (1000*freqtaus[i][:,0]).astype(np.int)[0:len(1000*freqtaus[i][:,0]):2] 
#    plt.legend(fontsize=14)
#    plt.xscale('log')
#    plt.yscale('log')
#    plt.xlim(xmin=(freqtaus[i][0,0])*950,xmax=(freqtaus[i][-1,0])*1050)
#    plt.xticks(ticksMHz,ticksMHz,fontsize=12)
#
#for i in range(2):
#    plt.subplot(2,2,i+3)
#    plt.xlabel(r'$\nu$ (MHz)',fontsize=16, labelpad=15.0)
#    plt.subplot(2,2,2*i+1)
#    plt.ylabel(r'$\tau$ (sec)',fontsize=16)
    

#fig = plt.figure(2,figsize=(10,8))
#for i in range(4,8):
#    plt.subplot(2,2,i-3)
#    fig.subplots_adjust(left = 0.08, right = 0.97, wspace=0.2, hspace = 0.3, bottom=0.1)
#    plt.rc('text', usetex=True)
#    plt.rc('font', family='serif')
#    plt.errorbar(1000*freqtaus[i][:,0],freqtaus[i][:,1],yerr=freqtaus[i][:,2],fmt='k*',markersize=9.0,capthick=2,linewidth=1.5,label=r'$\alpha = %.2f \pm %.2f$' %(alphas[i],alphaerrs[i]))
#    plt.errorbar(1000*freqtaus_1D[i][:,0],freqtaus_1D[i][:,1],yerr=freqtaus_1D[i][:,2],fmt='k^',markersize=9.0,alpha = 0.5, capthick=2,linewidth=1.5,label=r'$\alpha = %.2f \pm %.2f$' %(alphas1D[i],alphaerrs1D[i]))    
#    plt.plot(1000*freqtaus[i][:,0],fits[i],'k--',linewidth=1.5)
#    plt.plot(1000*freqtaus_1D[i][:,0],fits1D[i],'k--',linewidth=1.5)    
#    plt.title(pulsars[i])
#    plt.yticks(fontsize=12)
#    ticksMHz = (1000*freqtaus[i][:,0]).astype(np.int)[0:len(1000*freqtaus[i][:,0]):2] 
#    plt.legend(fontsize=14)
#    plt.xscale('log')
#    plt.yscale('log')
#    plt.xlim(xmin=(freqtaus[i][0,0])*950,xmax=(freqtaus[i][-1,0])*1050)
#    plt.xticks(ticksMHz,ticksMHz,fontsize=12)
#
#for i in range(2):
#    plt.subplot(2,2,i+3)
#    plt.xlabel(r'$\nu$ (MHz)',fontsize=16,labelpad=15.0)
#    plt.subplot(2,2,2*i+1)
#    plt.ylabel(r'$\tau$ (sec)',fontsize=16)
    
    
#fig = plt.figure(3,figsize=(10,4))
#for i in (8,9):
#    plt.subplot(1,2,i-7)
#    fig.subplots_adjust(left = 0.08, right = 0.97, wspace=0.2, hspace = 0.1, bottom=0.2)
#    plt.rc('text', usetex=True)
#    plt.rc('font', family='serif')
#    plt.errorbar(1000*freqtaus[i][:,0],freqtaus[i][:,1],yerr=freqtaus[i][:,2],fmt='k*',markersize=9.0,capthick=2,linewidth=1.5,label=r'$\alpha = %.2f \pm %.2f$' %(alphas[i],alphaerrs[i]))
#    plt.errorbar(1000*freqtaus_1D[i][:,0],freqtaus_1D[i][:,1],yerr=freqtaus_1D[i][:,2],fmt='k^',markersize=9.0,alpha = 0.5, capthick=2,linewidth=1.5,label=r'$\alpha = %.2f \pm %.2f$' %(alphas1D[i],alphaerrs1D[i]))    
#    plt.plot(1000*freqtaus[i][:,0],fits[i],'k--',linewidth=1.5)
#    plt.plot(1000*freqtaus_1D[i][:,0],fits1D[i],'k--',linewidth=1.5)
#    plt.title(pulsars[i])
#    plt.yticks(fontsize=12)
#    ticksMHz = (1000*freqtaus[i][:,0]).astype(np.int)[0:len(1000*freqtaus[i][:,0]):2] 
#    plt.legend(fontsize=14)
#    plt.xscale('log')
#    plt.yscale('log')
#    plt.xlim(xmin=(freqtaus[i][0,0])*950,xmax=(freqtaus[i][-1,0])*1050)
#    plt.xticks(ticksMHz,ticksMHz,fontsize=12)
#
#
#for i in range(2):
#    plt.subplot(1,2,i+1)
#    plt.xlabel(r'$\nu$ (MHz)',fontsize=16,labelpad=15.0)
#    plt.subplot(1,2,1+i*0)
#    plt.ylabel(r'$\tau$ (sec)',fontsize=16)
#
#
#for i in range(2):
#    plt.subplot(2,2,i+3)
#    plt.xlabel(r'$\nu$ (MHz)',fontsize=16,labelpad=15.0)
#    plt.subplot(2,2,2*i+1)
#    plt.ylabel(r'$\tau$ (sec)',fontsize=16)

#
#"""Plots associated with B1922 for proposal"""
#"""Calculate Tau at 230 MHz based on frequency close to 150 MHz"""
##indices = []
##for i in range(15):
##    indx = np.where((freqtaus[i][:,0] <0.156 ) & (freqtaus[i][:,0]>0.15))[0]
##    indices.append(indx)
#    
#indices = np.array([8,8,4,7,8,2,8,3,8,8,4,8,4,3,8])
#freqss = []
#for i in range(15):
#    freqs = freqtaus[i][:,0]
#    freqss.append(freqs)
#FF = []
#for i in range(15):
#    F = freqss[i][indices[i]]
#    FF.append(F)
#    
#newT = []
#for i in range(15):
#    newtaus1 = psr.tauatfreq(FF[i],freqtaus[i][indices[i],1],0.23,alphas[i])
#    newT.append(newtaus1)
#
#newT4 = []
#for i in range(15):
#    newtaus2 = psr.tauatfreq(FF[i],freqtaus[i][indices[i],1],0.23,3.4)
#    newT4.append(newtaus2)
#
#newT_340 = []
#for i in range(15):
#    newtaus3 = psr.tauatfreq(FF[i],freqtaus[i][indices[i],1],0.34,alphas[i])
#    newT_340.append(newtaus3)
#
#newT4_340 = []
#for i in range(15):
#    newtaus4 = psr.tauatfreq(FF[i],freqtaus[i][indices[i],1],0.34,4.0)
#    newT4_340.append(newtaus4)
#
#tt1 = newT[11]
#tt2 = newT4[11]
#tt3 =newT_340[11]
#tt4 = newT4_340[11]
#
#print tt1,tt2,tt3,tt4
#
#freq22 = freqss[11]
#taus22 = freqtaus[11][:,1]
#taus22err = freqtaus[11][:,2]
#erp = 0.1
#
#nF = np.append(freq22,0.23)
#nF_GBT = np.append(freq22,0.34)
#
#nT1 = np.append(taus22, tt1)
#nT2 = np.append(taus22, tt2)
#
#nTerr1 =np.append(taus22err, erp*tt1)
#nTerr2= np.append(taus22err, erp*tt2)
#
#powp1 = powmod.guess(nT1,x=nF)
#powo1 = powmod.fit(nT1,powp1,x=nF, weights=1/nTerr1**2)
#
#powp2 = powmod.guess(nT2[-5:17],x=nF[-5:17])
#powo2 = powmod.fit(nT2[-5:17],powp2,x=nF[-5:17], weights=1/nTerr2[-5:17]**2)
#
#
#
#plt.figure(10)
#i = 11
#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')
#plt.errorbar(1000*freqtaus[i][:,0],freqtaus[i][:,1],yerr=freqtaus[i][:,2],fmt='k*',markersize=10.0,capthick=2,linewidth=1.5,label=r'$\alpha = %.2f \pm %.2f$' %(alphas[i],alphaerrs[i]))
#plt.errorbar(230,nT1[-1],yerr=nTerr1[-1],fmt='go',markersize=9.0,capthick=2,linewidth=1.5,label=r'$\alpha = %.2f \pm %.2f$' %(-powo1.best_values['exponent'], powo1.params['exponent'].stderr))
#plt.errorbar(230,nT2[-1],yerr=nTerr2[-1],fmt='m^',markersize=9.0,capthick=2,linewidth=1.5,label=r'$\alpha = %.2f \pm %.2f$' %(-powo2.best_values['exponent'], powo2.params['exponent'].stderr))
#plt.plot(1000*nF,powo1.best_fit,'g-')
#plt.plot(1000*nF[-5:17],powo2.best_fit,'m:',lw=1.5)
#plt.plot(1000*freqtaus[i][:,0],fits[i],'k--',linewidth=1.5)
#plt.title('%s extrapolation' %pulsars[i])
#plt.yticks(fontsize=12)
##ticksMHz = (1000*freqtaus[i][:,0]).astype(np.int)[0:len(1000*freqtaus[i][:,0]):2] 
#ticksMHz = [110,120,130,140,150,160,170,180,190,200,210,220,230]
#plt.legend(fontsize=14)
#plt.xscale('log')
#plt.yscale('log')
#plt.xlim(xmin=(freqtaus[i][0,0])*950,xmax=nF[-1]*1050)
#plt.ylim(ymin=0.005)
#plt.xticks(ticksMHz,ticksMHz,fontsize=12)
##plt.xlabel(r'$\nu$ (MHz)',fontsize=16,labelpad=15.0)
#plt.xlabel(r'$\nu$ (MHz)',fontsize=16)
#plt.ylabel(r'$\tau$ (sec)',fontsize=16)

"""VLBI analysis"""
"""Using the measured values of tau at say 150 MHz (perhaps later change this to the modelled values)
I want to calculate the minimum size of a midway screen"""

indices = []
for i in range(15):
    indx = np.where((freqtaus[i][:,0] <0.160 ) & (freqtaus[i][:,0]>0.15))[0][0]
    indices.append(indx)

"""The frequencies close to 150 MHz for each dataset is:"""
freqss = []
for i in range(15):
    fre =  freqtaus[i][indices[i],0]
    freqss.append(fre)
#which is mainly 151, 154 MHz and for one 159 MHz
    
"""And the taus associated with the approx. 150 MHz:"""
taus150 = []
for i in range(15):
    tau150 =  freqtaus[i][indices[i],1]
    taus150.append(tau150)
taus150 = np.array(taus150)  


"""Redo using the modelfit taus at close to 150 MHz"""
"""Indices at which the freq are close to 150"""
indices = []
for i in range(15):
    indx = np.where((freqtaus[i][:,0] <0.160 ) & (freqtaus[i][:,0]>0.15))[0][0]
    indices.append(indx)

"""The frequencies close to 150 MHz for each dataset is:"""
freqss = []
for i in range(15):
    fre =  freqtaus[i][indices[i],0]
    freqss.append(fre)
#which is mainly 151, 154 MHz and for one 159 MHz


"""Taus at these freq. from the model fits"""
fittaus = []
for i in range(15):
    fittau =  fits[i][indices[i]]
    fittaus.append(fittau)
fittaus = np.array(fittaus)  
    
"""Taus at 150MHz, transferred using alpha"""
newT_150 = []
for i in range(15):
    newtaus150 = psr.tauatfreq(freqss[i],fittaus[i],0.15,alphas[i])
    newT_150.append(newtaus150)
newT_150 = np.array(newT_150) 
    
light = 9.71561189*10**-12   ##the speed of light in kpc/sec
mrad = 4.85*10**-9 ##convert mas to rad

GTab = np.loadtxt('/Users/marisa/Dropbox/Aris/LOFARCommCensus/DMDistAlphaGeyer.txt', skiprows = 1)
gdist = np.delete(GTab[:,1],4) ##Getting rid of the commissioning data for B0611+22, for which I only use the Cycle 5 data.


sigth_rad = np.sqrt(newT_150*light/gdist)  ##sigma_theta in rad
sigth_mas = np.sqrt(newT_150*light/gdist)/mrad


"""Assuming the screen size is defined by nn sigma, that means observing angle theta is"""

nn = 3
theta_rad = nn*sigth_rad
theta_mas = nn*sigth_mas

"""Alternatively using the Wucknitz2013 expression for tau when D = 2Ds"""

thetaW_rad = np.sqrt(2*newT_150*light/gdist) 
thetaW_mas = np.sqrt(2*newT_150*light/gdist)/mrad

"""Required baselines"""
wvlength = 2 #(@ 150MHz, in metres)

Bline = wvlength/theta_rad #in metres
Bline_km = Bline/1000. #in km

Bline_W_km = wvlength/(thetaW_rad*1000.)

print "Baselines:" 
print Bline_km
print "Baselines_W:" 
print Bline_W_km


