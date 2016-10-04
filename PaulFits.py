# -*- coding: utf-8 -*-
"""
Created on Thu Sep 15 17:14:12 2016

@author: marisa
"""

"""Code to fit profiles of B1937+21 and extract scattering"""


import argparse
import os, sys
#import pypsr_standalone as psr
import matplotlib.pyplot as plt
from lmfit import Model, conf_interval, printfuncs
from lmfit import minimize, Parameter, Parameters, fit_report
from lmfit.models import LinearModel, PowerLawModel, ExponentialModel, QuadraticModel
import numpy as np
import PaulFunctions as pf


from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)


"""Define options to the script"""
parser = argparse.ArgumentParser()
parser.add_argument('-p','--period',type=float,
                    help="Provide the period of the pulsar in seconds")
parser.add_argument('-r','--rawfile',
                    help="filepath to txt file, produced from shifting observing data")                   
parser.add_argument('-dc','--datacycle',
                    help="The type of data. Choose from comm., census, cycle5. Only used in creating filenames.")
parser.add_argument('-m','--method',
                    help="Choosing method to fit data or simulation. Choose between 'onedim', 'iso', 'aniso','postfold', multi, 'isoplusonedim'")                   


args = parser.parse_args()

"""Allocate variable names to the parsed options"""
pulseperiod = args.period
raw = args.rawfile
datac = args.datacycle
meth = args.method

args = parser.parse_args()

"""Create folder to save to"""
newpath = r'./SummaryPlots_Paul'
if not os.path.exists(newpath):
    os.makedirs(newpath)
    
"""I have split up the observations you have given me into separate observation files, this is not very elegant, but my previously elaborate code could not deal with yet another over loop"""
"""How I run over all is hopefully in some email"""


pulsar, nch, nbins, rawdata, freq = pf.read_Paul(raw)
    
print0 = "Pulsar name: %s" %pulsar
print1 = "Number of channels: %d" %nch
print2 = "Number of bins: %d" %nbins
for k in range(3):
    print eval('print{0}'.format(k))
    
data = rawdata
freqmsMHz = freq*1000

if meth in ('multi'): 
    """This first bit is historic, from when I propted the used to define the region of the Gaussians.
    Now because they are all aligned it's just hardcoded bw1, bw2 = np.linspace(660,720,5), np.linspace(725,750,5) for a free fit"""
    """To make sure you know what the fit is doing look at the function tau_fittermulti() in PaulFunctions - that's where you change it to have fixed fitting parameters if you like"""   
#    bw1 = raw_input("Provide bin-window for 1st the peak (e.g: [20,30]): ")
#    bw1 = eval(bw1)
#    bw2 = raw_input("Provide bin-window start and end 2nd peak (e.g: [35,45]): ")
#    bw2 = eval(bw2)
    bw1, bw2 = np.linspace(660,720,5), np.linspace(725,750,5)        
    peak1 = np.max(data[bw1[0]: bw1[-1]])
    peak2 = np.max(data[bw2[0]: bw2[-1]])
    maxmean = 2.0*bw2[-1]
    nms, bt, tstd, bp, bpstd, chis = [], [], [], [], [], []
    for k in range(0,len(bw1)):
        for j in range(0,len(bw2)):
			noiselessmodel, besttau, taustd, bestparams, bestparams_std, redchi = pf.tau_fittermulti(data,nbins,bw1[k],bw2[j], peak1, peak2, maxmean)
			nms.append(noiselessmodel)
			bt.append(besttau)
			tstd.append(taustd)
			bp.append(bestparams)
			bpstd.append(bestparams_std)
			chis.append(redchi)
    tstd = np.array(tstd)
    nonz = np.nonzero(tstd)
    tstd = np.array(tstd)
    nonz_min = np.argmin(tstd[np.nonzero(tstd)])
    tstd_ind = np.array(nonz).flatten()[nonz_min]     
    noiselessmodel, besttau, taustd, bestparams, bestparams_std, redchi = nms[tstd_ind], bt[tstd_ind], tstd[tstd_ind], bp[tstd_ind], bpstd[tstd_ind], chis[tstd_ind]               
          
else:
    print "Currently only set up for method 'multi' "
    
plt.close(1)
plt.figure(1,figsize=(6,4))
plt.rc('text', usetex = True)
plt.rc('font',family='serif')
plt.plot(data,'m',alpha=0.8, label="data")
plt.plot(noiselessmodel, 'c', alpha=0.7, label="model")
unscat = pf.TwoPeaksModel(np.linspace(1,nbins,nbins), 1024, bp[tstd_ind][4], bp[tstd_ind][5], bp[tstd_ind][2] , bp[tstd_ind][3], bp[tstd_ind][0], bp[tstd_ind][1], bp[tstd_ind][6])
plt.plot(unscat,'g--')    
plt.title('%s at %.1f MHz' %(pulsar, freqmsMHz))
plt.text(800,0.4,r'$\tau: %.4f \pm %.1e$ms' %(besttau*pulseperiod/(2*nbins)*1000, taustd*pulseperiod/(2*nbins)*1000),fontsize=11)
plt.legend(fontsize=12, loc='best')    
plt.xlim(xmin=500,xmax=1024)


#Kplot = '%s_%s.png'  % (pulsar,raw[-6:-4])
#picpathtau = newpath
#fileoutputtau = os.path.join(picpathtau,Kplot)
#plt.savefig(fileoutputtau, dpi=150), np.savetxt('%s_%s.txt' % (pulsar,raw[-6:-4]), np.column_stack((bt[tstd_ind], tstd[tstd_ind])).flatten())
#plt.savefig(fileoutputtau, dpi=150), np.savetxt('%s_%s_params.txt' % (pulsar,raw[-6:-4]), bp[tstd_ind])
#plt.savefig(fileoutputtau, dpi=150), np.savetxt('%s_%s_paramstd.txt' % (pulsar,raw[-6:-4]), bpstd[tstd_ind])
#print 'Saved %s in %s' %(Kplot,newpath) 


  
  

