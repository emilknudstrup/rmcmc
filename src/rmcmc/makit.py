#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
from statsmodels.nonparametric.kde import KDEUnivariate as KDE


## function to plot a KDE from the distributions
def plotKDE(results,colors=[],qs=[16,84],
            font=12,alpha=0.5,**kwargs):
    '''Plot a KDE from the distributions.'''
    if not len(colors):
        colors = ['C{}'.format(ii) for ii in range(9)]


    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(413,sharex=ax1)
    ax3 = fig.add_subplot(414)

    ax2.axhline(0,color='k',ls='--')

    lowerBounds = []
    upperBounds = []
    #subtract = None
    #unc_label = None
    
    for ii, key in enumerate(results.keys()):
        ldict = results[key]
        sols = ldict['distribution']
        kde = KDE(sols)
        kde.fit(kernel='gau', bw='scott', fft=True, gridsize=1000)
        ks = kde.support
        kd = kde.density/np.max(kde.density)
        
        nlam = np.median(sols)
        vals = np.percentile(sols,qs)
        low, up = vals[0],vals[1]
        lowerBounds.append(low)
        upperBounds.append(up)
            
        arr = ldict['observations']
        times = arr[:,0]
        rvs = arr[:,1]
        errs = arr[:,2]
        
        
        marr = ldict['model']
        xmin, xmax = min(marr[:,0]), max(marr[:,0])
        ax1.plot(marr[:,0],marr[:,1],color=colors[ii],lw=3)
        
        if not ii:
            ax1.errorbar(times,rvs,yerr=errs,fmt='.',color='k',zorder=10)
            subtract = ldict['subtract']

            lam_low_calc = nlam - low
            lam_up_calc = up - nlam
            unc_label = r'$\sigma (\lambda) = {:0.0f}^\circ$'.format(np.mean([lam_low_calc,lam_up_calc]))
            ax1.text(0.1,0.1,unc_label,fontsize=font*0.9,
            horizontalalignment='center',
            verticalalignment='center',
            transform = ax1.transAxes,
            bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2),zorder=15) 
            
        ax2.errorbar(times,rvs-subtract,yerr=errs,fmt='.',color=colors[ii],zorder=10)
        
        ax3.plot(ks,kd,label=key,**kwargs)
        shade = (ks > low) & (ks < up)
        ax3.fill_between(ks[shade],kd[shade],facecolor=colors[ii],interpolate=True,alpha=alpha)
        midx = np.argmin(abs(ks-nlam))
        ax3.plot((nlam,nlam),(0,kd[midx]),'k-',zorder=1)        
        lidx = np.argmin(abs(ks-low))
        ax3.plot((low,low),(0,kd[lidx]),'k--',zorder=1)
        uidx = np.argmin(abs(ks-up))
        ax3.plot((up,up),(0,kd[uidx]),'k--',zorder=1)
        
    ax1.set_ylabel(r'$\rm RV \ (m/s)$',fontsize=font)
    ax1.set_xlabel(r'$\rm Hours \ from \ midtransit$',fontsize=font)
    ax1.xaxis.tick_top()
    ax1.xaxis.set_label_position('top')
    ax2.xaxis.tick_top()
    ax2.set_ylabel(r'$\rm O-C$',fontsize=font*0.9)
    
    
    ax3.set_ylabel(r'$\rm KDE$',fontsize=font)
    ax3.set_xlabel(r'$\rm \lambda \, (^\circ)$',fontsize=font)
    
    ax1.set_xlim(xmin,xmax)
    ax3.set_ylim(0,1.1)
    ax3.set_xlim(min(lowerBounds)-5,max(upperBounds)+5)
    
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.subplots_adjust(hspace=0.08)
    #ax.hist(data,bins=100,density=True,alpha=0.5,label='data')
    #ax.legend()
    #return ax