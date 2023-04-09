#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
from statsmodels.nonparametric.kde import KDEUnivariate as KDE
import corner

## function to plot a KDE from the distributions
def plotKDE(results,colors=[],qs=[16,84],
			font=12,alpha=0.5,path='./name',save=1,
			usetex=False,figsize=(6.4,4.8),
			**kwargs):
	'''Plot a KDE from the distributions.'''
	if not len(colors):
		colors = ['C{}'.format(ii) for ii in range(9)]

	plt.rc('text',usetex=usetex)
	fig = plt.figure(figsize=figsize)
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
		sols = np.array(list(ldict['distribution']))
		#keep = np.percentile(sols,[1,99])
		#print(keep)
		#sols = sols[(sols > keep[0]) & (sols < keep[1])]
		#print(sols)
		#sols = sols[(sols < 15)]# & (sols < keep[1])]
		kde = KDE(sols)
		#kde.fit()#kernel='gau', bw='scott', fft=True, gridsize=1000)
		kde.fit(kernel='gau', bw='scott', fft=True, gridsize=1000)
		ks = kde.support
		kd = kde.density/np.max(kde.density)
		
		nlam = np.median(sols)
		vals = np.percentile(sols,qs)
		low, up = vals[0],vals[1]
		lowerBounds.append(low)
		upperBounds.append(up)
			
		
		
		marr = ldict['model']
		xmin, xmax = min(marr[:,0]), max(marr[:,0])
		ax1.plot(marr[:,0],marr[:,1],color=colors[ii],lw=3,label=r'$\lambda = {:0.1f}^\circ$'.format(key))
		#print(nlam,up,low)
		#if ii: ax3.plot(sols,np.zeros(len(sols)),marker='|',color='k',ls='none')
		
		if not ii:
			arr = ldict['observations']
			times = arr[:,0]
			rvs = arr[:,1]
			errs = arr[:,2]
			ax1.errorbar(times,rvs,yerr=errs,fmt='.',color='k',zorder=10)
			#subtract = ldict['subtract']

			lam_low_calc = nlam - low
			lam_up_calc = up - nlam
			unc_label = r'$\sigma (\lambda) = {:0.1f}^\circ$'.format(np.mean([lam_low_calc,lam_up_calc]))
			ax1.text(0.1,0.1,unc_label,fontsize=font*0.9,
			horizontalalignment='center',
			verticalalignment='center',
			transform = ax1.transAxes,
			bbox=dict(boxstyle="round", fc="white", ec="black", pad=0.2),zorder=15) 
			
		subtract = ldict['subtract']
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
	ax1.legend()
	plt.setp(ax2.get_xticklabels(), visible=False)
	plt.subplots_adjust(hspace=0.08)
	#plt.savefig(path+'_MC_KDE.pdf')
	#plt.savefig(path+'_MC_KDE.png',bbox_inches='tight',dpi=500)
	plt.savefig(path+'_MC_KDE.pdf',bbox_inches='tight')
	#ax.hist(data,bins=100,density=True,alpha=0.5,label='data')
	#ax.legend()
	#return ax

def plotCorner(filename='samples.npy',plotHist=False):
	arr = np.load(filename)
	keys = arr.dtype.names
	ndim = len(keys)
	data = np.empty(shape=(len(arr),ndim))
	for i,k in enumerate(keys):
		data[:,i] = arr[k]
	fig = corner.corner(data, labels=keys, quantiles=[0.16, 0.5, 0.84])

	if plotHist:
		fig, axes = plt.subplots(ndim, 1)
		for i in range(ndim):
			ax = axes[i]
			ax.hist(data[:,i],100,color='k',histtype='step')
			ax.set_yticks([])
			ax.set_xlabel(r"${}$".format(keys[i]))
			ax.set_ylabel(r"$p({})$".format(keys[i]))
			

# def plotHist(filename='samples.npy'):

# 	#def evince(self):
# 		#import matplotlib.pyplot as plt
	
# 		fig = plt.figure(figsize=(10,10))
# 		ax = fig.add_subplot(111)
# 		ax.hist(self.samples[:, 0], 100, color="k", histtype="step")
# 		plt.gca().set_yticks([]);