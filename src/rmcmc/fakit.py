#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import emcee
import multiprocessing
import scipy.optimize as sop
import scipy.stats as ss
from scipy.special import erf
import numpy as np
import pickle
from astroplan import FixedTarget, Observer
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.simbad import Simbad
from .dynamics import getRV
from .dynamics import totalDuration

## phasefold time onto a given period and reference
def time2phase(time,period,T0):
	'''Convert time to phase.

	Phase centered on the reference time, i.e., from -0.5 to 0.5.

	:param time: Times to convert.
	:type time: array
	:param per: Period.
	:type per: float
	:param T0: Reference time (mid-transit time).
	:type T0: float

	:rerturn: Phase.
	:rtype: array

	'''
	phase = ((time-T0)%period)/period
	for ii in range(len(phase)):
		if phase[ii] > 0.5:
			phase[ii] -= 1
	return phase

## read results from binary file
def readBinary(filename):
	with open(filename,'rb') as f:
		return pickle.load(f)
	

def flatPriorDis(r,x1,x2):
	'''Uniform prior distribution.
	
	:param r: Value to evaluate.
	:type r: float
	:param x1: Lower bound.
	:type x1: float
	:param x2: Upper bound.
	:type x2: float
	
	:return: Prior probability.
	
	
	'''
	return x1 + r*(x2 - x1)

def tgaussPriorDis(mu,sigma,xmin,xmax):
	'''Truncated Gaussian distribution.

	See :py:func:`tgauss_prior`.

	:param xmid: :math:`\mu`.
	:type xmid: float

	:param xwid: :math:`\sigma`.
	:type xwid: float
	
	:param xmin: :math:`a`.
	:type xmin: float

	:param xmax: :math:`b`.
	:type xmax: float

	:returns: :math:`f(x)`
	:rtype: float

	'''
	a = (xmin - mu)/sigma
	b = (xmax - mu)/sigma
	return ss.truncnorm.rvs(a,b,loc=mu,scale=sigma)

def flatPrior(val,xmin,xmax):
	'''Uniform prior.

	.. math:: 
		f(x) = \\frac{1}{b-a}, \, \ \mathrm{for} \ a \leq x \leq b,\ \mathrm{else} \, f(x) = 0 \, .
	
	:param val: :math:`x`.
	:type val: float
	
	:param xmin: :math:`a`.
	:type xmin: float

	:param xmax: :math:`b`.
	:type xmax: float

	:returns: :math:`f(x)`
	:rtype: float

	'''
	if val < xmin or val > xmax:
		ret = 0.0
	else:
		ret = 1/(xmax - xmin)
	return ret

def tgaussPrior(val,xmid,xwid,xmin,xmax):
	'''Truncated Gaussian prior.
	
	.. math:: 
		f (x; \mu, \sigma, a, b) = \\frac{1}{\sigma} \\frac{g(x)}{\Phi(\\frac{b-\mu}{\sigma}) - \Phi(\\frac{a-\mu}{\sigma})} \, ,

	where :math:`g(x)` is the Gaussian from :py:func:`gauss_prior` and :math:`\Phi` is the (Gauss) error function (`scipy.special.erf <https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.erf.html?highlight=erf#scipy.special.erf>`_ ).

	:param val: :math:`x`.
	:type val: float

	:param xmid: :math:`\mu`.
	:type xmid: float

	:param xwid: :math:`\sigma`.
	:type xwid: float
	
	:param xmin: :math:`a`.
	:type xmin: float

	:param xmax: :math:`b`.
	:type xmax: float	

	:returns: :math:`f(x)`
	:rtype: float

	'''
	if val < xmin or val > xmax:
		ret = 0.0
	else:
		nom = np.exp(-1.0*np.power(val-xmid,2) / (2*np.power(xwid,2))) / (np.sqrt(2*np.pi)*xwid)
		den1 = (1. + erf((xmax - xmid)/(np.sqrt(2)*xwid)))/2
		den2 = (1. + erf((xmin - xmid)/(np.sqrt(2)*xwid)))/2
		ret = nom/(den1 - den2)
	return ret

class Simulator(object):
	
	parameters = ['per','T0','ecc','w','K','RVsys','a',
				'rp','b','vsini','c1','c2','zeta','xi']
	
	def __init__(self,target,
				ndraws=500,nproc=1,
				times=np.array([]),
				precision=5,jitter=0.0,
				variables=['vsini','b']):
		self.target = target
		self.ndraws = ndraws
		self.times = times
		self.precision = precision
		self.jitter = jitter
		self.uncertainty = np.sqrt(self.precision**2 + self.jitter**2)
		self.variables = variables
		self.nproc = nproc
			
	#def model(self,target,RM=True,mpath=None,mtimes=np.array([])):
	def model(self,RM=True,mpath=None,mtimes=np.array([])):

		T0 = self.target.T0
		per = self.target.per
		ecc = self.target.ecc
		w = self.target.w
		K = self.target.K
		RVsys = self.target.RVsys
		a = self.target.a
		b = self.target.b

		Rp = self.target.rp
		
		## Projected stellar rotation
		vsini = self.target.vsini
		## Macroturbulence 
		zeta = self.target.zeta
		## Microturbulence 
		xi = self.target.xi

		## The all important obliquity
		lam = self.target.lam
		## LD coefficients
		c1, c2 = self.target.c1, self.target.c2
		
		if not len(mtimes): times = self.times
		else: times = mtimes
		
		RV = getRV(times,per,T0,ecc,w,K,RVsys,a,Rp,b,vsini,lam,c1,c2,zeta,xi,RM=RM,mpath=mpath)

		return RV

	## function to simulate time series with sampling of exposure time
	#def timeSimulator(self,target,exposure=900,start=1800,end=1800):
	def timeSimulator(self,exposure=900,start=1800,end=1800,overhead=50):
		#p, rp, ar, inc, ecc, ww = target.per, target.rp, target.a, target.inc, target.ecc, target.w
		p, rp, ar, b, ecc, ww = self.target.per, self.target.rp, self.target.a, self.target.b, self.target.ecc, self.target.w
		#inc *= np.pi/180.0
		inc = np.arccos(b/ar)
		ww = np.deg2rad(ww)
		duration = totalDuration(p,rp,ar,inc,ecc,ww)
		half = 0.5*duration
		T0 = self.target.T0

		exposure += overhead

		exposure /= 24*3600
		start /= 24*3600
		end /= 24*3600

		# number of exposures
		nexp = int((duration + start + end)/exposure)
		#self.times = np.linspace(T0-half-start,T0+half+end,nexp,endpoint=True)
		self.times = np.arange(T0-half-start,T0+half+end,exposure)

	def noiseSimulator(self):
		return np.random.normal(loc=0.0,scale=self.uncertainty,size=len(self.times))
	
	#def dataSimulator(self,target,**kwargs):
	def dataSimulator(self,**kwargs):
		if not len(self.times):
			self.timeSimulator(**kwargs)
		rvs = self.model()
		rvs += self.noiseSimulator()
		
		self.rvs = rvs
		self.error = np.ones(len(self.rvs))*self.uncertainty
		
	
	def modelSimulator(self):
		if not len(self.times):
			self.timeSimulator()
			
		self.mtimes = np.linspace(min(self.times)-0.005,max(self.times)+0.006,300)
		self.mRVs = self.model(mtimes=self.mtimes)
	


		

	# #def logLike(self,planet,simulated,lam):
	# def logLike(self,lam,per,T0,ecc,w,K,RVsys,a,Rp,b,vsini,c1,c2,zeta,xi,RM=True,mpath=None):
	# 	'''Log likelihood function.
		
	# 	Function to minimize.
		
	# 	'''
		
	# 	rvs = getRV(self.times,per,T0,ecc,w,K,RVsys,a,Rp,b,vsini,lam,c1,c2,zeta,xi,RM=RM,mpath=mpath)
	# 	rvs += self.noiseSimulator()
		
	# 	chi2 = (rvs - self.simulated)**2/self.precision**2
	# 	#den = np.log(1./(np.sqrt(2*np.pi)*self.precision))
	# 	#nom = -0.5*chi2
	# 	#loglike = np.sum(nom + den)
	# 	#print(np.sum(chi2),lam)
	# 	print(lam)
	# 	print(rvs ,self.simulated)
	# 	return np.sum(chi2)#-1*loglike

	#def logLike(self,planet,simulated,lam):
	def logLike(self,lam,per,T0,ecc,w,K,RVsys,a,Rp,b,vsini,c1,c2,zeta,xi,RM=True,mpath=None):
		'''Log likelihood function.
		
		Function to minimize.
		
		'''
		
		rvs = getRV(self.times,per,T0,ecc,w,K,RVsys,a,Rp,b,vsini,lam,c1,c2,zeta,xi,RM=RM,mpath=mpath)
		#rvs += self.noiseSimulator()
		
		#chi2 = (rvs - self.simulated)**2/self.precision**2
		chi2 = (rvs - self.simulated)**2/self.precision**2
		den = np.log(1./(np.sqrt(2*np.pi)*self.precision))
		nom = -0.5*chi2
		loglike = np.sum(nom + den)
		return -1*loglike#np.sum(chi2)#
	
	def flatPriorDis(self,r,x1,x2):
		'''Uniform prior distribution.
		
		:param r: Value to evaluate.
		:type r: float
		:param x1: Lower bound.
		:type x1: float
		:param x2: Upper bound.
		:type x2: float
		
		:return: Prior probability.
		
		
		'''
		return x1 + r*(x2 - x1)
	
	def calculateLambda(self,lam,per,T0,ecc,w,K,RVsys,a,Rp,b,vsini,c1,c2,zeta,xi,RM=True,mpath=None):
		'''Optimize lambda.
		
		Function to find the value of lambda that optimizes the chi-2.
		
		'''
		#self.target.lam = lam
		self.simulated = self.model() + self.noiseSimulator()
		#func = lambda ll : self.logLike(ll,per,T0,ecc,w,K,RVsys,a,Rp,b,vsini,c1,c2,zeta,xi,RM=RM,mpath=mpath)
		guess = self.flatPriorDis(np.random.uniform(),lam-30.,lam+30.)
		#res = sop.minimize(func,guess,method='Nelder-Mead',tol=1e-6,bounds=[(-180,180)])
		res = sop.minimize(self.logLike,guess,method='Nelder-Mead',tol=1e-6,bounds=[(-180,180)],
			args=(per,T0,ecc,w,K,RVsys,a,Rp,b,vsini,c1,c2,zeta,xi,RM,mpath))
		return res.x


	def createRun(self,lam):
		self.runs = []
		pars = self.target.__dict__
		for ii in range(self.ndraws):
			vals = [lam]
			for par in self.parameters:
				if par in self.variables:
					val = np.random.normal(loc=pars[par],scale=pars['s_'+par])
				else:
					val = pars[par]
				vals.append(val)
			self.runs.append(tuple(vals))
			
	
	def MC(self,lambdas,writeBinary=False,path='./',**kwargs):
		
		self.results = {}
		if not len(self.times): self.timeSimulator()
		for lam in lambdas:
			self.results[lam] = {}
			self.target.lam = lam
			#self.simulated = self.model() + self.noiseSimulator()
			self.createRun(lam)
			
			if self.nproc > 1:
				with multiprocessing.Pool(self.nproc) as pool:
					self.sols = pool.starmap(self.calculateLambda,self.runs)
			else:
				self.sols  = np.zeros(self.ndraws)
				for ii, run in enumerate(self.runs):
					self.sols[ii] = self.calculateLambda(*run)        
			
			self.dataSimulator()
			arr = np.zeros(shape=(len(self.times),3))
			arr[:,0] = time2phase(self.times,self.target.per,self.target.T0)*self.target.per*24
			arr[:,1] = self.rvs
			arr[:,2] = self.error
			
			self.modelSimulator()
			marr = np.zeros(shape=(len(self.mtimes),2))
			marr[:,0] = time2phase(self.mtimes,self.target.per,self.target.T0)*self.target.per*24
			marr[:,1] = self.mRVs
			
			self.results[lam] = {
				'distribution':self.sols,
				'observations':arr,
				'model':marr,
				'subtract':self.model(),
			}            
	
		## write results to binary file
		if writeBinary:
			## add trailing slash to path if not present
			if not path.endswith('/'): path += '/'
			
			filename = path + self.target.name + '.pkl'
			with open(filename,'wb') as f:
				pickle.dump(self.results,f)

class RVSimulator(object):
	
	parameters = ['per','T0','ecc','w','K','RVsys']
	fitParameters = ['K']
	
	def __init__(self,target,
				nproc=1,
				times=np.array([]),
				precision=5,jitter=0.0,
				variables=['K','ecc','w']):
		self.target = target
		#self.ndraws = ndraws
		#self.nwalkers = nwalkers
		self.times = times
		self.precision = precision
		self.jitter = jitter
		self.errors = np.sqrt(self.precision**2 + self.jitter**2)*np.ones(len(self.times))
		#self.uncertainty = np.sqrt(self.precision**2 + self.jitter**2)
		self.variables = variables
		self.nproc = nproc
	
	def priors(self):
		'''Prior probability.
		
		'''
		self.prior = {}
		self.prior['per'] = ['tgauss',self.target.per,self.target.s_per,0.5,self.target.per+0.5]
		self.prior['T0'] = ['tgauss',self.target.T0,self.target.s_T0,self.target.T0-1,self.target.T0+1]
		self.prior['ecc'] = ['tgauss',self.target.ecc,self.target.s_ecc,0.0,1.0]
		self.prior['w'] = ['tgauss',self.target.w,self.target.s_w,0.0,360.0]
		self.prior['K'] = ['uniform',self.target.K,self.target.s_K,0.0,100.0]
		self.prior['RVsys'] = ['tgauss',self.target.RVsys,self.target.s_RVsys,-100.0,100.0]

	def currentValues(self):
		'''Current values of parameters.
		
		'''
		self.current = {}
		self.current['per'] = self.target.per
		self.current['T0'] = self.target.T0
		self.current['ecc'] = self.target.ecc
		self.current['w'] = self.target.w
		self.current['K'] = self.target.K
		self.current['RVsys'] = self.target.RVsys

	def badWeather(self,fracDays=0.5):
		luck = np.random.uniform(0,1)
		if luck < fracDays:
			return True
		else:
			return False
	
	def telescopeSites(self,telescope):
		telescopes = {
			'SONG' : Observer(longitude = -16.509,#343.5033*u.deg,
				latitude = 28.2917*u.deg,
				elevation = 2395*u.m,
				timezone = 'Europe/London',
				name = 'Teide Observatory'), 
			'MUSCAT2' : Observer(longitude = -16.509,#343.5033*u.deg,
				latitude = 28.2917*u.deg,
				elevation = 2395*u.m,
				timezone = 'Europe/London',
				name = 'Teide Observatory'),
			'ORM' : Observer.at_site('Roque de los Muchachos'),
			'NOT' : Observer.at_site('Roque de los Muchachos'),
			'TNG' : Observer.at_site('Roque de los Muchachos'),
			'HARPS-N' : Observer.at_site('Roque de los Muchachos'),
			'DCT' : Observer.at_site('Discovery Channel Telescope'),
			'VLT' : Observer.at_site('Paranal Observatory'),
			'Paranal' : Observer.at_site('Paranal Observatory'),
			'LaSilla' : Observer.at_site('La Silla Observatory'),
			'HARPS' : Observer.at_site('La Silla Observatory'),
			'Keck' : Observer.at_site('Keck Observatory'),
			'Subaru' : Observer.at_site('Subaru'),
			'OHP' : Observer(longitude = 43.9308*u.deg,
				latitude = 5.7133*u.deg,
				elevation = 650*u.m,
				timezone = 'Europe/Paris',
				name = 'Haute-Provence Observatory'),
			'CASLEO' : Observer(longitude = -31.7986*u.deg,
				latitude = -69.2956*u.deg,
				elevation = 2483*u.m,
				timezone = 'Etc/GMT-3',
				name = 'Leoncito Astronomical Complex'),
			'BTA' : Observer(longitude = 43.6468*u.deg, 
				latitude = 41.4405*u.deg,
				elevation = 2070*u.m,
				timezone = 'Etc/GMT+5',
				name = 'Special Astrophysical Observatory'),
			'TLS' : Observer(longitude = 11.711167*u.deg,
				latitude = 50.98011*u.deg,
				elevation = 341*u.m,
				timezone = 'Europe/Berlin',
				name = 'Thueringer Landessternwarte Tautenbug')
		}
		return telescopes[telescope]
	
	def seeing(self,dist='uniform',low=0.2,high=2.0,mu=1.0,sigma=0.5):
		if dist == 'uniform':
			val = np.random.uniform(low,high)
		elif dist == 'normal':
			val = np.random.normal(mu,sigma)
			if val < low:
				val = low
		return val

	def timeSimulator(self,exposure=900,
					name='TIC 287156968',
					site='VLT',
					ra=None,dec=None,
					nrvs=50,
					start = Time('2023-10-01T12:00:00'),
					end = Time('2024-03-31T12:00:00'),
					sampling=1,
					seeingVals=['uniform',0.5,1.5,None,None],
					):
		
		if (ra is None) or (dec is None):
			result_table = Simbad.query_object(name)
			ra, dec = result_table['RA'][0], result_table['DEC'][0]		
		obs_target = FixedTarget(SkyCoord(ra=ra,dec=dec,unit='deg'))

		observatory = self.telescopeSites(site)

		semester = Time(np.arange(start.jd,end.jd,sampling),format='jd')

		self.times = np.array([])
		self.errors = np.array([])

		for day in semester:
			bw = self.badWeather()
			if bw:
				print('Bad weather on {}'.format(day.isot))
			else:
				print('Good weather on {}'.format(day.isot))
				#scope.scope_target.getVisPlot(tar.targets,tel.location,day.iso)

				#mnt = observatory.twilight_morning_nautical(time, which='next')
				eat = observatory.twilight_evening_astronomical(day, which='next')
				mat = observatory.twilight_morning_astronomical(day, which='next')
				nexps = (mat.jd - eat.jd)*24*3600/exposure

				grid = Time(np.linspace(eat.jd,mat.jd,int(nexps)),format='jd')
				up = observatory.target_is_up(grid, obs_target,horizon=25*u.deg)
				chance = len(grid[up])/len(grid)
				print('Chance of observing: {}'.format(chance))
				luck = np.random.uniform(0,1)
				if chance > luck:
					idx = np.random.randint(len(grid[up]))
					timestamp = grid[up][idx]
					self.times = np.append(self.times,timestamp.jd)
					#seeing = np.random.uniform(0.2,2.0)
					see = self.seeing(*seeingVals)
					err = np.sqrt((self.precision*(1+see))**2 + self.jitter**2)
					self.errors = np.append(self.errors,err)
					#noise = np.random.normal(0,err)
					#noises.append(noise)
					#rv = tracit.get_RV(timestamp.jd, op)
					#rvs = np.append(rvs,rv) + np.random.normal(0,err)
					altaz = observatory.altaz(timestamp, obs_target)
					air = 90 - 180.*np.arccos(1/altaz.alt.value)/np.pi
		


	## chi2 function
	def chi2(self):
		'''Chi-squared function.
		
		
		'''

		return np.sum(((self.obs-self.model)/self.errors)**2)

	def logLike(self):
		'''Log-likelihood function.
		
		
		'''
		nom = -0.5*self.chi2()
		den = np.log(1./(np.sqrt(2*np.pi)*self.errors))
		return np.sum(den+nom)

	## log prior
	def logProb(self,positions):
		'''Log-probability function.
		
		
		'''
		
		logp = 0.0

		for idx, key in enumerate(self.fitParameters):
			self.current[key] = positions[idx]
			dist, mu, sigma, low, high = self.prior[key]
			if dist == 'uniform':
				prob = flatPrior(positions[idx],low,high)
			elif dist == 'tgauss':
				prob = tgaussPrior(positions[idx],mu,sigma,low,high)
			
			if prob != 0.0:
				logp += np.log(prob)
			else:
				return -np.inf

		
		
		RV = getRV(self.times,
			self.current['per'],
			self.current['T0'],
			self.current['ecc'],
			self.current['w'],
			self.current['K'],
			self.current['RVsys'])

		self.model = RV

		return logp + self.logLike()

	
	def addParameter(self,name):
		'''Add a fitting parameter to the model.
		
		:param name: Name of the parameter. parameters: 'per','T0','ecc','w','K','RVsys'
		:type name: str
		
		'''
		self.fitParameters.append(name)

	def setPrior(self,parname,dist='uniform',mu=0.5,sigma=0.1,low=0.0,high=1.0):
		'''Set the prior distribution for a parameter.
		
		:param parname: Name of the parameter.
		:type parname: str
		:param mu: Mean for normal distribution.
		:type mu: float
		:param sigma: Standard deviation for normal distribution.
		:type sigma: float
		:param dist: Distribution to use. Default is uniform.
		:type dist: str
		:param low: Lower bound for uniform distribution.
		:type low: float
		:param high: Upper bound for uniform distribution.
		:type high: float
		
		'''
		self.prior[parname] = [dist,mu,sigma,low,high]
	
	def startValues(self):
		pars = self.fitParameters
		pos = []
		for ii in range(self.nwalkers):
			start = np.ndarray(self.ndim)
			for idx, par in enumerate(pars):
				dist, mu, sigma, low, high = self.prior[par]
				assert dist in ['tgauss','uniform'], print('{} is not a valid option for the starting distribution.'.format(dist))
				if dist == 'tgauss':
					start[idx] = tgaussPriorDis(mu,sigma,low,high)
				elif dist == 'uniform':
					start[idx] = flatPriorDis(np.random.uniform(),low,high)
			pos.append(start)
		return pos

	def setup(self,
		nwalkers=100,
		ndraws=10000,
		threads=1):

		print(
			'Setting up the model...'
			)
		print(
			'Parameters to fit: {}'.format(self.fitParameters)
			)
		print(
			"To add more parameters, use the addParameter('parname') method."
			)
		self.priors()
		for key in self.prior.keys():
			if key in self.fitParameters:
				pri = self.prior[key]
				print('Prior for {}: {}'.format(key,pri))
				print("Change it with the setPrior('par','distribution',mu,sigma,low,high) method.")
				print('Remember to call setup() again.')
		self.ndim = len(self.fitParameters)
		self.nwalkers = nwalkers
		self.ndraws = ndraws
		self.currentValues()
		if not len(self.times):
			self.times = np.arange(0,self.nrvs)
		RV = getRV(self.times,
			self.current['per'],
			self.current['T0'],
			self.current['ecc'],
			self.current['w'],
			self.current['K'],
			self.current['RVsys'])
		noise = [np.random.normal(0,err) for err in self.errors]		
		self.positions = self.startValues()
		self.obs = RV + noise

	def run(self,burnin=5000,filename='samples.npy'):
		print('Running the model...')
		sampler = emcee.EnsembleSampler(self.nwalkers,self.ndim,self.logProb)#,threads=self.threads)
		sampler.run_mcmc(self.positions,self.ndraws)
		#self.sampler = sampler
		print('Done.')
		print(
			'Discard the first {} steps as burn-in.'.format(burnin)
		)
		samples = sampler.get_chain(discard=burnin,flat=True)
		log_prob_samples = sampler.get_log_prob(discard=burnin, flat=True)
		all_samples = np.concatenate(
			(samples, log_prob_samples[:, None]),  axis=1)
		self.samples = all_samples
		labels = self.fitParameters + ['log_prob']
		dtypes = [(fp,float) for fp in labels]
		recarr = np.recarray((self.samples.shape[0],),dtype=dtypes)
		for ii, fp in enumerate(labels):
			recarr[fp] = self.samples[:,ii]
		np.save(filename,recarr)

