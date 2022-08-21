#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import multiprocessing
import scipy.optimize as sop
import numpy as np
import pickle
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
    def timeSimulator(self,exposure=900,start=1800,end=1800):
        #p, rp, ar, inc, ecc, ww = target.per, target.rp, target.a, target.inc, target.ecc, target.w
        p, rp, ar, b, ecc, ww = self.target.per, self.target.rp, self.target.a, self.target.b, self.target.ecc, self.target.w
        #inc *= np.pi/180.0
        inc = np.arccos(b/ar)
        ww *= np.pi/180.0
        duration = totalDuration(p,rp,ar,inc,ecc,ww)
        half = 0.5*duration
        T0 = self.target.T0

        exposure /= 24*3600
        start /= 24*3600
        end /= 24*3600

        # number of exposures
        nexp = int((duration + start + end)/exposure)
        
        self.times = np.linspace(T0-half-start,T0+half+end,nexp)

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
            
        self.mtimes = np.linspace(min(self.times)-0.005,max(self.times)+0.005,300)
        self.mRVs = self.model(mtimes=self.mtimes)
    


        

    #def logLike(self,planet,simulated,lam):
    def logLike(self,lam,per,T0,ecc,w,K,RVsys,a,Rp,b,vsini,c1,c2,zeta,xi,RM=True,mpath=None):
        '''Log likelihood function.
        
        Function to maximize.
        
        '''
        
        rvs = getRV(self.times,per,T0,ecc,w,K,RVsys,a,Rp,b,vsini,lam,c1,c2,zeta,xi,RM=RM,mpath=mpath)
        rvs += self.noiseSimulator()
        
        chi2 = (rvs - self.simulated)**2/self.precision**2
        den = np.log(1./(np.sqrt(2*np.pi)*self.precision))
        nom = -0.5*chi2
        loglike = np.sum(nom + den)
        return -1*loglike
    
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
        
        Function to find the minimum value of lambda.
        
        '''
        
        func = lambda ll : self.logLike(ll,per,T0,ecc,w,K,RVsys,a,Rp,b,vsini,c1,c2,zeta,xi,RM=RM,mpath=mpath)
        guess = self.flatPriorDis(np.random.uniform(),lam-30.,lam+30.)
        res = sop.minimize(func,guess,method='Nelder-Mead',tol=1e-6)
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
            self.simulated = self.model()
            self.createRun(lam)
        
            
            
            if self.nproc > 1:
                with multiprocessing.Pool(self.nproc) as pool:
                    self.sols  = pool.starmap(self.calculateLambda,self.runs)
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
                'subtract':self.simulated,
            }            
    
        ## write results to binary file
        if writeBinary:
            ## add trailing slash to path if not present
            if not path.endswith('/'): path += '/'
            
            filename = path + self.target.name + '.pkl'
            with open(filename,'wb') as f:
                pickle.dump(self.results,f)
            

    
## read results from binary file
def readBinary(filename):
    with open(filename,'rb') as f:
        return pickle.load(f)
    


# def model(time,target,RM=True,mpath=None):


#     RV = getRV(time,target,RM=RM,mpath=mpath)

#     return RV