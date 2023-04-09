#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd

class Target(object):
	'''Target class.

	Class to represent orbital/system parameters.

	Attributes:
		per (float): Period (days).
		T0 (float): Reference time (mid-transit time) (BJD).
		rp (float): Planet radius (Stellar radii).
		a (float): Semi-major axis (Stellar radii).
		b (float): Impact parameter (Stellar radii).
		ecc (float): Eccentricity.
		w (float): Argument of periastron (degrees).
		vsini (float): Rotational velocity (km/s).
		inc (float): Inclination (degrees).
		zeta (float): Macro-turbulence (km/s).
		xi (float): Micro-turbulence (km/s).
		lam (float): Projected obliquity (degrees).
		RVsys (float): Systemic radial velocity (m/s).
		K (float): Semi-amplitude (m/s).
		c1 (float): Linear limb-darkening coefficient.
		c2 (float): Quadratic limb-darkening coefficient.
		name (str): Name of the target.
		s_T0 (float): Uncertainty in T0.
		s_rp (float): Uncertainty in rp.
		s_a (float): Uncertainty in a.
		s_b (float): Uncertainty in b.
		s_ecc (float): Uncertainty in ecc.
		s_w (float): Uncertainty in w.
		s_vsini (float): Uncertainty in vsini.
		s_per (float): Uncertainty in per.
		
		fromCSV (function): Instantiate a target from a csv file.
		printColumn (function): Print the column of the csv file.
		printColumns (function): Print the column names of the csv file.
		
		

	'''

	def __init__(self,per=10.,T0=0.,rp=0.1,
				a=15.,b=0.4,ecc=0.,w=90.,vsini=5.,
				inc=None,zeta=2,xi=2,lam=0,
				K=0,s_K=1e-10,
				RVsys=0,s_RVsys=1e-10,
				c1=0.3,c2=0.2,
				name='Planet name',
				s_T0=1e-10,s_rp=0.01,
				s_a=1,s_b=0.1,
				s_ecc=0.1,s_w=10.,
				s_zeta=0.5,s_xi=0.5,
				s_vsini=0.5,s_per=1e-10):
		self.name = name
		self.per = per
		self.s_per = s_per
		self.T0 = T0
		self.s_T0 = s_T0
		self.rp = rp
		self.s_rp = s_rp
		self.a = a
		self.s_a = s_a
		self.b = b
		self.s_b = s_b
		self.inc = inc
		self.ecc = ecc
		self.s_ecc = s_ecc
		self.w = w
		self.s_w = s_w
		self.K = K
		self.s_K = s_K
		self.RVsys = RVsys
		self.s_RVsys = s_RVsys
		self.vsini = vsini
		self.s_vsini = s_vsini
		self.zeta = zeta
		self.s_zeta = s_zeta
		self.xi = xi
		self.s_xi = s_xi
		self.c1 = c1
		self.c2 = c2
		self.lam = lam
		
		
	## Instantiate a target from a csv file
	## set parameters in the csv file to be the attributes of the target
	def fromCSV(self,filename,name,skiprows=2):
		'''Instantiate a target from a csv file.
		
		:param filename: Name of the csv file.
		:param name: Name of the target.
		:param skiprows: Number of rows to skip.
		
		'''
		df = pd.read_csv(filename,skiprows=skiprows)
		df = df[df['name'] == name]
		
		keys = self.__dict__
		for key in keys:
			if key in list(df.keys()):
				self.__dict__[key] = df[key].values[0]
			#else:
			#    self.__dict__[key] = None

	def printColumn(self,filename,column='name',skiprows=2):
		'''Print the column of the csv file.'''
		df = pd.read_csv(filename,skiprows=skiprows)
		print(list(df[column]))
		
	def printColumns(self,filename,skiprows=2):
		'''Print the column names of the csv file.'''
		df = pd.read_csv(filename,skiprows=skiprows)
		print(list(df.keys()))