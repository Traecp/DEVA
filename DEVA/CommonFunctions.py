#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy as N
import operator
from scipy.fftpack import fft, fftfreq, fftshift
from lmfit import Parameters, minimize
	
def sort_table(table, col=0, reverse=False):
	#reverse = True: Z-->A
	#reverse = False: A-->Z
	return sorted(table, key=operator.itemgetter(col), reverse=reverse)

def list_to_table(main_store, sort_col=2):
	""" Convert the list of edf image into a table of 3 columns:
	EDF_name : EDF_prefix : EDF_number
	The EDF name must be in format: prefix_number.edf
	"""
	main_table = []
	if len(main_store)>0:
		for i in range(len(main_store)):
			try:
				if "_" in main_store[i]:
					e = []
					e.append(main_store[i])
					row=main_store[i].split("_")
					e.append(row[0])
					row = row[1]
					row = row.split(".")[0]
					if 'g' in row:
						row = row[:-1]
					e.append(int(row))
					main_table.append(e)
				else:
					continue
			except ValueError:
				continue
		store = sort_table(main_table,col=sort_col)
		return store
	else:
		return {}
	

def Fourier(X,vect):  
	Nb  = vect.size   #number of data points
	T  = X[1] - X[0] #sample spacing
	TF = fft(vect)
	
	xf = fftfreq(Nb,T)
	xf = fftshift(xf)
	yplot = fftshift(TF)
	yplot = N.abs(yplot)
	yplot = yplot[Nb/2:]
	xf    = xf[Nb/2:]
	return xf, yplot/yplot.max()

def flat_data(data,dynlow, dynhigh, log):
	""" Returns data where maximum superior than 10^dynhigh will be replaced by 10^dynhigh, inferior than 10^dynlow will be replaced by 10^dynlow"""
	if log:
		mi = 10**dynlow
		ma = 10**dynhigh
		data=N.minimum(N.maximum(data,mi),ma)
		data=N.log10(data)
	else:
		mi = dynlow
		ma = dynhigh
		data=N.minimum(N.maximum(data,mi),ma)
	return data

def psdVoigt(parameters,x):
	"""Define pseudovoigt function"""
	y0 = parameters['y0'].value
	xc = parameters['xc'].value
	A  = parameters['A'].value
	w  = parameters['w'].value
	mu = parameters['mu'].value

	y  = y0 + A * ( mu * (2/N.pi) * (w / (4*(x-xc)**2 + w**2)) + (1 - mu) * (N.sqrt(4*N.log(2)) / (N.sqrt(N.pi) * w)) * N.exp(-(4*N.log(2)/w**2)*(x-xc)**2) )

	return y

def objective(pars,y,x):
	#we will minimize this function
	err =  y - psdVoigt(pars,x)
	return err

def init(data_x,data_y,xc,arbitrary=False):
	""" param = [y0, xc, A, w, mu]
	Je veux que Xc soit la position que l'utilisateur pointe sur l'image pour tracer les profiles"""
	param = Parameters()
	if arbitrary:
		A  = data_y.max()
	else:
		idA=N.where(data_x==xc)[0][0]
		A  = data_y[idA]
	y0 = 20
	#xc = data_x[N.argmax(data_y)]
	w  = 0.5
	mu = 0.5
	param.add('y0', value=y0)
	param.add('xc', value=xc)
	param.add('A', value=A)
	param.add('w', value=w)
	param.add('mu', value=mu, min=0., max=1.)
	return param

def fit(data_x,data_y,xc, arbitrary=False):
	""" return: fitted data y, fitted parameters """
	param_init = init(data_x,data_y,xc,arbitrary)
	if data_x[0] > data_x[-1]:
		data_x = data_x[::-1]
	result = minimize(objective, param_init, args=(data_y,data_x))#changed from lmfit 0.9: One must use result.params for optimized parameters

	x = N.linspace(data_x.min(),data_x.max(),data_x.shape[0])
	y = psdVoigt(result.params,x)

	return result.params, y

def wpercentile(a, q, weights=None):
	"""Compute weighted percentiles *q* [%] of input 1D-array *a*."""

	a = N.asarray(a)
	if a.ndim > 1:
		raise NotImplementedError("implemented on 1D-arrays only")

	if weights is None:
		weights = N.ones_like(a)
	else:
		assert len(weights)==len(a), "incompatible weight and input arrays"
		assert (weights>0).all(), "weights are not always strictly positive"

	isorted = N.argsort(a)
	sa = a[isorted]
	sw = weights[isorted]
	sumw = N.cumsum(sw)                        # Strictly increasing
	scores = 1e2*(sumw - 0.5*sw)/sumw[-1]      # 0-100 score at center of bins

	def interpolate(q):
		i = scores.searchsorted(q)
		if i==0:                        # Below 1st score
			val = sa[0]
		elif i==len(a):                 # Above last score
			val = sa[-1]
		else:                           # Linear score interpolation
			val = (sa[i-1]*(scores[i] - q) + sa[i]*(q - scores[i-1])) / \
				(scores[i] - scores[i-1])
		return val

	out = N.array([ interpolate(qq) for qq in N.atleast_1d(q) ])

	return out.reshape(N.shape(q))      # Same shape as input q


def median_stats(a, weights=None, axis=None, scale=1.4826):
	"""Compute [weighted] median and :func:`nMAD` of array *a* along
	*axis*. Weighted computation is implemented for *axis* = None
	only."""

	if weights is not None:
		if axis is not None:
			raise NotImplementedError("implemented on 1D-arrays only")
		else:
			med  = wpercentile(a, 50., weights=weights)
			nmad = wpercentile(N.absolute(a - med), 50., weights=weights)*scale
	else:
		med = N.median(a, axis=axis)
		if axis is None:
			umed = med                      # Scalar
		else:
			umed = N.expand_dims(med, axis) # Same ndim as a
		nmad = N.median(N.absolute(a - umed), axis=axis) * scale

	return med,nmad

def get_index(arr, val):
	#Get index of val from a 1D array arr
	b = arr.min()
	a = (arr.max() - b)/(arr.size-1)
	x = (val - b)/a
	return int(x)

def get_motors(header):
	#Get a dictionnary like motors {name:position} from a header of an EDF image
	motor_pos  = header['motor_pos'].split()
	motor_mne  = header['motor_mne'].split()
	motor = {}
	for j in range(len(motor_mne)):
		motor[motor_mne[j]] = float(motor_pos[j])
	return motor
	
def get_counters(header):
	#Get a dictionnary like counters {name:position} from a header of an EDF image
	counter_name = header['counter_mne'].split()
	counter_value= header['counter_pos'].split()
	counters= {}
	for i in range(len(counter_name)):
		counters[counter_name[i]] = float(counter_value[i])
	return counters

def select_files_from_list(s, beg, end):
	""" select edf files limited between number 'beg' and 'end' from the list 's' """
	out=[]
	for i in range(len(s)):
		ss = s[i]
		if "-" in ss:
			spliter = "-"
		else:
			spliter = "_"
			
		l = ss.split(spliter)
		l = l[1]
		n = l.split(".")[0]
		if 'g' in n:
			n = n[:-1]
		n = int(n)
		if n in range(beg, end+1):
			out.append(ss)
	#out = sorted(out)
	return out

def get_img_list(edf_list):
	""" Get a list of image numbers from an image_name list """
	out = []
	for i in range(len(edf_list)):
		edf = edf_list[i]
		if edf != None:
			if "-" in edf:
				spliter = "-"
			else:
				spliter = "_"
			l = edf.split(spliter)
			l = l[1]
			n = l.split(".")[0]
			if 'g' in n:
				n = n[:-1]
			n = int(n)
			out.append(n)
		else:
			out.append(None)
	out = N.asarray(out)
	return N.asarray(out)

def get_column_from_table(table,col):
	out = []
	for i in range(len(table)):
		out.append(table[i][col])
	return out
