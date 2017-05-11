#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy as np
import operator
from scipy.fftpack import fft, fftfreq, fftshift
from lmfit import Parameters, minimize
	

def UB2B(UB):
	G_star=np.dot(UB.T,UB)
	G = np.linalg.inv(G_star)
	BB= np.sqrt(G)*2*np.pi
	BB=np.nan_to_num(BB)
	a=BB[0,0]
	b=BB[1,1]
	c=BB[2,2]
	gamma = np.arccos(BB[0,1]/a/b)
	gamma = np.degrees(gamma)
	beta = np.arccos(BB[0,2]/a/c)
	beta = np.degrees(beta)
	alpha = np.arccos(BB[1,2]/b/c)
	alpha = np.degrees(alpha)

	cosg = np.cos(np.radians(gamma))
	cosb = np.cos(np.radians(beta))
	cosa = np.cos(np.radians(alpha))
	sing = np.sin(np.radians(gamma))

	#Cell with a//Ox, b in Oxy, z perpendicular to Oxy
	A = np.zeros((3,3))
	A[0,0] = a
	A[0,1] = b * cosg
	A[1,1] = b * sing
	A[0,2] = c * cosb
	A[1,2] = -c*(cosb * cosg - cosa)/sing
	A[2,2] = c/sing * np.sqrt(sing**2 - cosb**2 +2*cosa*cosb*cosg - cosa**2)

	#Reciprocal cell
	B =2*np.pi*np.linalg.inv(A)
	B = B.transpose()
	return B

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
		for i in xrange(len(main_store)):
			try:
				if "_" in main_store[i]:
					e = []
					e.append(main_store[i])
					row=main_store[i].split("_")
					if len(row)>2:
						prefix=""
						for p in xrange(len(row)-1):
							prefix=prefix+row[p]+"_"
					else:
						prefix=row[0]+"_"
					e.append(prefix)
					row = row[-1]
					row = row.split(".")[0]
					if 'g' in row:
						row = row[:-1]
					e.append(int(row))#number
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
	yplot = np.abs(yplot)
	yplot = yplot[Nb/2:]
	xf    = xf[Nb/2:]
	return xf, yplot/yplot.max()

def flat_data(data,dynlow, dynhigh, log):
	""" Returns data where maximum superior than 10^dynhigh will be replaced by 10^dynhigh, inferior than 10^dynlow will be replaced by 10^dynlow"""
	if log:
		mi = 10**dynlow
		ma = 10**dynhigh
		data=np.minimum(np.maximum(data,mi),ma)
		data=np.log10(data)
	else:
		mi = dynlow
		ma = dynhigh
		data=np.minimum(np.maximum(data,mi),ma)
	return data

def psdVoigt(parameters,x):
	"""Define pseudovoigt function"""
	y0 = parameters['y0'].value
	xc = parameters['xc'].value
	A  = parameters['A'].value
	w  = parameters['w'].value
	mu = parameters['mu'].value

	y  = y0 + A * ( mu * (2/np.pi) * (w / (4*(x-xc)**2 + w**2)) + (1 - mu) * (np.sqrt(4*np.log(2)) / (np.sqrt(np.pi) * w)) * np.exp(-(4*np.log(2)/w**2)*(x-xc)**2) )

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
		idA=np.where(data_x==xc)[0][0]
		A  = data_y[idA]
	y0 = 20
	#xc = data_x[np.argmax(data_y)]
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

	x = np.linspace(data_x.min(),data_x.max(),data_x.shape[0])
	y = psdVoigt(result.params,x)

	return result.params, y

def wpercentile(a, q, weights=None):
	"""Compute weighted percentiles *q* [%] of input 1D-array *a*."""

	a = np.asarray(a)
	if a.ndim > 1:
		raise NotImplementedError("implemented on 1D-arrays only")

	if weights is None:
		weights = np.ones_like(a)
	else:
		assert len(weights)==len(a), "incompatible weight and input arrays"
		assert (weights>0).all(), "weights are not always strictly positive"

	isorted = np.argsort(a)
	sa = a[isorted]
	sw = weights[isorted]
	sumw = np.cumsum(sw)                        # Strictly increasing
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

	out = np.array([ interpolate(qq) for qq in np.atleast_1d(q) ])

	return out.reshape(np.shape(q))      # Same shape as input q


def median_stats(a, weights=None, axis=None, scale=1.4826):
	"""Compute [weighted] median and :func:`nMAD` of array *a* along
	*axis*. Weighted computation is implemented for *axis* = None
	only."""

	if weights is not None:
		if axis is not None:
			raise NotImplementedError("implemented on 1D-arrays only")
		else:
			med  = wpercentile(a, 50., weights=weights)
			nmad = wpercentile(np.absolute(a - med), 50., weights=weights)*scale
	else:
		med = np.median(a, axis=axis)
		if axis is None:
			umed = med                      # Scalar
		else:
			umed = np.expand_dims(med, axis) # Same ndim as a
		nmad = np.median(np.absolute(a - umed), axis=axis) * scale

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
	for j in xrange(len(motor_mne)):
		motor[motor_mne[j]] = float(motor_pos[j])
	return motor
	
def get_counters(header):
	#Get a dictionnary like counters {name:position} from a header of an EDF image
	counter_name = header['counter_mne'].split()
	counter_value= header['counter_pos'].split()
	counters= {}
	for i in xrange(len(counter_name)):
		counters[counter_name[i]] = float(counter_value[i])
	return counters

def select_files_from_list(s, beg, end):
	""" select edf files limited between number 'beg' and 'end' from the list 's' """
	out=[]
	for i in xrange(len(s)):
		ss = s[i]
		spliter = "_"
			
		l = ss.split(spliter)
		l = l[-1]
		n = l.split(".")[0]
		if 'g' in n:
			n = n[:-1]
		n = int(n)
		if n in xrange(beg, end+1):
			out.append(ss)
	#out = sorted(out)
	return out
def select_files_from_table(table, prefix, beg, end):
	#Get image name from table with the corresponding prefix, begin and end number
	out = []
	table = sort_table(table)
	N = len(table)
	targetrange = range(beg,end+1)
	for i in xrange(N):
		if table[i][1] == prefix:
			if table[i][2] in targetrange:
				out.append(table[i][0])
	return out
			
			

def get_img_list(edf_list):
	""" Get a list of image numbers from an image_name list """
	out = []
	for i in xrange(len(edf_list)):
		edf = edf_list[i]
		if edf != None:
			spliter = "_"
			l = edf.split(spliter)
			l = l[-1]
			n = l.split(".")[0]
			if 'g' in n:
				n = n[:-1]
			n = int(n)
			out.append(n)
		else:
			out.append(None)
	out = np.array(out)
	return out

def get_column_from_table(table,col):
	out = []
	for i in xrange(len(table)):
		out.append(table[i][col])
	return out

def group_image_from_table(table):
	#Get the prefix into groups
	#return: dict = {'prefix1': [img_num, img_num2, ...], 'prefix2':[img_num1,...]}
	out = {}
	table = sort_table(table)
	prefix = table[0][1]
	N = len(table)
	i = 0
	while i<N-1:
		img = []
		for j in xrange(i,N):
			if table[j][1]==prefix:
				img.append(table[j][2])
			else:
				prefix=table[j][1]
				break
		i=j
		img = np.asarray(img)
		out[prefix]=img
	return out
	



