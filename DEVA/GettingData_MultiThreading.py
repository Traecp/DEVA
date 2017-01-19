#!/usr/bin/python
# -*- coding: utf-8 -*-
import threading
from sys import stdout
import fabio
from os.path import join
from CommonFunctions import *
"""
Getting data from a scan - for 3D visualization
Using threading
"""
class get_scan_data_thread(threading.Thread):
	def __init__(self,i, edf_folder, edf_basename):
		threading.Thread.__init__(self)
		self.Data = ()
		self.order = i 
		self.edf_folder = edf_folder
		self.edf_basename = edf_basename
	def run(self):
		# MAIN_LOCK.acquire()
		edf    = join(self.edf_folder, self.edf_basename)
		data   = fabio.open(edf).data
		header = fabio.open(edf).header
		# stdout.write("Loading %s\n"%self.edf_basename)
		# stdout.flush()
		self.Data = (self.order,data,header)
		# MAIN_LOCK.release()

"""
Getting data from a scan - for 3D visualization
Using Process
"""		
def get_scan_data_process(queue, i, edf_folder, edf_basename):
	edf    = join(edf_folder, edf_basename)
	data   = fabio.open(edf).data
	header = fabio.open(edf).header
	stdout.write("Loading %s\n"%edf_basename)
	stdout.flush()
	DATA_out = (i,data,header)
	queue.put(DATA_out)
"""
Getting data for Pole Figure plotting
Using Threading
"""
class Pole_Figure_load_data(threading.Thread):
	def __init__(self,img_nbr, edf_folder, edf_basename, pole_2theta, azimuthalIntegration, plot_kphi, detector_type="D5"):
		threading.Thread.__init__(self)
		self.Data = []
		self.order = img_nbr
		self.edf_folder = edf_folder
		self.edf_basename = edf_basename
		self.edf = join(self.edf_folder, self.edf_basename)
		self.pole_2theta = pole_2theta
		self.detector_type = detector_type
		self.plot_kphi = plot_kphi
		self.azimuthalIntegration = azimuthalIntegration
		
	def run(self):
		# MAIN_LOCK.acquire()
		img = fabio.open(self.edf)
		this_motor = get_motors(img.header)
		if img.data.size%9600==0:
			data = correct_geometry(img.data, detector_type=self.detector_type)
			data = N.rot90(data)
			img.data = data 
		I,tth,chi = self.azimuthalIntegration.integrate2d(img.data,500,500,unit="2th_deg", method="cython")
		#2theta vs pixel: y = ax + b
		b = tth.min()
		a = (tth.max() - b)/500
		x = (self.pole_2theta - b)/a
		x = int(x)
		intensity = I[:,x-2:x+2].sum(axis=-1)
		intensity = intensity / 4.
		if self.plot_kphi:
			this_phi  = this_motor['kphi']
		else:
			this_phi  = this_motor['phi']
		self.chi_gonio= this_motor["chi"]
		self.omega    = this_motor["eta"]
		# stdout.write("\rLoading %s"%self.edf_basename)
		# stdout.flush()
		self.chi = chi
		self.Data = [self.order, intensity, this_phi]
		# MAIN_LOCK.release()

"""
Getting data for Pole Figure plotting
Using Process
"""
def Pole_Figure_load_data_process(queue, img_nbr, edf_folder, edf_basename, pole_2theta, azimuthalIntegration, plot_kphi, detector_type="D5"):
	edf = join(edf_folder, edf_basename)
	img = fabio.open(edf)
	this_motor = get_motors(img.header)
	if img.data.size%9600==0:
		data = correct_geometry(img.data, detector_type=detector_type)
		data = N.rot90(data)
		img.data = data 
	I,tth,chi = azimuthalIntegration.integrate2d(img.data,500,500,unit="2th_deg", method="cython")
	#2theta vs pixel: y = ax + b
	b = tth.min()
	a = (tth.max() - b)/500
	x = (pole_2theta - b)/a
	x = int(x)
	intensity = I[:,x-2:x+2].sum(axis=-1)
	intensity = intensity / 4.
	if plot_kphi:
		this_phi  = this_motor['kphi']
	else:
		this_phi  = this_motor['phi']
	chi_gonio     = this_motor["chi"]
	omega         = this_motor["eta"]
	# stdout.write("\rLoading %s"%self.edf_basename)
	# stdout.flush()
	DATA_out = [img_nbr, intensity, chi, this_phi, chi_gonio, omega]
	queue.put(DATA_out)