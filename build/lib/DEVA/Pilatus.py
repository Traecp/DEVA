#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy as N
import fabio

def readPilatus100kGeometry(header):
	geometry = {}
	crys_gonio_names = header[u"CRYSTAL_GONIO_NAMES"].split()
	crys_gonio_values= header[u"CRYSTAL_GONIO_VALUES"].split()
	for i in range(len(crys_gonio_names)):
		geometry[crys_gonio_names[i]] = float(crys_gonio_values[i])
	
	plt_gonio_names = header[u"PLT_GONIO_NAMES"].split()
	plt_gonio_values= header[u"PLT_GONIO_VALUES"].split()
	for i in range(len(plt_gonio_names)):
		geometry[plt_gonio_names[i]] = float(plt_gonio_values[i])
	source_wavelength = header[u"SOURCE_WAVELENGTH"].split()
	geometry["wavelength"] = float(source_wavelength[-1])*1e-10
	beam_position = header[u"PLT_SPATIAL_BEAM_POSITION"].split()
	geometry["beam_position"] = N.array([float(beam_position[0]), float(beam_position[1])])
	return geometry
	
class Pilatus100k_Img():
	def __init__(self, imgpath):
		self.path   = imgpath
		self.data   = None
		self.header = None
		self.prefix = None
		self.img_nbr= None
		self.geometry = {}
		self.readData()
		
	def readData(self):
		self.fabioImg = fabio.open(self.path)
		self.header   = self.fabioImg.header
		self.data     = self.fabioImg.data
		self.prefix   = self.fabioImg.filename.split("_")[0]
		self.img_nbr  = self.fabioImg.filenumber
		self.geometry = readPilatus100kGeometry(self.header)
		