#!/usr/bin/python
# -*- coding: utf-8 -*-
from DEVA.xpad import libXpad as libX

def correct_geometry(data,detector_type="D5"):
	"""Correct image geometry of XPAD D5 detector
	The default image shape is 960x560 """
	if detector_type == "D5":
		detector = libX.Detector()
	elif detector_type == "S70":
		detector = libX.Detector(nModules=1)
	else:
		print "Please specify the detector model."
		return None
	#Ajouter les gaps ## Code de Clement Buton @Soleil
	detector.set_data(array=data)
	detector.adjust_data()
	detector.set_physical_data()
	detector.reshape_pixels()
	data = detector.physical.data
	#data = N.rot90(data)
	return data