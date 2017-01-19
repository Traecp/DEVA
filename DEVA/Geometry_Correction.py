#!/usr/bin/python
# -*- coding: utf-8 -*-
import pyFAI
from pyFAI.distortion import Distortion

def correct_geometry(data,detector_type="D5"):
	"""Correct image geometry of XPAD detectors: S540 (or D5), S70 and S140
	The default image shape is 960x560 for D5
	The corrected shape:
		D5: 560x960 ==> 578x1153
		S70:560x120 ==> 578x120
	"""
	if detector_type == "D5":
		#8x7 modules, or S540
		detector = pyFAI.detector_factory("d5")
	elif detector_type == "S70":
		#1x7 modules
		detector = pyFAI.detector_factory("imxpad_s70")
	elif detector_type == "S140":
		#2x7 modules
		detector = pyFAI.detector_factory("imxpad_s140")
	else:
		print "Please specify the detector model."
		return None
	#geometry correction
	xpad_corr = Distortion(detector, resize=True)
	data = xpad_corr.correct(data)
	return data