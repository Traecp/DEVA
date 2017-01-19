#!/usr/bin/python
# -*- coding: utf-8 -*-
#Test somme images LaB6 en eta
import fabio
#import matplotlib.pyplot as plt
import numpy as np
from DEVA.Geometry_Correction import correct_geometry
from os.path import isfile,join

PIXEL_SIZE = 130./1000 #milimeter

def image_correction(data):
	#For D5 only, S70 to be added later if needed
	for col in range(1148):
		for row in range(81,84):
			data[row,col] = (data[row-4,col]+data[row+4,col])/2.0
		for row in range(164,167):
			data[row,col] = (data[row-4,col]+data[row+4,col])/2.0
		for row in range(247,250):
			data[row,col] = (data[row-4,col]+data[row+4,col])/2.0
		for row in range(330,333):
			data[row,col] = (data[row-4,col]+data[row+4,col])/2.0
		for row in range(413,416):
			data[row,col] = (data[row-4,col]+data[row+4,col])/2.0
		for row in range(496,499):
			data[row,col] = (data[row-4,col]+data[row+4,col])/2.0
	for col in range(477,560):
		data[482,col] = (data[481,col]+data[483,col])/2.0
	for col in range(882,1001):
		for row in [70,74,109,194,451,476]:
			data[row,col] = (data[row-1,col]+data[row+1,col])/2.0
	for col in range(1028,1147):
		for row in [6,552]:
			data[row,col] = (data[row-1,col]+data[row+1,col])/2.0
	for col in range(0,120):
		data[338,col] = (data[338-1,col]+data[338+1,col])/2.0
	for col in range(147,265):
		data[404,col] = (data[404-1,col]+data[404+1,col])/2.0
	for col in range(735,854):
		data[96,col] = (data[96-1,col]+data[96+1,col])/2.0
	for r in range(data.shape[0]):
		for c in range(1,data.shape[1]-1):
			if (data[r,c] > 20*data[r,c-1]) and (data[r,c] > 20*data[r,c+1]):
				data[r,c] = (data[r,c-1]+data[r,c+1])/2.0
	return data

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
	
def combine_edf(edf_1, edf_2, savename, detector_type="D5"):
	
	if detector_type == "D5":
		DET_SIZE_X = 1153
		DET_SIZE_Y = 578
	elif detector_type == "S70":
		DET_SIZE_X = 120
		DET_SIZE_Y = 578
		
	img_1 = fabio.open(edf_1)
	img_2 = fabio.open(edf_2)
	
	if img_1.data.shape==(960,560) or img_1.data.shape==(120,560):
		d = correct_geometry(img_1.data, detector_type = detector_type)
		img_1.data = d
	if img_2.data.shape==(960,560) or img_2.data.shape==(120,560):
		d = correct_geometry(img_2.data, detector_type = detector_type)
		img_2.data = d
	motor_1 = get_motors(img_1.header)
	motor_2 = get_motors(img_2.header)
	X1 = motor_1['Xdet']
	X2 = motor_2['Xdet']
	Z1 = motor_1['Zdet']
	Z2 = motor_2['Zdet']
	
	deltaX = X2-X1
	deltaZ = Z2-Z1
	pixX   = int(abs(deltaX)/PIXEL_SIZE)
	pixZ   = int(abs(deltaZ)/PIXEL_SIZE)
	mat_decal_Z = np.zeros(shape=(pixZ, DET_SIZE_Y))
	mat_decal_X = np.zeros(shape=(DET_SIZE_X+pixZ,pixX))
	
	if detector_type == "D5":
		data_1 = image_correction(img_1.data)
	data_1 = np.rot90(data_1)
	if deltaZ < 0:
		data_1 = np.vstack((mat_decal_Z,data_1))
	else:
		data_1 = np.vstack((data_1, mat_decal_Z))
	if deltaX >0:
		data_1 = np.hstack((mat_decal_X,data_1))
	else:
		data_1 = np.hstack((data_1, mat_decal_X))
	norm_1 = data_1 == 0.

	if detector_type == "D5":
		data_2 = image_correction(img_2.data)
	data_2 = np.rot90(data_2)
	if deltaZ < 0:
		data_2 = np.vstack((data_2,mat_decal_Z))
	else:
		data_2 = np.vstack((mat_decal_Z, data_2))
	if deltaX > 0:
		data_2 = np.hstack((data_2,mat_decal_X))
	else:
		data_2 = np.hstack((mat_decal_X, data_2))
	norm_2 = data_2 == 0.

	
	data = (data_1+data_2)/2.0
	norm = norm_1 + norm_2
	data[norm] = data[norm] * 2.0
	
	#Saving the final image
	img_1.data = data
	img_1.write(savename, force_type=np.float32)
	print "Image %s saved"%savename
	
	#plt.imshow(np.log10(data_1+1e-6), origin="lower",vmin=0,vmax=3.3)
	#plt.savefig("I1.png")
	#plt.show()

	#plt.imshow(np.log10(data_2+1e-6), origin="lower",vmin=0,vmax=3.3)
	#plt.savefig("I2.png")
	#plt.show()
	
	#plt.imshow(np.log10(data+1e-6),origin="lower",vmin=0,vmax=3.3)
	#plt.savefig("I3.png")
	#plt.show()
	
#list_1 = [92]
#source_folder = "../raw"
#target_folder = "../Somme_EDF"
#edf_prefix = "14Nov26_"
#monitor_ref = 1e6 #counts per second
#monitor_col = "vct3"
#for i in range(len(list_1)):
	#suivante = list_1[i]+2
	#edf_1_basename = edf_prefix+"%04d.edf"%list_1[i]
	#edf_2_basename = edf_prefix+"%04d.edf"%suivante
	#edf_1 = join(source_folder, edf_1_basename)
	#edf_2 = join(source_folder, edf_2_basename)
	#savename = edf_prefix+"Somme-%04d-%04d.edf"%(list_1[i], suivante)
	#savename = join(target_folder, savename)
	#combine_edf(edf_1, edf_2, savename)

