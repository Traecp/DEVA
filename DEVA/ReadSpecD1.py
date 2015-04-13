#!/usr/bin/python
# -*- coding: utf-8 -*-
import os
import numpy as N

class Scan_D1(object):
	""" Object containing one scan
	Reading from old format dataspec recorded by Xpad D1 detector
	"""
	def __init__(self,scanid,command,colnames,header, scan_status, data):
		self.nr = int(scanid)
		self.command = command
		self.colnames= colnames
		self.header = header
		self.motors  = {}
		self.scan_status = scan_status
		self.data = data
		self.count_time = 0
	def ReadData(self):
		for line in self.header:
			if line.startswith("#P0"):
				line = line.split()
				line = line[1:]
				self.motors['del'] = float(line[0])
				self.motors['eta'] = float(line[1])
				self.motors['chi'] = float(line[2])
				self.motors['phi'] = float(line[3])
				self.motors['nu'] = float(line[4])
				self.motors['mu'] = float(line[5])
				#break
			elif line.startswith("#T"):
				line = line.split()
				self.count_time = float(line[1])
			else:
				continue
		
	def __str__(self):
		return self.command
		

class Read_Spec_D1(object):
	""" Read a scan from dataspec file, i.e. kappapsic.DATE - old format used by D1 xpad detector
	Use: scan = Read_Spec_D1(specfile)
	Where specfile is the kappapsic.DATE file
	"""
	def __init__(self, specfile):
		self.full_filename = specfile
		self.filename = os.path.basename(self.full_filename)
		self.scan_list = []
		self.Parse()
		
	def Parse(self):
		
		f=open(self.full_filename,'r+')
		for line in f:
			if line.startswith("#S "):
				print "Parsing ",line
				header = []
				data   = {}
				cmd    = line
				cmd    = cmd.split()
				scanid = int(cmd[1])
				cm     = cmd[2:]
				command= ""
				for i in range(len(cm)):
					command+=cm[i]
					if i<len(cm)-1:
						command+=" "
				#The next 14 lines is the header of this scan
				for i in range(14):
					header.append(f.next().split("\n")[0])
				colnames = f.next().split("\n")[0]
				colnames = colnames.split()
				colnames = colnames[1:]
				for c in range(len(colnames)):
					data[colnames[c]] = []
				
				stop = False
				while not stop:
					l=f.next()
					if l.startswith("#R %d"%scanid) or l.startswith("#C"):
						stop = True
						if l.startswith("#R %d"%scanid):
							scan_status = "OK"
						else:
							scan_status = "ABORTED"
					else:
						stop = False
						if l.startswith("#"):
							continue
						else:
							item = l.split()
							for i in range(len(item)):
								data[colnames[i]].append(float(item[i]))
				
				for k in data.keys():
					data[k] = N.asarray(data[k])
				s = Scan_D1(scanid,command,colnames,header, scan_status, data)
				self.scan_list.append(s)
			else:
				continue
		f.close()
		
	def __str__(self):
		out=""
		out+="Spec file: %s"%self.full_filename
		out+= "\nThere are %d scans recorded"%len(self.scan_list)
		return out
		
