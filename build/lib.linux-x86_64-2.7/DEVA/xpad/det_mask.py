#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# This file is licensed under the CeCILL licence
# the french "GPL" licence agreed by CNRS
# see  http://www.cecilll.info

"""det_mask.py  --  Bm2py.Detector modules  ....
   mask defective or unwished pixels
"""
version =" $Id: det_mask.py,v 1.2 2011/02/01 13:11:56 berar Exp $"

# If the packages has been installed correctly, this should work:
import re
import numpy as np
#------------------------------------------------------------
#from   Bm2py import  colored
#from   Bm2py.Detector import model
#import Bm2py
#------------------------------------------------------------
#   1.00 2011/01/19 new det_mask.py file with det_mask class previously included
#       in det_base.py. Typical usage is an instance for detector and instances for
#       each image needing specific masks
#------------------------------------------------------------
class det_mask :
    """mask defective or unwished pixels
    add mask to detector data, expected entry like
        Detector.detmask = [['r',(x1,y1),(x2,y2),'comments'],
                ]
        Detector.imgmask = [
            ['f','%FILE%/edf/mask'],['r',(x1,y1),(x2,y2),'comments'],
                ]
        to be done : allow invert mask and circles
        """
    dataShape = None
    dataRoiShape = None
    def __init__(self, mask_list, shape, name='DET_MASK', previous=None):
        self.name=name
        if previous == None :
            self.unmasked = np.ones(shape, dtype=np.bool)
        else :
            self.unmasked = previous.get_unmasked().copy()
        if name == 'DET_MASK':
            self.add_masks(mask_list)
        else :
            try :
                fn = open(name)
            except :
                print "'%s' mask file not found, ignored !"%(name)
            else :
                for line in fn.readlines() :
                    if line.startswith('[') :
                        self.add_masks(eval('['+line+']'))
#------------------------------------------------------------
    def add_masks(self, mask_list):
        #print "add_masks  mask_list", mask_list, type(mask_list)
        for mask in mask_list :
            #print "add_masks mask", mask, type(mask), mask[0]
            if mask[0]=='#'   :
                continue
            if mask[0]=='t'   :
                self.add_triangle(mask)
            elif mask[0]=='r' :
                self.add_rectangle(mask)
            elif mask[0]=='p' :
                self.add_rectangle(mask,size=1)
            elif mask[0]=='f' :
                maskFileName = mask[1]
                if maskFileName.startswith('%FILE%') :
                    split = maskFileName.split('/',2)
                    maskFileName = re.sub(split[1],split[2],self.name)
                try :
                    fn = open(maskFileName)
                except :
                    print "'%s' mask file not found, ignored !"%(maskFileName)
                    continue
                for line in fn.readlines() :
                    if line.startswith('[') :
                        self.add_masks(eval('['+line+']'))
            else : raise Bm2_Error, "invalid mask key not(t,r,p,f) : "+ mask[0]
#------------------------------------------------------------
    def add_triangle(self, mask):
        print "add triangular mask : ", mask
        if len(mask) < 4 : raise Bm2_Error, "invalid triangular mask "+repr(mask)
        xmin=self.unmasked.shape[0]-1
        xmax=0
        ymin=self.unmasked.shape[1]-1
        ymax=0
        removed = 0
        for i in range(1,4) :
            if int(mask[i][0]) < xmin : xmin = int(mask[i][0])
            if int(mask[i][0]) > xmax : xmax = int(mask[i][0])
            if int(mask[i][1]) < ymin : ymin = int(mask[i][1])
            if int(mask[i][1]) > ymax : ymax = int(mask[i][1])
        if xmax > self.unmasked.shape[0] :  xmax = self.unmasked.shape[0]-1
        if ymax > self.unmasked.shape[1] :  ymax = self.unmasked.shape[1]-1
        surface=abs((mask[2][0]-mask[1][0])*(mask[3][1]-mask[1][1])-(mask[2][1]-mask[1][1])*(mask[3][0]-mask[1][0]))
        for x in range(xmin,xmax) :
            for y in range(ymin,ymax) :
                try :
                    surf=abs((mask[1][0]-x)*(mask[2][1]-y)-(mask[1][1]-y)*(mask[2][0]-x))+\
                        abs((mask[2][0]-x)*(mask[3][1]-y)-(mask[2][1]-y)*(mask[3][0]-x))+\
                        abs((mask[3][0]-x)*(mask[1][1]-y)-(mask[3][1]-y)*(mask[1][0]-x))
                    if surf <= surface :
                        self.unmasked[x,y]=0
                except : pass
        print "add triangular mask : ", mask, "removing", removed, "pixels, remaining", len(np.nonzero(self.unmasked.flatten())[0])
#------------------------------------------------------------
    def add_rectangle(self, mask, size=2):
        if size != 2 and size != 1 :
            raise Bm2_Error, "invalid size "+repr(size)+"for rectangular mask : "+repr(mask)
        if len(mask) < size+1 : raise Bm2_Error, "invalid rectangular(point) mask "+repr(mask)
        xmin=self.unmasked.shape[0]-1
        xmax=0
        ymin=self.unmasked.shape[1]-1
        ymax=0
        for i in range(1,size+1) :
            if int(mask[i][0]) < xmin : xmin = int(mask[i][0])
            if int(mask[i][0]) > xmax : xmax = int(mask[i][0])
            if int(mask[i][1]) < ymin : ymin = int(mask[i][1])
            if int(mask[i][1]) > ymax : ymax = int(mask[i][1])
        if xmax > self.unmasked.shape[0] :  xmax = self.unmasked.shape[0]-1
        if ymax > self.unmasked.shape[1] :  ymax = self.unmasked.shape[1]-1
        removed = 0
        for x in range(xmin,xmax+1) :
            for y in range(ymin,ymax+1) :
                if self.unmasked[x,y] != 0 :
                    self.unmasked[x,y]=0
                    removed += 1
        if size == 2   : print self.name, "add rectangular mask : ",
        elif size == 1 : print self.name, "add point mask : ",
        print mask, "removing", removed, "pixels, remaining", len(np.nonzero(self.unmasked.flatten())[0])
#------------------------------------------------------------
    def get_unmasked(self):
        return self.unmasked
#------------------------------------------------------------
    def apply_mask(self, data, roi=None):
        #print "apply_mask  data.shape", data.shape, "class", det_mask.dataShape, det_mask.dataRoiShape, "roi", roi
        if det_mask.dataShape == None :
            det_mask.dataShape = data.shape
        if roi != None and  det_mask.dataRoiShape == None :
            det_mask.dataRoiShape = roi.take(data).shape
        if data.shape != det_mask.dataShape and data.shape != det_mask.dataRoiShape :
            print colored('det_mask.apply_mask : dataShape has changed !'), data.shape, det_mask.dataShape, det_mask.dataRoiShape, "roi", roi
            det_mask.dataShape = data.shape
            if roi != None :
                det_mask.dataRoiShape = roi.take(data).shape
        if roi != None :
            data = np.where(roi.take(self.unmasked), roi.take(data), 0)
        else :
            data = np.where(self.unmasked, data, 0)
        return data.copy()
#------------------------------------------------------------
    def hot_mask(self, data, **kwargs) :
        #sigma=5, ratio=5
        for key in kwargs :
            print ('hot_mask %s = %s;\n'%(key, repr(kwargs[key])))
        if 'window' in kwargs :
            window = kwargs['window']
        else :
            window = (3,3)
        if 'minHot' in kwargs :
            minHot = kwargs['minHot']
        else :
            MinHot = 1000
            
#------------------------------------------------------------
