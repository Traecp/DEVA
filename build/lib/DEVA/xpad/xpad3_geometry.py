#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This file is licensed under the CeCILL licence
# the french "GPL" licence agreed by CNRS
# see  http://www.cecilll.info
#

"""xpad3_geometry.py [options]  args ...
    map of xpad3 pixel coordinates taking into account tiled geometry
    map of associated flat detector
    ? 1.0x  (07fev12) xpad3_geometry.py handle tiled or untiled geometry
    1.10  9/10/12 allow direct use of y or z rotations
    1.09  multiple geometry ... xpad3_geometry.py
    1.08  release before multiple geometry : xpad3_tiled_geometry.py
    1.07 efficient treatment of noisy pixels to be done to avoid pollution of neighbouring ones
    1.06 17/01/11 overlappedMask activated
    1.05 12/01/11 no more negative surface, self.surfYZ0 restricted in 0...1
    """
__author__ = "J.-F. Berar (CNRS/NEEL)"
__version__ = "1.10"
__date__ = "11/10/2010"

import os, os.path, sys
import time, math, glob, re
import numpy as np

#try :
    #import pylab
#except :
    #import matplotlib as pylab
#import numpy.masked_array as ma
default_float_precision = 3
np.set_printoptions(precision=default_float_precision, suppress=True)

try :
    import DEVA.xpad.det_mask
except :
    sys.path.append('/home/xpad/PROGS/BM2PY/Bm2py/Detectors/XPAD3_GEOMETRY')
    import det_mask

# prgName is set to the sys.argv to recognise geometry 
prgName = os.path.basename(os.path.splitext(sys.argv[0])[0])
refName = os.path.basename(os.path.splitext(__file__)[0])+"_"+__version__+" ("+__date__+")"

print prgName+" : "+refName

#--------------------------------------------------------------
# parameters
#--------------------------------------------------------------
#class Xpad3Geometry :
    #chipsPerModule = 7

#class D1Geometry(Xpad3Geometry)  :
    #modulesNumber = 8
    
geometries = {
    "common" : {
        "chipsPerModule" : 7,
        "pixelYSize"  : 0.130,    # mm
        "pixelZSize"  : 0.130,    # mm
        "xpad3Rows"   : 120,
        "xpad3Cols"   : 80,
        "xpad3DiodeGap" : 0,
        "zDiodeGardRing"  : 0,    # mm overlap apres le dernier pixel
        "diodeThickness"  : 0.500, # mm needed in overlap
        "moduleAngle" : 0,        # degre
        "modulesGeometry": { }    #  { moduleNum : { 'o' : (y,z), 'e' : (y,z)},  }    o : ecart du pixel(0,0), e : pixel(0,79)}
        },
    "D1" : {
        "modulesNumber" : 8,
        "widePixelYSize"  : 2.4615, # les pixels de bord mesurent 320 microns   5*130:650
        "zDiodeGardRing"  : 0.20,   # mm overlap apres le dernier pixel
        "xpad3DiodeGap" : 3,
        "modulePeriodicity" : 14.88,  # mm metrologie
        "moduleAngle" : -7.0,         # degre
        },
    "D5" : {
        "modulesNumber" : 8,
        "widePixelYSize"  : 2.50 ,   # ???
        "zDiodeGardRing"  : 0.20,   # mm overlap apres le dernier pixel
        "xpad3DiodeGap"   : 3,
        "modulePeriodicity" : 14.88,  # mm metrologie
        },
    "imXpad2" : {
        "modulesNumber" : 2,
        "widePixelYSize"  : 2.50 ,   # ???
        "zDiodeGardRing"  : 0.20,   # mm overlap apres le dernier pixel
        "xpad3DiodeGap"   : 3,
        "modulePeriodicity" : 14.88,  # mm metrologie
        },
    "S70" : {
        "modulesNumber" : 1,
        "widePixelYSize"  : 2.50 ,   # ???
        "zDiodeGardRing"  : 0.20,   # mm overlap apres le dernier pixel
        "xpad3DiodeGap"   : 3,
        "modulePeriodicity" : 14.88,  # mm metrologie
        }
    }


#--------------------------------------------------------------
class Xpad3Coordinates ( object ) :
    """ coordinates transformation for tiled xpad3 detector """
    validMaskArray = None
    diodeMaskArray = None
    correctedMaskArray = None
    numberUnMasked = 0
    strategy = { 'min' : 0.50 }

    def __init__ (self, geometryName, rectangular=False, rawCenter=None, distance=0, yRotation=0,
        zRotation=0, unMaskOverlapped=False, maskStrategy='auto', options={}, **kwargs ) :
        detectorGeometry = geometries["common"]
        self.detectorGeometry = geometries["common"]
        if geometryName in geometries.keys() :
            for k in geometries[geometryName].keys() :
                detectorGeometry[k] = geometries[geometryName][k]
        else :
            raise "eval file geometry "+geometryName+".geo to be implemented"
            
        self.angle = detectorGeometry["moduleAngle"]
        #self["moduleAngle"] = detectorGeometry["moduleAngle"]
        #print dir(self)
        #print self.angle
        #print self.moduleAngle
        #sys.exit()
        angleRadian = math.radians(self.angle)
        cosa = math.cos(angleRadian)
        sina = math.sin(angleRadian)
        self.periodicity = detectorGeometry["modulePeriodicity"]
        self.overlap = detectorGeometry["xpad3Rows"] - self.periodicity*cosa/detectorGeometry["pixelZSize"] 
        self.gap =  self.periodicity * sina
        print "\nXpad3Coordinates, geometry : periodicity=%.2f mm angle=%.1f deg gap=%.2f mm ovelap=%.2f pixels"%(self.periodicity, self.angle, self.gap, self.overlap)
        self.mappedPixelZSize = detectorGeometry["pixelZSize"]
        if rectangular : self.mappedPixelZSize *= cosa

        self.geometry = detectorGeometry["modulesGeometry"]
        self.modulesNumber = detectorGeometry["modulesNumber"]
        self.chipsPerModule = detectorGeometry["chipsPerModule"]
        
        self.rawShape = (self.chipsPerModule*detectorGeometry["xpad3Cols"], self.modulesNumber*detectorGeometry["xpad3Rows"])
        self.diodeCorrectedShape = (self.rawShape[0] + detectorGeometry["xpad3DiodeGap"] * (self.chipsPerModule-1), self.rawShape[1] )
        #self.diodeShape = (self.rawShape[0], self.rawShape[1] + xpad3DiodeGap*(self.chipsPerModule-1))
        self.maskStrategy = maskStrategy

        self.distance = distance
# A METTRE AU POINT !!!
        self.options = options
        print options, type(options), "self.options", type(self.options)
        testPrint = False
        #testPrint = True
        #self.options={}
        #self.options['debug'] = True
        
        self.flatShape = None
        if rawCenter == None :
            self.rawCenter = (self.rawShape[0]/2, self.rawShape[1]/2)
        else :
            self.rawCenter = rawCenter
        self.yRotation = yRotation
        self.zRotation = zRotation
        print " image %d*%d shape"%(self.chipsPerModule, self.modulesNumber), self.rawShape, "seen at %.1fmm, rawCenter pixel coords %.1f %.1f"%(self.distance,self.rawCenter[0],self.rawCenter[1])," rotation errors (Y %.2f, Z%.2f deg)"%(self.yRotation,self.zRotation)
        self.moduleCoordinates = None
        self.flatCoordinates = None
        self.geometricArrays = None
        self.unMaskOverlapped = unMaskOverlapped
        self.overlappedRowMask = [] 
        # tous les pixels meme partiellement recouverts
        overlapPixels = int(self.overlap + detectorGeometry["zDiodeGardRing"]/detectorGeometry["pixelZSize"] - detectorGeometry["diodeThickness"]*sina/detectorGeometry["pixelZSize"]) + 1  
        # elargissement utile (JF 9/10/12)
        overlapPixels = overlapPixels + 1
        # build basic geometry with 6 points per module : pixel[0,0], pixel[80,0], pixel[0,119],
        oneModuleCoordinates = np.zeros((3,6), dtype=float)  #coordinates of right/left chip corners and guardring
        #oneModuleCoordinates[1,1:6:2] = (self.diodeCorrectedShape[0]-1) * pixelYSize 
        oneModuleCoordinates[1,1:6:2] = (self.diodeCorrectedShape[0]-1) * detectorGeometry["pixelYSize"] + (self.chipsPerModule-1) *  (2*detectorGeometry["widePixelYSize"] - (2+detectorGeometry["xpad3DiodeGap"])) * detectorGeometry["pixelYSize"]
        oneModuleCoordinates[2,2:] = (detectorGeometry["xpad3Rows"]-1) * detectorGeometry["pixelZSize"]
        oneModuleCoordinates[2,4:] += detectorGeometry["zDiodeGardRing"]
        if testPrint :
            print  oneModuleCoordinates.shape, "oneModuleCoordinates : [Y,Z]  [x,y,z] mm,  0<=Y<7*80+6*3=578  0<=Z<120"
            print "    pixel[0,0]", oneModuleCoordinates[:,0], "          pixel[577,0]", oneModuleCoordinates[:,1]
            print "    pixel[0,119]", oneModuleCoordinates[:,2], "    pixel[577,119]", oneModuleCoordinates[:,3]
            print "    guardRing", oneModuleCoordinates[:,4] ,"       ", oneModuleCoordinates[:,5]
        self.moduleOrigPixel = np.empty((2, self.modulesNumber), dtype=float)
        self.flatModuleOrigPixel = np.zeros((2, self.modulesNumber), dtype=float)
        moduleCoordinates = np.empty((3,6, self.modulesNumber), dtype=float)
        for imodule in range(self.modulesNumber) :
            coordinates = oneModuleCoordinates.copy()
            if self.geometry.has_key(imodule) :
                eym = (self.geometry[imodule]['o'][0] + self.geometry[imodule]['e'][0]) / 2.0
                ezm = (self.geometry[imodule]['o'][1] + self.geometry[imodule]['e'][1]) / 2.0
                dy = (self.geometry[imodule]['e'][0] - self.geometry[imodule]['o'][0]) / 2.0
                dz = (self.geometry[imodule]['e'][1] - self.geometry[imodule]['o'][1]) / 2.0
                eyo = eym - dy
                ezo = ezm - dz
                errAngle = 2.0 * dz / detectorGeometry["xpad3Cols"]
                Y = coordinates[1,:,:].copy()
                Z = coordinates[2,:,:].copy()
                coordinates[1,:] = Y + Z * errAngle + eyo
                coordinates[2,:] = -Y * errAngle + Z + ezo
            moduleCoordinates[0,:,imodule] = coordinates[0,:] * cosa - coordinates[2,:] * sina
            moduleCoordinates[1,:,imodule] = coordinates[1,:]
            moduleCoordinates[2,:,imodule] = coordinates[0,:] * sina + coordinates[2,:] * cosa + imodule * (detectorGeometry["xpad3Rows"] - self.overlap) * detectorGeometry["pixelZSize"]
            self.moduleOrigPixel[0,imodule] = 0
            self.moduleOrigPixel[1,imodule] = imodule * detectorGeometry["xpad3Rows"]
            if self.angle < 0 and imodule != self.modulesNumber-1 :
                ie=int(self.moduleOrigPixel[1,imodule]+detectorGeometry["xpad3Rows"])
                self.overlappedRowMask.extend(range(ie-overlapPixels, ie))
                if testPrint :
                    print imodule,  "self.overlappedRowMask", self.overlappedRowMask, ie-overlapPixels, ie
            elif self.angle > 0 and not self.unMaskOverlapped :
                raise "overlappedRowMask not implemented for alpha >0 "
        self.moduleCoordinates = np.empty_like(moduleCoordinates)
        cosy = math.cos(math.radians(self.yRotation))
        siny = math.sin(math.radians(self.yRotation))
        cosz = math.cos(math.radians(self.zRotation))
        sinz = math.sin(math.radians(self.zRotation))
        self.moduleCoordinates[0,:,:] = (moduleCoordinates[0,:,:] * cosy - moduleCoordinates[2,:,:] * siny) * cosz - moduleCoordinates[1,:,:] * sinz 
        self.moduleCoordinates[1,:,:] =  (moduleCoordinates[0,:,:] * cosy - moduleCoordinates[2,:,:] * siny)  * sinz + moduleCoordinates[1,:,:] * cosz
        self.moduleCoordinates[2,:,:] = moduleCoordinates[0,:,:] * siny + moduleCoordinates[2,:,:] * cosy
        self.moduleDxyzDYZ = np.empty((3,2,self.modulesNumber))
        self.moduleDxyzDYZ[:,0,:] =  (self.moduleCoordinates[:,1,:]-self.moduleCoordinates[:,0,:]) / self.diodeCorrectedShape[0]
        self.moduleDxyzDYZ[:,1,:] =  (self.moduleCoordinates[:,2,:]-self.moduleCoordinates[:,0,:]) / detectorGeometry["xpad3Rows"]
        delta = self.moduleDxyzDYZ[1,0,:] * self.moduleDxyzDYZ[2,1,:] - self.moduleDxyzDYZ[2,0,:] * self.moduleDxyzDYZ[1,1,:]
        #print "delta", delta.shape, delta
        self.moduleDYZDyz = np.empty((2,2,self.modulesNumber))
        self.moduleDYZDyz[0,0,:] = self.moduleDxyzDYZ[2,1,:]/delta
        self.moduleDYZDyz[1,0,:] = self.moduleDxyzDYZ[1,1,:]/delta
        self.moduleDYZDyz[0,1,:] = self.moduleDxyzDYZ[2,0,:]/delta
        self.moduleDYZDyz[1,1,:] = self.moduleDxyzDYZ[1,0,:]/delta
        if testPrint :
            print  "moduleCoordinates[xyz,pos,module]", self.moduleCoordinates.shape
            for i in range(self.modulesNumber) :
                print  "module", i, 'pixel[Y,Z] [x,y,z]'
                print "    pixel[0,0]", self.moduleCoordinates[:,0,i], "    pixel[577,0]", self.moduleCoordinates[:,1,i]
                print "    pixel[0,119]", self.moduleCoordinates[:,2,i], "  pixel[577,119]", self.moduleCoordinates[:,3,i]
                print "    guardRing", self.moduleCoordinates[:,4,i] ,"     ", self.moduleCoordinates[:,5,i]
                print "    d_xyz_dYZ", self.moduleDxyzDYZ[:,0,i], self.moduleDxyzDYZ[:,1,i]
                print "    d_YZ_d_yz", self.moduleDYZDyz[:,0,i], self.moduleDYZDyz[:,1,i]
        #self.surfaceRatio = 1.0/cosa/cosz/cosy
        #print "self.surfaceRatio", self.surfaceRatio, 1/cosa, 1/cosa/cosa
        
    def tiled_module_coordinates(self, module=None) :
        "return modules corners coordinates"
        if module == None :
            return self.moduleCoordinates[:,:4,:]
        else :
            return self.moduleCoordinates[:,:4,module]

    def get_diode_pixel_coordinates(self, raw, exactGap=False) :
        """ return the diode corrected coordinates of a pixel located in the raw image"""
        y = int(raw[0])
        ny = y / self.detectorGeometry["xpad3Cols"]
        if not exactGap :
            return raw[0] + ny * self.detectorGeometry["xpad3DiodeGap"], raw[1]
        return raw[0] + ny * (2*self.detectorGeometry["widePixelYSize"]-2 ) , raw[1]  
            

    def get_flat_pixel_coordinates(self, raw, diodeCorrect=True, exactGap=False) :
        """ return the new coordinates of a pixel located in the raw image"""
        if diodeCorrect :
            raw = self.get_diode_pixel_coordinates(raw)
        y = int(raw[0])
        z = int(raw[1])
        if self.flatCoordinates == None :
            self.flat_coordinates()
        if y < 0 or y > self.rawShape[0] or z < 0 or z > self.rawShape[1] :
            raise NotImplementedError("pixel out of raw shape, but it remains possible")
        ny = (self.flatCoordinates[0,y,z] + (self.flatCoordinates[0,y+1,z]-self.flatCoordinates[0,y,z])*(raw[0]-y) - self.flatCoordinates[0,0,0]) / self.mappedPixelZSize
        nz = (self.flatCoordinates[1,y,z] + (self.flatCoordinates[1,y,z+1]-self.flatCoordinates[1,y,z])*(raw[1]-z) - self.flatCoordinates[1,0,0]) / self.detectorGeometry["pixelYSize"]
        return (ny,nz)

    def flat_coordinates(self, exactGap=True) :
        """ return the coordinates of pixels in external reference
        if distance is 0, simple projection like from infinite
        else conic projection from source set at distance on pixel(0,0) plane"""
        print "flat_coordinates", self.distance
        diodeCorrectedCenter = self.get_diode_pixel_coordinates(self.rawCenter)
        print  "diodeCorrectedCenter", diodeCorrectedCenter,  diodeCorrectedCenter[0]
        self.flatCoordinates = np.empty((2, self.diodeCorrectedShape[0], self.diodeCorrectedShape[1]), dtype=float)
        xyzCoordinates = np.empty((3, self.diodeCorrectedShape[0], self.diodeCorrectedShape[1]), dtype=float)
        if not exactGap or exactGap :
            YZ = np.mgrid[0:self.diodeCorrectedShape[0],0:self.detectorGeometry["xpad3Rows"]]
        for module in range(self.modulesNumber) :
            xyzCoordinates[0,:,module*self.detectorGeometry["xpad3Rows"]:(module+1)*self.detectorGeometry["xpad3Rows"]] = self.moduleCoordinates[0,0,module] + self.moduleDxyzDYZ[0,0,module] * YZ[0,:,:] + self.moduleDxyzDYZ[0,1,module] * YZ[1,:,:]
            xyzCoordinates[1,:,module*self.detectorGeometry["xpad3Rows"]:(module+1)*self.detectorGeometry["xpad3Rows"]] = self.moduleCoordinates[1,0,module] + self.moduleDxyzDYZ[1,0,module] * YZ[0,:,:] + self.moduleDxyzDYZ[1,1,module] * YZ[1,:,:]
            xyzCoordinates[2,:,module*self.detectorGeometry["xpad3Rows"]:(module+1)*self.detectorGeometry["xpad3Rows"]] = self.moduleCoordinates[2,0,module] + self.moduleDxyzDYZ[2,0,module] * YZ[0,:,:] + self.moduleDxyzDYZ[2,1,module] * YZ[1,:,:]
        print "xyzCoordinates", xyzCoordinates.shape, xyzCoordinates[:,0,0], xyzCoordinates[:,-1,0], xyzCoordinates[:,0,-1], xyzCoordinates[:,-1,-1]
        if self.distance == 0 or self.distance > 5000 :
            self.flatCoordinates = xyzCoordinates[1:,...]
        else :
            xyzCoordinates[0,...] += self.distance
            xyzCoordinates[1,...] -= diodeCorrectedCenter[0] * self.detectorGeometry["pixelYSize"]
            xyzCoordinates[2,...] -= diodeCorrectedCenter[1] * self.mappedPixelZSize
            print "xyzCoordinates", xyzCoordinates[:,0,0], xyzCoordinates[:,-1,0], xyzCoordinates[:,0,-1], xyzCoordinates[:,-1,-1]
            self.flatCoordinates[0,:,:] = self.distance * xyzCoordinates[1,:,:] / xyzCoordinates[0,:,:] + diodeCorrectedCenter[0] * self.detectorGeometry["pixelYSize"]
            self.flatCoordinates[1,:,:] = self.distance * xyzCoordinates[2,:,:] / xyzCoordinates[0,:,:]+ diodeCorrectedCenter[1] * self.mappedPixelZSize
        np.save("xyzCoordinates.npy",xyzCoordinates)
        try :
            print  "flatCoordinates : center Pixel (y=%.0f, z=%.0f)"%self.rawCenter ," , corners :",self.flatCoordinates.shape,  self.flatCoordinates[:,0,0], self.flatCoordinates[:,-1,0], self.flatCoordinates[:,0,-1], self.flatCoordinates[:,-1,-1]
        except :
            print  "EXCEPT flatCoordinates : center Pixel (y=?.0f, z=?.0f)",self.rawCenter ," , corners :",self.flatCoordinates.shape,  self.flatCoordinates[:,0,0], self.flatCoordinates[:,-1,0], self.flatCoordinates[:,0,-1], self.flatCoordinates[:,-1,-1]
            
        ySize = np.amax(self.flatCoordinates[0,:,:])- np.amin(self.flatCoordinates[0,:,:])
        zSize = np.amax(self.flatCoordinates[1,:,:])- np.amin(self.flatCoordinates[1,:,:])
        self.flatShape = (int(ySize/self.detectorGeometry["pixelYSize"]+2), int(zSize/self.mappedPixelZSize+2))
        self.flatCenter = self.get_flat_pixel_coordinates(diodeCorrectedCenter, diodeCorrect=False)
        print "mapped_geometry", self.flatShape,"distance %.1fmm,"%(self.distance), "new center pixel (y=%.2f, z=%.2f)"%self.flatCenter, "new center pixel exactGap(y=%.2f, z=%.2f)"%self.get_flat_pixel_coordinates(diodeCorrectedCenter, diodeCorrect=False,exactGap=True)
        return self.flatCoordinates

    def map_diodes(self, rawValues) :
        " return new count array mapped for wide pixel between diodes"
        mapValues = np.zeros(self.diodeCorrectedShape, dtype = float)
        print "map_diodes", rawValues.shape, mapValues.shape
        mapValues[0,:] = rawValues[0,:]
        mapValues[-1,:] = rawValues[-1,:]
        for i in range(0,self.chipsPerModule) :
            io = i*self.detectorGeometry["xpad3Cols"]
            jo = i*(self.detectorGeometry["xpad3Cols"]+self.detectorGeometry["xpad3DiodeGap"])
            mapValues[jo+1:jo+self.detectorGeometry["xpad3Cols"]-1,:] = rawValues[io+1:io+self.detectorGeometry["xpad3Cols"]-1,:]
            if i != 0 :
                v = rawValues[io-1,:]/self.detectorGeometry["widePixelYSize"]
                w = rawValues[io,:]/self.detectorGeometry["widePixelYSize"]
                vv = (v + w)/2.0
                mapValues[jo-4,:] = v
                mapValues[jo-3,:] = v
                mapValues[jo-1,:] = w
                mapValues[jo,:] = w
                mapValues[jo-2,:] = vv    # STRATEGY
                if self.diodeMaskArray != None :
                    ind = np.where(self.diodeMaskArray[jo-2,:] == True)[0]
                    if len(ind) :
                        v_ind = np.where(self.validMaskArray[io-1,:]==False)
                        if len(v_ind) :
                            mapValues[jo-2,:v_ind] = w
                        w_ind = np.where(self.validMaskArray[io,:]==False)
                        if len(w_ind) :
                            mapValues[jo-2,:w_ind] = v
                    print v_ind, w_ind, "multi pixels corrected for mask effect", i,    
        return  mapValues

    def mapped_geometric_arrays(self) :
        " geometrical arrays for projection"
        if self.flatCoordinates == None :
            self.flat_coordinates()
        self.geometricArrays = True
        #print "mapped_geometric_arrays", self.flatCoordinates.shape
        # each pixel can be split on four, due to the projection they are smaller
        yPixel = (self.flatCoordinates[0,:,:] - self.flatCoordinates[0,0,0]) / self.detectorGeometry["pixelYSize"]
        zPixel = (self.flatCoordinates[1,:,:] - self.flatCoordinates[1,0,0]) / self.mappedPixelZSize
        deltaY = np.empty(self.diodeCorrectedShape)
        deltaY[0:-1,:] = yPixel[0:-1,:] - yPixel[1:,:]
        deltaY[-1,:] = deltaY[-2,:]
        try:
            if self.options.debug :
                print "deltaY", deltaY.mean(),  deltaY.std(),  deltaY.min(),  deltaY.max()
        except : pass
        self.surfY0 = np.where(yPixel != np.ceil(yPixel), (yPixel - np.ceil(yPixel)) / deltaY, 1.0)
        self.surfY0 = np.where(self.surfY0 <= 1.0, self.surfY0, 1.0)
        deltaZ = np.empty(self.diodeCorrectedShape)
        deltaZ[:,0:-1] = zPixel[:,0:-1] - zPixel[:,1:]
        deltaZ[:,-1] = deltaZ[:,-2]
        try:
            if self.options.debug :
                print "deltaZ", deltaZ.mean(),  deltaZ.std(),  deltaZ.min(),  deltaZ.max()
        except : pass
        self.surfZ0 = np.where(zPixel != np.ceil(zPixel), (zPixel - np.ceil(zPixel)) / deltaZ, 1.0)
        self.surfZ0 = np.where(self.surfZ0 <= 1.0, self.surfZ0, 1.0)
        self.yIndex = yPixel.astype(int)
        self.zIndex = zPixel.astype(int)
        self.newCountsShape = (self.yIndex.max()+2, self.zIndex.max()+2)
        print "mapped_geometric_arrays yIndex",self.yIndex.min() , self.yIndex.max(), "zIndex", self.zIndex.min(), self.zIndex.max(), "newCountsShape", self.newCountsShape 

    def map_geometry(self, rawValues) :
        """ return new count array mapped for tiled geometry
        each count is mapped into the four mapping pixels"""
        if self.geometricArrays == None :
            self.mapped_geometric_arrays()
        if self.diodeMaskArray == None :
            self.diode_corrected_mask()
        newCounts = np.zeros(self.newCountsShape, dtype=float).flatten()
        newSurfYZ = np.zeros(self.newCountsShape, dtype=float).flatten()
        newMask   = np.zeros(self.newCountsShape, dtype=float).flatten()
        print "map_geometry", rawValues.shape, self.masked_stats(rawValues, string=''), self.newCountsShape
        surfY1 = 1 - self.surfY0
        surfZ1 = 1 - self.surfZ0

        print "mapped_geometric_arrays yIndex",self.yIndex.min() , self.yIndex.max(), "zIndex", self.zIndex.min(), self.zIndex.max(), "newCountsShape", self.newCountsShape

        SY=[self.surfY0, 1.0 - self.surfY0]
        SZ=[self.surfZ0, 1.0 - self.surfZ0]
        
        # as array operation do not allow reentrance in +=
        # loop has to be performed to work on non close pixels which may overlap
        nstep = 5  # to ensure that same pos in never encountered

        for sy in range(2) :
            for sz in range(2) :
                indices = (self.yIndex + sy) * self.newCountsShape[1] + (self.zIndex + sz)               
                surfYZ = SY[sy] * SZ[sz]
                for y in range(nstep) :
                    for z in range(nstep) :
                        ind = indices[y::nstep,z::nstep].flatten()
                        surf = surfYZ[y::nstep,z::nstep]
                        mask = np.where(self.diodeMaskArray[y::nstep,z::nstep]==True, surf, 0.0)
                        maskedCnt = mask * rawValues[y::nstep,z::nstep]
                        newCounts[ind] +=  maskedCnt.flatten() 
                        newSurfYZ[ind] +=  surf.flatten()
                        newMask[ind] +=  mask.flatten()
                        
        #indexDBG = '[103:105,333:344]', '[103:105,445:456]'
        #newSurfYZ.shape = self.newCountsShape
        #newMask.shape = self.newCountsShape
        #newCounts.shape = self.newCountsShape
        #for DBG in indexDBG :
            #for KBG in 'newMask', 'newSurfYZ', 'newCounts' :
                ##eval(KBG+'.shape = self.newCountsShape')
                #print KBG+DBG, "  ", eval(KBG+DBG)*1e3
                ##eval(KBG+'.shape = '+KBG+'.size')
        #newSurfYZ.shape = (newSurfYZ.size)
        #newMask.shape = (newMask.size)
        #newCounts.shape = (newCounts.size)
            
        #apply mask and related corrections ....
        self.correctedMaskArray = np.ones(self.newCountsShape, dtype=bool).flatten()
        ind = np.where(newMask < self.strategy['min'])[0]
        if len(ind) > 0 :
            self.correctedMaskArray[ind] = False
            newCounts[ind] = 0.0
        # surf supposée ok entre 0.95, 1.1 valeurs empiriques
        # correct for surface, pixels with masked neighbouring but avoid in case of overlapped
        ind = np.where(np.logical_and(newMask > 0, newMask < newSurfYZ))[0]
        if len(ind) > 0 :
            # newSurfYZ can be false due to overlap
            # then we replace it by the mean of surrounding ones if valuable
            addedSurf = np.zeros((len(ind)), dtype=float)
            summedSurf = np.zeros((len(ind)), dtype=float)
            indY = ind / rawValues.shape[1]
            indZ = ind % rawValues.shape[1]
            for i in range(-10,11) :
                indZZ = indZ + i
                jnd = np.where(np.logical_and(indZZ>=0, indZZ<self.newCountsShape[1]))[0]
                knd = indY[jnd]*self.newCountsShape[1]+indZZ[jnd]
                lnd = np.where(np.logical_and(newSurfYZ[knd]>0.95, newSurfYZ[knd]<1.1))[0]
                knd_ok = knd[lnd]
                addedSurf[lnd] += 1.0
                summedSurf[lnd] += newSurfYZ[knd_ok]
            #print "addedSurf", addedSurf , addedSurf.size , addedSurf.mean(), len(np.where(addedSurf>0)[0])
            jnd = np.where(addedSurf>0)[0]
            correctedSurf = np.zeros((len(ind)), dtype=float)
            correctedSurf[jnd] = summedSurf[jnd]/addedSurf[jnd]
            newCounts[ind] = newCounts[ind]*correctedSurf/newMask[ind]
            #newSurfYZ[ind] = correctedSurf
        self.correctedMaskArray.shape = self.newCountsShape
        newSurfYZ.shape = self.newCountsShape
        newMask.shape = self.newCountsShape
        newCounts.shape = self.newCountsShape
        #for DBG in indexDBG :
            #for KBG in 'newMask', 'newSurfYZ', 'newCounts' :
                #print KBG+DBG, "  ", eval(KBG+DBG)*1e3
        if options.debug :
            print "newCounts  %10.2f  %10.2f  %10.2f "%( newCounts.mean(), newCounts.min(), newCounts.max())            
            np.savetxt('xpad3_geometry_hsurf.dat', newSurfYZ.sum(axis=-1) , fmt="%12.6G")
            np.savetxt('xpad3_geometry_hmask.dat', newMask.sum(axis=-1) , fmt="%12.6G")
            np.savetxt('xpad3_geometry_vsurf.dat', newSurfYZ.sum(axis=0) , fmt="%12.6G")
            np.savetxt('xpad3_geometry_vmask.dat', newMask.sum(axis=0) , fmt="%12.6G")
        if options.show :
            fig = pylab.figure(figsize=(10,6))
            ax = pylab.subplot(121)  #row, col, num
            pylab.title('newSurfYZ '+ repr(self.newCountsShape) + '?') # repr(values.shape))
            mean = newSurfYZ.mean()
            img = ax.imshow(newSurfYZ.T, vmin=0.92, vmax=1.05, origin='lower')
            ax = pylab.subplot(122)  #left, bottom, width, height
            img = ax.imshow(newMask.T, vmin=0.92, vmax=1.05, origin='lower')
            pylab.colorbar(img, ax=ax, orientation='horizontal', fraction=0.03, pad=0.04, format='%.1f' )
        newCounts.shape = self.newCountsShape
        return newCounts[0:self.flatShape[0],0:self.flatShape[1]]

    def flatten_counts(self, countsArray, diodeCorrection=True, geometryCorrection=True, normalize = "") :
        """ return new count array mapped for wide pixel between diodes and tiled geometry,
        it call if needed the coordinates calculation"""
        if countsArray.shape == self.rawShape and   diodeCorrection :
            countsArray = self.map_diodes(countsArray)
            if self.unMaskOverlapped :
                print "self.unMaskOverlapped", self.unMaskOverlapped
            else :
                print "flatten_counts : unMaskOverlapped -> pixels set to zero \noverlappedRowMask", self.overlappedRowMask
                countsArray[:,self.overlappedRowMask] = 0
                
        if normalize != "" :
            coef = XpadSpec(normalize, monitor=options.monitor).monitor_value(os.path.basename(self.imageName), show=options.debug)
            if type(coef) != type({}) :
                coef = float(coef)
                if coef != 0 :
                    countsArray /= coef
                    print "COEF '%s'"%os.path.basename(self.imageName), "image NORMALIZE by", coef, "from file", normalize
                else :
                    print "COEF '%s'"%os.path.basename(self.imageName), "image not NORMALIZE,coef is Null ", coef
            else :
                print "COEF '%s'"%os.path.basename(self.imageName), "image not NORMALIZE", coef
        if not geometryCorrection :
            return countsArray
        return self.map_geometry(countsArray)


    def load(self, arg, noisyMask=None) :
        "ensure that image are read in the good order"        
        from StringIO import StringIO   # StringIO behaves like a file object
        filename, extension = os.path.splitext(arg)
        if extension == ".gz" or ( not os.path.exists(arg) and os.path.exists(arg+'.gz')) :
            import gzip
            if extension == ".gz" :
                content = StringIO(gzip.open(arg, 'rb').read())
                filename, extension = os.path.splitext(filename)
            else :
                content = StringIO(gzip.open(arg+'.gz', 'rb').read())
        else :
            content = StringIO(open(arg, 'rb').read())
        self.imageName = os.path.basename(filename)
        if extension == ".dat" or extension == ".txt" :
            values = np.loadtxt(content, dtype=int)[::-1,:].T
        elif  arg.endswith(".edf")  :
            splitter=re.compile("^\s*(?P<key>\S+)\s*=\s*(?P<val>.+)$")
            endspace=re.compile("\s+$")
            semi=re.compile(";$")
            line = content.readline()
            while line[0] != "}" :
                ok=splitter.search(line)
                if ok :
                    val=ok.group('val')
                    if endspace.search(val) : val=endspace.sub('',val)
                    if semi.search(val)      : val=semi.sub('',val)
                    if endspace.search(val) : val=endspace.sub('',val)
                    if ok.group('key') == 'ByteOrder' :
                        if not re.match('LowByteFirst',val) :
                            raise prgName, 'can only handle LowByteFirst argument(s) not "'+line+'"'
                    elif ok.group('key') == 'DataType' :
                        if not (re.match('SignedInteger',val) or re.match('FloatValue',val)) :
                            raise prgName, 'can only handle SignedInteger or FloatValue argument(s) not "'+line+'"'
                        else : datatype=val
                    elif ok.group('key') == 'Dim_1' :
                        dimy=int(val)
                    elif ok.group('key') == 'Dim_2' :
                        dimx=int(val)
                    #elif ok.group('key') == 'Exposure_sec' :
                        #self.exposure=float(val)
                    #elif ok.group('key') == 'Image' :
                        #self.pict_number=int(val)
                line = content.readline()
            print "fit2d, dimx, dimy", dimx, dimy
            if datatype == 'SignedInteger'  :
                values=np.fromstring(content.read(4*dimx*dimy), dtype=np.int32, count=dimx*dimy).reshape((dimx,-1)).T
            elif datatype == 'FloatValue' :
                values=np.fromstring(content.read(4*dimx*dimy), dtype=np.float32, count=dimx*dimy).reshape((dimx,-1)).T
            else :
                raise prgName, "No datatype specified in EDF file"
            # edf_write(values, arg, ext='EDF.edf' )
        elif  extension == ".npy" :
            values=np.load(content)
        else :
            values = np.loadtxt(content, dtype=int)[::-1,:].T
        if noisyMask :
            values = values & noisyMask
            print arg, 'noisymask 0x%-0X applied to values'%(options.noisymask)
        return values


    def raw_pixel_to_xyz_coordinates(self, yPixel, zPixel=None, withModule=False) :
        " return xyz coordinates of a pixel, no approximation "
        if zPixel == None :
            zPixel = yPixel[1]
            yPixel = yPixel[0]
        pixelModule = int(zPixel/xpad3Rows)
        zPixelModule = math.fmod(zPixel, xpad3Rows)
        pixelChip = int(yPixel/xpad3Cols)
        yPixelChip =  yPixel + pixelChip * xpad3DiodeGap
        print "self.moduleCoordinates[:,:,pixelModule]", self.moduleCoordinates[:,:,pixelModule]
        xyzCoordinates = self.moduleCoordinates[:,0,pixelModule] + self.moduleDxyzDYZ[:,0,pixelModule] * yPixelChip + self.moduleDxyzDYZ[:,1,pixelModule] * zPixelModule
        if testPrint :
            print "raw_pixel_to_xyz_coordinates(y=%.2f, z=%.2f)"%(yPixel, zPixel),"in module %d, zp=%.2f"%(pixelModule, zPixelModule), xyzCoordinates
        if not withModule :
            return xyzCoordinates
        else :
            return xyzCoordinates, pixelModule

    def fictious(self, rad=True) :
        "return a fictious value build for distance and center"
        intensity = 1.0e3
        diodesValues = np.zeros(self.diodeCorrectedShape, dtype=float) + 1.11
        #xyzRawCenter,centerModule = self.raw_pixel_to_xyz_coordinates(self.rawCenter, withModule=True)
        print "build fictious data for distance %.1f and center %.1f, %.1f"%(self.distance, self.rawCenter[0], self.rawCenter[1])
        if self.geometricArrays == None :
            self.mapped_geometric_arrays()
        xyzRawCenter = self.raw_pixel_to_xyz_coordinates(self.rawCenter)
        e = self.distance / pixelZSize
        if rad :
            for radius in range(100,self.rawShape[1],100) :
                "circle intersect with modules"
                #radius = pixelRadius * pixelYSize
                for module in  range(0,self.modulesNumber) :
                    "(y-yc)2+(z-zc)2=R2*(e+x)2/e2, z=x/a, x>0 repere du module"
                    xm = self.moduleCoordinates[0,2,module]
                    zo = self.moduleCoordinates[2,0,module]
                    a = xm / (self.moduleCoordinates[2,2,module] - zo)
                    ycp = (xyzRawCenter[1]- self.moduleCoordinates[1,0,module]) / pixelYSize
                    zcp = (xyzRawCenter[2]-zo) / pixelZSize
                    zop = zo / pixelZSize
                    ad = (radius*a)**2-e**2
                    bd = radius**2*e*a+e**2*zcp
                    delta = bd**2-ad*((radius*e)**2-(e*zcp)**2)
                    if delta <= 0 :
                        continue
                    delta = math.sqrt(delta)
                    if ad >= 0 :
                        continue
                    z1 = (-bd + delta)/ad
                    z2 = (-bd - delta)/ad
                    if z1 >= 120.0 or z2 <= 0 :
                        continue
                    if z1 < 0.0 : z1 = 0.0
                    if z2 > 120.0 : z2 = 120.0
                    z = np.arange(math.floor(z1), math.ceil(z2)+1)
                    y = np.sqrt((radius*(e+a*z)/e)**2 - (z-zcp)**2)

                    def drawIntersectionAlongZ(z, yi, maxY) :
                        for iz in range(z.shape[0]) :
                            zz = int(z[iz])
                            if zz >= 120 : continue
                            if not math.isnan(yi[iz]) and 0 <= yi[iz] < maxY :
                                yy = int(yi[iz])
                                ni = intensity * (yy+1-yi[iz])
                                diodesValues[yy, zz+120*module] = ni
                                if zz < 120-1 :
                                    diodesValues[yy, zz+120*module+1] = intensity - ni
                    drawIntersectionAlongZ(z, ycp - y, diodesValues.shape[0]-1)
                    drawIntersectionAlongZ(z, ycp + y, diodesValues.shape[0]-1)

                    zzo =  - bd/ad
                    y0 = max(math.sqrt(max(0, (radius*(e+a*z1)/e)**2 - (z1-zcp)**2)), math.sqrt(max(0, (radius*(e+a*z2)/e)**2 - (z2-zcp)**2)))
                    y1 = math.floor(max(0,ycp-y0))
                    y2 = math.ceil(min(ycp+y0+1,578))
                    y = np.arange(y1, y2)
                    delta = bd**2-ad*((radius*e)**2-(e*zcp)**2-(e*(y-ycp))**2)
                    delta = np.sqrt(delta)/ad

                    def drawIntersectionAlongY(y, zi, maxZ) :
                        for iy in range(y.shape[0]) :
                            diode_corrected_maskyy = int(y[iy])
                            if not math.isnan(zi[iy]) and 0 <= zi[iy] < maxZ :
                                zz = int(zi[iy])
                                ni = intensity * (zz+1-zi[iy])
                                diodesValues[yy, zz+120*module] = ni
                                if yy < 578-1 :
                                    diodesValues[yy+1, zz+120*module] = intensity - ni
                    drawIntersectionAlongY(y, zzo - delta, 119)
                    drawIntersectionAlongY(y, zzo + delta, 119)
        rawValues = np.zeros(self.rawShape, dtype=float)
        for i in range(0,self.chipsPerModule) :
            io = i*xpad3Cols
            jo = i*(xpad3Cols+xpad3DiodeGap)
            rawValues[io:io+xpad3Cols,:] = diodesValues[jo:jo+xpad3Cols,:]
            if i != 0 :
                rawValues[io-1,:] += diodesValues[jo-3,:]+ diodesValues[jo-2,:]/2
                rawValues[io,:] += diodesValues[jo-1,:]+ diodesValues[jo-2,:]/2
        return rawValues.astype(int)


    # cleaning related procedures
    def masked_stats(self, values, string=False, histo=False) :
        values = values.astype(float)
        valuesSum = np.sum(values)        
        try :
            std = math.sqrt((np.sum(np.square(values))-valuesSum*valuesSum/ self.numberUnMasked)/self.numberUnMasked)
        except :
            print "masked_stats", self.numberUnMasked, valuesSum
            print "sum(x2)", np.sum(np.square(values)),"-(sx)2/n)", valuesSum*valuesSum/ self.numberUnMasked
            std = 0
        mean = valuesSum / self.numberUnMasked
        maxi = np.max(values)
        mean0 = valuesSum/values.size
        #print "masked_stats", mean, valuesSum, self.numberUnMasked
        if string == False:
            return mean, std, maxi, mean0
        string =string+" mean=%.2f(%.2lf) max=%.0f mean_0=%.2f "%(mean, std, maxi, mean0)
        if not histo :
            return string
        return string + " histo :" + str(np.histogram(values)[0])

    def diode_corrected_mask(self) :
        if self.validMaskArray == None :
            self.build_mask()
        self.diodeMaskArray = np.ones(self.diodeCorrectedShape, dtype=bool)
        self.diodeMaskArray[0,:] = self.validMaskArray[0,:]
        self.diodeMaskArray[-1,:] = self.validMaskArray[-1,:]
        for i in range(0,self.chipsPerModule) :
            io = i*self.detectorGeometry["xpad3Cols"]
            jo = i*(self.detectorGeometry["xpad3Cols"]+self.detectorGeometry["xpad3DiodeGap"])
            self.diodeMaskArray[jo+1:jo+self.detectorGeometry["xpad3Cols"]-1,:] = self.validMaskArray[io+1:io+self.detectorGeometry["xpad3Cols"]-1,:]
            if i != 0 :
                v = self.validMaskArray[io-1,:]
                w = self.validMaskArray[io,:]
                vv = np.logical_or(v, w)      # STRATEGY ?
                self.diodeMaskArray[jo-4,:] = v
                self.diodeMaskArray[jo-3,:] = v
                self.diodeMaskArray[jo-2,:] = vv
                self.diodeMaskArray[jo-1,:] = w
                self.diodeMaskArray[jo,:] = w

    def build_mask(self, geometrical=True, fromFile=None, noisyAbove=None, values=None, save=None) :
        if self.validMaskArray == None :
            self.validMaskArray = np.ones(self.rawShape, dtype=bool)
        if geometrical :
            self.validMaskArray[:, self.overlappedRowMask] = False
            self.numberUnMasked = len(np.where(self.validMaskArray==True)[0])
            print "build_mask geometrical %d / %d"%(self.numberUnMasked, self.validMaskArray.size,)
        if fromFile != None :
            if fromFile == True :
                fromFile = "xpad3_tiled_geometry_raw.mask"
            if os.path.exists(fromFile) :
                mask = det_mask.det_mask(None, self.rawShape, name=fromFile).get_unmasked()
                self.validMaskArray = np.logical_and(self.validMaskArray, mask)
            self.numberUnMasked = len(np.where(self.validMaskArray==True)[0])
            print "build_mask fromFile '%s' %d / %d"%(fromFile, self.numberUnMasked, self.validMaskArray.size,)
        cleanValues = np.where( self.validMaskArray, values, 0)
        print arg, self.masked_stats(cleanValues, string='base_mask', histo=True)
        if values!=None and noisyAbove != None :
            noisy = np.where( cleanValues > noisyAbove)
            cleanValues[noisy] = 0
            self.validMaskArray[noisy] = False
            self.numberUnMasked = len(np.where(self.validMaskArray==True)[0])
            mean, std, maxi, mean0 = self.masked_stats(cleanValues)
            print arg, self.masked_stats(cleanValues, string='%d noisy > %d removed, %d remaining,'%(noisy[0].size, noisyAbove, self.numberUnMasked), histo=True)
        if save !=None :
            self.save_mask(format='ALL', name=save)
        return cleanValues


    def load_mask(self, name=None) :
        if name.endswith("FMaskI.npy"):
            self.validMaskArray = np.ones(self.rawShape, dtype=bool)
            flatMaskIndices = np.load(name)
            self.validMaskArray.flat[flatMaskIndices] = False
            print name , " has been loaded"
	    return True
	else : 
	    #print name , " mask not found"
	    return False
        pass

    def save_mask(self, format='FMI', name=None):
        if name == None :
            name = 'xpad3_tiled_geometry_raw'
        if format == 'FMI' or 'ALL' :
            flatMaskIndices = np.where(self.validMaskArray.flatten() == False)[0]
            np.savetxt(name+'_FMaskI.dat', flatMaskIndices, fmt='%d')
        if format == 'BIN' or 'ALL' :
            array = self.validMaskArray.T.astype(np.int8)
            array.tofile(name+'_Mask.bin')
        if format == 'NPY' or 'ALL' :
            flatMaskIndices = np.where(self.validMaskArray.flatten() == False)[0]
            np.save(name+'_FMaskI.npy', flatMaskIndices)
        if format == 'EDF' or 'ALL' :
            edf_write(self.validMaskArray.astype(np.int), name, ext='_Mask.edf' )

    def apply_mask(self, values, mask=None, geometrical=True):
        if mask == None  :
            mask = self.validMaskArray        
        if geometrical :
            mask[:, self.overlappedRowMask] = False
        flatMaskIndices = np.where(mask.flatten() == False)[0]
        self.numberUnMasked = self.validMaskArray.size - len(flatMaskIndices)
        values.flat[flatMaskIndices] = 0
        print self.masked_stats(values, string='apply_mask', histo=True)
        return values

    def get_corrected_mask(self):
        return self.correctedMaskArray

    def get_raw_mask(self):
        return self.validMaskArray

    def get_mask_strategy(self):
        pass

    def set_mask_strategy(self, strategy=None):
        pass


################################################################################
def radial_distribution0(center, data, maxi, removeZero=True) :
    print "radial_distribution center=%.2f, %.2f "%center
    yzList = np.indices(data.shape)
    distances = np.hypot(yzList[0,...]-center[0], yzList[1,...]-center[1]).flatten()
    hash_1d = distances.argsort(0)
    pixel_bin = np.arange(0, distances[hash_1d[-1]], 1)
    bins = np.searchsorted(distances[hash_1d], pixel_bin)
    # bin will contain the index of the 1st element > pixel_bin.
    # Thus, the slice distances[bin[i]:bin[i+1]] will give all elements w/in bin
    means = np.empty(bins.shape, dtype=float)
    radii = np.empty(bins.shape, dtype=float)
    for ibin in range(bins.size):
        pixels = hash_1d[bins[ibin]:bins[ibin+1]]
        intensities = data.flatten()[pixels]    #, dtype='ushort')
        u_distances = distances[pixels]   #, dtype='ushort')
        if removeZero :
            nz = np.where(intensities != 0)
            intensities = intensities[nz]
            u_distances = u_distances[nz]
        radii[ibin] = u_distances.mean()
        means[ibin] = intensities.mean()
        if radii[ibin] > maxi :
            return radii[:ibin], means[:ibin]
    return radii, means
################################################################################
def radial_distributions(center, data, maxi, removeZero=True, sectors=4) :
    print "radial_distributions with %d sectors (1st sym around X)"%(sectors)+" center=%.2f, %.2f "%center
    yzList = np.indices(data.shape)
    #print "yzList", yzList.shape 
    distances = np.hypot(yzList[0,...]-center[0], yzList[1,...]-center[1]).flatten()
    hash_1d = distances.argsort(0)
    pixel_bin = np.arange(0, distances[hash_1d[-1]], 1)
    bins = np.searchsorted(distances[hash_1d], pixel_bin)
    # bins will contain the index of the 1st element > pixel_bin.
    # Thus, the slice distances[bin[i]:bin[i+1]] will give all elements w/in bin
    yzIndexed = np.reshape(yzList,(2,-1))
    #print "yzIndexed", yzIndexed.shape, yzIndexed
    if bins.size > maxi :
        size = maxi
    else :
        size = bins.size
    radii = np.arange(0.5,size+0.5,1.0, dtype=float)
    results = np.zeros((size,sectors), dtype=float)
    for ibin in range(size):
        pixels = hash_1d[bins[ibin]:bins[ibin+1]]
        angles = np.arctan2(yzIndexed[0,pixels]-center[0], yzIndexed[1,pixels]-center[1]).flatten()
        sectorRange = 2*np.pi/sectors
        angles = np.where(angles >= -sectorRange/2.0, angles, angles + 2*np.pi)
        for sector in range(sectors) :
            start = (sector-0.5)*sectorRange
            selected = np.where(np.logical_and(angles >= start , angles < start+sectorRange))
            selectedPixels = pixels[selected]    
            intensities = data.flatten()[selectedPixels]    # dtype='ushort')
            u_distances = distances[selectedPixels]         # dtype='ushort')
            if removeZero :
                nz = np.where(intensities != 0)
                intensities = intensities[nz]
                u_distances = u_distances[nz]
            if u_distances.size > 0 :
                results[ibin, sector] = intensities.mean()
    return radii, results


################################################################################
def edf_write(data, name, ext='g.edf', verbose=False,**kwargs) :
    fname = os.path.basename(os.path.splitext(name)[0]+ext)
    out = open(fname, 'w')
    out.write("""{
 HeaderID  = EH:000001:000000:000000;
 ByteOrder = LowByteFirst ;
 Image     = 1 ;
""" )
    out.write(' Dim_1 = %d;\n' % data.shape[0])
    out.write(' Dim_2 = %d;\n' % data.shape[1])
    if data.dtype == np.float64 :
        out.write(' DataType = FloatValue ;\n')
        data=data.astype(np.float32)
    elif data.dtype ==  np.float32 :
        out.write(' DataType = FloatValue ;\n')
    #elif repr(type(data[0,0])) == "<type 'numpy.int64'>" :
    elif data.dtype == np.int64 :
        out.write(' DataType = SignedInteger ;\n')
        data=data.astype(np.int32)
    else :
        out.write(' DataType = SignedInteger ;\n')
    for key in kwargs :
        out.write(' %s = %s;\n'%(key, repr(kwargs[key])))
    for i in range(out.tell(),509) : out.write(' ')
    out.write('\n}\n')
    out.write(data.T.tostring())
    out.close()
    if verbose : print "edf image",fname, "has been writen"
################################################################################
def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.
    output:
        the smoothed signal
    TODO: the window parameter could be the window itself if an array instead of a string
    """
    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."
    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."
    if window_len<3:
        return x
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"

    s=np.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
    if window == 'flat': #moving average
        w=numpy.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')
    y=np.convolve(w/w.sum(),s,mode='same')
    return y[window_len:-window_len+1]
    
################################################################################
def smoothedValues(values) :
    smouthValues = 4.0*values.copy()
    smouthValues[1:,:] += values[:-1,:]
    smouthValues[:-1,:] += values[1:,:]
    smouthValues[:,1:] += values[:,:-1]
    smouthValues[:,:-1] += values[:,1:]
    smouthValues /= 8.0
    return smouthValues
    
def process(detector, arg, options) :
    global pylab
    print "xpad3_tiled_geometry.process", arg
    if options.fictious :
        values = detector.fictious(rad=options.rad)
    else :
        values = detector.load(arg, noisyMask=options.noisymask)
        if arg.endswith(".gz") :
            arg = arg[:-3]
        if arg.endswith(".dat") :
            arg = arg[:-4]
    mean = values.mean()
    print arg, values.shape, 'raw max= %.1f mean=%.1f'%(np.amax(values.flat), mean)
    if options.allEdf or options.rawOnly :
        edf_write(values, arg, ext='raw.edf' )
    if options.rawOnly :
        return
        
    if options.mask == 'search' :
        for maskName in glob.glob(options.maskLoad) :
            if os.path.exists(maskName) :
                if detector.load_mask(maskName) :
                    break
        if detector.validMaskArray != None :
            values = detector.apply_mask(values)
        else :
            print "no valid mask has been found"
    elif options.mask == 'build' :
        values = detector.build_mask(geometrical=True, fromFile=True,
            noisyAbove=options.noisyabove, values=values, save=options.maskSave)
    elif options.mask == 'None' and options.noisyabove != None :
        noisy = np.where( values > options.noisyabove)
        values[noisy] = 0
        detector.numberUnMasked=values.size-noisy[0].size
        #mean, std, maxi, mean0 = self.masked_stats(values)
        print arg, detector.masked_stats(values, string='%d noisy > %d removed, %d remaining,'%(noisy[0].size, options.noisyabove, detector.numberUnMasked), histo=True)
    elif options.mask != 'None' :
        raise NotImplementedError("options.mask != None, build, search")
    if options.allEdf or options.maskedOnly :
        edf_write(values, arg, ext='maskedRaw.edf' )
    if options.maskedOnly :
        return

        values = values.astype(float)

    final = detector.flatten_counts(values, normalize=options.normalize)
    if options.multiply != 1.0 :
    	final = float(options.multiply)*final
    	print "image multiplied by",  options.multiply
    if options.center != None:
        center = detector.get_flat_pixel_coordinates(options.center)
    else :
        center = (final.shape[0]/2.0,final.shape[1]/2.0)
    if options.save == '.edf' or options.edf :
        if options.center != None:
            edf_write(final, arg, center=center, distance=options.distance)
        else :
            edf_write(final, arg, distance=options.distance)
    if  options.save != None and options.save != 'edf':
        print 'save option to be implemented'
    if options.edf :
        try :
            open("f2d.mac","w").write("I ACCEPT\nIMAGE PROCESSING (GENERAL)\nINPUT\n"+os.path.basename(os.path.splitext(arg)[0])+"g.edf\n")
            os.system("fit2d -dim1000x1000 -macf2d.mac &")
        except : pass
    if options.rad or options.showrad :
        nRadii=448
        radii, radial0 = radial_distribution0(center, final, nRadii)
        rad = open(os.path.splitext(arg)[0]+"_cen_%.1f_%.1f.dat"%center,"w")
        for i in range(nRadii) :
            rad.write("%7.2f %8.2f\n"%(radii[i],radial0[i]))
        rad.close()
        if options.sectors > 1 :
            radii, radial = radial_distributions(center, final, nRadii, sectors=options.sectors)
            rad = open(os.path.splitext(arg)[0]+"_%d"%(options.sectors)+"cen_%.1f_%.1f.dat"%center,"w")
            for i in range(nRadii) :
                rad.write("%7.2f "%(radii[i]))
                for s in range(options.sectors) :
                    rad.write(" %8.2f"%(radial[i,s]))
                rad.write("\n")
            rad.close()

    if options.test != None :
        fig = pylab.figure(figsize=(10,6))
        #ax = pylab.subplot(131)  #row, col, num
        ax = fig.add_axes((0.04,0.05,0.28,0.8), label='axes1')  #left, bottom, width, height
        pylab.title('raw data '+repr(values.shape)+'\n center %.1f %.1f'%(options.center))
        smouthValues = smoothedValues(values)
        mean = final.mean()
        img = ax.imshow(np.log(smouthValues.T), vmin=0, vmax=math.log(3*mean), origin='lower')
#            ax = pylab.subplot(132)  #left, bottom, width, height
        ax = fig.add_axes((0.37,0.05,0.28,0.8), label='axes2')  #left, bottom, width, height
        pylab.title(detector.masked_stats(final,string="")+'\ngeo. cor. '+repr(final.shape)+'\n center %.1f %.1f'%center)
        #img = ax.imshow(final.T, vmin=0, vmax=2*mean, origin='lower')
        img = ax.imshow(np.log(final.T), vmin=0, vmax=math.log(3*mean), origin='lower')
        pylab.colorbar(img, ax=ax, orientation='horizontal', fraction=0.03, pad=0.04, format='%.1f' )

        #ax = pylab.subplot(236)
        ax = fig.add_axes((0.72,0.05,0.26,0.40))  #left, bottom, width, height
        if options.rad :
            nRadii=448
            rawradii, rawradial = radial_distribution(options.center, values, nRadii, "raw")
            #radial = np.log(radial)
            pylab.plot(rawradii, rawradial, 'g-')
            #radial = np.log(radial)
            pylab.plot(radii, radial, 'r-')
            mean = radial[0:nRadii].mean()
            ax.set_yscale('log')
            ax.set_ylim(0.5*mean , 5*mean )
        elif options.test :
            raw_mask = detector.get_raw_mask().sum(axis=0)
            rawVProj = values.sum(axis=0)/raw_mask
            rawVProj = np.where(raw_mask != 0, rawVProj, 0.0 )
            try :
                cor_mask = detector.get_corrected_mask().sum(axis=0)
                corVProj = final.sum(axis=0)/cor_mask
                corVProj = np.where(cor_mask != 0, corVProj, 0.0 )
            except :
                print "except : detector.get_corrected_mask"
                corVProj = final.sum(axis=0)/final.shape[0]
            pylab.savetxt('xpad3_geometry_hproj_cor.dat', final.sum(axis=-1) , fmt="%12.6G")
            pylab.savetxt('xpad3_geometry_vproj_cor.dat', corVProj, fmt="%12.6G")
            pylab.savetxt('xpad3_geometry_vproj_raw.dat', rawVProj, fmt="%12.6G")
            mean = corVProj.mean()
            pylab.plot(rawVProj, 'g-'  )
            pylab.plot(corVProj, 'r-' )
            w=5
            rawSmooth=np.append(np.zeros((w),dtype=np.float),smooth(rawVProj, window_len=2*w+1)+0.1*mean)
            pylab.plot(rawSmooth, 'b-'  )
            corSmooth=np.append(np.zeros((w),dtype=np.float),smooth(corVProj, window_len=2*w+1)+0.1*mean)
            pylab.plot(corSmooth, 'r-'  )
            ax.set_yscale('linear')
            ax.set_ylim( 0.9*mean, 1.2*mean)
            print    "rawVProj",  rawVProj.shape, rawVProj.min(), rawVProj.max(), rawVProj.mean(), "corVProj",  corVProj.shape, corVProj.min(), corVProj.max(), mean


        #pylab.subplot(233)
        ax = fig.add_axes((0.71,0.55,0.27,0.35))  #left, bottom, width, height
        d = -options.distance
        if d < -30 or d > -10 : d= -30
        pixelZSize = detector.detectorGeometry["pixelZSize"]
        pylab.plot([d,0], [center[1]*pixelZSize,center[1]*pixelZSize])
        pylab.plot([0,0], [0, final.shape[1]*pixelZSize])
        for m in range(detector.detectorGeometry["modulesNumber"]) :
            detectorCoordinates = detector.tiled_module_coordinates(module=m)
            pylab.plot(detectorCoordinates[0,:], detectorCoordinates[2,:], 'ro-')
            if m == 0 :
                pylab.plot([d,0], [center[1]*pixelZSize+(detectorCoordinates[2,0]-center[1]*pixelZSize)*(options.distance+d)/options.distance, detectorCoordinates[2,0]])
            if m == 7 :
                pylab.plot([d,detectorCoordinates[0,3]], [center[1]*pixelZSize+(detectorCoordinates[2,3]-center[1]*pixelZSize)*(options.distance+d)/options.distance, detectorCoordinates[2,3]])
        pylab.xlabel('Y')
        pylab.ylabel('Z')
        pylab.title('Tiled geometry\ndistance %.0fmm'%(options.distance))
        pylab.grid(True)
        pylab.ylim( -10, detectorCoordinates[2,3]+10 )
        if options.fictious :
            pylab.savefig('xpad3_geometry_%d.png'%(options.distance))
        else :
            pylab.savefig('xpad3_geometry.png')

    elif options.show :
        fig = pylab.figure(figsize=(10,6))
        ax = pylab.subplot(111)  #left, bottom, width, height
        pylab.title('geo. cor. '+repr(final.shape)+'\n center %.1f %.1f'%center)
        img = ax.imshow(final.T, vmin=0, vmax=2*mean, origin='lower')

    elif options.showrad :
        fig = pylab.figure(figsize=(10,6))
        ax = pylab.subplot(121)  #left, bottom, width, height
        pylab.title('geo. cor. '+repr(final.shape)+'\n center %.1f %.1f'%center)
        smoothed= np.log(smoothedValues(final))
        mean = smoothed.mean()
        print mean
        if mean > 0 :
            img = ax.imshow(smoothed.T, vmin=0, vmax=3*mean, origin='lower')
        else :
            img = ax.imshow(smoothed.T, origin='lower')            
        pylab.colorbar(img, ax=ax, orientation='horizontal', fraction=0.03, pad=0.04, format='%.1f' )
        ax = pylab.subplot(122)  #left, bottom, width, height
        ax.plot(radii, radial0, '-', label="all")
        if options.sectors > 1 :
            for s in range(options.sectors) :
                ax.plot(radii, radial[:,s], '-', label=repr(s))
        ax.legend( loc=1, ncol=1,  borderaxespad=0.)
       #mean = radial0[0:nRadii].mean()
        #mini = radial0[0:nRadii].min()
        ax.set_yscale('log')
        #ax.set_ylim(0 , 10*mean )

    

################################################################################
if __name__ == '__main__':


    import optparse
    parser = optparse.OptionParser(usage=__doc__)
    parser.add_option("--geometry=", dest="geometry", action="store", default = None, help="name of detector : D1, D5, imXpad2 or filename name.geo to be exec")
    parser.add_option("--distance=", dest="distance", action="store", default = 0, help="distance (mm) from scatering center")
    parser.add_option("--center=", dest="center", action="store", default = None, help="projection of the scattering center")
    parser.add_option("--save=", dest="save", action="store", default = '.edf', help="save axtension or None")
    parser.add_option("--show", dest="show", action="store_true", help="show result image using matplotlib")
    parser.add_option("--edf", dest="edf", action="store_true", help="fork edf image to fit2d")
    parser.add_option("--rad", dest="rad", action="store_true", help="build a radial distribution")
    parser.add_option("--showrad", dest="showrad", action="store_true", help="build a radial distribution and show it")
    parser.add_option("--sectors", dest="sectors", action="store", help="sectors used in radial distribution ", default = 4)
    parser.add_option("--unMaskOverlapped", dest="unMaskOverlapped", action="store_true", help="do not mask pixels overlapped by another module for source at infinite, often badly calibrated and ignored by default")
    parser.add_option("--yRotation=", dest="yRotation", action="store", default = 0.0, help="rotation error around Y (deg)")
    parser.add_option("--zRotation=", dest="zRotation", action="store", default = 0.0, help="rotation error around Z (deg)")
    #parser.add_option("--modules", dest="modules", action="store", default=modulesNumber, help="specify number of modules/detector")
    #parser.add_option("--chips", dest="chips", action="store", default=chipsPerModule, help="specify number of chips/module")
    parser.add_option("--badNonZero", dest="badNonZero", action="store_true", help="bad are consodered as unknow and non zero (TO BE IMPLELENTED)")
    parser.add_option("--rectangular=", dest="rectangular", action="store_true", help="output on squared pixel and not on rectangular ones")
    parser.add_option("--noTiling", dest="noTiling", action="store_true", help="do not apply correction for tiled detector")
    parser.add_option("--test", dest="test", action="store", default = None, help="use default test data = 1(rad), 2 ,3")
    parser.add_option("--fictious", dest="fictious", action="store_true", help="build a default test data")
    parser.add_option("--noisyabove", dest="noisyabove", action="store", default = None, help="remove spot above")
    parser.add_option("--noisymask", dest="noisymask", action="store", default = None, help="remove spot using bit mask")
    parser.add_option("--mask", dest="mask", action="store", default = 'search', help="mask procedure to be applied (search=default, build, None)")
    parser.add_option("--maskLoad", dest="maskLoad", action="store", default = "*Mask*", help="Name of the mask file to be loaded if exists")
    parser.add_option("--maskSave", dest="maskSave", action="store", default = "xpad3_tiled_geometry_raw", help="Name of the mask file to be saved ")
    parser.add_option("--allEdf", dest="allEdf", action="store_true", help="save all edf not only last one")
    parser.add_option("--rawOnly", dest="rawOnly", action="store_true", help="do not correct for geometry, just build edf")
    parser.add_option("--maskedOnly", dest="maskedOnly", action="store_true", help="do not correct for geometry, just build masked edf")
    parser.add_option("--normalize=", dest="normalize", action="store", default = "", help="normalize from spec file")
    parser.add_option("--monitor=", dest="monitor", action="store", default = "vct1_3", help="monitor for normalization, default vct1_3")
    parser.add_option("--multiply=", dest="multiply", action="store", default = 1.0, help="multiply counts by, default 1.0")
    parser.add_option("--debug", dest="debug", action="store_true")

    (options, args) = parser.parse_args()
    print options

    testData = { 1 : {'distance' : 1300, 'center' : (298.2, 356.4), 'args' : "05Oct10_48.dat", 'rad' : True },
        2 : {'distance' : 735.5, 'center' : (272, 354), 'args' : "21Jul10_1710.dat.gz", 'rad' : False, 'noisyabove' : 100},
        3 : {'distance' : 319.0, 'center' : (275, 467), 'args' : "28Jan11_298.dat.gz", 'rad' : False}
        }
    
    if options.fictious :
        options.test = 1
    #print "options.test", options.test, testData[options.test]
    if options.test != None :
        print "options.test", options.test
        options.test = eval(options.test)
        print "options.test", options.test, "testData", testData[options.test]
        options.show = True
        options.allEdf = True
        if testData[options.test].has_key('rad') :
            options.rad = testData[options.test]['rad'] 
        if options.distance == 0 :
            options.distance = testData[options.test]['distance']
        else :
            options.distance = eval(options.distance)
        if options.center == None :
            options.center = testData[options.test]['center']
        else :
            options.center = eval(options.center)
        if not options.fictious :
            if not options.noisyabove :
                if testData[options.test].has_key('noisyabove') :
                    options.noisyabove = testData[options.test]['noisyabove']
                else :
                    options.noisyabove = 20000
            else :
                options.noisyabove = eval(options.noisyabove)
            if not options.noisymask :
                options.noisymask = 0x0003FFF
            else :
                options.noisymask = eval(options.noisymask)
            from DEVA.xpad.xpad3_spec import XpadSpec
        #options.edf = True
        #options.rad = True
        if args == [] :
            if options.fictious :
                args.append("fictious_%d"%(options.distance))
            else :
                args.append(testData[options.test]['args'])
    else :
        if options.center != None :
            options.center = eval(options.center)
        if options.distance != 0 :
            options.distance = eval(options.distance)
        if options.noisyabove :
                options.noisyabove = eval(options.noisyabove)
        if options.noisymask :
                options.noisymask = eval(options.noisymask)

    if options.show or options.showrad :
        import pylab
        import warnings
        warnings.simplefilter('ignore')

    if options.normalize :
        from DEVA.xpad.xpad3_spec import XpadSpec

    if options.geometry == None :
        if prgName == "xpad3_tiled_geometry" :
            options.geometry = "D1"
        else :
            c = prgName.split('_')
            if len(c) >1 and c[0] == "xpad3":
                options.geometry = c[1]
            if options.geometry == None or options.geometry == "geometry" :
                options.geometry = "S70"
        print "options.geometry set to :", options.geometry
    if options.multiply :
        print "options.multiply a revoir nb 30/04/12, tj actuel ? jf 9/10/12 "
        #options.multiply = eval(options.multiply)

        
    # create the detector with its properties
    detector = Xpad3Coordinates( options.geometry,
                    distance=options.distance,
                    rawCenter=options.center,
                    unMaskOverlapped=options.unMaskOverlapped,
                    yRotation=options.yRotation,
                    zRotation=options.zRotation,
                    options=options)

    print "args", args

    if args == [] :
        parser.print_help()
    for arg in args :
        arg=arg.strip()
        process(detector, arg, options)

    if options.show or options.showrad :
        #pylab.ion()
        pylab.show()
