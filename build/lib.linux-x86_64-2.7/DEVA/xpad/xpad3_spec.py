#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This file is licensed under the CeCILL licence
# the french "GPL" licence agreed by CNRS
# see  http://www.cecilll.info
#
""" xpad3_spec.py [options]  args ...
    1.07 semble ok
    """
__author__ = "J.-F. Berar (CNRS/NEEL)"
__version__ = "1.07"
__date__ = "26/09/2011"

import os, os.path, re
default_float_precision = 3

#--------------------------------------------------------------
# test sur fourc.26Oct11 ~maret
#--------------------------------------------------------------
class XpadSpec ( object ) :
    """ get Xpad 3 image normaization in spec file """
    def __init__ (self, fileName = "", image= None, monitor=None) :
        lines = []
        self.tocValues = {}
        self.tocHeads = {}
        self.keys = []
        self.monitor = monitor
        try :
            lines = open(fileName).readlines()
        except :
            raise "unable to open or read " + fileName
        inScan = False
        scanImageName = None
        #print "XpadSpec(filename=",fileName,"image=",image,"monitor=",monitor,")"
        for line in lines :
            line = line.strip()
            if line == "" or line.startswith("#C ")  :
                continue
            elif line.startswith("#R") and inScan :
                inScan = False
            elif line.startswith("#L") :
                countersName = line[1+len("#L"):]
            elif line.startswith("#XPAD_IMAGE_NUMBER") :
                h,name = line.split()
                if inScan :
                    scanImageName = name
                #    try :
                #        self.tocValues[name] = countersValue
                #        self.tocHeads[name] = countersName
                #    except :
                #        print "error line=", line
                #        print "error name=", name
                #        print "error countersName=", countersName
                #        print "error countersValue=", countersValue
                #if name.startswith("sl47e9") :
                   #try : print "XPAD_IMAGE_NUMBER", name, self.tocValues[name]
                   #except : pass
            elif line.startswith("#XPADSCAN") :
                inScan = True
                scanImageName = None
            elif line.startswith('#COUNTERS') :
                #print name, line
                self.tocValues[name] = line[1+len("#COUNTERS"):]
                self.tocHeads[name] = countersName
                inScan = False
                #if name.startswith("sl47e9") :
                   #try : print "COUNTERS", name, self.tocValues[name]
                   #except : pass
            elif inScan and scanImageName != None and not line.startswith('#') :
                self.tocValues[name] = line
                self.tocHeads[name] = countersName
					 #if name.startswith("sl47e9") :
                   #try : print "INSCAN", name, self.tocValues[name]
                   #except : pass
        self.keys = sorted(self.tocValues.keys())
        #for k in self.keys :
            #if k.startswith("sl47e9") : print k, self.tocValues[k], self.tocHeads[k]
        if image != None :
            self.get_monitor(image, monitor=monitor)

    def get_monitor(self, image, monitor=None, show=False) :
        if monitor != None :
            self.monitor = monitor
        for k in self.keys :
            if image == None or re.search(image, k) :
                if show :
                    print k, self.tocValues[k], self.tocHeads[k]
                monitors = self.tocHeads[k].split()
                for m in range(len(monitors)) :
                    if monitors[m]  == self.monitor :
                        value = self.tocValues[k].split()[m]
                print k, self.monitor, value

    def monitor_value(self, image, monitor=None, show=False) :
        returnValues = {}
        show=True
        #print "XpadSpec.monitor_value(image=",image," monitor=", monitor,")"
        if monitor != None :
            self.monitor = monitor
        for k in self.keys :
#            if image == None or re.search(image, k) :
            if image == None or re.match(image+"[^0-9]", k) :
                if show :
                    print k, self.tocValues[k], self.tocHeads[k]
                monitors = self.tocHeads[k].split()
                for m in range(len(monitors)) :
                    if monitors[m]  == self.monitor :
                        returnValues[k] =self.tocValues[k].split()[m]
        #print "DEBUG", returnValues, self.monitor
        if len(returnValues.keys()) == 1 :
            return returnValues[returnValues.keys()[0]]
        else :
            return returnValues
        
################################################################################


if __name__ == '__main__':


    import sys, optparse
    parser = optparse.OptionParser(usage=__doc__)
    parser.add_option("--fileName=", dest="fileName", action="store", default = "specfile", help="specfile.Name")
    parser.add_option("--image=", dest="image", action="store", default = None, help="name of image")
    parser.add_option("--monitor=", dest="monitor", action="store", default = "vct1_3", help="monitor to use")
    parser.add_option("--show=", dest="show", action="store_true", default = False, help="show all counters")

    (options, args) = parser.parse_args()
    print sys.argv[0],__version__,__date__, options, sys.argv[1:]
    
        
    if len(args) == 0 :
        xpadSpec = XpadSpec(fileName=options.fileName, monitor=options.monitor)
        print xpadSpec.get_monitor(None, show=options.show)
    elif len(args) == 1 :
        print args[0], XpadSpec(fileName=options.fileName, monitor=options.monitor).monitor_value(args[0], show=options.show)
    else :
        xpadSpec = XpadSpec(fileName=options.fileName, image=options.image, monitor=options.monitor)
        for arg in args :
            print arg, xpadSpec.monitor_value(arg, show=options.show)
        
