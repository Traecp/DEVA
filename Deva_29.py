#!/usr/bin/python
# -*- coding: utf-8 -*-
################ DEVA software: D2AM Edf images Visualisation and Analysis ##############
###### Dependencies: numpy, scipy, matplotlib, lmfit (sudo easy_install -U lmfit), pyFAI, fabio, Basemap
#from gi.repository import GObject
import gtk,gobject,threading, time
import os
import gc
from os import listdir
from os.path import isfile,join
import math
import numpy as N
from numpy import unravel_index
from scipy import ndimage
from lmfit import Parameters, minimize
from DEVA.xpad import libXpad as libX
##############
## Graphic library ##
##############
import matplotlib as mpl
mpl.use('GtkAgg')
from mpl_toolkits.basemap import Basemap
from matplotlib.figure import Figure
#from matplotlib.axes import Subplot
from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas
from matplotlib.backends.backend_gtkagg import NavigationToolbar2GTKAgg as NavigationToolbar
from matplotlib.cm import jet # colormap
from matplotlib.widgets import Cursor
from matplotlib.patches import Rectangle
from matplotlib.ticker import MaxNLocator
#calculation of two theta,chi
import fabio
try:
	import pyFAI
except:
	from DEVA import pyFAI
try:
	import xrayutilities
except:
	from DEVA import xrayutilities

__author__="Tra NGUYEN THANH"
__version__ = "1.2.9"
__date__="26/11/2014"

#mpl.rcParams['font.size'] = 18.0
mpl.rcParams['axes.labelsize'] = 'large'
mpl.rcParams['legend.fancybox'] = True
mpl.rcParams['legend.handletextpad'] = 0.5
mpl.rcParams['legend.fontsize'] = 'medium'
mpl.rcParams['figure.subplot.bottom'] = 0.13
mpl.rcParams['figure.subplot.top'] = 0.93
mpl.rcParams['figure.subplot.left'] = 0.14
mpl.rcParams['figure.subplot.right'] = 0.915
#mpl.rcParams['image.cmap'] = jet
mpl.rcParams['savefig.dpi'] = 300

#Global variables
_PIXEL_SIZE = 0.130 #mm
_SPEC_IMG_COL = "img" #column containing image number in spec file

def flat_data(data,dynlow, dynhigh):
	""" Returns data where maximum superior than 10^dynhigh will be replaced by 10^dynhigh, inferior than 10^dynlow will be replaced by 10^dynlow"""
	mi = 10**dynlow
	ma = 10**dynhigh
	return N.minimum(N.maximum(data,mi),ma)

def psdVoigt(parameters,x):
	"""Define pseudovoigt function"""
	y0 = parameters['y0'].value
	xc = parameters['xc'].value
	A  = parameters['A'].value
	w  = parameters['w'].value
	mu = parameters['mu'].value

	y  = y0 + A * ( mu * (2/N.pi) * (w / (4*(x-xc)**2 + w**2)) + (1 - mu) * (N.sqrt(4*N.log(2)) / (N.sqrt(N.pi) * w)) * N.exp(-(4*N.log(2)/w**2)*(x-xc)**2) )

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
		idA=N.where(data_x==xc)[0][0]
		A  = data_y[idA]
	y0 = 20
	#xc = data_x[N.argmax(data_y)]
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
	result = minimize(objective, param_init, args=(data_y,data_x))

	x = N.linspace(data_x.min(),data_x.max(),data_x.shape[0])
	y = psdVoigt(param_init,x)

	return param_init, y

def wpercentile(a, q, weights=None):
	"""Compute weighted percentiles *q* [%] of input 1D-array *a*."""

	a = N.asarray(a)
	if a.ndim > 1:
		raise NotImplementedError("implemented on 1D-arrays only")

	if weights is None:
		weights = N.ones_like(a)
	else:
		assert len(weights)==len(a), "incompatible weight and input arrays"
		assert (weights>0).all(), "weights are not always strictly positive"

	isorted = N.argsort(a)
	sa = a[isorted]
	sw = weights[isorted]
	sumw = N.cumsum(sw)                        # Strictly increasing
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

	out = N.array([ interpolate(qq) for qq in N.atleast_1d(q) ])

	return out.reshape(N.shape(q))      # Same shape as input q


def median_stats(a, weights=None, axis=None, scale=1.4826):
	"""Compute [weighted] median and :func:`nMAD` of array *a* along
	*axis*. Weighted computation is implemented for *axis* = None
	only."""

	if weights is not None:
		if axis is not None:
			raise NotImplementedError("implemented on 1D-arrays only")
		else:
			med  = wpercentile(a, 50., weights=weights)
			nmad = wpercentile(N.absolute(a - med), 50., weights=weights)*scale
	else:
		med = N.median(a, axis=axis)
		if axis is None:
			umed = med                      # Scalar
		else:
			umed = N.expand_dims(med, axis) # Same ndim as a
		nmad = N.median(N.absolute(a - umed), axis=axis) * scale

	return med,nmad

def select_files_from_list(s, beg, end):
	""" select edf files limited between number 'beg' and 'end' from the list 's' """
	out=[]
	for i in range(len(s)):
		ss = s[i]
		l = ss.split("_")
		l = l[1]
		n = l.split(".")[0]
		n = int(n)
		if n in range(beg, end+1):
			out.append(ss)
	out = sorted(out)
	return out

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

def get_img_list(edf_list):
	""" Get a list of image numbers from an image_name list """
	out = []
	for i in range(len(edf_list)):
		edf = edf_list[i]
		if edf != None:
			l = edf.split("_")
			l = l[1]
			n = l.split(".")[0]
			n = int(n)
			out.append(n)
		else:
			out.append(None)
	#out = N.asarray(out)
	return N.asarray(out)
	
class MyMainWindow(gtk.Window):

	def __init__(self):
		super(MyMainWindow, self).__init__()
		self.set_title("DEVA - D2AM EDF Visualisation and Analysis - version %s - Last update: %s"%(__version__, __date__))
		self.set_size_request(1200, 900)
		#self.modify_bg(gtk.STATE_NORMAL, gtk.gdk.Color(6400, 6400, 6440))
		self.set_position(gtk.WIN_POS_CENTER)
		self.set_border_width(10)
		#self.set_icon_from_file('DEVA/D2AM.jpg')
		##################### TOOL BAR ####################################################
		self.detector_type = "D5" # My detector by default is XPAD D5
		self.menubar = gtk.MenuBar()

		self.detector_menu = gtk.Menu()
		self.detm = gtk.MenuItem("Detector")
		self.detm.set_submenu(self.detector_menu)

		self.detector_D5 = gtk.MenuItem("XPAD D5")
		self.detector_D5.connect("activate", self.set_detector, "D5") #a definir fonction set_detector
		self.detector_menu.append(self.detector_D5)

		self.detector_D1 = gtk.MenuItem("XPAD D1")
		self.detector_D1.connect("activate", self.set_detector, "D1") #a definir fonction set_detector
		self.detector_menu.append(self.detector_D1)

		self.detector_S70 = gtk.MenuItem("XPAD S70")
		self.detector_S70.connect("activate", self.set_detector, "S70")
		self.detector_menu.append(self.detector_S70)

		self.menubar.append(self.detm)

		self.toolbar = gtk.Toolbar()
		self.toolbar.set_style(gtk.TOOLBAR_ICONS)

		self.refreshtb = gtk.ToolButton(gtk.STOCK_REFRESH)
		self.opentb = gtk.ToolButton(gtk.STOCK_OPEN)
		self.savetb = gtk.ToolButton(gtk.STOCK_SAVE)
		self.sep = gtk.SeparatorToolItem()
		self.quittb = gtk.ToolButton(gtk.STOCK_QUIT)
		self.sep2 = gtk.SeparatorToolItem()
		self.zoomtb = gtk.ToggleToolButton(gtk.STOCK_ZOOM_IN)
		self.hometb = gtk.ToolButton(gtk.STOCK_HOME)
		self.aspecttb = gtk.ToolButton(gtk.STOCK_PAGE_SETUP)
		self.loadcalibtb = gtk.ToolButton(gtk.STOCK_CONVERT)

		self.toolbar.insert(self.opentb, 0)
		self.toolbar.insert(self.refreshtb, 1)

		self.toolbar.insert(self.sep, 2)
		self.toolbar.insert(self.savetb, 3)
		self.toolbar.insert(self.zoomtb, 4)
		self.toolbar.insert(self.hometb, 5)
		self.toolbar.insert(self.aspecttb, 6)
		self.toolbar.insert(self.loadcalibtb, 7)

		self.toolbar.insert(self.sep2, 8)
		self.toolbar.insert(self.quittb, 9)

		self.tooltips = gtk.Tooltips()
		self.tooltips.set_tip(self.refreshtb,"Reload data files")
		self.tooltips.set_tip(self.opentb,"Open a folder containing EDF images")
		self.tooltips.set_tip(self.savetb,"Save image")
		self.tooltips.set_tip(self.quittb,"Quit the program")
		self.tooltips.set_tip(self.zoomtb,"Zoom in")
		self.tooltips.set_tip(self.hometb,"Reset image")
		self.tooltips.set_tip(self.aspecttb,"Change the graph's aspect ratio")
		self.tooltips.set_tip(self.loadcalibtb,"Load a calibration file (PONI file)")

		#self.newtb.set_sensitive(False)
		#self.aspecttb.set_sensitive(False)
		#self.savetb.set_sensitive(False)
		self.opentb.connect("clicked", self.choose_folder)
		self.refreshtb.connect("clicked",self.folder_update)
		self.savetb.connect("clicked", self.save_image)
		self.quittb.connect("clicked", gtk.main_quit)
		self.zoomtb.connect("toggled", self.zoom_on)
		self.hometb.connect("clicked", self.reset_image)
		self.aspecttb.connect("clicked", self.change_aspect_ratio)
		self.loadcalibtb.connect("clicked", self.load_calibration)
		self.graph_aspect = False

		############################# BOXES ###############################################
		vbox = gtk.VBox()
		vbox.pack_start(self.menubar,False,False,0)
		vbox.pack_start(self.toolbar,False,False,0)
		hbox=gtk.HBox()
		######################### TREE VIEW #############################################
		self.sw = gtk.ScrolledWindow()
		self.sw.set_shadow_type(gtk.SHADOW_ETCHED_IN)
		self.sw.set_policy(gtk.POLICY_NEVER, gtk.POLICY_AUTOMATIC)

		hbox.pack_start(self.sw, False, False, 0)
		self.store={}
		self.MODEL = gtk.TreeStore(str)
		#self.list_store = gtk.ListStore(str)
		#self.treeView = gtk.TreeView(self.list_store)
		self.treeView = gtk.TreeView(self.MODEL)
		self.treeView.connect("row-activated",self.on_changed_edf)

		rendererText = gtk.CellRendererText()
		self.TVcolumn = gtk.TreeViewColumn("EDF images", rendererText, text=0)
		self.TVcolumn.set_sort_column_id(0)    
		self.treeView.append_column(self.TVcolumn)

		self.sw.add(self.treeView)
		self.threads = list()  # I made this be part of the full application
		""" That allows us to wait for all threads on close-down, don't know how
		necessary it is."""

		self.current_folder = os.getcwd()
		self.edf_folder     = self.current_folder
		#if self.current_folder is not os.getcwd():
		#	glib.timeout_add_seconds(5, self.folder_update)
		############################# NOTEBOOK ################################################
		self.notebook = gtk.Notebook()
		self.page_single_figure = gtk.HBox()
		self.midle_panel = gtk.VBox()
		self.page_pole_figure   = gtk.HBox()
		self.page_batch_correction = gtk.HBox() #For batch correction of EDF
		
		#**************************** Geometry setup ******************************************
		self.geometry_setup_tbl = gtk.Table(3,2,False)
		self.geometry_manual    = gtk.ToggleButton("CLICK HERE to use the above setup for 2Theta-Chi calculation")
		self.geometry_manual.connect("toggled", self.manual_calibration)
		
		self.geometry_distance_txt = gtk.Label("Sample-Detector distance (m):")
		self.geometry_distance     = gtk.Entry()
		self.geometry_distance.set_text("")
		self.geometry_direct_beam_txt = gtk.Label("Direct beam position X,Y (separated by comma):")
		self.geometry_direct_beam     = gtk.Entry()
		self.geometry_direct_beam.set_text("")
		self.geometry_energy_txt      = gtk.Label("Energy (eV):")
		self.geometry_energy          = gtk.Entry()
		self.geometry_energy.set_text("")
		self.geometry_distance_txt.set_alignment(0,0.5)
		self.geometry_direct_beam_txt.set_alignment(0,0.5)
		self.geometry_energy_txt.set_alignment(0,0.5)
		
		self.geometry_setup_tbl.attach(self.geometry_energy_txt, 0,1,0,1)
		self.geometry_setup_tbl.attach(self.geometry_energy, 1,2,0,1)
		self.geometry_setup_tbl.attach(self.geometry_distance_txt, 0,1,1,2)
		self.geometry_setup_tbl.attach(self.geometry_distance, 1,2,1,2)
		self.geometry_setup_tbl.attach(self.geometry_direct_beam_txt, 0,1,2,3)
		self.geometry_setup_tbl.attach(self.geometry_direct_beam, 1,2,2,3)
		
		############################# PAGE1: FIGURES ##########################################
		self.edf = ""
		self.edf_choosen = ""
		self.my_notes = []
		self.lines = []
		self.points=[]
		self.arb_lines_X=[]
		self.arb_lines_Y=[]
		self.arb_line_points=0
		self.SELECTED_IMG_NUM = None
		self.IMG_INIT = False

		self.header = {}
		self.calibrated = False #The detector need to be calibrated using pyFAI program
		self.calibrated_quantitative = False
		self.fig=Figure(dpi=100)
		#self.ax  = self.fig.add_subplot(111)
		self.ax  = self.fig.add_axes([0.12,0.2,0.7,0.7])
		self.MAIN_XLABEL = self.ax.set_xlabel("X (pixel)")
		self.MAIN_YLABEL = self.ax.set_ylabel("Y (pixel)")
		#self.fig.text(0.5, 0.92, self.edf, ha='center', fontsize=26)
		self.xDim0 = 0
		self.xDim1 = 560
		self.yDim0 = 0
		self.yDim1 = 960
		self.MAIN_EXTENT = (self.xDim0,self.xDim1,self.yDim0,self.yDim1)
		self.fig.subplots_adjust(left=0.1,bottom=0.20, top=0.90)
		self.data = N.zeros(shape=(self.yDim1,self.xDim1))
		self.vmin = 0
		self.vmax = 1000
		self.vmax_range = self.vmax

		#self.init_image()
		self.img = self.ax.imshow(self.data,origin='lower',vmin=self.vmin, vmax=self.vmax, cmap=jet, interpolation='nearest',aspect='auto')

		self.canvas  = FigureCanvas(self.fig)
		#self.main_figure_navBar = NavigationToolbar(self.canvas, self)
		self.cursor = Cursor(self.ax, color='k', linewidth=1, useblit=True)
		#Global color bar
		self.cax = self.fig.add_axes([0.85, 0.20, 0.03, 0.70])#left,bottom,width,height
		self.cb  = self.fig.colorbar(self.img, cax = self.cax, format="%.f")#format=fm
		self.cb.set_label(r'$Intensity\ (Counts\ per\ second)$', fontsize=18)
		self.cb.locator = MaxNLocator(nbins=6)
		#Horizontal colorbar
		self.cax2 = self.fig.add_axes([0.12, 0.1, 0.7, 0.03])#left,bottom,width,height
		self.cb2  = self.fig.colorbar(self.img, cax = self.cax2, orientation='horizontal', format="%.f")#format=fm
		self.cb2.set_label(r'$Intensity\ (Counts\ per\ second)$', fontsize=18)
		self.cb2.locator = MaxNLocator(nbins=6)
		self.cb2.ax.set_visible(False)
		########### ZOOM ACTION #################################
		self.zoom = False #Press the Zoom button on the tool bar
		self.zoom_press = False #Button pressed in the axe to draw a rectangle
		self.rect = Rectangle((0,0),0,0)
		self.x0 =0
		self.y0 =0
		self.x1 =0
		self.y1 =0
		self.ax.add_patch(self.rect)
		self.canvas.mpl_connect("motion_notify_event",self.on_motion)
		self.canvas.mpl_connect("button_press_event",self.on_press)
		self.canvas.mpl_connect("button_release_event",self.on_release)
		self.mouse_moved = False #If click without move: donot zoom the image

		#self.midle_panel.pack_start(self.main_figure_navBar, False,False, 2)
		#*********************************** SCAN SLIDER *****************************************
		#Variables for scan slider:
		self.SCAN_IMG = []
		self.SPEC_IMG = []
		self.SPEC_FILE = ""
		self.SELECTED_IMG_NUM = None
		self.SPEC_IS_LOADED = False
		self.SPEC_DATA = None #SPEC object from xrayutilities.io.SPECFile(specfile)
		self.SPEC_SCAN_LIST = None #list of all spec scan, each spec scan is an object from spec_data.scan_list 
		self.SPEC_ACTUAL_SCAN = None #actual scan number
		self.SPEC_ACTUAL_SCAN_IMG = 0
		self.DATA_IS_LOADED = False
		self.SPEC_ACTUAL_SCAN_DATA = []
		
		self.scan_slider_frame = gtk.Frame()
		self.scan_slider_table_align = gtk.Alignment(0,0.5,1,1)
		self.scan_slider_table_align.set_padding(10,5,5,5)
		
		self.scan_slider_frame.set_label("Scan Slider")
		self.scan_slider_frame.set_label_align(0.5,0.5)
		self.scan_slider_table = gtk.Table(2,3,False)
		self.scan_slider_specfile_txt = gtk.Label("Spec file")
		self.scan_slider_specfile_txt.set_alignment(0,0.5)
		self.scan_slider_scanNumber_txt = gtk.Label("Scan #")
		self.scan_slider_scanNumber_txt.set_alignment(0,0.5)
		self.scan_slider_path = gtk.Entry()
		self.scan_slider_browseSpec = gtk.Button("Browse spec file")
		#self.scan_slider_browseSpec.set_size_request(80,-1)
		self.scan_slider_browseSpec.connect("clicked",self.load_specFile)
		scan_slider_spin_adj = gtk.Adjustment(1,0,1,1,10,0)#actual, min, max, step increment, page increment, page size
		
		self.scan_slider_spinButton = gtk.SpinButton(scan_slider_spin_adj, 1,0)
		self.scan_slider_imgSlider  = gtk.HScale()
		self.scan_slider_imgSlider.set_range(0,1)
		self.scan_slider_imgSlider.set_value(1)
		self.scan_slider_imgSlider.set_digits(0)
		self.scan_slider_imgSlider.set_increments(1,1)
		self.scan_slider_spinButton.set_numeric(True)
		self.scan_slider_spinButton.connect("value-changed",self.update_scan_slider)
		self.scan_slider_imgSlider.connect("value-changed", self.slider_plot_scan)
		
		self.scan_slider_table.attach(self.scan_slider_specfile_txt, 0,1,0,1)
		self.scan_slider_table.attach(self.scan_slider_path, 1,2,0,1)
		self.scan_slider_table.attach(self.scan_slider_browseSpec, 2,3,0,1)
		self.scan_slider_table.attach(self.scan_slider_scanNumber_txt, 0,1,1,2)
		self.scan_slider_table.attach(self.scan_slider_spinButton, 1,2,1,2)
		self.scan_slider_table.attach(self.scan_slider_imgSlider, 2,3,1,2)
		
		self.scan_slider_table_align.add(self.scan_slider_table)
		self.scan_slider_frame.add(self.scan_slider_table_align)
		
		self.midle_panel.pack_start(self.geometry_setup_tbl, False,False, 0)
		self.midle_panel.pack_start(self.geometry_manual, False,False, 2)
		self.midle_panel.pack_start(self.canvas, True,True, 2)
		self.midle_panel.pack_start(self.scan_slider_frame, False,False,0)
		self.page_single_figure.pack_start(self.midle_panel, True,True, 0)

		########################################## Check Buttons RIGHT PANEL ###################

		self.right_panel = gtk.VBox(False,0)

		self.detector_disposition_horizontal = gtk.ToggleButton("Vertical detector")
		self.detector_disposition_horizontal.connect("toggled", self.detector_disposition)
		self.horizontal_detector = False #By default, the detector is in the vertical position, i.e. 960 rows x 560 cols

		self.linear_scale_btn = gtk.ToggleButton("Linear scale")
		self.linear_scale_btn.connect("toggled",self.log_update)
		#self.linear_scale_btn.set_size_request(50,30)

		self.log_scale=0

		self.adj_btn = gtk.CheckButton("Geometry correction")
		#self.adj_btn.set_size_request(600,30)
		self.adj_btn.connect("toggled", self.plot_update)

		self.cln_btn = gtk.CheckButton("Clean data")
		#self.cln_btn.set_size_request(130,30)
		self.cln_btn.connect("toggled", self.plot_update)

		self.tth_chi_space_btn = gtk.CheckButton("2Theta-Chi space")
		#self.tth_chi_space_btn.set_size_request(130,30)
		#self.tth_chi_space_btn.connect("toggled", self.change_space,"TTH")
		self.tth_chi_space_btn.connect("toggled", self.plot_update)

		self.save_adj_btn = gtk.Button("Save Corrected EDF")
		self.save_adj_btn.connect("clicked",self.save_adjust)

		self.separator = gtk.HSeparator()
		#self.plotXYprofiles_btn = gtk.CheckButton("Plot X,Y profiles") #Plot a cross profile of X and Y data
		self.plotXYprofiles_btn = gtk.RadioButton(None,"Plot X,Y profiles")
		self.plotXYprofiles_btn.set_active(True)
		self.arbitrary_profiles_btn = gtk.RadioButton(self.plotXYprofiles_btn,"Arbitrary profiles")


		self.q_space_btn = gtk.CheckButton("Q space") #Plot the map in the reciprocal space Qx Qz
		self.q_space_btn.set_sensitive(False)
		#self.q_space_btn.connect("toggled",self.change_space,"Q")

		self.show_chi_delta_btn = gtk.CheckButton("Show 2Theta, Chi, d")
		self.show_chi_delta_btn.connect("toggled",self.show_chi_delta)
		self.show_chi_delta_flag = False
		self.show_chi_txt   = gtk.Label()
		self.show_chi_txt.set_alignment(0,0)
		self.show_delta_txt = gtk.Label()
		self.show_delta_txt.set_alignment(0,0)
		self.show_d_txt = gtk.Label()
		self.show_d_txt.set_alignment(0,0)
		#### Pack these options in a table
		self.option_table = gtk.Table(5,3,False) #5 rows, 3 cols, homogeneous
		self.option_table.attach(self.detector_disposition_horizontal, 0,1,0,1)
		self.option_table.attach(self.linear_scale_btn, 1,2,0,1)
		self.option_table.attach(self.save_adj_btn, 2,3,0,1)
		self.option_table.attach(self.adj_btn, 1,2,1,2)
		self.option_table.attach(self.cln_btn, 1,2,2,3)
		self.option_table.attach(self.tth_chi_space_btn,0,1,1,2)
		self.option_table.attach(self.q_space_btn,0,1,2,3)
		self.option_table.attach(self.plotXYprofiles_btn,0,1,3,4)
		self.option_table.attach(self.arbitrary_profiles_btn,1,2,3,4)
		self.option_table.attach(self.show_chi_delta_btn,2,3,1,2)
		self.option_table.attach(self.show_delta_txt,2,3,2,3)
		self.option_table.attach(self.show_chi_txt, 2,3,3,4)
		self.option_table.attach(self.show_d_txt, 2,3,4,5)


		### Options for profile plots
		self.profiles_log_btn = gtk.ToggleButton("Y-Log")
		self.profiles_log_btn.connect("toggled",self.profiles_update)
		self.profiles_export_data_btn = gtk.Button("Export data")
		self.profiles_export_data_btn.connect("clicked",self.profiles_export)

		self.profiles_option_box = gtk.HBox(False,0)
		self.profiles_option_box.pack_start(self.profiles_log_btn, False, False, 0)
		self.profiles_option_box.pack_start(self.profiles_export_data_btn, False, False, 0)
		### Figure of profiles plot
		self.fig_profiles = Figure()
		self.profiles_ax1 = self.fig_profiles.add_subplot(211)
		self.profiles_ax1.set_title("Y profile", size=12)
		self.profiles_ax2 = self.fig_profiles.add_subplot(212)
		self.profiles_ax2.set_xlabel("X profile", size=12)
		self.profiles_canvas = FigureCanvas(self.fig_profiles)
		self.profiles_canvas.set_size_request(450,50)
		self.profiles_navBar = NavigationToolbar(self.profiles_canvas, self)
		self.cursor_pro1 = Cursor(self.profiles_ax1, color='k', linewidth=1, useblit=True)
		self.cursor_pro2 = Cursor(self.profiles_ax2, color='k', linewidth=1, useblit=True)


		#### Results of fitted curves
		self.fit_results_table = gtk.Table(7,3, False)
		title = gtk.Label("Fitted results:")
		self.chi_title = gtk.Label("Y")
		self.tth_title = gtk.Label("X")
		y0 = gtk.Label("y0:")
		xc = gtk.Label("xc:")
		A = gtk.Label("A:")
		w = gtk.Label("FWHM:")
		mu = gtk.Label("mu:")
		y0.set_alignment(0,0.5)
		xc.set_alignment(0,0.5)
		A.set_alignment(0,0.5)
		w.set_alignment(0,0.5)
		mu.set_alignment(0,0.5)

		self.chi_fitted_y0 = gtk.Label()
		self.chi_fitted_xc = gtk.Label()
		self.chi_fitted_A = gtk.Label()
		self.chi_fitted_w = gtk.Label()
		self.chi_fitted_mu = gtk.Label()

		self.tth_fitted_y0 = gtk.Label()
		self.tth_fitted_xc = gtk.Label()
		self.tth_fitted_A = gtk.Label()
		self.tth_fitted_w = gtk.Label()
		self.tth_fitted_mu = gtk.Label()

		self.fit_results_table.attach(title,0,3,0,1)
		self.fit_results_table.attach(self.chi_title,1,2,1,2)
		self.fit_results_table.attach(self.tth_title,2,3,1,2)
		self.fit_results_table.attach(y0,0,1,2,3)
		self.fit_results_table.attach(xc,0,1,3,4)
		self.fit_results_table.attach(A,0,1,4,5)
		self.fit_results_table.attach(w,0,1,5,6)
		self.fit_results_table.attach(mu,0,1,6,7)

		self.fit_results_table.attach(self.chi_fitted_y0,1,2,2,3)
		self.fit_results_table.attach(self.chi_fitted_xc,1,2,3,4)
		self.fit_results_table.attach(self.chi_fitted_A,1,2,4,5)
		self.fit_results_table.attach(self.chi_fitted_w,1,2,5,6)
		self.fit_results_table.attach(self.chi_fitted_mu,1,2,6,7)

		self.fit_results_table.attach(self.tth_fitted_y0,2,3,2,3)
		self.fit_results_table.attach(self.tth_fitted_xc,2,3,3,4)
		self.fit_results_table.attach(self.tth_fitted_A,2,3,4,5)
		self.fit_results_table.attach(self.tth_fitted_w,2,3,5,6)
		self.fit_results_table.attach(self.tth_fitted_mu,2,3,6,7)

		#### PACK the right panel
		self.right_panel.pack_start(self.option_table, False, False, 0)
		self.right_panel.pack_start(self.profiles_option_box,False,False,0)
		self.right_panel.pack_start(self.profiles_navBar,False,False,0)
		self.right_panel.pack_start(self.profiles_canvas,True,True,0)
		self.right_panel.pack_start(self.fit_results_table, False, False, 0)

		#hbox.pack_end(fixed,False,False,5)
		self.page_single_figure.pack_end(self.right_panel,False, False,5)
		self.notebook.append_page(self.page_single_figure,gtk.Label("EDF viewer"))
		self.notebook.append_page(self.page_pole_figure,gtk.Label("Pole figure"))
		self.notebook.append_page(self.page_batch_correction,gtk.Label("EDF batch correction"))
		###################### PAGE 2: POLE FIGURES ###################################
		self.left_panel_polefigure = gtk.VBox() #This box is to contain options for pole figure construction
		self.right_panel_polefigure = gtk.VBox()

		self.images_used = gtk.Label("Use images")
		self.images_used.set_alignment(0,0.5)
		self.images_from = gtk.Label("From #: ")
		self.images_from_nb = gtk.Entry()
		self.images_from_nb.set_usize(30,0)
		self.images_from.set_alignment(0,0.5)

		self.images_to   = gtk.Label("To #: ")
		self.images_to_nb = gtk.Entry()
		self.images_to_nb.set_usize(30,0)
		self.images_to.set_alignment(0,0.5)

		self.pole_2theta_txt  = gtk.Label("2 Theta: ")
		self.pole_2theta_field = gtk.Entry()
		self.pole_2theta_field.set_usize(30,0)
		self.pole_2theta_txt.set_alignment(0,0.5)
		
		self.pole_chi_min_txt  = gtk.Label("Chi minimum: ")
		self.PF_chi_min = gtk.Entry()
		self.PF_chi_min.set_usize(30,0)
		self.pole_chi_min_txt.set_alignment(0,0.5)

		self.select_phi = gtk.RadioButton(None, "Plot PHI")
		self.select_kphi= gtk.RadioButton(self.select_phi, "Plot KPHI")

		#self.check_chi_positive = gtk.CheckButton("Chi positive only")
		#self.check_chi_positive.set_size_request(180,30)
		#self.check_chi_positive.connect("toggled", self.chi_positive_only)

		self.plot_pole_figure_btn = gtk.Button("Plot")
		self.plot_pole_figure_btn.set_size_request(60,30)
		self.plot_pole_figure_btn.connect("clicked",self.plot_pole_figure)
		self.log_pole_figure_btn = gtk.ToggleButton("Log scale")
		self.log_pole_figure_btn.set_size_request(60,30)
		self.log_pole_figure_btn.connect("toggled", self.replot_PF)
		self.data_loading = gtk.ProgressBar()

		self.input_table = gtk.Table(5,2,True) #rows, cols, homogeneous
		self.input_table.attach(self.images_from, 0,1,0,1)
		self.input_table.attach(self.images_from_nb, 1,2,0,1)
		self.input_table.attach(self.images_to, 0,1,1,2)
		self.input_table.attach(self.images_to_nb,1,2,1,2)
		self.input_table.attach(self.pole_2theta_txt, 0,1,2,3)
		self.input_table.attach(self.pole_2theta_field, 1,2,2,3)
		#self.input_table.attach(self.pole_chi_min_txt, 0,1,3,4)
		#self.input_table.attach(self.PF_chi_min, 1,2,3,4)
		self.input_table.attach(self.select_phi, 0,1,4,5)
		self.input_table.attach(self.select_kphi, 1,2,4,5)

		self.left_panel_polefigure.pack_start(self.images_used, False, False, 5)
		self.left_panel_polefigure.pack_start(self.input_table, False, False, 2)
		#self.left_panel_polefigure.pack_start(self.check_chi_positive, False, False, 2)
		self.left_panel_polefigure.pack_start(self.plot_pole_figure_btn, False, False, 2)
		#self.left_panel_polefigure.pack_start(self.log_pole_figure_btn, False, False, 2)
		self.left_panel_polefigure.pack_start(self.data_loading, False, False, 2)

		self.page_pole_figure.pack_start(self.left_panel_polefigure,False,False,10)

		self.polefig=Figure(figsize=(10,8),dpi=100)
		self.polar_ax  = self.polefig.add_axes([0.1,0.1,0.75,0.8], polar='True')
		#self.polar_ax  = self.polefig.add_axes([0.1,0.1,0.75,0.8])
		#self.polar_ax.set_frame_on(False)
		self.polar_cax = self.polefig.add_axes([0.85,0.1,0.03,0.8])
		self.pole_canvas  = FigureCanvas(self.polefig)
		#self.PF_cursor = Cursor(self.polar_ax, color='k', linewidth=1, useblit=True)
		self.PF_navBar = NavigationToolbar(self.pole_canvas, self)

		self.right_panel_polefigure.pack_start(self.PF_navBar,False,False,0)
		self.right_panel_polefigure.pack_start(self.pole_canvas,True,True,0)
		self.page_pole_figure.pack_start(self.right_panel_polefigure,True, True,10)
		
		#*****************************************************************************
		#                  PAGE 3: CORRECTION OF A BATCH EDF
		#*****************************************************************************
		self.table_1 = gtk.Table(2,3, False)
		self.table_2 = gtk.Table(6,5, False)
		self.t1_src_folder_txt= gtk.Label("Source folder:")
		self.t1_src_folder_txt.set_alignment(0,0.5)
		self.t1_src_path  = gtk.Entry()
		self.t1_src_path.set_usize(100,0)
		
		self.src_dialog = gtk.FileChooserDialog(title="Select a sources folder",action=gtk.FILE_CHOOSER_ACTION_SELECT_FOLDER, buttons = (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL, gtk.STOCK_OPEN, gtk.RESPONSE_OK))
		self.src_dialog.set_current_folder(self.current_folder)
		self.t1_src_button= gtk.Button("Browse")
		#self.t1_src_button.set_title("Browse")
		#self.t1_src_button.connect("clicked", self.select_folder,self.src_dialog,self.t1_src_path, "S")
		self.t1_src_button.connect("clicked", self.select_source_folder)
		
		self.t1_des_folder_txt = gtk.Label("Destination folder:")
		self.t1_des_folder_txt.set_alignment(0,0.5)
		self.t1_des_path = gtk.Entry()
		self.t1_des_path.set_usize(100,0)
		self.des_dialog = gtk.FileChooserDialog(title="Select a destination folder",action=gtk.FILE_CHOOSER_ACTION_SELECT_FOLDER, buttons = (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL, gtk.STOCK_OPEN, gtk.RESPONSE_OK))
		self.des_dialog.set_current_folder(self.current_folder)
		self.t1_des_button = gtk.Button("Browse")
		#self.t1_des_button.set_title("Browse")
		#self.t1_des_button.connect("clicked", self.select_folder, self.des_dialog, self.t1_des_path, "D")
		self.t1_des_button.connect("clicked", self.select_destination_folder)
		
		self.t2_detector_txt = gtk.Label("Detector: ")
		self.t2_detector_txt.set_alignment(0,0.5)
		self.t2_det_combobox = gtk.combo_box_new_text()
		self.t2_det_combobox.append_text("D1")
		self.t2_det_combobox.append_text("D5")
		self.t2_det_combobox.append_text("S70")
		self.t2_det_combobox.set_active(1)
		self.t2_det_combobox.connect("changed", self.batch_change_detector)
		
		self.d1_specfile_txt = gtk.Label("Spec file:")
		self.d1_specfile_txt.set_alignment(1,0.5)
		self.d1_specfile_path = gtk.Entry()
		self.d1_specfile_browse = gtk.Button("Browse")
		self.d1_specfile_browse.connect("clicked", self.select_normalisation_file,self.d1_specfile_path)
		
		self.d1_center_txt = gtk.Label("Center, separated by comma (X,Y):")
		self.d1_center_txt.set_alignment(1,0.5)
		self.d1_center = gtk.Entry()
		
		self.d1_distance_txt = gtk.Label("Distance (mm):")
		self.d1_distance_txt.set_alignment(1,0.5)
		self.d1_distance = gtk.Entry()
		
		self.t2_monitor_txt = gtk.Label("Monitor counter: ")
		self.t2_monitor_txt.set_alignment(0,0.5)
		self.t2_mon_combobox = gtk.combo_box_new_text()
		self.t2_mon_combobox.append_text("vct1")
		self.t2_mon_combobox.append_text("vct2")
		self.t2_mon_combobox.append_text("vct3")
		self.t2_mon_combobox.append_text("d0_cps")
		self.t2_mon_combobox.append_text("pm01")
		self.t2_mon_combobox.append_text("pm02")
		self.t2_mon_combobox.set_active(0)
		
		self.t2_ref_mon_txt = gtk.Label("Reference monitor (?):")
		self.t2_ref_mon_txt.set_alignment(0,0.5)
		self.t2_ref_mon_entry = gtk.Entry()
		self.t2_img_start_txt = gtk.Label("Image begin (?):")
		self.t2_img_start_txt.set_alignment(0,0.5)
		self.t2_img_start_entry = gtk.Entry()
		self.t2_img_end_txt = gtk.Label("Image end (?):")
		self.t2_img_end_txt.set_alignment(0,0.5)
		self.t2_img_end_entry = gtk.Entry()
		self.t2_ascii_out = gtk.Label("Save data as ASCII file")
		self.t2_ascii_out.set_alignment(0,0.5)
		self.ascii_out = gtk.CheckButton()
		
		self.t2_tooltip = gtk.Tooltips()
		self.t2_tooltip.set_tip(self.t2_ref_mon_txt, "Reference value of the monitor which corresponds to the highest value of the current machine. This is for the data normalization. If empty, the data will not be normalized")
		self.t2_tooltip.set_tip(self.t2_img_start_txt, "Starting image to be corrected. If empty the whole folder will be proceeded.")
		self.t2_tooltip.set_tip(self.t2_img_end_txt, "Ending image to be corrected. If empty the whole folder will be proceeded.")
		
		self.run_button = gtk.Button("RUN")
		self.run_button.connect("clicked", self.batch_transform)
		self.run_button.set_usize(30,0)
		self.run_button.set_alignment(0.5,0.5)
		separator = gtk.HSeparator()
		self.show_process = gtk.Label("Processing:")
		self.show_process.set_alignment(0,0.5)
		self.show_process_info = gtk.Label()
		self.progressbar = gtk.ProgressBar()
		
		
		#**** Packaging *************************************************************
		self.table_1.attach(self.t1_src_folder_txt, 0,1,0,1)
		self.table_1.attach(self.t1_src_path, 1,2,0,1)
		self.table_1.attach(self.t1_src_button, 2,3,0,1)
		self.table_1.attach(self.t1_des_folder_txt, 0,1,1,2)
		self.table_1.attach(self.t1_des_path, 1,2,1,2)
		self.table_1.attach(self.t1_des_button, 2,3,1,2)
		
		self.table_2.attach(self.t2_detector_txt, 0,1,0,1)
		self.table_2.attach(self.t2_det_combobox, 1,2,0,1)
		self.table_2.attach(self.t2_monitor_txt, 0,1,1,2)
		self.table_2.attach(self.t2_mon_combobox, 1,2,1,2)
		self.table_2.attach(self.t2_ref_mon_txt, 0,1,2,3)
		self.table_2.attach(self.t2_ref_mon_entry, 1,2,2,3)
		self.table_2.attach(self.t2_img_start_txt, 0,1,3,4)
		self.table_2.attach(self.t2_img_start_entry, 1,2,3,4)
		self.table_2.attach(self.t2_img_end_txt, 0,1,4,5)
		self.table_2.attach(self.t2_img_end_entry, 1,2,4,5)
		self.table_2.attach(self.t2_ascii_out, 0,1,5,6)
		self.table_2.attach(self.ascii_out, 1,2,5,6)
		
		self.table_2.attach(self.d1_specfile_txt,2,3,0,1)
		self.table_2.attach(self.d1_specfile_path,3,4,0,1)
		self.table_2.attach(self.d1_specfile_browse,4,5,0,1)
		self.table_2.attach(self.d1_center_txt,2,3,1,2)
		self.table_2.attach(self.d1_center,3,4,1,2)
		self.table_2.attach(self.d1_distance_txt,2,3,2,3)
		self.table_2.attach(self.d1_distance,3,4,2,3)
		
		
		self.page_batch_correction_container = gtk.VBox()
		self.page_batch_correction_container.pack_start(self.table_1, False, False, 10)
		self.page_batch_correction_container.pack_start(self.table_2, False, False, 10)
		self.page_batch_correction_container.pack_start(self.run_button, False, False, 10)
		#self.page_batch_correction_container.pack_start(separator, False, False, 10)
		self.table_3 = gtk.Table(1,2,False)
		self.table_3.attach(self.show_process, 0,1,0,1)
		self.table_3.attach(self.show_process_info,1,2,0,1)
		self.page_batch_correction_container.pack_start(self.table_3, False, False, 10)
		self.page_batch_correction_container.pack_start(self.progressbar, False, False, 10)
		self.page_batch_correction.pack_start(self.page_batch_correction_container, False, False, 20)
		
		###################### PACK THE NOTEBOOK #####################################
		hbox.pack_start(self.notebook)
		vbox.pack_start(hbox,True,True,0)
		############################### Sliders ######################################
		#sld_box = gtk.Fixed()
		sld_box = gtk.HBox(False,2)

		self.vmin_txt = gtk.Label("Vmin")
		self.vmin_txt.set_alignment(0,0.5)
		#self.vmin_txt.set_justify(gtk.JUSTIFY_CENTER)
		self.vmax_txt = gtk.Label("Vmax")
		self.vmax_txt.set_alignment(0,0.5)
		#self.vmax_txt.set_justify(gtk.JUSTIFY_CENTER)
		self.sld_vmin = gtk.HScale()
		self.sld_vmax = gtk.HScale()

		self.sld_vmin.set_size_request(200,25)
		self.sld_vmax.set_size_request(200,25)
		self.sld_vmin.set_range(0,self.vmax)
		self.sld_vmax.set_range(0,self.vmax)
		self.sld_vmax.set_value(self.vmax)
		self.sld_vmin.set_value(0)
		self.sld_vmin.connect('value-changed',self.scale_update)
		self.sld_vmax.connect('value-changed',self.scale_update)

		vmax_spin_adj         = gtk.Adjustment(self.vmax, 0, self.vmax_range, 0.5, 10.0, 0.0)
		self.vmax_spin_btn    = gtk.SpinButton(vmax_spin_adj,1,1)
		self.vmax_spin_btn.set_numeric(True)
		self.vmax_spin_btn.set_wrap(True)
		self.vmax_spin_btn.set_size_request(80,-1)
		#self.vmax_spin_btn.set_alignment(0,0.5)
		self.vmax_spin_btn.connect('value-changed',self.scale_update_spin)

		vmin_spin_adj         = gtk.Adjustment(self.vmin, 0, self.vmax_range, 0.5, 10.0, 0.0)
		self.vmin_spin_btn    = gtk.SpinButton(vmin_spin_adj,1,1)
		self.vmin_spin_btn.set_numeric(True)
		self.vmin_spin_btn.set_wrap(True)
		self.vmin_spin_btn.set_size_request(80,-1)
		#self.vmax_spin_btn.set_alignment(0,0.5)
		self.vmin_spin_btn.connect('value-changed',self.scale_update_spin)

		self.slider_reset_btn = gtk.Button("Reset")
		self.slider_reset_btn.set_size_request(50,30)
		self.slider_reset_btn.set_alignment(1,0.5)
		self.slider_reset_btn.connect('clicked',self.reset_scale)

		sld_box.pack_start(self.vmin_txt,False,False,0)
		sld_box.pack_start(self.sld_vmin,False,False,0)
		sld_box.pack_start(self.vmin_spin_btn,False,False,0)
		sld_box.pack_start(self.vmax_txt,False,False,0)
		sld_box.pack_start(self.sld_vmax,False,False,0)
		sld_box.pack_start(self.vmax_spin_btn,False,False,0)
		sld_box.pack_start(self.slider_reset_btn,False,False,0)

		vbox.pack_start(sld_box,False,False,3)

		################# Status bar #################################################
		#self.status_bar = gtk.EventBox()
		self.status_bar = gtk.HBox(False,2)
		#self.status_bar.modify_bg(gtk.STATE_NORMAL, self.status_bar.get_colormap().alloc_color("white"))
		self.stt = gtk.Fixed()
		self.x_pos = gtk.Label("X =")
		self.y_pos = gtk.Label("Y =")
		self.z_pos = gtk.Label("Z =")
		#self.edf_pos = gtk.Label("EDF choosen: ")
		self.del_pos = gtk.Label()
		self.eta_pos = gtk.Label()
		self.phi_pos = gtk.Label()
		self.kphi_pos = gtk.Label()
		self.chi_pos = gtk.Label()
		self.nu_pos  = gtk.Label()
		self.mu_pos  = gtk.Label()
		self.time_pos  = gtk.Label()
		self.stt.put(self.x_pos,10,5)
		self.stt.put(self.y_pos,70,5)
		self.stt.put(self.z_pos,130,5)
		#self.stt.put(self.edf_pos,250,5)
		self.stt.put(self.del_pos,220,5)
		self.stt.put(self.eta_pos,325,5)
		self.stt.put(self.phi_pos,425,5)
		self.stt.put(self.kphi_pos,525,5)
		self.stt.put(self.chi_pos,625,5)
		self.stt.put(self.nu_pos,710,5)
		self.stt.put(self.mu_pos,790,5)
		self.stt.put(self.time_pos,880, 5)
		self.status_bar.add(self.stt)
		vbox.pack_start(self.status_bar,False,False,0)


		self.add(vbox)
		self.connect("destroy", gtk.main_quit)
		self.show_all()
		self.progressbar.hide()
		self.data_loading.hide()
		self.d1_center.hide()
		self.d1_center_txt.hide()
		self.d1_distance.hide()
		self.d1_distance_txt.hide()
		self.d1_specfile_browse.hide()
		self.d1_specfile_path.hide()
		self.d1_specfile_txt.hide()

######################################## Definitions ########################################################################
	def pro_format_coord(self,x,y):
		return 'x=%.4f, y=%.1f'%(x,y)

	def init_image(self):
		self.ax.clear()
		self.cax.clear()
		self.cax2.clear()
		self.ax.add_patch(self.rect)
		self.img = self.ax.imshow(self.data,origin='lower',vmin=self.vmin, vmax=self.vmax, cmap=jet, interpolation='nearest',aspect='auto')
		self.img.set_extent(self.MAIN_EXTENT)
		if self.log_scale == 0:
			clabel = r'$Intensity\ (Counts\ per\ second)$'
		else:
			clabel = r'$Log_{10}\ (Counts\ per\ second)\ [arb.\ units]$'
		self.cb  = self.fig.colorbar(self.img, cax = self.cax, format="%.f")#format=fm
		self.cb.set_label(clabel, fontsize=18)
		self.cb.locator = MaxNLocator(nbins=6)
		#Horizontal colorbar
		self.cb2  = self.fig.colorbar(self.img, cax = self.cax2, orientation='horizontal', format="%.f")#format=fm
		self.cb2.set_label(clabel, fontsize=18)
		self.cb2.locator = MaxNLocator(nbins=6)
		self.cb2.ax.set_visible(False)
		
		if self.tth_chi_space_btn.get_active():
			self.MAIN_XLABEL.set_text("2Theta (deg.)")
			self.MAIN_YLABEL.set_text("Chi (deg.)")
		else:
			self.MAIN_XLABEL.set_text("X (pixel)")
			self.MAIN_YLABEL.set_text("Y (pixel)")
		
		if self.graph_aspect == True:
			self.ax.set_aspect('equal')
		else:
			self.ax.set_aspect('auto')
			
		self.IMG_INIT = True
		#self.cursor = Cursor(self.ax, color='k', linewidth=1, useblit=True)
		
	def change_aspect_ratio(self,w):
		self.graph_aspect = not (self.graph_aspect)
		if self.graph_aspect == True:
			self.ax.set_aspect('equal')
		else:
			self.ax.set_aspect('auto')
		self.canvas.draw()

	def popup_info(self,info_type,text):
		""" info_type = WARNING, INFO, QUESTION, ERROR """
		if info_type.upper() == "WARNING":
			mess_type = gtk.MESSAGE_WARNING
		elif info_type.upper() == "INFO":
			mess_type = gtk.MESSAGE_INFO
		elif info_type.upper() == "ERROR":
			mess_type = gtk.MESSAGE_ERROR
		elif info_type.upper() == "QUESTION":
			mess_type = gtk.MESSAGE_QUESTION

		self.warning=gtk.MessageDialog(self, gtk.DIALOG_DESTROY_WITH_PARENT, mess_type, gtk.BUTTONS_CLOSE,text)
		self.warning.run()
		self.warning.destroy()

	def calibration(self, widget):
		""" Choose an edf image to do calibration
		The calibration uses the pyFAI-calib program writen by Jerome Kieffer
		This program takes an EDF image as calibrant, and gives a calibration file (or poni file) xxxx.poni
		The poni file contains: sample-detector distance, direct beam center, rotation angles in radians around the 3 axis"""
		cmd = "pyFAI-calib -p 130 -S LaB6.D -w 1.512 "+self.edf
		os.system(cmd)
    
	def load_calibration(self, widget):
		""" load a pre-calib file , PONI file """
		dialog = gtk.FileChooserDialog("Select a PONI file",None,gtk.FILE_CHOOSER_ACTION_OPEN,(gtk.STOCK_CANCEL,gtk.RESPONSE_CANCEL, gtk.STOCK_OPEN, gtk.RESPONSE_OK))
		filtre = gtk.FileFilter()
		filtre.set_name("PONI")
		filtre.add_pattern("*.poni")
		dialog.add_filter(filtre)
		response = dialog.run()
		if response == gtk.RESPONSE_OK:
			self.ponifile = dialog.get_filename()
			print "Calibration file is: ",self.ponifile
			self.azimuthalIntegration = pyFAI.load(self.ponifile)
			#if self.detector_type=="D5":
				#self.tableTwoTheta = self.azimuthalIntegration.twoThetaArray((578,1148))
				#self.tableChi      = self.azimuthalIntegration.chiArray((578,1148))
				#self.tableChi      = N.degrees(self.tableChi)-90
			
			#elif self.detector_type=="S70":
				##self.azimuthalIntegration.setChiDiscAtZero()
				#self.tableTwoTheta = self.azimuthalIntegration.twoThetaArray((120,578))
				#self.tableChi      = self.azimuthalIntegration.chiArray((120,578))
				#self.tableChi      = N.degrees(self.tableChi)+90
			#elif self.detector_type=="D1":
				#self.tableTwoTheta = self.azimuthalIntegration.twoThetaArray((577,913))
				#self.tableChi      = self.azimuthalIntegration.chiArray((577,913))
				#self.tableChi      = N.degrees(self.tableChi)-90

				
			#self.table_dSpace = self.azimuthalIntegration.wavelength / (2*N.sin(self.tableTwoTheta/2.0)) * 1e10 #d in Angstrom
			#self.tableTwoTheta = N.degrees(self.tableTwoTheta)
			##self.tableChi      = self.azimuthalIntegration.chiArray((578,1148))
			##self.tableChi      = N.degrees(self.tableChi)-90
			##self.tableQ        = self.azimuthalIntegration.qArray((578,1148))

			self.calibrated = True
			self.calibrated_quantitative = True
			print "is calibrated? ",self.calibrated
			self.geometry_manual.set_active(False)
			self.canvas.draw()
			s = self.ponifile.split("/")[-1]
			self.popup_info("info","This detector is calibrated with the PONI file %s!!!"%s)
		else:
			pass
		dialog.destroy()
	
	def manual_calibration(self,widget):
		#***** Check input data
		if self.geometry_manual.get_active():
			distance = self.geometry_distance.get_text()
			energy   = self.geometry_energy.get_text()
			direct_beam = self.geometry_direct_beam.get_text()
			distance = float(distance)
			energy   = float(energy)
			direct_beam = direct_beam.split(",")
			direct_beam = [float(direct_beam[0]),float(direct_beam[1])]
			from scipy.constants import h,c,e
			poni1 = direct_beam[1]*_PIXEL_SIZE/1000.
			poni2 = direct_beam[0]*_PIXEL_SIZE/1000.
			wavelength = h*c/e/energy
			
			self.azimuthalIntegration = pyFAI.azimuthalIntegrator.AzimuthalIntegrator(dist=distance,
																					  poni1=poni1,
																					  poni2=poni2,
																					  rot1=None,
																					  rot2=None,
																					  rot3=None,
																					  pixel1=_PIXEL_SIZE/1000.,
																					  pixel2=_PIXEL_SIZE/1000.,
																					  wavelength=wavelength)
			self.calibrated=True
			self.calibrated_quantitative = False
			#self.calculation_angular_coordinates()
		else:
			self.calibrated=False
			
	def calculation_angular_coordinates(self):
		if self.detector_type=="D5":
			self.tableTwoTheta = self.azimuthalIntegration.twoThetaArray((578,1148))
			self.tableChi      = self.azimuthalIntegration.chiArray((578,1148))
			self.tableChi      = N.degrees(self.tableChi)-90
		
		elif self.detector_type=="S70":
			#self.azimuthalIntegration.setChiDiscAtZero()
			self.tableTwoTheta = self.azimuthalIntegration.twoThetaArray((120,578))
			self.tableChi      = self.azimuthalIntegration.chiArray((120,578))
			self.tableChi      = N.degrees(self.tableChi)+90
		elif self.detector_type=="D1":
			self.tableTwoTheta = self.azimuthalIntegration.twoThetaArray((577,913))
			self.tableChi      = self.azimuthalIntegration.chiArray((577,913))
			self.tableChi      = N.degrees(self.tableChi)-90

			
		self.table_dSpace = self.azimuthalIntegration.wavelength / (2*N.sin(self.tableTwoTheta/2.0)) * 1e10 #d in Angstrom
		self.tableTwoTheta = N.degrees(self.tableTwoTheta)
		#self.tableChi      = self.azimuthalIntegration.chiArray((578,1148))
		#self.tableChi      = N.degrees(self.tableChi)-90
		#self.tableQ        = self.azimuthalIntegration.qArray((578,1148))
	def set_detector(self,widget, det_type):
		if det_type.upper()=="D5":
			print "D5 selected"
			self.detm.set_label("XPAD D5")
			self.detector_type = "D5"

		elif det_type.upper()=="S70":
			print "S70 selected"
			self.detm.set_label("XPAD S70")
			self.detector_type = "S70"
			self.xDim1 = 560
			self.yDim1 = 120

		elif det_type.upper()=="D1":
			print "D1 selected"
			self.detm.set_label("XPAD D1")
			self.detector_type = "D1"
			self.xDim1 = 577
			self.yDim1 = 913

		self.init_image()
		self.canvas.draw()
		return 0

	def check_azimuthal_integrator(self):
		if not self.calibrated_quantitative:
			rot1 = N.radians(self.nu)*(-1.)
			rot2 = N.radians(self.delta)*(-1.)
			rot3 = N.radians(90-self.chi)
			self.azimuthalIntegration.rot1 = rot1
			self.azimuthalIntegration.rot2 = rot2
			self.azimuthalIntegration.rot3 = rot3
		self.calculation_angular_coordinates()
		
	def calc_qxqz(self):
		""" Calculate Qx Qz """
		if self.calibrated:
			self.check_azimuthal_integrator()
			psi = N.zeros(shape=self.tableTwoTheta.shape)
			Qmod= self.tableQ
			Qx  = Qz = psi

			eta = N.ones(shape=self.tableTwoTheta.shape)
			eta = eta * self.eta
			psi = eta - self.tableTwoTheta/2.0
			Qx  = Qmod * N.sin(N.radians(psi))
			Qz  = Qmod * N.cos(N.radians(psi))

			#self.Qx = Qx
			#self.Qz = Qz
			return Qx, Qz
		else:
			self.popup_info('warning','Impossible to calculate Qx Qz because the detector is not calibrated')
			return None,None

	def on_changed_edf(self,widget,row,col):

		self.clear_notes()
		self.init_image()
		model = widget.get_model()
		#print "TRA CHOOSE: ",model
		#print "ROW: ",row,len(row)
		#print "MODEL[ROW]: ",model[row]
		self.edf_choosen = model[row][0]
		self.edf_folder = self.current_folder
		if len(row)==2:
			#print "Selected Directory: ",model[row[0]][0]
			self.edf_folder = join(self.current_folder,model[row[0]][0])
			#print "edf_folder: ",self.edf_folder
		self.edf = join(self.edf_folder,self.edf_choosen)
		#self.edf_pos.set_text("EDF choosen:  %s"%self.edf_choosen)
		try:
			self.MAIN_TITLE.set_text(self.edf_choosen,fontsize=18)
		except:
			self.MAIN_TITLE = self.ax.set_title(self.edf_choosen,fontsize=18)
		#self.MAIN_TITLE = self.fig.suptitle(self.edf_choosen)
		### Data Loading #########
		self.fabioIMG = fabio.open(self.edf)
		self.header = self.fabioIMG.header
		
		#print self.header
		if self.detector_type != "D1":
			self.counter = get_counters(self.header)
			self.motor = get_motors(self.header)
			motor_mne = self.motor.keys()
			if len(motor_mne)>1:
				if 'xsamp' in motor_mne:
					manip = "gisaxs"
				else:
					manip = "kappapsic"
			else:
				manip = "Unknown"
			if manip == "kappapsic":
				self.delta = self.motor['del']
				self.eta   = self.motor['eta']
				self.chi   = self.motor['chi']
				self.phi   = self.motor['phi']
				self.nu    = self.motor['nu']
				self.mu    = self.motor['mu']
				self.kphi  = self.motor['kphi']
			elif manip == "gisaxs":
				self.delta = 0
				self.eta = 0
				self.phi = self.kphi = self.chi = self.nu = self.mu = 0			
			self.count_time = self.counter['sec']
			#
			if manip == "kappapsic":
				self.del_pos.set_text("%s: %.2f"%("Del",self.delta))
				self.eta_pos.set_text("%s: %.2f"%("Eta",self.eta))
				self.phi_pos.set_text("%s: %.2f"%("Phi",self.phi))
				self.kphi_pos.set_text("%s: %.2f"%("Kphi",self.kphi))
				self.chi_pos.set_text("%s: %.2f"%("Chi",self.chi))
				self.nu_pos.set_text("%s: %.2f"%("Nu",self.nu))
				self.mu_pos.set_text("%s: %.2f"%("Mu",self.mu))
			elif manip == "gisaxs":
				self.del_pos.set_text("%s: %.2f"%("Xsamp",self.motor['xsamp']))
				self.eta_pos.set_text("%s: %.2f"%("Zsamp",self.motor['zsamp']))
				self.phi_pos.set_text("%s: %.2f"%("Rsamp",self.motor['rsamp']))
				self.kphi_pos.set_text("%s: %.2f"%("XstoP",self.motor['XstoP']))
				self.chi_pos.set_text("%s: %.2f"%("ZstoP",self.motor['ZstoP']))
				self.nu_pos.set_text("%s: %.2f"%("Xdet",self.motor['Xdet']))
				self.mu_pos.set_text("%s: %.2f"%("Zdet",self.motor['Zdet']))
			self.time_pos.set_text("Seconds:  %d"%self.count_time)
		gc.collect() # Clear unused variables
		#self.data = self.fabioIMG.data
		self.plot_data()
		num = self.edf_choosen.split("_")[1]
		num = num.split(".")[0]
		self.SELECTED_IMG_NUM = int(num)
		if self.SPEC_IS_LOADED and len(self.SPEC_IMG) != 0:
			self.check_and_update_scan_slider()#Check the scan number and image number --> set the spin button for scan num and the slider for img num
		return

	def load_specFile(self,widget):
		dialog = gtk.FileChooserDialog("Select spec file",None,gtk.FILE_CHOOSER_ACTION_OPEN,(gtk.STOCK_CANCEL,gtk.RESPONSE_CANCEL, gtk.STOCK_OPEN, gtk.RESPONSE_OK))
		dialog.set_current_folder(self.current_folder)
		response = dialog.run()
		if response == gtk.RESPONSE_OK:
			file_choosen = dialog.get_filename()
			self.scan_slider_path.set_text(file_choosen)
			self.SPEC_FILE = file_choosen
		else:
			pass
		dialog.destroy()
		while gtk.events_pending():
			gtk.main_iteration()
		self.SPEC_DATA = xrayutilities.io.SPECFile(self.SPEC_FILE)
		self.update_spec_data()
		
		return
	
	def update_spec_data(self):
		self.SPEC_IMG = []
		self.SPEC_SCAN_LIST = self.SPEC_DATA.scan_list
		if _SPEC_IMG_COL not in self.SPEC_SCAN_LIST[0].colnames:
			self.popup_info("error","Spec file does not containt image field. The default image coloumn is 'img'")
			return
		else:
			first_scan_num = self.SPEC_SCAN_LIST[0].nr
			last_scan_num = self.SPEC_SCAN_LIST[-1].nr
			self.SPEC_SCAN_RANGE = (first_scan_num, last_scan_num)
			#print first_scan_num, last_scan_num
			self.SPEC_SCAN_NUM_LIST = []
			for i in range(len(self.SPEC_SCAN_LIST)):
				#item = self.SPEC_SCAN_LIST[i]
				self.SPEC_SCAN_NUM_LIST.append(self.SPEC_SCAN_LIST[i].nr)
				if self.SPEC_SCAN_LIST[i].scan_status == 'OK':
					self.SPEC_SCAN_LIST[i].ReadData()
					this_img_list = self.SPEC_SCAN_LIST[i].data[_SPEC_IMG_COL]
					this_img_list = this_img_list.astype('int')
				else:
					this_img_list = []
				self.SPEC_IMG.append(this_img_list)
				
			self.SPEC_IS_LOADED = True
			self.scan_slider_spinButton.set_range(self.SPEC_SCAN_RANGE[0], self.SPEC_SCAN_RANGE[1])
			
			#print "SPEC_IMG len: ",len(self.SPEC_IMG)
			#print "SPEC SCAN NUM LIST: ",self.SPEC_SCAN_NUM_LIST
			self.check_and_update_scan_slider()
		#return
	
	def update_scan_slider(self,widget):
		self.update_scan_slider_now()
		return
	
	def update_scan_slider_now(self):
		actual_scan_num = self.scan_slider_spinButton.get_value()
		actual_scan_num = int(actual_scan_num)
		#print "Actual scan number: ",actual_scan_num
		for i in range(len(self.SPEC_SCAN_NUM_LIST)):
			if actual_scan_num == self.SPEC_SCAN_NUM_LIST[i]:
				self.SPEC_ACTUAL_SCAN = self.SPEC_SCAN_LIST[i]
				break
		#print "Actual scan number: ", self.SPEC_ACTUAL_SCAN.nr
		
		#self.SPEC_ACTUAL_SCAN.ReadData()#All scan are Data Ready
		self.SPEC_ACTUAL_SCAN_IMG = []
		self.SPEC_ACTUAL_SCAN_IMG = self.SPEC_ACTUAL_SCAN.data[_SPEC_IMG_COL]
		self.SPEC_ACTUAL_SCAN_IMG = self.SPEC_ACTUAL_SCAN_IMG.astype('int')
		#print "Actual scan images: ",self.SPEC_ACTUAL_SCAN_IMG
		actual_img_num  = self.SELECTED_IMG_NUM
		#actual_img_num = self.scan_slider_imgSlider.get_value()
		#actual_img_num = int(actual_img_num)
		if (actual_img_num == None) or (actual_img_num not in self.SPEC_ACTUAL_SCAN_IMG):
			actual_img_num = self.SPEC_ACTUAL_SCAN_IMG[0]
			for k in self.store.keys():
				if actual_img_num in self.store_img[k]:
					self.edf_folder = k
					break
				
		#while gtk.events_pending():
			#gtk.main_iteration()
		if self.SPEC_ACTUAL_SCAN_IMG[0] == self.SPEC_ACTUAL_SCAN_IMG[-1]:
			self.scan_slider_imgSlider.set_range(self.SPEC_ACTUAL_SCAN_IMG[0], self.SPEC_ACTUAL_SCAN_IMG[-1]+1)
		else:
			self.scan_slider_imgSlider.set_range(self.SPEC_ACTUAL_SCAN_IMG[0], self.SPEC_ACTUAL_SCAN_IMG[-1])
		
		self.scan_slider_imgSlider.set_value(actual_img_num)#This will call the plot_scan too
		#print "Actual image number: ",actual_img_num
		#print self.SPEC_ACTUAL_SCAN_IMG
		
		self.SPEC_ACTUAL_SCAN_DATA = []
		try:
			self.SPEC_ACTUAL_SCAN_IMG_NAMES = select_files_from_list(self.store[self.edf_folder], self.SPEC_ACTUAL_SCAN_IMG[0], self.SPEC_ACTUAL_SCAN_IMG[-1])
			
			img_list = self.SPEC_ACTUAL_SCAN_IMG_NAMES
			#print "Actual images: ",img_list
			
			for j in range(len(img_list)):
				edf = join(self.edf_folder, img_list[j])
				data= fabio.open(edf).data
				self.SPEC_ACTUAL_SCAN_DATA.append(data)
			self.SPEC_ACTUAL_SCAN_DATA = N.asarray(self.SPEC_ACTUAL_SCAN_DATA)
			#print "Test if SCAN DATA is loaded: ",self.SPEC_ACTUAL_SCAN_DATA.shape
		except:
			pass
			#self.popup_info("warning","Attention: Data not found for this scan.")
		
		self.SPEC_SCAN_MOTOR_NAME = self.SPEC_ACTUAL_SCAN.colnames[0]
		self.SPEC_SCAN_MOTOR_DATA = self.SPEC_ACTUAL_SCAN.data[self.SPEC_SCAN_MOTOR_NAME]
		self.plot_scan()
		
		if not self.IMG_INIT:
			self.init_image()
		return
	
	def check_and_update_scan_slider(self):
		"""Get the actual scan object corresponding to the scan number and image number selected"""
		if self.DATA_IS_LOADED and self.SELECTED_IMG_NUM != None:
			for i in range(len(self.SPEC_IMG)):
				if self.SELECTED_IMG_NUM in self.SPEC_IMG[i]:
					self.SPEC_ACTUAL_SCAN = self.SPEC_SCAN_LIST[i]
					break
			
		else:
			self.SPEC_ACTUAL_SCAN = self.SPEC_SCAN_LIST[-1]#if no image is selected, the last scan will be displayed
			
		self.scan_slider_spinButton.set_value(self.SPEC_ACTUAL_SCAN.nr)#This will call the update scan slider too
		self.update_scan_slider_now()
	
	def slider_plot_scan(self, widget):
		self.plot_scan()
		
	def plot_scan(self):
		try:
			img_num = self.scan_slider_imgSlider.get_value()
			img_num = int(img_num)
			self.SELECTED_IMG_NUM = img_num
			#print "Image number: ",img_num
			img_index = N.where(self.SPEC_ACTUAL_SCAN_IMG == img_num)
			img_index = img_index[0][0]
			self.data = self.SPEC_ACTUAL_SCAN_DATA[img_index]
			if self.adj_btn.get_active():
				if self.data.shape == (960,560) or self.data.shape == (120,560):
					self.data = self.correct_geometry(self.data)
			if self.horizontal_detector:
				self.data = N.rot90(self.data)
			this_title = self.SPEC_ACTUAL_SCAN_IMG_NAMES[img_index]
			scan_motor = self.SPEC_SCAN_MOTOR_NAME
			this_motor_value = self.SPEC_SCAN_MOTOR_DATA[img_index]
			this_title = this_title +" - %s = %s"%(scan_motor, this_motor_value)
			self.MAIN_TITLE.set_text(this_title)
			self.scale_plot()
			self.canvas.draw()
			if isfile(join(self.edf_folder, self.SPEC_ACTUAL_SCAN_IMG_NAMES[img_index])):
				self.fabioIMG = fabio.open(join(self.edf_folder, self.SPEC_ACTUAL_SCAN_IMG_NAMES[img_index]))
			else:
				pass
		except:
			pass
		return
	
	def change_scale(self,button, data):
		if button.get_active():
			button.set_label("Log scale")
			data = N.log10(data+1e-6)
			
			actual_vmin = self.sld_vmin.get_value()
			actual_vmax = self.sld_vmax.get_value()
			#print data
			#print "Log scale called. Actual VMAX: ",actual_vmax
			self.vmax = N.log10(actual_vmax) if self.log_scale == 0 else actual_vmax
			if actual_vmin == 0:
				self.vmin=0
			elif actual_vmin >0:
				self.vmin = N.log10(actual_vmin) if self.log_scale == 0 else actual_vmin

			#vmax_range = N.max(data[N.isfinite(data)])
			self.vmax_range = data.max()
			#self.vmax_range = self.vmax_range if self.vmax_range < 8*med_log else 8*med_log
			self.log_scale = 1

		else:
			button.set_label("Linear scale")
			actual_vmin = self.sld_vmin.get_value()
			actual_vmax = self.sld_vmax.get_value()
			#print "Linear scale called, Actual vmax: ",actual_vmax
			if self.log_scale == 0:
				self.vmax = actual_vmax
			elif self.log_scale == 1:
				self.vmax = N.power(10.,actual_vmax)
			#print data
			self.vmax_range = data.max()
			#vmax_range = N.max(data[N.isfinite(data)])
			#if self.vmax_range > 100*med_lin:
			#       self.vmax_range = 100*med_lin

			if actual_vmin ==0:
				self.vmin = 0
			elif actual_vmin>0:
				if self.log_scale == 0:
					self.vmin = actual_vmin
				elif self.log_scale == 1:
					self.vmin = N.power(10.,actual_vmin)
			self.log_scale = 0

		self.sld_vmax.set_range(0,self.vmax_range)
		self.sld_vmin.set_range(0,self.vmax_range)
		return data
	
	
	def scale_plot(self):
		data = self.data.copy()
		data = self.change_scale(self.linear_scale_btn, data)
		self.img.set_data(data)
		#self.img.set_extent(self.MAIN_EXTENT)
		#self.canvas.draw()
    
	def log_update(self,widget):
		self.scale_plot()
		if self.log_scale==1:
			self.cb.set_label(r'$Log_{10}\ (Counts\ per\ second)\ [arb.\ units]$',fontsize=18)
			self.cb2.set_label(r'$Log_{10}\ (Counts\ per\ second)\ [arb.\ units]$',fontsize=18)
		else:
			self.cb.set_label(r'$Intensity\ (Counts\ per\ second)$', fontsize=18)
			self.cb2.set_label(r'$Intensity\ (Counts\ per\ second)$', fontsize=18)
		self.slider_update()

	def scale_update(self,widget):
		self.vmin = self.sld_vmin.get_value()
		self.vmax = self.sld_vmax.get_value()
		self.vmin_spin_btn.set_value(self.vmin)
		self.vmax_spin_btn.set_value(self.vmax)
		self.slider_update()

	def scale_update_spin(self,widget):
		self.vmin = self.vmin_spin_btn.get_value()
		self.vmax = self.vmax_spin_btn.get_value()
		self.slider_update()

	def slider_update(self):
		if self.notebook.get_current_page() == 0:
			self.img.set_clim(self.vmin, self.vmax)
			self.sld_vmax.set_value(self.vmax)
			self.sld_vmin.set_value(self.vmin)
			if self.linear_scale_btn.get_active():
				self.vmin_spin_btn.set_adjustment(gtk.Adjustment(self.vmin, 0, self.vmax_range, 0.1, 1.0, 0))
				self.vmax_spin_btn.set_adjustment(gtk.Adjustment(self.vmax, 0, self.vmax_range, 0.1, 1.0, 0))
			else:
				self.vmin_spin_btn.set_adjustment(gtk.Adjustment(self.vmin, 0, self.vmax_range, 10, 100, 0))
				self.vmax_spin_btn.set_adjustment(gtk.Adjustment(self.vmax, 0, self.vmax_range, 10, 100, 0))
			self.vmax_spin_btn.update()
			
			self.ax.relim()
			self.img.set_extent(self.MAIN_EXTENT)
			self.ax.set_xlim(self.MAIN_EXTENT[0], self.MAIN_EXTENT[1])
			self.ax.set_ylim(self.MAIN_EXTENT[2], self.MAIN_EXTENT[3])
			self.canvas.draw()
			#self.ax.figure.canvas.draw_idle()
		elif self.notebook.get_current_page() == 1:
			self.polar_img.set_clim(self.vmin, self.vmax)
			self.sld_vmax.set_value(self.vmax)
			self.sld_vmin.set_value(self.vmin)
			self.vmin_spin_btn.set_adjustment(gtk.Adjustment(self.vmin, 0, self.vmax_range, 0.5, 10.0, 0))
			self.vmax_spin_btn.set_adjustment(gtk.Adjustment(self.vmax, 0, self.vmax_range, 0.5, 10.0, 0))
			self.vmax_spin_btn.update()
			#self.polar_ax.relim()
			self.plot_PF_2()

			
			#self.pole_canvas.draw()
    
	def detector_disposition(self,widget):
		#data = self.data.copy()
		if self.detector_disposition_horizontal.get_active():
			self.horizontal_detector = True
			self.detector_disposition_horizontal.set_label("Horizontal detector")
			if self.data.shape[0]>self.data.shape[1]:#self.data.shape != (578,1148):
				self.data = N.rot90(self.data)
			self.img.set_array(self.data)
			#self.cb.ax.set_visible(False)
			#self.cb2.ax.set_visible(True)
		#self.img.set_extent((0,self.xDim,0,self.yDim))
		else:
			self.horizontal_detector = False
			self.detector_disposition_horizontal.set_label("Vertical detector")

			self.data = N.rot90(self.data,3)
			self.img.set_array(self.data)
			#self.cb.ax.set_visible(True)
			#self.cb2.ax.set_visible(False)

		imshape = self.data.shape
		self.xDim1 = imshape[1]
		self.yDim1 = imshape[0]
		if self.xDim1 > self.yDim1:
			self.cb.ax.set_visible(False)
			self.cb2.ax.set_visible(True)
		else:
			self.cb.ax.set_visible(True)
			self.cb2.ax.set_visible(False)
		self.img.set_extent((0,self.xDim1,0,self.yDim1))
		self.ax.set_xlim(0,self.xDim1)
		self.ax.set_ylim(0,self.yDim1)
		self.slider_update()

	def calcul_chi_2theta_d(self,event):
		""" EVENT Calculate chi, 2theta of the corresponding points on the image """
		x = event.xdata
		y = event.ydata
		self.check_azimuthal_integrator()
		chi = self.tableChi[y,x] -90.#+ self.chi
		tth = self.tableTwoTheta[y,x]
		d   = self.table_dSpace[y,x]
		return chi, tth, d


	def show_chi_delta(self,widget):
		""" If the show_chi_delta_btn is checked and the calibration file is loaded """
		if self.show_chi_delta_btn.get_active():
			if self.calibrated:
				self.show_chi_delta_flag = True
				self.show_chi_txt.set_visible(True)
				self.show_delta_txt.set_visible(True)
				self.show_d_txt.set_visible(True)
			else:
				self.show_chi_delta_btn.set_active(False)
				self.show_chi_delta_flag = False
				self.popup_info("warning","Please calibrate the detector before checking this box!")

		else:
			self.show_chi_delta_flag = False
			self.show_chi_txt.set_visible(False)
			self.show_delta_txt.set_visible(False)
			self.show_d_txt.set_visible(False)

	def status_update(self,event):
		if event.inaxes==self.ax:
			self.x_pos.set_text("X = %d"%event.xdata)
			self.y_pos.set_text("Y = %d"%event.ydata)
			self.z_pos.set_text("Z = %d"%self.data[event.ydata,event.xdata])
			if self.show_chi_delta_flag == True:
				chi,tth,d = self.calcul_chi_2theta_d(event)
				self.show_chi_txt.set_text("Chi = %.2f"%chi)
				self.show_delta_txt.set_text("2Theta = %.2f"%tth)
				self.show_d_txt.set_text("d = %.4f A"%d)

	def plot_update(self,widget):
		if self.tth_chi_space_btn.get_active():
			self.MAIN_XLABEL.set_text("2Theta (deg.)")
			self.MAIN_YLABEL.set_text("Chi (deg.)")
		else:
			self.MAIN_XLABEL.set_text("X (pixel)")
			self.MAIN_YLABEL.set_text("Y (pixel)")
		self.plot_data()

	def zoom_on(self,widget):
		"""For the Zoom button"""
		if self.zoomtb.get_active():
			self.zoom = True
			self.canvas.window.set_cursor(gtk.gdk.Cursor(gtk.gdk.CROSS))
			self.cursor.visible = False
		else:
			self.zoom = False
			self.canvas.window.set_cursor(None)
			self.cursor.visible = True

	def zoom_image(self):
		tmp_x = [self.x0, self.x1]
		tmp_y = [self.y0, self.y1]
		tmp_x = N.sort(tmp_x)
		tmp_y = N.sort(tmp_y)
		self.xDim0, self.xDim1 = tmp_x[0], tmp_x[1]
		self.yDim0, self.yDim1 = tmp_y[0], tmp_y[1]

		self.ax.set_xlim(self.xDim0, self.xDim1)
		if self.detector_type not in ["D5", "D1"]:
			self.ax.set_ylim(self.yDim1, self.yDim0)
		else:
			self.ax.set_ylim(self.yDim0, self.yDim1)
		self.cb.ax.set_visible(False)
		self.canvas.draw()

	def reset_scale(self,widget):
		if self.linear_scale_btn.get_active():
			self.vmin = 0
			self.vmax = N.log10(10*self.med)
		else:
			self.vmin = 0
			self.vmax = 10*self.med
		self.slider_update()

	def reset_image(self,widget):
		"""For the Home button"""
		self.xDim0 = 0
		self.xDim1 = self.data.shape[1]
		self.yDim0 = 0
		self.yDim1 = self.data.shape[0]
		if self.tth_chi_space_btn.get_active():
			self.xDim0 = self.tth_pyFAI.min()
			self.xDim1 = self.tth_pyFAI.max()
			self.yDim0 = self.chi_pyFAI.min()
			self.yDim1 = self.chi_pyFAI.max()
			
		self.MAIN_EXTENT = (self.xDim0,self.xDim1,self.yDim0, self.yDim1)
		self.ax.set_xlim(self.xDim0,self.xDim1)
		if self.detector_type=="S70":
			self.ax.set_ylim(self.yDim1,self.yDim0)
		else:
			self.ax.set_ylim(self.yDim0,self.yDim1)
		self.img.set_extent(self.MAIN_EXTENT)
		if self.data.shape[0] < self.data.shape[1]:
			self.cb.ax.set_visible(False)
			self.cb2.ax.set_visible(True)
		else:
			self.cb.ax.set_visible(True)
			self.cb2.ax.set_visible(False)

		self.canvas.draw()
	
	def get_list_dir(self,this_dir):
		""" Get the list of directories inside this_dir"""
		self.MODEL.clear()
		main_dir = this_dir
		list_dir = listdir(main_dir)
		no_of_threads = len(list_dir)
		for i in range(no_of_threads):
			#t = threading.Thread(target=self._thread_scanning,args=(main_dir,list_dir[i],))
			#self.threads.append(t)
			#t.start()
			self._thread_scanning(main_dir,list_dir[i])

		#gobject.timeout_add(200, self._callback)  # This will cause the main app to
		#check every 200 ms if the threads are done.
	def _callback(self):
		if threading.active_count() == 1:  # If only one left, scanning is done
			return False  # False make callback stop
		#print threading.active_count()
		return True
	
	def _thread_scanning(self,main_d,list_d):
		path = os.sep.join((main_d, list_d))  # Made use of os's sep instead...
		if os.path.isdir(path):
			list_subd = os.listdir(path)
			grand_parent = self.MODEL.append(None,[list_d])
			#for sub in list_subd:
				#parent=self.MODEL.append(grand_parent,[sub])
			self.scan_EDF_files(grand_parent,path)

		#time.sleep(3)  # Useless other than to delay finish of thread.
	
	def scan_EDF_files(self,parent, this_dir):
		#print "Scanning this directory: ",this_dir
		main_store= [i for i in listdir(this_dir) if isfile(join(this_dir,i)) and i.endswith(".edf") or i.endswith(".edf.gz")]
		main_store = sorted(main_store)
		if len(main_store)>0:
			self.store[str(this_dir)] = main_store
			for f in main_store:
				self.MODEL.append(parent,[f])
		else:
			self.store[str(this_dir)] = [None]
		#print "Scanning... ",self.store.keys()
		#print this_dir
		
	
	def choose_folder(self,w):
		dialog = gtk.FileChooserDialog(title="Select an EDF folder",action=gtk.FILE_CHOOSER_ACTION_SELECT_FOLDER, buttons = (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL, gtk.STOCK_OPEN, gtk.RESPONSE_OK))
		dialog.set_current_folder(self.current_folder)
		response=dialog.run()

		if response==gtk.RESPONSE_OK:
			folder=dialog.get_filename()
			folder_basename = folder.split("/")[-1]
			#print folder
			main_store= [i for i in listdir(folder) if isfile(join(folder,i)) and i.endswith(".edf") or i.endswith(".edf.gz")]
			self.store = {}
			if len(main_store)>0:
				main_store = sorted(main_store)
				self.store[str(folder)] = main_store

			self.current_folder = folder
			#print self.store
			self.get_list_dir(self.current_folder)
			if len(main_store)>0:
				#self.list_store.clear()
				#self.MODEL.clear()
				for i in main_store:
					self.MODEL.append(None,[i])
				self.TVcolumn.set_title(folder_basename)
				
			else:
				pass
			self.DATA_IS_LOADED = True
			
		else:
			pass
		dialog.destroy()
		if self.DATA_IS_LOADED:
			self.store_img = {}
			for k in self.store.keys():
				self.store_img[k] = get_img_list(self.store[k])
			#print "### ",self.store.keys()
			if self.SPEC_DATA:
				self.SPEC_DATA.Update()
				self.update_spec_data()

	def folder_update(self,widget):
		folder = self.current_folder
		if folder is not os.getcwd():
			main_store= [i for i in listdir(folder) if isfile(join(folder,i)) and i.endswith(".edf") or i.endswith(".edf.gz")]
			main_store = sorted(main_store)
			self.store={}
			#self.list_store.clear()
			self.store[self.current_folder] = main_store
			self.get_list_dir(self.current_folder)
			self.store_img = {}
			for k in self.store.keys():
				self.store_img[k] = get_img_list(self.store[k])
			
			for i in main_store:
				self.MODEL.append(None,[i])
			
			self.DATA_IS_LOADED = True
			if self.SPEC_DATA:
				self.SPEC_DATA.Update()
				self.update_spec_data()
		#return 1

	def save_image(self,widget):
		dialog = gtk.FileChooserDialog(title="Save image", action=gtk.FILE_CHOOSER_ACTION_SAVE, buttons = (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL, gtk.STOCK_SAVE, gtk.RESPONSE_OK))
		filename = self.edf_choosen.split(".")[0] if self.edf_choosen != "" else ""
		dialog.set_current_name(filename+".pdf")
		#dialog.set_filename(filename)
		dialog.set_current_folder(self.current_folder)
		filtre = gtk.FileFilter()
		filtre.set_name("images")
		filtre.add_pattern("*.png")
		filtre.add_pattern("*.jpg")
		filtre.add_pattern("*.pdf")
		filtre.add_pattern("*.ps")
		filtre.add_pattern("*.eps")
		dialog.add_filter(filtre)
		filtre = gtk.FileFilter()
		filtre.set_name("Other")
		filtre.add_pattern("*")
		dialog.add_filter(filtre)
		response = dialog.run()

		if response==gtk.RESPONSE_OK:
			self.fig.savefig(dialog.get_filename())
		dialog.destroy()

	def save_adjust(self,widget):
		""" Save the current EDF image, the file name will be name+adjusted """
		self.fabioIMG.data = self.data
		name = self.edf.split(".")[0]+"_adjusted"
		ext  = self.edf.split(".")[1]
		filename = name+"."+ext
		self.fabioIMG.write(filename)
		self.popup_info("info", "Image %s is successfully saved !"%filename)

	def change_space(self,widget,space):
		if space.upper()=="TTH":
			self.q_space_btn.set_active(False)
			self.tth_chi_space_btn.set_active(True)
		elif space.upper()=="Q":
			self.q_space_btn.set_active(True)
			self.tth_chi_space_btn.set_active(False)
		self.plot_data()

	def plot_data(self):
		"""plot the selected edf image"""
		self.data = self.fabioIMG.data
		if self.adj_btn.get_active():
			adjust = True
		else:
			adjust = False

		if self.cln_btn.get_active():
			clean = True
		else:
			clean = False

		#self.detector = libX.Detector()


		### Calculate the median and the deviation of this data ###
		self.med, self.nMAD = median_stats(self.data)
		#If the loaded EDF image is already adjusted, shape = (578,1148)
		if self.detector_type=="D5":
			self.detector = libX.Detector()
			if self.data.shape != (960,560):
				self.adj_btn.set_sensitive(False)
				adjust = False
			else:
				self.adj_btn.set_sensitive(True)
				adjust = self.adj_btn.get_active()

		elif self.detector_type == "D1":
			self.adj_btn.set_sensitive(False)
			adjust = False

		elif self.detector_type=="S70":
			self.detector = libX.Detector(nModules=1)

			if self.data.shape != (120,560):
				self.adj_btn.set_sensitive(False)
				adjust = False
			else:
				self.adj_btn.set_sensitive(True)
				adjust = self.adj_btn.get_active()

		if clean:
			nBad = 0
			bad = self.data > (self.med + 100*self.nMAD)
			nBad+=len(self.data[bad])
			self.data[bad]=0
		
		if adjust:

			self.detector.set_data(array=self.data)
			self.detector.adjust_data()
			self.detector.set_physical_data()
			self.detector.reshape_pixels()
			self.data = self.detector.physical.data

		imshape = self.data.shape
		self.xDim0=0
		self.yDim0=0
		if self.horizontal_detector == True and imshape != (578,1148):
			self.xDim1 = imshape[0]
			self.yDim1 = imshape[1]
			self.data = N.rot90(self.data) #rotation in the clock-wise direction - right rotation
		else:
			self.xDim1 = imshape[1]
			self.yDim1 = imshape[0]
		self.MAIN_EXTENT = (self.xDim0, self.xDim1, self.yDim0, self.yDim1)
		
		if self.tth_chi_space_btn.get_active():
			if self.calibrated == False:
				self.popup_info('warning','Please calibrate the detector before checking this!')
				self.tth_chi_space_btn.set_active(False)
			elif self.calibrated==True:
				self.show_chi_delta_btn.set_sensitive(False)
				self.show_chi_delta_flag=False
				self.check_azimuthal_integrator()
				if self.data.shape == (578,1148) and self.detector_type=="D5":
					self.data,self.tth_pyFAI,self.chi_pyFAI = self.azimuthalIntegration.integrate2d(self.data,578,1148,unit="2th_deg")
					self.chi_pyFAI = self.chi_pyFAI - 90#180 + self.chi
				elif self.data.shape==(1148,578) and self.detector_type=="D5":
					self.data = N.rot90(self.data)
					self.data,self.tth_pyFAI,self.chi_pyFAI = self.azimuthalIntegration.integrate2d(self.data,578,1148,unit="2th_deg")
					self.chi_pyFAI = self.chi_pyFAI - 90# 180 + self.chi
				elif self.data.shape==(120,578) and self.detector_type=="S70":
					#self.azimuthalIntegration.setChiDiscAtZero()
					self.data,self.tth_pyFAI,self.chi_pyFAI = self.azimuthalIntegration.integrate2d(self.data,120,578,unit="2th_deg")
					self.chi_pyFAI = self.chi_pyFAI + 90
				elif self.data.shape == (577,913) and self.detector_type=="D1":
					self.data,self.tth_pyFAI,self.chi_pyFAI = self.azimuthalIntegration.integrate2d(self.data,577,913,unit="2th_deg")
					self.chi_pyFAI = self.chi_pyFAI - 90.#A corriger avec chi gonio
				else:
					self.popup_info('warning','Please adjust the image to proceed this operation!')
					self.tth_chi_space_btn.set_active(False)
				
				self.MAIN_EXTENT = (self.tth_pyFAI.min(), self.tth_pyFAI.max(), self.chi_pyFAI.min(), self.chi_pyFAI.max())
				#print self.MAIN_EXTENT
		else:
			self.show_chi_delta_btn.set_sensitive(True)
			self.show_chi_delta_flag = self.show_chi_delta_btn.get_active()

		
		self.scale_plot()

		if self.data.shape[0] < self.data.shape[1]:
			self.cb.ax.set_visible(False)
			self.cb2.ax.set_visible(True)
		else:
			self.cb.ax.set_visible(True)
			self.cb2.ax.set_visible(False)
		#self.img.set_extent(self.MAIN_EXTENT)
		#self.ax.set_xlim(0,self.xDim1)
		#if self.detector_type in ["D5", "D1"]:
			#self.ax.set_ylim(0,self.yDim1)
		#elif self.detector_type=="S70":
			#self.ax.set_ylim(self.yDim1,0)
		self.slider_update()

	def plot_profiles(self, x, y, cross_line=True):
		"""Line x = 2theta profile, Column y = Chi profile, if not transform in 2theta,chi space, the X,Y profiles are plotted in detector coordinates"""
		if cross_line:
			x=x[0]
			y=y[0]
			#Draw the cross-lines in the main figure:
			hline = self.ax.axhline(y, color='k', ls='--', lw=1)
			self.lines.append(hline)
			vline = self.ax.axvline(x, color='k', ls='--', lw=1)
			self.lines.append(vline)
			
			if self.tth_chi_space_btn.get_active():
				x = get_index(self.tth_pyFAI,x)
				y = get_index(self.chi_pyFAI,y)
			x= int(x)
			y= int(y)        

			self.profiles_data_X = self.data[y-5:y+5,:].sum(axis=0)
			self.profiles_data_X = self.profiles_data_X / 10.

			self.profiles_data_Y = self.data[:, x-5:x+5].sum(axis=-1)
			self.profiles_data_Y = self.profiles_data_Y / 10.

			if self.tth_chi_space_btn.get_active():
				X_label = "2 Theta (deg)"
				Y_label = "Chi (deg)"
				yc = self.chi_pyFAI[y]
				xc = self.tth_pyFAI[x]
				coor_X = self.tth_pyFAI
				coor_Y = self.chi_pyFAI
				self.chi_title.set_text("Chi")
				self.tth_title.set_text("2 Theta")
			else:
				X_label = "X (pixels)"
				Y_label = "Y (pixels)"
				xc = x
				yc = y
				coor_X = N.arange(self.profiles_data_X.shape[0])
				coor_Y = N.arange(self.profiles_data_Y.shape[0])
				self.chi_title.set_text("Y")
				self.tth_title.set_text("X")
		else:
			#print "Data shape: ",self.data.shape
			num = int(N.hypot(x[1]-x[0], y[1]-y[0])) #Number of points to be taken along the line
			print "Number of points selected: ",num
			xi, yi = N.linspace(x[0], x[1], num), N.linspace(y[0], y[1], num)
			#self.profiles_data_X = self.profiles_data_Y = self.data[yi.astype(N.int), xi.astype(N.int)]
			self.profiles_data_X = self.profiles_data_Y = ndimage.map_coordinates(self.data, N.vstack((yi,xi)))
			if self.tth_chi_space_btn.get_active():
				X_label = "2 Theta (deg)"
				Y_label = "Chi (deg)"
				coor_X = N.linspace(self.tth_pyFAI[int(x[0])], self.tth_pyFAI[int(x[1])], num)
				coor_Y = N.linspace(self.chi_pyFAI[int(y[0])], self.chi_pyFAI[int(y[1])], num)

				self.chi_title.set_text("Chi")
				self.tth_title.set_text("2 Theta")
			else:
				X_label = "X (pixels)"
				Y_label = "Y (pixels)"
				coor_X = xi
				coor_Y = yi
				self.chi_title.set_text("Y")
				self.tth_title.set_text("X")
			tmpX = N.sort(coor_X)
			tmpY = N.sort(coor_Y)
			xc = tmpX[self.profiles_data_X.argmax()]
			yc = tmpY[self.profiles_data_Y.argmax()]
		self.coor_X_export = coor_X
		self.coor_Y_export = coor_Y
		self.profiles_ax1.cla()
		self.profiles_ax2.cla()
		self.profiles_ax1.format_coord = self.pro_format_coord
		self.profiles_ax2.format_coord = self.pro_format_coord
		# The CHI or Vertical (Y) profile (ax1):

		self.profiles_ax1.plot(coor_Y, self.profiles_data_Y, color='blue', lw=1.5)
		Y_fitted_params, Y_fitted_data = fit(coor_Y, self.profiles_data_Y, yc, arbitrary= not cross_line)
		self.profiles_ax1.plot(coor_Y, Y_fitted_data, color='red', lw=1, alpha=0.8)
		self.profiles_ax1.set_title(Y_label, size=12)

		# The TTH or Horizontal (X) profile (ax2):
		self.profiles_ax2.plot(coor_X, self.profiles_data_X, color='blue', lw=1.5)
		X_fitted_params, X_fitted_data = fit(coor_X, self.profiles_data_X, xc, arbitrary= not cross_line)
		self.profiles_ax2.plot(coor_X, X_fitted_data, color='red', lw=1, alpha=0.8)
		self.profiles_ax2.set_xlabel(X_label, size=12)
		self.profiles_canvas.draw()
		# Show the fitted results
		self.chi_fitted_y0.set_text("%.2f"%Y_fitted_params['y0'].value)
		self.chi_fitted_xc.set_text("%.2f"%Y_fitted_params['xc'].value)
		self.chi_fitted_A.set_text("%.2f"%Y_fitted_params['A'].value)
		self.chi_fitted_w.set_text("%.2f"%Y_fitted_params['w'].value)
		self.chi_fitted_mu.set_text("%.2f"%Y_fitted_params['mu'].value)

		self.tth_fitted_y0.set_text("%.2f"%X_fitted_params['y0'].value)
		self.tth_fitted_xc.set_text("%.2f"%X_fitted_params['xc'].value)
		self.tth_fitted_A.set_text("%.2f"%X_fitted_params['A'].value)
		self.tth_fitted_w.set_text("%.2f"%X_fitted_params['w'].value)
		self.tth_fitted_mu.set_text("%.2f"%X_fitted_params['mu'].value)

		self.profiles_refresh()
		self.canvas.draw()

	def draw_pointed(self, x, y, finished=False):
		p=self.ax.plot(x,y,'ro')
		self.points.append(p[0])
		if finished:
			l=self.ax.plot(self.arb_lines_X, self.arb_lines_Y, '--',linewidth=1.5, color='white')
			self.lines.append(l[0])
		else:
			pass

		if self.detector_type=="S70":
			self.ax.axis([self.xDim0,self.xDim1, self.yDim1,self.yDim0])
		else:
			self.ax.axis([self.xDim0, self.xDim1, self.yDim0, self.yDim1])

		self.canvas.draw()


	def profiles_refresh(self):
		""" """
		if self.profiles_log_btn.get_active():
			self.profiles_ax1.set_yscale('log')
			self.profiles_ax2.set_yscale('log')

		else:
			self.profiles_ax1.set_yscale('linear')
			self.profiles_ax2.set_yscale('linear')

		self.profiles_canvas.draw()
		#return

	def profiles_update(self, widget):
		self.profiles_refresh()

	def profiles_export(self,widget):
		""" Export X,Y profiles data in the same folder as the EDF image """
		proX_fname = self.edf.split(".")[0]+"_X_profile.dat"
		proY_fname = self.edf.split(".")[0]+"_Y_profile.dat"
		proX_export= N.vstack([self.coor_X_export, self.profiles_data_X])
		proX_export=proX_export.T
		proY_export= N.vstack([self.coor_Y_export, self.profiles_data_Y])
		proY_export=proY_export.T
		try:
			N.savetxt(proX_fname, proX_export)
			N.savetxt(proY_fname, proY_export)
			self.popup_info('info','Data are successfully exported in %s and %s!'%(proX_fname, proY_fname))

		except:
			self.popup_info('error','ERROR! Data not exported!')

	def draw_rect(self):
		self.rect.set_width(self.x1 - self.x0)
		self.rect.set_height(self.y1 - self.y0)
		self.rect.set_xy((self.x0, self.y0))
		self.rect.set_linestyle('solid')
		self.rect.set_facecolor("white")
		self.rect.set_alpha(0.3)
		self.rect.set_edgecolor("black")
		#self.rect.set_visible(self.zoom_press)
		self.rect.set_visible(self.zoom_press)
		self.canvas.draw()

	def peak_max_old(self, event):
		x = int(event.xdata)
		y = int(event.ydata)
		img = self.data.copy()
		#### range problem
		if y-25<0:
			down = 0
			j0 = y
		else:
			down = y-25
			j0 = 25
		if y+25>self.yDim1:
			up = self.yDim1-1
		else:
			up = y+25
      
		if x-25<0:
			left = 0
			i0 = x
		else:
			left = x-25
			i0 = 25
		if x+25>self.xDim1:
			right = self.xDim1-1
		else:
			right = x+25

		#img = N.asarray(img[y-25:y+25,x-25:x+25])
		img = N.asarray(img[down:up,left:right])
		img = ndimage.median_filter(img,5)
		j,i = N.unravel_index(img.argmax(), img.shape)
		x = x + i - i0
		y = y + j - j0
		chi = self.tableChi[y,x] - 90 + self.chi
		tth = self.tableTwoTheta[y,x]
		d   = self.table_dSpace[y,x]
		chi_text = r'$ \chi \ = \ %.2f$'%chi
		tth_text = r'$ 2\theta \ = \ %.2f$'%tth
		d_text = r'$ d \ = \ %.4f \ \AA$'%d
		txt1 = self.ax.text(x,y+50,tth_text, fontsize=14, color="white")
		txt2 = self.ax.text(x,y,chi_text, fontsize=14, color="white")
		txt3 = self.ax.text(x,y-50,d_text, fontsize=14, color="white")
		self.my_notes.append(txt1)
		self.my_notes.append(txt2)
		self.my_notes.append(txt3)
		self.canvas.draw()

	def peak_max(self, event):
		x = int(event.xdata)
		y = int(event.ydata)
		chi,tth,d=self.calcul_chi_2theta_d(event)
		chi_text = r'$ \chi \ = \ %.2f$'%chi
		tth_text = r'$ 2\theta \ = \ %.2f$'%tth
		d_text = r'$ d \ = \ %.4f \ \AA$'%d
		txt1 = self.ax.text(x,y+50,tth_text, fontsize=14, color="black")
		txt2 = self.ax.text(x,y,chi_text, fontsize=14, color="black")
		txt3 = self.ax.text(x,y-50,d_text, fontsize=14, color="black")
		self.my_notes.append(txt1)
		self.my_notes.append(txt2)
		self.my_notes.append(txt3)
		self.canvas.draw()

	def on_press(self, event):
		#********** Zoom action ***************************************
		if (self.zoom) and (event.inaxes == self.ax) and (event.button==1):
			#print 'press'
			self.zoom_press = True
			self.x0 = event.xdata
			self.y0 = event.ydata
		#******** Plot cross profiles *********************************
		elif (event.inaxes == self.ax) and (event.button==3) and self.plotXYprofiles_btn.get_active():
			x = event.xdata
			y = event.ydata
			xx=[]
			yy=[]
			xx.append(x)
			yy.append(y)
			try:
				self.clear_notes()
			except:
				pass
			self.plot_profiles(xx,yy,cross_line=True)
		#******** Plot arbitrary profiles  ****************************
		elif (event.inaxes == self.ax) and (event.button==1) and self.arbitrary_profiles_btn.get_active():
			self.arb_line_points +=1
			#print "Number of points clicked: ",self.arb_line_points
			if self.arb_line_points>2:
				self.clear_notes()
				self.arb_line_points=1

			x = event.xdata
			y = event.ydata
			self.arb_lines_X.append(x)
			self.arb_lines_Y.append(y)
			if len(self.arb_lines_X)<2:
				finished=False
			elif len(self.arb_lines_X)==2:
				finished = True

			self.draw_pointed(x,y,finished)#If finished clicking, connect the two points by a line
			if finished:
				self.plot_profiles(self.arb_lines_X, self.arb_lines_Y, cross_line=False)
				self.arb_lines_X=[]
				self.arb_lines_Y=[]
		#******** Mark the peak information ****************************
		elif (event.inaxes == self.ax) and (event.button==3) and (self.show_chi_delta_flag == True):
			self.peak_max(event)

		elif event.button==2:
			self.clear_notes()
		else:
			return

	def clear_notes(self):
		if len(self.my_notes)>0:
			for txt in self.my_notes:
				txt.remove()
		if len(self.lines)>0:
			for line in self.lines:
				line.remove()
		if len(self.points)>0:
			for p in self.points:
				p.remove()

		self.canvas.draw()
		self.my_notes = []
		self.lines=[]
		self.points=[]
		self.arb_lines_X=[]
		self.arb_lines_Y=[]
		self.arb_line_points = 0

	def on_motion(self,event):
		if event.inaxes == self.ax:
			self.status_update(event)
			if self.zoom_press:
				self.mouse_moved = True
				self.x1 = event.xdata
				self.y1 = event.ydata
				self.draw_rect()
		else:
			return

	def on_release(self, event):
		if (self.zoom) and (event.inaxes == self.ax):
			#print 'release'

			self.zoom_press = False
			self.x1 = event.xdata
			self.y1 = event.ydata
			self.draw_rect()
			if self.mouse_moved==True:
				self.zoom_image()
				self.mouse_moved = False


	def check_input_data(self):
		""" Check if the user has entered correctly the data needed to construct the pole figure
		We need the start and end number of images, and the 2 theta for which we want to plot """
		try:
			self.img_deb = int(self.images_from_nb.get_text())
			self.img_fin = int(self.images_to_nb.get_text())
			self.pole_2theta = float(self.pole_2theta_field.get_text())
			#self.plot_chi_min = float(self.PF_chi_min.get_text())
			return True
		except:
			self.popup_info("error","Values are not entered correctly. Please check again!")
			return False

	def correct_geometry(self, data):
		"""Correct image geometry of XPAD D5 detector
		The default image shape is 960x560 """
		if self.detector_type == "D5":
			detector = libX.Detector()
		elif self.detector_type == "S70":
			detector = libX.Detector(nModules=1)
		else:
			self.popup_info("warning","Please specify the detector model.")
			return None
		#Ajouter les gaps ## Code de Clement Buton @Soleil
		detector.set_data(array=data)
		detector.adjust_data()
		detector.set_physical_data()
		detector.reshape_pixels()
		data = detector.physical.data
		#data = N.rot90(data)
		return data
		   
	def load_data(self, img_deb, img_fin, pole_2theta, plot_chi_min, logscale=0):
		"""
		folder: directory where edf images are found
		img_deb, img_fin: number of beginning and ending images for the maps
		pole_2theta: 2theta (degrees) to construct the pole figure
		logscale: (1 or 0) to present the pole figure in log scale or linear scale
		"""
		self.data_loading.show()
		phi_table = []
		intensity = []
		img_list  = select_files_from_list(self.store, img_deb, img_fin)
		if self.detector_type == "D1":
			first_dim = 577
			second_dim= 913
		elif self.detector_type == "D5":
			first_dim = 578
			second_dim= 1148
		elif self.detector_type == "S70":
			first_dim = 120
			second_dim= 578
		#if not self.horizontal_detector:
			#temp=first_dim
			#first_dim = second_dim
			#second_dim = temp
			
		total = len(img_list)
		processed = 0.
		for j in range(len(img_list)):
			gc.collect()
			edf_basename = img_list[j]
			edf = join(self.current_folder, edf_basename)
			try:
				img = fabio.open(edf)
				this_motor = get_motors(img.header)
				rot1 = N.radians(this_motor['nu'])*(-1.)
				rot2 = N.radians(this_motor['del'])*(-1.)
				rot3 = N.radians(90-this_motor['chi'])
				self.azimuthalIntegration.rot1 = rot1
				self.azimuthalIntegration.rot2 = rot2
				self.azimuthalIntegration.rot3 = rot3
				if img.data.shape == (960,560):
					data = self.correct_geometry(img.data)
					data = N.rot90(data)
					img.data = data 
					
				I,tth,chi = self.azimuthalIntegration.integrate2d(img.data,first_dim,second_dim,unit="2th_deg")
				#2theta vs pixel: y = ax + b
				b = tth.min()
				a = (tth.max() - b)/first_dim
				x = (pole_2theta - b)/a
				x = int(x)
				#print "Center pixel used: ",x
				i = I[:,x-5:x+5].sum(axis=-1)
				i = i / 10.
				intensity.append(i)
				
				this_kphi = this_motor['kphi']
				this_phi  = this_motor['phi']
				chi_gonio = this_motor['chi']
				self.sample_tilt = chi_gonio
				self.omega = this_motor['eta']
				nu_gonio  = this_motor['nu']
				if self.select_phi.get_active():			
					kphi = this_phi
				elif self.select_kphi.get_active():
					kphi = this_kphi
				#print "kPhi: ",kphi
				phi_table.append(kphi)
				#phi rows and chi colunms
				#print "Image %s is successully loaded"%edf_basename
				
			except:
				#print "Image %s does not exist!"%edf_basename
				continue
			processed +=1.
			fr = processed/total
			self.data_loading.set_fraction(fr)
			self.data_loading.set_text("Data loading: "+str(int(fr*100))+"%")
			self.data_loading.set_show_text(True)
			while gtk.events_pending():
				gtk.main_iteration()

		intensity = N.asarray(intensity)
		if logscale==1:
			intensity = intensity+1e-6
			intensity = N.log10(intensity)
		if self.detector_type=="S70":
			chi = chi + 90
		else:
			chi = chi-90
		phi_table = N.asarray(phi_table)
		
		#Prendre la partie negative
		if chi_gonio > 70:
			b = chi.min()
			a = (chi.max() - b)/(chi.shape[0])
			x = (plot_chi_min - b)/a
			x = int(x)
			chi = chi[:x]
			intensity = intensity[:,:x]
		
		return chi, phi_table, intensity

	def coordinates_transform(self, psi, phi, tth):
		#Transformer les coordonnees du detecteur en coordonnees polaires, avec phi et psi sont des matrices 2D, omega=angle incidente, tth=2theta bragg
		#Ref: Marie-Ingrid Richard, J. Appl. Cryst. 46,(2013),1842
		#tilt = self.sample_tilt
		#tilt = 90-tilt
		tilt = 0
		omega= self.omega
		phi = N.radians(phi)
		psi = N.radians(psi)
		
		cosPhi = N.cos(phi)
		sinPhi = N.sin(phi)
		cosPsi = N.cos(psi)
		sinPsi = N.sin(psi)
		cosTil = N.cos(N.radians(tilt))
		sinTil = N.sin(N.radians(tilt))
		cosTta = N.cos(N.radians(tth/2.))
		sinTta = N.sin(N.radians(tth/2.))
		cosOme = N.cos(N.radians(omega))
		sinOme = N.sin(N.radians(omega))
		
		#Phi: rotate around z
		#Omega: Rotate around y
		#Chi: rotate around -x
		Rx = [[1,0,0], [0,cosTil,sinTil], [0,-sinTil,cosTil]]
		Ry = [[cosOme,0,sinOme], [0,1,0], [-sinOme,0,cosOme]]
		Rz = [[cosPhi,-sinPhi,0], [sinPhi, cosPhi,0], [0,0,1]]
		qLab = [-sinTta, -cosTta*sinPsi, cosTta*cosPsi]
		
		Rx = N.asarray(Rx)
		Ry = N.asarray(Ry)
		Rz = N.asarray(Rz)
		#qLab=N.asarray(qLab)
		
		Rot = N.dot(N.dot(Rz,Rx),Ry) #Be careful of the order
		#Rot = N.dot(N.dot(Rz,Ry),Rx)
		qSample = N.dot(Rot, qLab)
		psi_pole = N.arccos(abs(qSample[2]))
		cosF    = qSample[0] / N.sqrt(qSample[0]**2 + qSample[1]**2)
		sinF    = qSample[1] / N.sqrt(qSample[0]**2 + qSample[1]**2)
		phi_pole = N.arccos(cosF)
		neg = N.sign(sinF)<0
		phi_pole[neg] = 2.0*N.pi-phi_pole[neg]
		
		return N.degrees(psi_pole), N.degrees(phi)
	
	def plot_pole_figure(self,widget):
		""" Click the PLOT button to plot the pole figure """
		check = self.check_input_data()
		if check:
			negative_chi = 0
			if self.detector_type is not "D1":
				self.PF_psi, self.PF_phi, self.PF_intensity = self.load_data(self.img_deb, self.img_fin, self.pole_2theta, 0, logscale=self.log_scale)
				self.plot_PF_2()
		else:
			self.popup_info("error","Please make sure that you have correctly entered all of the fields!")
		self.data_loading.hide()

	def replot_PF(self,widget):
		""" Replot the pole figure with logscale option"""
		self.PF_intensity = self.change_scale(widget,self.PF_intensity)
		
	def plot_PF(self):
		#********** Plot the pole figure ***************
		gc.collect()
		self.polar_ax.clear()
		self.polar_cax.clear()
		PF_data = self.PF_intensity.copy()
		if self.log_scale == 1:
			PF_data = flat_data(PF_data, self.vmin, self.vmax)
			PF_data = N.log10(PF_data+1e-6)
		else:
			dynhigh = N.where(PF_data > self.vmax)
			dynlow  = N.where(PF_data < self.vmin)
			PF_data[dynhigh] = self.vmax
			PF_data[dynlow]  = self.vmin
		psi,phi   = N.meshgrid(self.PF_psi,self.PF_phi)
		psi_neg   = N.where(psi<0.)
		psi[psi_neg] = 0. #Correction of the overplotted area
		psi     = 90.-psi
		m = Basemap(boundinglat=psi.min()-1, lon_0=270., resolution=None, projection = 'npstere', ax=self.polar_ax)#resolution= 'l'ow, 'h'igh, 'f'ull, None
		X,Y = m(phi, psi)
		#clevels = N.linspace(self.vmin, self.vmax,100)
		self.polar_img = self.polar_ax.contourf(X,Y,PF_data, 100)
		#self.polar_img = self.polar_ax.pcolormesh(X,Y,PF_data, vmin=self.vmin, vmax=self.vmax)

		nbcircles = N.arange(int(psi.min()), int(psi.max()), int(psi.max()-psi.min())/7)
		m.drawparallels(nbcircles, labels=[0,0,0,0], color='black', labelstyle='+/-', dashes = [2,2])
		m.drawmeridians(N.arange(0,360,90), labels=[1,1,1,1], color='black', labelstyle='+/-', dashes=[2,2])
		def PF_format_coord(x,y):
			pp,cc = m(x,y,inverse=True)
			cc=90-cc
			return 'Phi = %.2f, Psi = %.2f'%(pp,cc)
		if self.log_scale == 0:
			clabel = r'$Intensity\ (Counts\ per\ second)$'
		else:
			clabel = r'$Log_{10}\ (Counts\ per\ second)\ [arb.\ units]$'
		self.polar_cb  = self.polefig.colorbar(self.polar_img,cax=self.polar_cax, format="%d")
		self.polar_cb.set_label(clabel, fontsize=18)
		self.polar_cb.locator = MaxNLocator(nbins=6)
		for i in range(1,len(nbcircles)):
			p = nbcircles[i]
			xx,yy = m(45,p)
			self.polar_ax.text(xx,yy, str(90-p), color='white')
		title   = "2 Theta = %.2f Deg."%self.pole_2theta
		self.polar_ax.text(0.5, 1.08, title, horizontalalignment='center', transform = self.polar_ax.transAxes, fontsize=20)
		self.polar_ax.format_coord = PF_format_coord
		self.pole_canvas.draw()
		
	def plot_PF_2(self):
		#********** Plot the pole figure ***************
		gc.collect()
		self.polar_ax.clear()
		self.polar_cax.clear()
		PF_data = self.PF_intensity.copy()
		if self.log_scale == 1:
			PF_data = flat_data(PF_data, self.vmin, self.vmax)
			PF_data = N.log10(PF_data+1e-6)
			clabel = r'$Log_{10}\ (Counts\ per\ second)\ [arb.\ units]$'
			fmt = "%.1f"
		else:
			clabel = r'$Intensity\ (Counts\ per\ second)$'
			fmt = "%d"
			dynhigh = N.where(PF_data > self.vmax)
			dynlow  = N.where(PF_data < self.vmin)
			PF_data[dynhigh] = self.vmax
			PF_data[dynlow]  = self.vmin
		psi,phi   = N.meshgrid(self.PF_psi,self.PF_phi)
		
		psi,phi = self.coordinates_transform(psi, phi, self.pole_2theta)
		
		phi = N.radians(phi)
		clevels = N.linspace(self.vmin, self.vmax,30)
		self.polar_img = self.polar_ax.contourf(phi, psi, PF_data, clevels, vmin=self.vmin, vmax=self.vmax)
		self.polar_ax.set_rmin(psi.min())
		#self.polar_ax.set_rmin(0)
		#self.polar_ax.set_ylim(psi.min(), psi.max())
		self.polar_cb  = self.polefig.colorbar(self.polar_img,cax=self.polar_cax, format=fmt)
		self.polar_cb.set_label(clabel, fontsize=18)
		self.polar_cb.locator = MaxNLocator(nbins=6)
		title   = "2 Theta = %.2f Deg."%self.pole_2theta
		self.polar_ax.text(0.5, 1.08, title, horizontalalignment='center', transform = self.polar_ax.transAxes, fontsize=20)
		self.polar_ax.relim()
		self.pole_canvas.draw()
	#***************************************************************************
	#     For batch images correction 
	#***************************************************************************
	def select_source_folder(self, widget):
		dialog = gtk.FileChooserDialog(title="Select a sources folder",action=gtk.FILE_CHOOSER_ACTION_SELECT_FOLDER, buttons = (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL, gtk.STOCK_OPEN, gtk.RESPONSE_OK))
		dialog.set_current_folder(self.current_folder)
		response = dialog.run()
		if response==gtk.RESPONSE_OK:
			folder=dialog.get_filename()
			self.t1_src_path.set_text(folder)
			self.current_folder = folder
			self.src_folder = folder
		else:
			pass
		dialog.destroy()
	
	def select_destination_folder(self, widget):
		dialog = gtk.FileChooserDialog(title="Select the destination folder",action=gtk.FILE_CHOOSER_ACTION_SELECT_FOLDER, buttons = (gtk.STOCK_CANCEL, gtk.RESPONSE_CANCEL, gtk.STOCK_OPEN, gtk.RESPONSE_OK))
		dialog.set_current_folder(self.current_folder)
		response = dialog.run()
		if response==gtk.RESPONSE_OK:
			folder=dialog.get_filename()
			self.t1_des_path.set_text(folder)
			self.current_folder = folder
			self.des_folder = folder
		else:
			pass
		dialog.destroy()
	
	def show_proc(self,out_info):
		self.show_process_info.set_label(out_info)
	
	def batch_change_detector(self, widget):
		detector = self.t2_det_combobox.get_active_text()
		if detector == "D1":
			self.d1_center.show()
			self.d1_center_txt.show()
			self.d1_distance.show()
			self.d1_distance_txt.show()
			self.d1_specfile_browse.show()
			self.d1_specfile_path.show()
			self.d1_specfile_txt.show()
		else:
			self.d1_center.hide()
			self.d1_center_txt.hide()
			self.d1_distance.hide()
			self.d1_distance_txt.hide()
			self.d1_specfile_browse.hide()
			self.d1_specfile_path.hide()
			self.d1_specfile_txt.hide()
		
	def select_normalisation_file(self,widget,path):
		dialog = gtk.FileChooserDialog("Select spec file",None,gtk.FILE_CHOOSER_ACTION_OPEN,(gtk.STOCK_CANCEL,gtk.RESPONSE_CANCEL, gtk.STOCK_OPEN, gtk.RESPONSE_OK))
		dialog.set_current_folder(self.current_folder)
		response = dialog.run()
		if response == gtk.RESPONSE_OK:
			file_choosen = dialog.get_filename()
			path.set_text(file_choosen)
			self.current_folder = os.path.dirname(file_choosen)
			#if label == "A":
				#self.attenuation_file = file_choosen.decode('utf8')
			#elif label == "S":
			self.spec_file = file_choosen.decode('utf8')
			#elif label == "M":
				#self.mca_file = file_choosen.decode('utf8')
		else:
			pass
		dialog.destroy()
		
	def geometric_correction_D1(self,edf,adjusted_folder):
		distance = self.d1_distance.get_text()
		center   = self.d1_center.get_text()
		normalize= self.spec_file
		program  = 'xpad3_geometry.py'
		cmd = program+" --geometry=D1 --mask=search --distance="+distance+" --center="+center+" --normalize="+normalize+" --monitor=d0_cps --multiply="+str(self.monitor_ref)+" "+edf
		print cmd
		os.system(cmd)
		menage = "mv *g.edf %s"%self.des_folder
		os.system(menage)
		
		
	def transform(self,edf,adjusted_folder,normalisation):
		""" transformer edf en edf adjusted (en rajoutant les gaps)
		Corriger les lignes mortes, nettoyer les pixels chauds
		Normaliser avec le Imachine = 200 mA: Icorr = Imesure*monitor_ref(@200 mA)/monitor_mesure, cela donne les intensites en CPS
		"""
		if self.t2_det_combobox.get_active_text() == "D5":
			detector = libX.Detector()
		elif self.t2_det_combobox.get_active_text() == "S70":
			detector = libX.Detector(nModules=1)
		else:
			self.popup_info("warning","Please specify the detector model.")
		edf_img = fabio.open(edf)
		data = edf_img.data#edfFile.data
		header = edf_img.header#edfFile.header
		counter_name = header['counter_mne'].split()
		counter_value= header['counter_pos'].split()
		counters= {}
		for i in range(len(counter_name)):
			counters[counter_name[i]] = float(counter_value[i])
		monitor_col = self.t2_mon_combobox.get_active_text()
		monitor = counters[monitor_col]
		#Ajouter les gaps ## Code de Clement Buton @Soleil
		detector.set_data(array=data)
		detector.adjust_data()
		detector.set_physical_data()
		detector.reshape_pixels()
		data = detector.physical.data
		data = N.rot90(data)
		
		#Normalisation
		if normalisation:
			monitor_ref = self.monitor_ref
			norm_factor = monitor_ref/monitor
			data = data * norm_factor
		#sauver EDF apres correction
		name = edf.split("/")[-1]
		name_adjusted = name.split(".")[0]+"_corrected"
		if self.ascii_out.get_active():
			ext = "dat"
		else:
			ext  = name.split(".")[1]
		filename = adjusted_folder+"/"+name_adjusted+"."+ext
		if self.ascii_out.get_active():
			N.savetxt(filename, data, header=str(header))
		else:
			edf_img.data = data
			edf_img.write(filename)
		
	def processing(self):
		self.progressbar.show()
		try:
			src_folder = self.src_folder
			des_folder = self.des_folder
			try:
				img_beg = int(self.t2_img_start_entry.get_text())
				img_end = int(self.t2_img_end_entry.get_text())
				beg_end = True
			except:
				beg_end = False
			try:
				self.monitor_ref = float(self.t2_ref_mon_entry.get_text())
				normalisation = True
			except:
				self.monitor_ref = 1
				normalisation = False
			
			filelist=[i for i in listdir(src_folder) if isfile(join(src_folder,i)) and i.endswith(".edf") or i.endswith(".edf.gz") or i.endswith(".dat") or i.endswith(".dat.gz")]
			if beg_end:
				edf_list = select_files_from_list(filelist, img_beg, img_end)
			else:
				edf_list = filelist
			total = len(edf_list)
			processed = 0
			for edf in edf_list:
				try:
					edf_base = edf
					edf = join(src_folder, edf_base)
					if self.t2_det_combobox.get_active_text() == "D1":
						self.geometric_correction_D1(edf, des_folder)
					else:
						self.transform(edf, des_folder, normalisation)
					out_info = "Image %s saved successfully!"%edf_base.split(".")[0]
				except:
					out_info = "Image %s does not existe or cannot be corrected!"%edf
				processed +=1.
				fraction = processed/total *100.
				fraction = int(fraction)
				self.progressbar.set_fraction(fraction/100.)
				self.progressbar.set_text(str(fraction)+"%")
				self.progressbar.set_show_text(True)
				#while gtk.events_pending():
					#gtk.main_iteration()
				self.show_process_info.set_text(out_info)
				yield True	
		except:
			self.popup_info("warning","Please check that you have carefully entered all of the fields.")
		#self.progressbar.hide()
		yield False
		
	def batch_transform(self,widget):
		task = self.processing()
		gobject.idle_add(task.next)
#***********************************************************************************************
#                                       FINISHED
#***********************************************************************************************
if __name__ == "__main__":
	gobject.threads_init()
	MyMainWindow()
	gtk.main()
