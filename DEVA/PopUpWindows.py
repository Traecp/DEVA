#!/usr/bin/python
# -*- coding: utf-8 -*-
import gtk
import numpy as N
from scipy import ndimage, stats

import matplotlib as mpl
from matplotlib.figure import Figure
from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas
from matplotlib.backends.backend_gtkagg import NavigationToolbar2GTKAgg as NavigationToolbar
from matplotlib.widgets import Cursor
from matplotlib.ticker import MaxNLocator
#mpl.rcParams['font.size'] = 18.0
mpl.rcParams['axes.labelsize'] = 'medium'
mpl.rcParams['legend.fancybox'] = True
mpl.rcParams['legend.handletextpad'] = 0.5
mpl.rcParams['legend.fontsize'] = 'medium'
mpl.rcParams['figure.subplot.bottom'] = 0.25
#mpl.rcParams['figure.subplot.top'] = 0.93
#mpl.rcParams['figure.subplot.left'] = 0.14
#mpl.rcParams['figure.subplot.right'] = 0.915
mpl.rcParams['image.cmap'] = 'jet'


class PopUpFringes(object):
	def __init__(self, xdata, xlabel, ylabel, title):
		self.popupwin=gtk.Window()
		self.popupwin.set_size_request(600,550)
		self.popupwin.set_position(gtk.WIN_POS_CENTER)
		self.popupwin.set_border_width(10)
		self.xdata = xdata
		vbox = gtk.VBox()
		self.fig=Figure(dpi=100)
		self.ax  = self.fig.add_subplot(111)
		self.canvas  = FigureCanvas(self.fig)
		self.main_figure_navBar = NavigationToolbar(self.canvas, self)
		self.cursor = Cursor(self.ax, color='k', linewidth=1, useblit=True)
		self.ax.set_xlabel(xlabel, fontsize = 18)
		self.ax.set_ylabel(ylabel, fontsize = 18)
		self.ax.set_title(title, fontsize = 18)
		
		xi = N.arange(len(self.xdata))		
		slope, intercept, r_value, p_value, std_err = stats.linregress(self.xdata,xi)
		fitline = slope*self.xdata+intercept
		
		self.ax.plot(self.xdata, fitline, 'r-',self.xdata,xi, 'bo')
		self.ax.axis([self.xdata.min(),self.xdata.max(),xi.min()-1, xi.max()+1])
		
		self.ax.text(0.3, 0.9,'Slope = %.4f +- %.4f' % (slope, std_err),
								horizontalalignment='center',
								verticalalignment='center',
								transform = self.ax.transAxes,
								color='red')
		vbox.pack_start(self.main_figure_navBar, False, False, 0)
		vbox.pack_start(self.canvas, True, True, 2)
		self.popupwin.add(vbox)
		self.popupwin.connect("destroy", self.dest)
		self.popupwin.show_all()
	
	def dest(self,widget):
		self.popupwin.destroy()
	
class PopUpImage(object):
	def __init__(self, xdata, ydata, xlabel, ylabel, title):
		self.popupwin=gtk.Window()
		self.popupwin.set_size_request(600,550)
		self.popupwin.set_position(gtk.WIN_POS_CENTER)
		self.popupwin.set_border_width(10)
		self.xdata = xdata
		self.ydata = ydata
		vbox = gtk.VBox()
		self.fig=Figure(dpi=100)
		self.ax  = self.fig.add_subplot(111)
		self.canvas  = FigureCanvas(self.fig)
		self.main_figure_navBar = NavigationToolbar(self.canvas, self)
		self.cursor = Cursor(self.ax, color='k', linewidth=1, useblit=True)
		self.canvas.mpl_connect("button_press_event",self.on_press)
		self.ax.set_xlabel(xlabel, fontsize = 18)
		self.ax.set_ylabel(ylabel, fontsize = 18)
		self.ax.set_title(title, fontsize = 18)
		self.ax.plot(self.xdata, self.ydata, 'b-', lw=2)
		
		self.textes = []
		self.plots  = []
		vbox.pack_start(self.main_figure_navBar, False, False, 0)
		vbox.pack_start(self.canvas, True, True, 2)
		self.popupwin.add(vbox)
		self.popupwin.connect("destroy", self.dest)
		self.popupwin.show_all()
	
	def dest(self,widget):
		self.popupwin.destroy()
	
	def on_press(self, event):
		if event.inaxes == self.ax and event.button==3:
			self.clear_notes()
			xc = event.xdata
			#***** Find the closest x value *****
			residuel = self.xdata - xc
			residuel = N.abs(residuel)
			j = N.argmin(residuel)
			#y = self.ydata[i-1:i+1]
			#yc= y.max()
			#j = N.where(self.ydata == yc)
			#j = j[0][0]
			xc= self.xdata[j]
			x_fit = self.xdata[j-3:j+3]
			y_fit = self.ydata[j-3:j+3]
			fitted_param, fitted_data = fit(x_fit, y_fit, xc, True)
			x_fit = N.linspace(x_fit.min(), x_fit.max(), 200)
			y_fit = psdVoigt(fitted_param, x_fit)
			period = fitted_param['xc'].value
			std_err= fitted_param['xc'].stderr
			
			p = self.ax.plot(x_fit, y_fit,'r-')
			p2 = self.ax.axvline(period,color='green',lw=2)
			
			txt=self.ax.text(0.05, 0.9, 'Period = %.4f +- %.4f (nm)'%(period, std_err), transform = self.ax.transAxes, color='red')
			self.textes.append(txt)
			self.plots.append(p[0])
			self.plots.append(p2)
		elif event.inaxes == self.ax and event.button==2:
			dif = N.diff(self.ydata)
			dif = dif/dif.max()
			p3  = self.ax.plot(dif,'r-')
			self.plots.append(p3[0])
		self.canvas.draw()
	
	def clear_notes(self):
		if len(self.textes)>0:
			for t in self.textes:
				t.remove()
		if len(self.plots)>0:
			for p in self.plots:
				p.remove()
		self.textes = []
		self.plots  = []

class PopUpRSM(object):
	def __init__(self, mapData,H,K,L):
#		super(PopUpRSM, self).__init__()
		self.window = gtk.Window()
		self.window.set_size_request(1100,850)
		self.window.set_position(gtk.WIN_POS_CENTER)
		self.window.set_border_width(10)
		self.mapData = mapData
		self.H = H
		self.K = K
		self.L = L
		self.X_pts, self.Y_pts, self.Z_pts = self.mapData.shape
		self.center_X = self.X_pts/2
		self.center_Y = self.Y_pts/2
		self.center_Z = self.Z_pts/2
		

		self.slicing_frame = gtk.Frame()
		self.slicing_frame.set_label("Slicing control. Nx,Ny,Nz = %d, %d, %d"%(self.X_pts, self.Y_pts, self.Z_pts))
		self.slicing_frame.set_label_align(0.5,0.5)
		self.scale_frame = gtk.Frame()
		self.scale_frame.set_label("Color scale control")
		self.scale_frame.set_label_align(0.5,0.5)

		#Slicing sliders
		self.slicing_box = gtk.HBox(False,2)

		self.slicing_reset_btn = gtk.Button("Reset")
		self.slicing_reset_btn.connect("clicked", self.reset_image)

		self.activate_slicing_btn = gtk.ToggleButton("Activate slicing")
		self.activate_slicing_btn.connect("toggled", self.init_image)

		self.hk_txt = gtk.Label("HK slice")
		self.hk_txt.set_alignment(0,0.5)
		self.kl_txt = gtk.Label("KL slice")
		self.kl_txt.set_alignment(0,0.5)
		self.hl_txt = gtk.Label("HL slice")
		self.hl_txt.set_alignment(0,0.5)
		self.slider_hk = gtk.HScale()
		self.slider_kl = gtk.HScale()
		self.slider_hl = gtk.HScale()

		self.slider_hk.set_size_request(200,25)
		self.slider_kl.set_size_request(200,25)
		self.slider_hl.set_size_request(200,25)

		self.slider_hk.set_range(0,self.Z_pts-1)
		self.slider_kl.set_range(0,self.X_pts-1)
		self.slider_hl.set_range(0,self.Y_pts-1)

		self.slider_hk.set_value(self.center_Z)
		self.slider_kl.set_value(self.center_X)
		self.slider_hl.set_value(self.center_Y)
		
		self.slider_hk.set_digits(0)
		self.slider_hk.set_increments(1,1)
		self.slider_kl.set_digits(0)
		self.slider_kl.set_increments(1,1)
		self.slider_hl.set_digits(0)
		self.slider_hl.set_increments(1,1)

		self.slider_hk.connect('value-changed',self.slice_update, "HK")
		self.slider_kl.connect('value-changed',self.slice_update, "KL")
		self.slider_hl.connect('value-changed',self.slice_update, "HL")

		self.slicing_box.pack_start(self.slicing_reset_btn,False,False,0)
		self.slicing_box.pack_start(self.activate_slicing_btn,False,False,0)
		self.slicing_box.pack_start(self.hk_txt,False,False,0)
		self.slicing_box.pack_start(self.slider_hk,False,False,10)
		self.slicing_box.pack_start(self.kl_txt,False,False,0)
		self.slicing_box.pack_start(self.slider_kl,False,False,10)
		self.slicing_box.pack_start(self.hl_txt,False,False,0)
		self.slicing_box.pack_start(self.slider_hl,False,False,10)

		self.vmin = self.mapData.min()
		self.vmax = self.mapData.max()

		self.scale_box = gtk.HBox(False, 2)
		self.vmin_txt = gtk.Label("Vmin")
		self.vmin_txt.set_alignment(0,0.5)
		self.vmax_txt = gtk.Label("Vmax")
		self.vmax_txt.set_alignment(0,0.5)
		self.sld_vmin = gtk.HScale()
		self.sld_vmax = gtk.HScale()
		self.slicing_activated = False

		self.sld_vmin.set_size_request(200,25)
		self.sld_vmax.set_size_request(200,25)
		self.sld_vmin.set_range(self.vmin-0.5*abs(self.vmin),self.vmax)
		self.sld_vmax.set_range(self.vmin,self.vmax*2.0)
		self.sld_vmax.set_value(self.vmax)
		self.sld_vmin.set_value(self.vmin)
		self.sld_vmin.connect('value-changed',self.scale_update)
		self.sld_vmax.connect('value-changed',self.scale_update)

		self.scale_box.pack_start(self.vmin_txt,False,False,0)
		self.scale_box.pack_start(self.sld_vmin,False,False,0)
		self.scale_box.pack_start(self.vmax_txt,False,False,0)
		self.scale_box.pack_start(self.sld_vmax,False,False,0)

		self.slicing_frame_align = gtk.Alignment()
		self.slicing_frame_align.set_padding(15,10,10,10)
		self.slicing_frame_align.set(0.5, 0.5, 1.0, 1.0)
		self.slicing_frame_align.add(self.slicing_box)
		self.slicing_frame.add(self.slicing_frame_align)


		self.scale_frame_align = gtk.Alignment()
		self.scale_frame_align.set_padding(15,10,10,10)
		self.scale_frame_align.set(0.5, 0.5, 1.0, 1.0)
		self.scale_frame_align.add(self.scale_box)
		self.scale_frame.add(self.scale_frame_align)


		self.fig=Figure(dpi=100)
		self.ax1  = self.fig.add_axes([0.1, 0.25, 0.2, 0.7])
		self.ax2  = self.fig.add_axes([0.4, 0.25, 0.2, 0.7])
		self.ax3  = self.fig.add_axes([0.7, 0.25, 0.2, 0.7])
		self.canvas  = FigureCanvas(self.fig)
		self.navBar = NavigationToolbar(self.canvas, self)
		#self.canvas.mpl_connect("button_press_event",self.on_press)
		self.ax1.set_xlabel("H", fontsize = 18)
		self.ax1.set_ylabel("K", fontsize = 18)
		self.ax2.set_xlabel("K", fontsize = 18)
		self.ax2.set_ylabel("L", fontsize = 18)
		self.ax3.set_xlabel("H", fontsize = 18)
		self.ax3.set_ylabel("L", fontsize = 18)
		#self.ax.set_title(title, fontsize = 18)
		self.img1 = self.ax1.contourf(self.H, self.K, self.mapData.sum(axis=2), 50)
		self.img2 = self.ax2.contourf(self.K, self.L, self.mapData.sum(axis=0), 50)
		self.img3 = self.ax3.contourf(self.H, self.L, self.mapData.sum(axis=1), 50)
		self.cax = self.fig.add_axes([0.15, 0.1, 0.7, 0.03])#left,bottom,width,height
		self.cb  = self.fig.colorbar(self.img1, cax = self.cax, orientation='horizontal', format="%.2f")#format=fm
		self.cb.set_label("Log10_intensity (arb. unit)", fontsize=18)
		self.cb.locator = MaxNLocator(nbins=6)
		self.vbox = gtk.VBox()
		self.vbox.pack_start(self.slicing_frame, False, False, 0)
		self.vbox.pack_start(self.navBar, False, False, 0)
		self.vbox.pack_start(self.canvas, True, True, 2)
		self.vbox.pack_start(self.scale_frame, False, False, 0)

		self.window.add(self.vbox)
		self.window.connect("destroy", self.destroy_me)
		self.window.show_all()
		self.scale_frame.hide()
	def destroy_me(self,widget):
		self.window.destroy()

	def init_image(self,w):
		if self.activate_slicing_btn.get_active():
			self.activate_slicing_btn.set_label("Deactivate slicing")
			self.slicing_activated = True
			self.scale_frame.show()
		
			self.ax1.clear()
			self.ax2.clear()
			self.ax3.clear()
			self.cax.clear()
			self.ax1.set_xlabel("H", fontsize = 18)
			self.ax1.set_ylabel("K", fontsize = 18)
			self.ax2.set_xlabel("K", fontsize = 18)
			self.ax2.set_ylabel("L", fontsize = 18)
			self.ax3.set_xlabel("H", fontsize = 18)
			self.ax3.set_ylabel("L", fontsize = 18)
			self.img1 = self.ax1.imshow(self.mapData[:,:,self.center_Z-10:self.center_Z+10].sum(axis=2)/20.,origin='lower',vmin=self.vmin, vmax=self.vmax*1.3, interpolation='nearest',aspect='auto', extent=(self.H.min(), self.H.max(), self.K.min(), self.K.max()))
			self.img2 = self.ax2.imshow(self.mapData[self.center_X-10:self.center_X+10,:,:].sum(axis=0)/20.,origin='lower',vmin=self.vmin, vmax=self.vmax*1.3, interpolation='nearest',aspect='auto', extent=(self.K.min(), self.K.max(), self.L.min(), self.L.max()))
			self.img3 = self.ax3.imshow(self.mapData[:,self.center_Y-10:self.center_Y+10,:].sum(axis=1)/20.,origin='lower',vmin=self.vmin, vmax=self.vmax*1.3, interpolation='nearest',aspect='auto', extent=(self.H.min(), self.H.max(), self.L.min(), self.L.max()))

			self.cb  = self.fig.colorbar(self.img1, cax = self.cax, orientation='horizontal', format="%.2f")#format=fm
			self.cb.set_label("Log10_intensity (arb. unit)", fontsize=12)
			self.cb.locator = MaxNLocator(nbins=6)
			self.canvas.draw()
		else:
			self.activate_slicing_btn.set_label("Activate slicing")
			self.reset_image_now()


	def reset_image(self, w):
		self.reset_image_now()
	def reset_image_now(self):
		self.slicing_activated = False
		self.scale_frame.hide()
		self.activate_slicing_btn.set_active(False)
		self.activate_slicing_btn.set_label("Activate slicing")
		self.ax1.clear()
		self.ax2.clear()
		self.ax3.clear()
		self.cax.clear()
		self.ax1.set_xlabel("H", fontsize = 18)
		self.ax1.set_ylabel("K", fontsize = 18)
		self.ax2.set_xlabel("K", fontsize = 18)
		self.ax2.set_ylabel("L", fontsize = 18)
		self.ax3.set_xlabel("H", fontsize = 18)
		self.ax3.set_ylabel("L", fontsize = 18)
		#self.ax.set_title(title, fontsize = 18)
		self.img1 = self.ax1.contourf(self.H, self.K, self.mapData.sum(axis=2), 50)
		self.img2 = self.ax2.contourf(self.K, self.L, self.mapData.sum(axis=0), 50)
		self.img3 = self.ax3.contourf(self.H, self.L, self.mapData.sum(axis=1), 50)
		self.cb  = self.fig.colorbar(self.img1, cax = self.cax, orientation='horizontal', format="%.2f")#format=fm
		self.cb.set_label("Log10_intensity (arb. unit)", fontsize=12)
		self.cb.locator = MaxNLocator(nbins=6)
		self.canvas.draw()

	def scale_update(self,widget):
		self.vmin = self.sld_vmin.get_value()
		self.vmax = self.sld_vmax.get_value()
		self.img1.set_clim(self.vmin, self.vmax)
		self.ax1.relim()
		self.img2.set_clim(self.vmin, self.vmax)
		self.ax2.relim()
		self.img3.set_clim(self.vmin, self.vmax)
		self.ax3.relim()
		self.canvas.draw()

	def slice_update(self, widget, target):
		index = widget.get_value()
		index = int(index)
#		print index
		if target=="HK":
			self.img1.set_data(self.mapData[:,:,index-10:index+10].sum(axis=2)/20.)
		elif target=="KL":
			self.img2.set_data(self.mapData[index-10:index+10,:,:].sum(axis=0)/20.)
		elif target=="HL":
			self.img3.set_data(self.mapData[:,index-10:index+10,:].sum(axis=1)/20.)
		self.canvas.draw()
			
