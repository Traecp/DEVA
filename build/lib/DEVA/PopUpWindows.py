#!/usr/bin/python
# -*- coding: utf-8 -*-
import gtk
import numpy as N
from scipy import ndimage, stats

from matplotlib.figure import Figure
from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas
from matplotlib.backends.backend_gtkagg import NavigationToolbar2GTKAgg as NavigationToolbar
from matplotlib.widgets import Cursor


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
