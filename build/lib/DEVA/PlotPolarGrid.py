#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy as N
"""
Plot the polar grid on a specific axe
"""
def plotCircles(ax, min_r, max_r, nbr):
	phi=N.linspace(0,2*N.pi,50)
	r = N.linspace(min_r,max_r*0.95,nbr)
	for i in r:
		ax.plot(i*N.sin(phi),i*N.cos(phi),'k:')
	decorChi(ax, min_r, max_r, nbr)
	
def plotLines(ax, r, nbr):
	phi = N.arange(0,360,360/nbr)
	phi = N.radians(phi)
	v   = N.linspace(-r,r,50)
	for p in phi:
		x = v*N.cos(p)
		y = v*N.sin(p)
		ax.plot(x,y,"k:")
	decorPhi(ax, r, nbr)
	
def decorChi(ax, r_min, r_max, nbr):
	Tx=N.linspace(r_min, r_max*0.95, nbr)
	Ty=Tx*N.tan(N.radians(15.))
	for j in range(1,len(Tx)):
		ax.text(Tx[j],Ty[j],"%.1f"%(N.degrees(Tx[j])),color="k",horizontalalignment="center", verticalalignment="center",fontsize=14)

def decorPhi(ax, r, nbr):
	phi_d = N.arange(0,360,360/nbr)
	phi = N.radians(phi_d)
	x   = r*N.cos(phi)*1.06
	y   = r*N.sin(phi)*1.06
	for i in range(len(x)):
		ax.text(x[i],y[i],"%d"%phi_d[i],horizontalalignment="center", verticalalignment="center",fontsize=14)
		
def plotPolarGrid(ax, min_r, max_r, nbr_circles, nbr_lines):
	plotCircles(ax, min_r, max_r, nbr_circles)
	plotLines(ax, max_r, nbr_lines)