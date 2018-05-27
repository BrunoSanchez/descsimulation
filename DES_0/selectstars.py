import math
import numpy as np
from pylab import *
import os
import sys

def selectstars(allstars,plot):

	nstars=allstars.shape[0]
	ncol=allstars.shape[1]
	goodstars=np.zeros((nstars,ncol),float) # good stars
    
		
	#extract nobjameters
	e1=allstars[:,18]
	erre1=allstars[:,34]
	erre2=allstars[:,35]
	ab=allstars[:,10]
	e2=allstars[:,19]
	posx= allstars[:,2]+allstars[:,6]
	posy= allstars[:,3]+allstars[:,7]
	e=(e1**2+e2**2)**0.5
	theta=(np.arctan2(e2,e1)/2.0)+(np.pi/2.0)
	
	#print '------------------------------------------------------------------------'
	#print '                        START THE STAR SELECTION                        '
	#print '========================================================================'
	
	#discard stars with e greater than 0.2
	
	stars=allstars[(e<0.2),:]
	i=stars.shape[0]
	#~ print 'primer selec', i
	
    
	e1=stars[:,18]
	e2=stars[:,19]
	posx= stars[:,2]+stars[:,6]
	posy= stars[:,3]+stars[:,7]
	ab=stars[:,10]
	e=(e1**2+e2**2)**0.5
	theta=(np.arctan2(e2,e1)/2.0)+(np.pi/2.0)
	a=(ab*((1.0+e)/(1.0-e)))**0.5
	a1=a*np.cos(theta)
	a2=a*np.sin(theta)
	if plot in ('s', 'S', 'si', 'Si', 'SI'):
		posygraf=insert(posy,0,9500)
		posxgraf=insert(posx,0,10.0)
		a1graf=insert(a1,0,3.0)
		a2graf=insert(a2,0,0.0)
		
		quiver(posxgraf,posygraf,a1graf,a2graf,headlength=10)
		plt.title('Mayor axis map after PSF correction')
	
	
	#select the stars for the psf map
	
	ngood=0 
	vecinas=10 #this is something you may be will want to change
	
	for j in range(i):
		x=posx[j]
		y=posy[j]
		difpos=((posx-x)**2+(posy-y)**2)**0.5
		cercanas=argsort(difpos)[1:vecinas+1]
		e1close=(e1[cercanas])
		e2close=(e2[cercanas])
		thetaprom=((np.arctan2(e2close,e1close)/2.0)+(np.pi/2.0)).mean()
		sigma=(np.sqrt(e1close**2+e2close**2)).std()	
		difej=np.sqrt(e1[j]**2.+e2[j]**2.)-(np.sqrt(e1close**2+e2close**2)).mean()

		if difej < 2.0*sigma and theta[j] < thetaprom + (np.pi/4.0) and theta[j] > thetaprom - (np.pi/4.0) :
			goodstars[ngood,:]=stars[j,:]
			ngood=ngood+1
	
	print 'cantidad de estrellas seleccionadas',ngood
	
	goodstars=goodstars[:ngood,:]
	
	if plot in ('s', 'S', 'si', 'Si', 'SI'):
		stars=goodstars
		e1=stars[:,18]
		e2=stars[:,19]
		posx= stars[:,2]+stars[:,6]
		posy= stars[:,3]+stars[:,7]
		ab=stars[:,10]
		e=(e1**2+e2**2)**0.5
		theta=(np.arctan2(e2,e1)/2.0)+(np.pi/2.0)
		a=(ab*((1.0+e)/(1.0-e)))**0.5
		a1=a*np.cos(theta)
		a2=a*np.sin(theta)
		
		posygraf=insert(posy,0,9500)
		posxgraf=insert(posx,0,10.0)
		a1graf=insert(a1,0,3.0)
		a2graf=insert(a2,0,0.0)
		print 'cantidad de estrellas seleccionadas',ngood
		quiver(posxgraf,posygraf,a1graf,a2graf,headlength=10,color='r')
		plt.title('Mayor axis map after PSF correction')
		plt.show()
    
	
	return goodstars
			
		

