import math
import numpy as np
import os
import sys

def rot(posx,posy,e1,e2,xy0):

	#print '------------------------------------------------------------------------'
	#print '                     CHANGING THE CORDINATE AXIS                        '
	#print '========================================================================'
	

	e=(e1**2+e2**2)**0.5
	theta=(np.arctan2(e2,e1)/2.0)+(np.pi/2.0)
			
	u1 = posx-xy0[0] #Xran[j]
	u2 = posy-xy0[1]
	fi=np.arctan2(u2,u1)
	mask_u2_neg = u2 < 0.0
	fi[mask_u2_neg] += 2.0*np.pi
	beta2=theta-fi#+(np.pi/4.0)
	et=(-1.0*e*np.cos(2.0*beta2))
	ex=(e*np.sin(2.0*beta2))
				
	return et,ex
