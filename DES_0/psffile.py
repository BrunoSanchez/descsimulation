import math
import numpy as np
import os
import sys

def psffile(goodstars,objcat,vecinas,psfout):

		idgood=goodstars[:,1]
		e1=goodstars[:,18]
		e2=goodstars[:,19]
		amp=goodstars[:,11]
		ab=goodstars[:,10]
		posx= goodstars[:,2]+goodstars[:,6]
		posy= goodstars[:,3]+goodstars[:,7]
		e=(e1**2+e2**2)**0.5
	
		nobj=objcat.shape[0]
		psfobj=np.zeros((nobj,8),float)
		
		#print 'n gx',nobjcat	
		
		if nobj == 1:
			gxposx= np.array([objcat[1]])
			gxposy= np.array([objcat[2]])
		else:
			gxposx= objcat[:,1]
			gxposy= objcat[:,2]

		for n in range(nobj):
			x=gxposx[n]
			y=gxposy[n]
			difpos=((posx-x)**2+(posy-y)**2)**0.5
			cercanas=np.argsort(difpos)[0:vecinas]
			prome1=0.0	
			prome2=0.0
			promab=0.0
			for k in range(len(cercanas)):
				prome1=prome1+e1[cercanas[k]]
				prome2=prome2+e2[cercanas[k]]	
				promab=promab+ab[cercanas[k]]
		
			prome1=prome1/float(vecinas)
			prome2=prome2/float(vecinas)
			promab=promab/float(vecinas)
			thetaprom=((np.arctan2(prome2,prome1)/2.0))
			eprom=(prome1**2+prome2**2)**0.5
			psfobj[n,:]=[x,y,0.0,0.0,eprom,thetaprom,promab,1.0]
		
		#print 'Writing the psf file for galaxies in ',psfobjfile
		f1=open(psfout,'w')
		f1.write('1\n')
		f1.write(str(nobj))
		f1.write('\n')
		f1.write('6\n')
		np.savetxt(f1, psfobj, fmt=['%15.10f']*4+['%13.10f']*4)
		f1.close()
		
		return psfout
