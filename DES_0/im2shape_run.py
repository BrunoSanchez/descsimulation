import math
import numpy as np
import os
import sys


def im2shape_run(imagen,catfilein,catfileout,psffile,num,chains,niter):
		print 'corriendo im2shape'
		print imagen,catfilein,catfileout,psffile,num,chains,niter
		#print '----------------------------------------------------------------'
		#print '             RUNNING im2shape ON THE STARS                      '
		#print '================================================================'
			
		#print'-------------------------------------------------------------------------------------'
		#print'------------------------IMAGEN'+str(nobj)+'--------------------------------------------'	
		os.system('rm '+catfileout)		
		f1=open('conf'+str(num)+'.in','w')
		f1.write('0 \n')
		f1.write(catfileout+' \n')
		f1.write(imagen+'\n')
		f1.write(catfilein+'\n')
		f1.write('16\n')
		f1.write('2 0.000000 20.000000\n')
		f1.write('1 0.000000 1000.000000\n')
		f1.write('1\n')
		f1.write('1 0.000000 16.000000\n')
		f1.write('1 0.000000 16.000000\n')
		f1.write('3 0.000000 1.000000\n')
		f1.write('1 0.000000 3.141600\n')
		f1.write('2 -2.000000 8.000000\n')
		f1.write('2 0.000000 20.000000\n')
		f1.write(psffile+'\n')
		f1.write(chains+'\n')
		f1.write(niter+'\n')
		f1.write('model.txt\n')
		f1.write('2.000000\n')
		f1.write('0\n')
		f1.close()	
		
		os.system('./im2shape conf'+str(num)+'.in > im2shape.output')
		
		
		f = open(catfileout, 'r')
		linecat=f.readline()
		ncolcat = int(linecat[:2])
		nobj = int(linecat[3:-1])
		f.close()
		#print ncolcat,ngxcat
		cat = np.loadtxt(catfileout, comments='#',skiprows=1) #sextractor catalogue for galaxies
		
		return cat, nobj, ncolcat
	
		


