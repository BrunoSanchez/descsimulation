import math
import numpy as np
from pylab import *
import cosmolopy.distance as cd
from multiprocessing import Pool
from multiprocessing import Process
import os
import pyfits
import sys
from pyraf import iraf

cosmo = {'omega_M_0':0.3, 'omega_lambda_0':0.7, 'omega_k_0':0.0, 'h':0.7}



def seeing(imagen,pixsize,corrida,filtro,magmax,magmin,fwhmmax,plot):

	

	
	
	#run sextrator

	sexfile='first'+filtro+str(corrida)+'.sex' #name of configuration file of SExtractor
	salida_sex='run1-'+filtro+str(corrida)+'.cat' #exit catalogue from SExtractor

	#print '----------------------------------------------------------'
	#print '           WRITING SExtractor CONFIGURATION FILE          '
	#print '=========================================================='
	
	os.system('rm '+sexfile)
	os.system('rm '+salida_sex)	
	f1=open(sexfile,'w')
	f1.write('#-------------------------------- Catalog ------------------------------------\n')
	f1.write('CATALOG_NAME     '+salida_sex+'   # name of the output catalog\n')
	f1.write('CATALOG_TYPE     ASCII_HEAD     # NONE,ASCII,ASCII_HEAD, ASCII_SKYCAT,\n')
	f1.write('                                # ASCII_VOTABLE, FITS_1.0 or FITS_LDAC\n')
	f1.write('PARAMETERS_NAME  first.param     # name of the file containing catalog contents \n')
	f1.write('#------------------------------- Extraction ----------------------------------\n')
	f1.write('DETECT_TYPE		CCD			# "CCD" or "PHOTO" (*)\n')
	f1.write('FLAG_IMAGE		flag.fits		# filename for an input FLAG-image\n')
	f1.write('DETECT_MINAREA	5			# minimum number of pixels above threshold\n')
	f1.write('DETECT_THRESH		5. 			# <sigmas> or <threshold>,<ZP> in mag.arcsec-2\n')
	f1.write('ANALYSIS_THRESH	2.			# <sigmas> or <threshold>,<ZP> in mag.arcsec-2\n')
	f1.write('FILTER		Y			# apply filter for detection ("Y" or "N")?\n')
	f1.write('FILTER_NAME	  	default.conv	# name of the file containing the filter\n')
	f1.write('DEBLEND_NTHRESH	32			# Number of deblending sub-thresholds\n')
	f1.write('DEBLEND_MINCONT	0.005			# Minimum contrast parameter for deblending\n')
	f1.write('CLEAN			Y			# Clean spurious detections? (Y or N)?\n')
	f1.write('CLEAN_PARAM		1.0			# Cleaning efficiency\n')
	f1.write('MASK_TYPE		CORRECT		# type of detection MASKing: can be one of\n')
	f1.write('					# "NONE", "BLANK" or "CORRECT"\n')
	f1.write('#------------------------------ Photometry -----------------------------------\n')
	f1.write('PHOT_APERTURES	10			# MAG_APER aperture diameter(s) in pixels\n')
	f1.write('PHOT_AUTOPARAMS	2.5, 3.5		# MAG_AUTO parameters: <Kron_fact>,<min_radius>\n')
	f1.write('SATUR_LEVEL		65000.0		# level (in ADUs) at which arises saturation\n')
	f1.write('MAG_ZEROPOINT		26.73			# magnitude zero-point\n')
	f1.write('MAG_GAMMA		4.0			# gamma of emulsion (for photographic scans)\n')
	f1.write('GAIN			2.00			# detector gain in e-/ADU.\n')
	f1.write('PIXEL_SCALE		'+str(pixsize)+'	# size of pixel in arcsec (0=use FITS WCS info).\n')
	f1.write('#------------------------- Star/Galaxy Separation ----------------------------\n')
	f1.write('SEEING_FWHM		0.498			# stellar FWHM in arcsec\n')
	f1.write('STARNNW_NAME	default.nnw		# Neural-Network_Weight table filename\n')
	f1.write('#------------------------------ Background -----------------------------------\n')
	f1.write('BACK_SIZE		64			# Background mesh: <size> or <width>,<height>\n')
	f1.write('BACK_FILTERSIZE	3			# Background filter: <size> or <width>,<height>\n')
	f1.write('BACKPHOTO_TYPE	GLOBAL		# can be "GLOBAL" or "LOCAL" (*)\n')
	f1.write('BACKPHOTO_THICK	24			# thickness of the background LOCAL annulus (*)\n')
	f1.write('#------------------------------ Check Image ----------------------------------\n')
	f1.write('CHECKIMAGE_TYPE	NONE			# can be one of "NONE", "BACKGROUND",\n')
	f1.write('						# "MINIBACKGROUND", "-BACKGROUND", "OBJECTS",\n')
	f1.write('						# "-OBJECTS", "SEGMENTATION", "APERTURES",\n')
	f1.write('						# or "FILTERED" (*)\n')
	f1.write('CHECKIMAGE_NAME	apertures.fits	# Filename for the check-image (*)\n')
	f1.write('#--------------------- Memory (change with caution!) -------------------------\n')
	f1.write('MEMORY_OBJSTACK	8000			# number of objects in stack\n')
	f1.write('MEMORY_PIXSTACK	400000		# number of pixels in stack\n')
	f1.write('MEMORY_BUFSIZE	1024			# number of lines in buffer\n')
	f1.write('#----------------------------- Miscellaneous ---------------------------------\n')
	f1.write('VERBOSE_TYPE	NORMAL		# can be "QUIET", "NORMAL" or "FULL" (*)\n')
	f1.close()


	callsex='sex '+imagen+' -c '+sexfile+' > sex_output'
		
	os.system(callsex)

	#print imagen
	#wait=raw_input('chequear y presionar ENTER')
	cat = np.loadtxt(salida_sex, comments='#') #reads the catalogue of the first run of sextractor
	
	#print '----------------------------------------------------------'
	#print '              COMPUTE SATURATION LEVEL                    '
	#print '=========================================================='


	#compute the saturation level
	FLUXMAX = cat[:,8]
	SATUR=0.8*FLUXMAX.max()
	#print 'SATUR_LEVEL =',SATUR
	
	#print '----------------------------------------------------------'
	#print '                COMPUTE THE SEEING                        '
	#print '=========================================================='


	#now, starts to compute the seeing
	
	FWHM = cat[:,5]
	FLAG = cat[:,11]
	MAGBEST = cat[:,6]
	nobj=len(MAGBEST)
	fwhm=np.zeros((nobj),float)
	mag=np.zeros((nobj),float)
		
	if plot in ('s', 'S', 'si', 'Si', 'SI'):
		plt.plot(FWHM,MAGBEST, 'b.')
		plt.ylabel('MAG_BEST')
		plt.xlabel('FWHM')
		plt.show()
	
	
	j=0
	


	#make the cuts in magnitude to get mostly the stars
	mask2=(MAGBEST < magmax)*(MAGBEST > magmin)*(FWHM < fwhmmax)
	fwhm=FWHM[mask2]
	mag=MAGBEST[mask2]
	j=len(fwhm)
	
	
	
	
	
	
	
	
	#and get the maximun value of the FWHM distribution
	
	
	step=0.05
	intervalos=int((fwhm.max()-fwhm.min())/step)
	interval1=fwhm.min()
	interval2=interval1+step
	maximo=0
	moda=90.
	for j in range(intervalos):
		contar=0
		for x in fwhm:
			if x > interval1 and x < interval2:
				contar=contar+1
			if contar > maximo:
				maximo=contar
				moda=(interval1+interval2)/2.
		interval1=interval1+step
		interval2=interval2+step
	
		
	if plot in ('s', 'S', 'si', 'Si', 'SI'):	
		plt.plot(FWHM,MAGBEST,'k.')
		plt.plot(fwhm,mag,'r.')
		plt.show()
		print 'SEEING in pix',moda
		plt.hist(fwhm,15)
		plt.show()
	
	#print 'SEEING in arcsec',moda*pixsize
	seeing=moda*pixsize
	print ' '
	print ' ------------ seeing: ',seeing
	print ' ------------ imagen: ',imagen
	print ' '
	puntox=np.array([moda],float)
	puntoy=np.array([18.],float)

	if plot in ('s', 'S', 'si', 'Si', 'SI'):
		plt.plot(FWHM,MAGBEST, 'b.',puntox,puntoy,'ro')
		plt.ylabel('MAG_BEST')
		plt.xlabel('FWHM')
		#~ plt.axis([0,30,16,26])
		plt.show()
		
	
	
	#~ print '----------------------------------------------------------'
	#~ print '           WRITING SExtractor CONFIGURATION FILE          '
	#~ print '=========================================================='
	
	sexfile='second'+filtro+str(corrida)+'.sex' #name of configuration file of SExtractor
	salida_sex='run2-'+filtro+str(corrida)+'.cat' #exit catalogue from SExtractor

	
	os.system('rm '+salida_sex)
	os.system('rm '+sexfile)
	#~ print corrida
	
	f1=open(sexfile,'w')
	f1.write('#-------------------------------- Catalog ------------------------------------\n')
	f1.write('CATALOG_NAME     '+salida_sex+'   # name of the output catalog\n')
	f1.write('CATALOG_TYPE     ASCII_HEAD     # NONE,ASCII,ASCII_HEAD, ASCII_SKYCAT,\n')
	f1.write('                                # ASCII_VOTABLE, FITS_1.0 or FITS_LDAC\n')
	f1.write('PARAMETERS_NAME  second.param     # name of the file containing catalog contents \n')
	f1.write('#------------------------------- Extraction ----------------------------------\n')
	f1.write('DETECT_TYPE		CCD			# "CCD" or "PHOTO" (*)\n')
	f1.write('FLAG_IMAGE		flag.fits		# filename for an input FLAG-image\n')
	f1.write('DETECT_MINAREA	5			# minimum number of pixels above threshold\n')
	f1.write('DETECT_THRESH	1.5 			# <sigmas> or <threshold>,<ZP> in mag.arcsec-2\n')
	f1.write('ANALYSIS_THRESH	2.			# <sigmas> or <threshold>,<ZP> in mag.arcsec-2\n')
	f1.write('FILTER		Y			# apply filter for detection ("Y" or "N")?\n')
	f1.write('FILTER_NAME	  	default.conv	# name of the file containing the filter\n')
	f1.write('DEBLEND_NTHRESH	32			# Number of deblending sub-thresholds\n')
	f1.write('DEBLEND_MINCONT	0.005			# Minimum contrast parameter for deblending\n')
	f1.write('CLEAN			Y			# Clean spurious detections? (Y or N)?\n')
	f1.write('CLEAN_PARAM		1.0			# Cleaning efficiency\n')
	f1.write('MASK_TYPE		CORRECT		# type of detection MASKing: can be one of\n')
	f1.write('						# "NONE", "BLANK" or "CORRECT"\n')
	f1.write('#------------------------------ Photometry -----------------------------------\n')
	f1.write('PHOT_APERTURES	10			# MAG_APER aperture diameter(s) in pixels\n')
	f1.write('PHOT_AUTOPARAMS	2.5, 3.5		# MAG_AUTO parameters: <Kron_fact>,<min_radius>\n')
	f1.write('SATUR_LEVEL      '+str(int(SATUR))+'        # level (in ADUs) at which arises saturation\n')
	f1.write('MAG_ZEROPOINT	26.73			# magnitude zero-point\n')
	f1.write('MAG_GAMMA		4.0			# gamma of emulsion (for photographic scans)\n')
	f1.write('GAIN			2.00			# detector gain in e-/ADU.\n')
	f1.write('PIXEL_SCALE     '+str(pixsize)+'        # size of pixel in arcsec (0=use FITS WCS info)\n')
	f1.write('#------------------------- Star/Galaxy Separation ----------------------------\n')
	f1.write('SEEING_FWHM      '+str('%.1f' % seeing)+'            # stellar FWHM in arcsec\n')
	f1.write('STARNNW_NAME	default.nnw		# Neural-Network_Weight table filename\n')
	f1.write('#------------------------------ Background -----------------------------------\n')
	f1.write('BACK_SIZE		64			# Background mesh: <size> or <width>,<height>\n')
	f1.write('BACK_FILTERSIZE	3			# Background filter: <size> or <width>,<height>\n')
	f1.write('BACKPHOTO_TYPE	GLOBAL		# can be "GLOBAL" or "LOCAL" (*)\n')
	f1.write('BACKPHOTO_THICK	24			# thickness of the background LOCAL annulus (*)\n')
	f1.write('#------------------------------ Check Image ----------------------------------\n')
	f1.write('CHECKIMAGE_TYPE	NONE			# can be one of "NONE", "BACKGROUND",\n')
	f1.write('						# "MINIBACKGROUND", "-BACKGROUND", "OBJECTS",\n')
	f1.write('						# "-OBJECTS", "SEGMENTATION", "APERTURES",\n')
	f1.write('						# or "FILTERED" (*)\n')
	f1.write('CHECKIMAGE_NAME	apertures.fits	# Filename for the check-image (*)\n')
	f1.write('#--------------------- Memory (change with caution!) -------------------------\n')
	f1.write('MEMORY_OBJSTACK	8000			# number of objects in stack\n')
	f1.write('MEMORY_PIXSTACK	400000		# number of pixels in stack\n')
	f1.write('MEMORY_BUFSIZE	1024			# number of lines in buffer\n')
	f1.write('#----------------------------- Miscellaneous ---------------------------------\n')
	f1.write('VERBOSE_TYPE	NORMAL		# can be "QUIET", "NORMAL" or "FULL" (*)\n')
	f1.write('#------------------------------- New Stuff -----------------------------------\n')
	f1.close()
	
	del(cat)
	del(FWHM)
	del(FLAG)
	del(MAGBEST)
	del(nobj)
	del(fwhm)
	del(mag)
	
	
	
	return seeing
