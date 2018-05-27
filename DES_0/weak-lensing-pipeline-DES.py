import math
import numpy as np
from pylab import *
import cosmolopy.distance as cd
from multiprocessing import Pool
from multiprocessing import Process
import os
import sys
from star_gx_des import *
from rot import *
from psffile import *
from selectstars import *
from im2shape_run import *
from seeing_des import *
from astropy.io import fits
from fit_gaussian import *

cosmo = {'omega_M_0':0.3, 'omega_lambda_0':0.7, 'omega_k_0':0.0, 'h':0.7}



ok = raw_input('Desea continuar ejecucion (S/N) ')
PSF = raw_input('Chequear PSF (S/N) ')
plot=raw_input('Graficar (S/N) ')
pulls=1#int(raw_input('Enter number of pulls '))
nobj=1#int(raw_input('Emepzar en imagen (0 a n) '))
proc='LN'#raw_input('Maquina ')
path='/home/elizabeth/Documentos/PhD/Analisis-Lentes/DES_cluster/individualimagefiles/Imagen-tipo-6/imagenes-LN-PSF/'
#~ path='./'
pix_size=0.27
pixsize=pix_size

MAGMIN=14.
MAGMAX=35.
z_nobj=0.33 
clase= 0
DELTA0=-32.081192
ALFA0=336.75
extincion=0.




salidafile='GX_DES_'+proc+'.cat'

registro='reg_DES'+proc

print 'Archivos de salida ', salidafile, registro



if ok in ('s', 'S', 'si', 'Si', 'SI'):

	f122=open(registro,'a')

	f8=open(salidafile,'a')

if PSF in ('s', 'S', 'si', 'Si', 'SI'):
	f12=open('stars_LN.out','a')
	f13=open('stars_correc_LN.out','a')

if ok in ('n', 'N', 'no', 'No', 'NO'):
	
	print 'escribiendo archivos'
	f122=open(registro,'w')

	f8=open(salidafile,'w')


	print 'registro'

	f122.write('# nobj id clase(1 shell) z \n')
	f122.write('#1 nombre imagen \n')
	f122.write('#2 seeing in arcsec \n')
	f122.write('#3 n stars \n')
	f122.write('#4 n gx \n')
	f122.write('#5 cantidad de galaxias background seleccionadas\n')
	
	print 'gx_back'
	

	f8.write('#1   idfainti \n')
	f8.write('#2   idfaints \n')
	f8.write('#3   xfaint \n')
	f8.write('#4   yfaint \n')
	f8.write('#5   alfafain \n')
	f8.write('#6   deltafai \n')
	f8.write('#7   r \n')
	f8.write('#8   fi \n')
	f8.write('#9   faintmag \n')
	f8.write('#10  faintcol \n')
	f8.write('#11  etfaint \n')
	f8.write('#12  exfaint \n')
	f8.write('#13  sigmae \n')
	f8.write('#14  thetafai \n')
	f8.write('#15  e1faint \n')
	f8.write('#16  e2faint \n')
	f8.write('#17  peso \n')
	f8.write('#18  r_kpc \n')
	f8.write('#19  alfa \n')
	f8.write('#20  abfaint \n')
	f8.write('#21  faintfwh \n')
	f8.write('#22  iden \n')
	f8.write('#23  dife \n')
	f8.write('#24  z_nobj \n')
	f8.write('#25  ALFA0 \n')
	f8.write('#26  DELTA0 \n')
	f8.write('#27  beta \n')
	f8.write('#28  radio \n')
	f8.write('#29  seeing \n')
	f8.write('#30  MAGMIN \n')
	f8.write('#31  zback \n')
	f8.write('#32  xy0[0] \n')
	f8.write('#33  xy0[1] \n')
	f8.write('#34  xy1[0] \n')
	f8.write('#35  xy1[1] \n')
	f8.write('#36  xy2[0] \n')
	f8.write('#37  xy2[1] \n')
TOTAL_nobjes=0
TOTAL_gxback=0

def radianes(x):
	y=(x*np.pi)/180.
	return y
	
def grados(x):
	y=(x*180.)/np.pi
	return y





def analisis(entrada):
	print 'entrada para la funcion analisis'
	print entrada
	print '-------------------------------'
	

	iden=entrada[0]
	beta=entrada[1]
	z_nobj=entrada[2]
	ALFA0=entrada[3]
	DELTA0=entrada[4]
	clase=entrada[5]
	
	corrida=int(entrada[6])
	
	seeing=entrada[7]
	
	extincion=entrada[8]
	
	MAGMIN=entrada[9]
	MAGMAX=entrada[10]
	
	
	xy0=[entrada[11],entrada[12]]

	zback=entrada[13]
	
	#corrida=int(corrida)

	#LEE EL ARCHIVO CON LAS IMAGENES QUE VA A USAR
	
	image=path+'im_06_'+str(int(iden))+'.fits'
	
	TOTAL_gxback=0
	
	D_ang=cd.angular_diameter_distance(z_nobj, z0=0, **cosmo)
	kpcscale=D_ang*(((1.0/3600.0)*np.pi)/180.0)*1000.0
	

	sexfile='secondr'+str(corrida)+'.sex' #name of configuration file of SExtractor


	salida_sex='run2-r'+str(corrida)+'.cat' #exit catalogue from SExtractor
	fwhm=seeing/pix_size #en pix
	


	#run sextrator
	#print '----------------------------------------------------------------'
	#print '                      RUNNING SExtractor                                     '
	#print '================================================================'

	imagen='sex '+image+' -c '+sexfile+' > sex_output'

		
	print imagen
	print '--------------------'	
	os.system(imagen)

	
	

	#some nobjameters
    
	#print '----------------------------------------------------------------'
	#print '                    START THE SELECTION                         '
	#print '================================================================'
    

	fwhm=seeing/pix_size #en pix
    
	stars, gx = star_gx(salida_sex,plot,fwhm)
	
	idenstars=stars[:,0]
	starsmag=stars[:,11]
	starsmu=stars[:,23]
	starsfw=stars[:,20]
	
	gxmag=gx[:,11]
	gxmu=gx[:,23]
	gxfw=gx[:,20]
	alfa=gx[:,3]
	delta=gx[:,4]
	posx=gx[:,1]
	posy=gx[:,2]
	galaxias=len(gx)
	estrellas=len(stars)
    
	ngxcat=galaxias
	
	
	#print '------------------------------------------------------------------------'
	#print '                     SELECTING BACKGROUND GALAXIES                      '
	#print '========================================================================'
    
    


    
    
    
	MAGBESTr=gxmag
   
	
		
	#BACKGROUND GALAXY SELECTION
    
	#print 'selecting galaxies between mag r ',MAGMIN,' and ',MAGMAX
    

	mask=(MAGBESTr<MAGMAX)*(MAGBESTr>MAGMIN)*(gxfw>4.)
	faintgx=gx[mask,:]
	faintmag=MAGBESTr[mask]
	faintcolor=MAGBESTr[mask]
	
	faintfwhm=gxfw[mask]
	alfafaint=alfa[mask]
	deltafaint=delta[mask]
			
	nf=faintgx.shape[0]		
			
	
	TOTAL_gxback=TOTAL_gxback+nf
		
	if nf > 0:
		
				#print '----------------------------------------------------------------'
		#print '                    ESCRIBIENDO ARCHIVOS                        '
		#print '================================================================'
		os.system('rm stars'+str(corrida)+'.cat')
		f1=open('stars'+str(corrida)+'.cat','w')
		
		f1.write('26 ')
		f1.write(str(len(stars)))
		f1.write('\n')
		np.savetxt(f1, stars, fmt='%6i   %8.3f   %7.3f   %10.7f	%10.7f	%7.3f	%7.3f	%5.1f	 %7.4f	 %7.4f  %5.1f   %8.4f	 %8.4f	 %8.4f	%8.4f	%8.4f	%8.4f	%8.4f	%8.4f 	%8i	%8.2f	%8.3f	%10.4f	  %8.4f	  %5.2f %2i')
		f1.close()
			
		os.system('rm gx'+str(corrida)+'.cat')		
		f2=open('gx'+str(corrida)+'.cat','w')
		f2.write('26 ')
		f2.write(str(nf))
		f2.write('\n')	
		np.savetxt(f2, faintgx, fmt='%6i   %8.3f   %7.3f   %10.7f	%10.7f	%7.3f	%7.3f	%5.1f	 %7.4f	 %7.4f  %5.1f   %8.4f	 %8.4f	 %8.4f	%8.4f	%8.4f	%8.4f	%8.4f	%8.4f 	%8i	%8.2f	%8.3f	%10.4f	  %8.4f	  %5.2f %2i')
		f2.close()
		
		
		
		   
		#Corro im2shape en el catalogo de estrellas
		print 'CORRIDA DE im2shape... estrellas '
		
		
		x_star=stars[:,1]
		y_star=stars[:,2]
		fitshape=(16,16)
		allstars2=fit_gaussian_parche(idenstars,image,x_star,y_star,fitshape)
		
		allstars,nstars,ncol=im2shape_run(image,'stars'+str(corrida)+'.cat',
							 'stars'+str(corrida)+'.out','example_psf.dat',corrida,
							 '4','1000')	

		
		import ipdb; ipdb.set_trace()
		
		
		#print '----------------------------------------------------------------'
		#print '             SELECCION DE ESTRELLAS PARA PSF_FILE               '
		#print '================================================================'
	
		
		#arranges where we are going to save information
		
		goodstars=selectstars(allstars,plot)		
		
		


		
		if PSF in ('s', 'S', 'si', 'Si', 'SI'):
			#print '------------------------------------------------------------------------'
			#print '                        CHECK THE PSF CORRECTION                        '
			#print '========================================================================'
			
			print 'CORRIDA DE im2shape... estrellas  con PSF'
			
			psfstarfile='psfstars'+str(corrida)+'.dat'
			psffile_stars=psffile(goodstars,stars,5,psfstarfile)
			allstars2,nstars,ncol=im2shape_run(image,'stars'+str(corrida)+'.cat',
								 'stars_correc'+str(corrida)+'.out',psfstarfile,corrida,
								 '4','1000')	
	    #~ 
		        
		#print '------------------------------------------------------------------------'
		#print '               COMPUTING THE PSF FOR THE GALAXIES                       '
		#print '========================================================================'
	
		psffile_gx=psffile(goodstars,faintgx,5,'psfgx'+str(corrida)+'.dat')
		
		#print '----------------------------------------------------------------'
		#print '             RUNNING im2shape ON THE GALAXIES                   '
		#print '================================================================'
							 
		print 'CORRIDA DE im2shape en gx... OBJ ',iden
		gxim2,ngxcat,ncolcat=im2shape_run(image,'gx'+str(corrida)+'.cat',
							 'gx'+str(corrida)+'.out',psffile_gx,corrida,
							 '4','1000')	

		
		
				#-----------------------EXTRACT AND COMPUTE SOME PARAMETERS--------------------------------
		
		# corte en dife
	

		gxcat = faintgx

		nf=gxcat.shape[0]	
		#-----------------------EXTRACT AND COMPUTE SOME PARAMETERS--------------------------------
		
		
		if nf == 1:
					#nobjameters from sextractor's catalogues
			idsex=np.array([gxcat[0]])
			x=np.array([gxcat[1]])
			y=np.array([gxcat[2]])
			alfa=np.array([gxcat[3]])
			delta=np.array([gxcat[4]])
			MAG1=np.array([gxcat[11]])
			fwhm=np.array([gxcat[20]])
			
			MUMAX=np.array([gxcat[23]])
				
			#nobjameters from im2shape catalogue
			idsex2=np.array([gxim2[1]])
			idim2=np.array([gxim2[0]])
			posx= np.array([gxim2[2]+gxim2[6]])
			posy= np.array([gxim2[3]+gxim2[7]])
			e1=np.array([gxim2[18]])
			e2=np.array([gxim2[19]])
			erre1=np.array([gxim2[34]])
			erre2=np.array([gxim2[35]])
			ab=np.array([gxim2[10]])
		else:
			#nobjameters from sextractor's catalogues
			idsex=gxcat[:,0]
			x=gxcat[:,1]
			y=gxcat[:,2]
			alfa=gxcat[:,3]
			delta=gxcat[:,4]
			MAG1=gxcat[:,11]
			fwhm=gxcat[:,20]
			
			MUMAX=gxcat[:,23]
				
			#nobjameters from im2shape catalogue
			idsex2=gxim2[:,1]
			idim2=gxim2[:,0]
			posx= gxim2[:,2]+gxim2[:,6]
			posy= gxim2[:,3]+gxim2[:,7]
			e1=gxim2[:,18]
			e2=gxim2[:,19]
			erre1=gxim2[:,34]
			erre2=gxim2[:,35]
			ab=gxim2[:,10]
		
		e=(e1**2+e2**2)**0.5
		theta=(np.arctan2(e2,e1)/2.0)+(np.pi/2.0)
	
		et,ex=rot(posx,posy,e1,e2,xy0)
			
			
			
		#----------------------------MAKE THE CATALOGUE FOR BACKGROUND GALAXIES---------------------------------
		columnas=33
		faint=np.zeros((nf,columnas),float)
		theta=np.arctan2(e2,e1)/2.0+np.pi
		idfaintsex=gxcat[:,0]
		alfafaint=gxcat[:,3]
		deltafaint=gxcat[:,4]
		r=np.zeros(nf,float)
		fi=np.zeros(nf,float)
		mufaint=gxcat[:,23]
		peso=np.zeros(nf,float)


			
		xfaint=posx
		yfaint=posy
		etfaint=et
		exfaint=ex
		e1faint=e1
		e2faint=e2
		abfaint=ab
		thetafaint=theta
		sigmae=(erre1**2.0+erre2**2.0)**0.5
		idfaintim2=idim2
			
		r=grados(np.arccos(np.sin(radianes(delta))*np.sin(radianes(DELTA0))+np.cos(radianes(delta))*np.cos(radianes(DELTA0))*np.cos(radianes(alfa-ALFA0))))
		peso=1.0
		rsec=r*3600.0
		r_kpc=rsec*kpcscale
    
		# Asign a weigh to each galaxy given the probability that it's a background galaxy
    

						
		
		#print '------------------------------------------------------------------------'
		#print '          MAKING THE CATALOGUE FOR ALL BACKGROUND GALAXIES              '
		#print '========================================================================'
		
		print 'checkkkkk'
		print faint.shape
		print nf
		print idfaintim2.shape
		print idfaintsex.shape
		faint[:,0]= idfaintim2
		faint[:,1]= idfaintsex
		faint[:,2]= xfaint
		faint[:,3]= yfaint
		faint[:,4]= alfafaint
		faint[:,5]= deltafaint
		faint[:,6]= r
		faint[:,7]= fi
		faint[:,8]= faintmag
		faint[:,9]= faintcolor
		faint[:,10]=etfaint
		faint[:,11]=exfaint
		faint[:,12]=sigmae
		faint[:,13]=thetafaint
		faint[:,14]=e1faint
		faint[:,15]=e2faint
		faint[:,16]=peso
		faint[:,17]=r_kpc
		faint[:,18]=alfa
		faint[:,19]=abfaint
		faint[:,20]=faintfwhm
		faint[:,21]=iden
		faint[:,22]=sigmae
		faint[:,23]=z_nobj
		faint[:,24]=ALFA0
		faint[:,25]=DELTA0
		faint[:,26]=beta
		faint[:,27]=beta
		faint[:,28]=seeing
		faint[:,29]=MAGMIN
		faint[:,30]=zback
		faint[:,31]=xy0[0]
		faint[:,32]=xy0[1]
		
		os.system('rm gx_back_1_'+str(corrida)+'.cat')
		f1=open('gx_back_1_'+str(corrida)+'.cat','w')
		np.savetxt(f1, faint, fmt=['%6i']*2+['%12.5f']*19+['%6i']+['%11.5f']*11)
		f1.close()	
		
		
			
	print 'accaaaaaaaaaaaaaaaaaaaaaaaaaa',corrida,galaxias,estrellas,nf

	salida=[iden,corrida,faint,allstars,allstars2,nf]
	print salida
	return salida




N_TOT=62


while nobj <= N_TOT:
	entrada=np.zeros((pulls,14),float)
	contar=0
	while contar < pulls:
		iden = nobj
		f122.write(str(iden)+'  ')
		f122.write(str(ALFA0)+'  ')
		f122.write(str(DELTA0)+'  ')

		print 'VOY POR LA IMAGEN ',iden
		
		
		print '----------------------------------------------------------'
		print '                    OBTAINING IMAGES                      '
		print '=========================================================='
		
		print 'ENTRADA nobjA LA FUNC SEEING'
		
		IMAGE=path+'im_06_'+str(int(nobj))+'.fits'
		SEEING=seeing(IMAGE,0.27,contar,'r',15.,12.,5.,plot)
		xy0=np.array([1052.,2048.])
		f122.write(str(SEEING)+'  ')

		print 'salida best...',IMAGE,SEEING		

		f122.write('\n')			
		if SEEING < 1.3:

			print 'nobj numero',iden

			
			
			#--------------------------------CHANGE THE COORDINATE AXIS-------------------------------
			#Now is going to change the cordinates axis of e1 and e2 to get et ex
			
			beta=0.53
			zback=0.8
			if MAGMIN in ('completar','calcular','compute'):
				MAGMIN_obj=salida_beta[2]
			else:
				MAGMIN_obj=MAGMIN
				
			entrada[contar,0]=iden
			entrada[contar,1]=beta
			entrada[contar,2]=z_nobj
			entrada[contar,3]=ALFA0
			entrada[contar,4]=DELTA0
			entrada[contar,5]=clase
			entrada[contar,6]=contar
			entrada[contar,7]=SEEING
			entrada[contar,8]=extincion
			entrada[contar,9]=MAGMIN_obj
			entrada[contar,10]=MAGMAX
			entrada[contar,11]=xy0[0]
			entrada[contar,12]=xy0[1]
			entrada[contar,13]=zback


			contar=contar+1
	
		if nobj == N_TOT:
			contar=pulls	
		nobj=nobj+1
		#wait=raw_input('chequear imagen, presionar ENTER')
		
		
	if contar > 0:
		
		
		entrada=entrada[:contar,:]
		
		print 'entrada',entrada
		#~ wait=raw_input('chequear imag	en, presionar ENTER')
		f122.write('ultima '+IMAGE+'\n')
		if pulls>1:
			pool = Pool(processes=(contar))              # start 4 worker processes
			salida=pool.map(analisis, entrada)
			pool.terminate()
		else:
			entrada=entrada[0,:]
			salida=analisis(entrada)
		print 'corrida',salida[0]
	
		for j in range(contar):
			if pulls > 1:
				corrida=salida[j]
			else:
				corrida=salida
			idim=corrida[0]
			f8.write('#'+str(int(corrida[0]))+'\n')
			print 'nobja la corrrida ',j,' obtuvo ',corrida[5],'galaxias'
			
			if corrida[5] > 0:
				faint1= corrida[2]
				if corrida[5] == 1:
					faint1=np.array([faint1])
					np.savetxt(f8, faint1, fmt=['%6i']*2+['%12.5f']*19+['%6i']+['%11.5f']*11)
				else:                                                                       
					np.savetxt(f8, faint1, fmt=['%6i']*2+['%12.5f']*19+['%6i']+['%11.5f']*11)
			if PSF in ('s', 'S', 'si', 'Si', 'SI'):

				idarray=np.zeros(corrida[3].shape[0])
				idarray.fill(idim)

				stars0=np.column_stack((corrida[3],idarray))
				stars1=np.column_stack((corrida[4],idarray))

				np.savetxt(f12, stars0, fmt=['%6i']*4+['%12.5f']*39+['%6i']*2)
				np.savetxt(f13, stars1, fmt=['%6i']*4+['%12.5f']*39+['%6i']*2)
			#wait=raw_input('termino una vez')
	
	


f7.close()
f8.close()

print 'END OF THE PROGRAM :)'
