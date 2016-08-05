import math
import numpy as np
from pylab import *
import os
import pyfits
import sys



def star_gx(salida_sex,plot,fwhm):
    #~ print '######## CLASIFICANDO LOS OBJETOs #######'

    cat = np.loadtxt(salida_sex, comments='#') #reads the catalogue of the first run of sextractor

    #extract some nobjameters from the catalogue
    FWHM = cat[:,20]
    MAGBEST = cat[:,11]
    MUMAX = cat[:,23]
    CLASS=cat[:,24]
    FLAG=cat[:,25]
    nobj=len(MAGBEST)
    ALFA=cat[:,3]
    DELTA=cat[:,4]
    control2=np.zeros(4,float)
    comnobjar=0.0

    # STARS CANDIDATES
    candidatas=np.zeros((nobj,26),float)
    i=0
    for j in range(nobj):
        if CLASS[j] > 0.85  and (fwhm-0.5 < FWHM[j] < fwhm+0.8):
            candidatas[i,:]=cat[j,:]
            i=i+1


    if plot in ('s', 'S', 'si', 'Si', 'SI'):
        print 'candidatas ',i#'m,n',m, n
        print 'seeing', fwhm
        plt.plot(MAGBEST,MUMAX, 'k.')
        plt.xlabel('MAG_BEST')
        plt.ylabel('MU_MAX')
        plt.show()

        plt.plot(MAGBEST,CLASS, 'k.')
        plt.xlabel('MAG_BEST')
        plt.ylabel('MU_MAX')
        plt.show()

        plt.plot(FWHM,MAGBEST, 'k.')
        plt.ylabel('MAGBEST')
        plt.xlabel('FWHM')
        plt.show()

    X=candidatas[:i,11]
    Y=candidatas[:i,23]
    nrot=i
    X2=X**2
    Y2=Y**2
    XY=X*Y
    Sx=X.sum()
    Sy=Y.sum()
    Sxy=XY.sum()
    Sxx=X2.sum()
    Syy=Y2.sum()
    m=(nrot*Sxy-Sx*Sy)/(nrot*Sxx-Sx*Sx)
    n=(Sxx*Sy-Sx*Sxy)/(nrot*Sxx-Sx*Sx)
    muaj=n+20.0
    varx=np.zeros(2,float)
    vary=np.zeros(2,float)
    varx[0]=X.min()
    varx[1]=X.max()

    for j in range(2):
        vary[j]=m*varx[j]+n

    # x y of the line to do the selection
    #x=range(int(MAGBEST.max()-MAGBEST.min()))+MAGBEST.min()
    zero=MAGBEST.min()
    x=range(15)+zero
    y=m*x+(muaj-20.)


    anchomag=0.4  #width in magnitudes
    ancho=anchomag

    #print 'ancho',anchomag, ancho

    mumin=7.8#float(raw_input('Ingrese mu minimo '))
    mumax=20.5
    stars=np.zeros((nobj,26),float)
    gx=np.zeros((nobj,26),float)
    g=0
    i=0


    for j in range(nobj):
        mu=m*MAGBEST[j]+(muaj-20.)
        if MUMAX[j] > mumin and FWHM[j] > fwhm-0.5 and FLAG[j] < 4.0:
            if MUMAX[j] < mu+ancho and MUMAX[j] < mumax and FWHM[j] < fwhm+0.8  and MUMAX[j] > mu-ancho:
                stars[i,:]=cat[j,:]
                i=i+1
            elif CLASS[j] < 0.8:
                if MUMAX[j] > mumax:
                    gx[g,:]=cat[j,:]
                    g=g+1
                elif MUMAX[j] > mu-ancho:
                    gx[g,:]=cat[j,:]
                    g=g+1


    stars=stars[:i,:]
    gx=gx[:g,:]
    #print 'fwhm', fwhm
    #print 'n obj', nobj
    print 'cantidad de estrellas seleccionadas: ', len(stars)
    print 'cantidad de galaxias seleccionadas: ', len(gx)

    starsmag=stars[:,11]
    starsmu=stars[:,23]
    starsfw=stars[:,20]

    gxmag=gx[:,11]
    gxmu=gx[:,23]
    gxfw=gx[:,20]

    if plot in ('s', 'S', 'si', 'Si', 'SI'):
        print 'm,n',m, n
        plt.plot(MAGBEST,MUMAX, 'k.',starsmag,starsmu,'b.',gxmag, gxmu, 'g.',x,y,'r',varx,vary,'b')
        plt.xlabel('MAG_BEST')
        plt.ylabel('MU_MAX')
        plt.show()

        plt.plot(FWHM,MAGBEST, 'k.',starsfw,starsmag,'b*',gxfw,gxmag,'g.')
        plt.ylabel('MAGBEST')
        plt.xlabel('FWHM')
        plt.show()


    return stars, gx
