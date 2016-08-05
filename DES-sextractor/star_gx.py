import os

import numpy as np
import matplotlib.pyplot as plt

from astropy.io import ascii


def star_gx(salida_sex, fwhm, plot=None):
    #~ print '######## CLASIFICANDO LOS OBJETOs #######'
    #reads the catalogue of the first run of sextractor
    cat = ascii.read(salida_sex, format='sextractor')

    #extract some nobjameters from the catalogue
    FWHM = cat['FWHM_IMAGE']
    MAGBEST = cat['MAG_BEST']
    MUMAX = cat['MU_MAX']
    CLASS=cat['CLASS_STAR']
    FLAG=cat['FLAGS']
    nobj=len(cat)
    ALFA=cat['ALPHA_J2000']
    DELTA=cat['DELTA_J2000']
    control2=np.zeros(4,float)
    comnobjar=0.0

    # STARS CANDIDATES
    hi_class = cat['CLASS_STAR'] > 0.85
    hifwhm = cat['FWHM_IMAGE'] - fwhm > -0.5
    lofwhm = cat['FWHM_IMAGE'] - fwhm < 0.8

    candidatas = cat[hi_class & hifwhm & lofwhm

    if plot is not None:
        print 'Total candidatas = {}'.format(len(candidatas))#'m,n',m, n
        print 'Valor de seeing = {}'.format(fwhm)
        plt.plot(MAGBEST, MUMAX, 'k.')
        plt.xlabel('MAG_BEST')
        plt.ylabel('MU_MAX')
        plt.show()

        plt.plot(MAGBEST, CLASS, 'k.')
        plt.xlabel('MAG_BEST')
        plt.ylabel('MU_MAX')
        plt.show()

        plt.plot(FWHM, MAGBEST, 'k.')
        plt.ylabel('MAGBEST')
        plt.xlabel('FWHM')
        plt.show()

    X = candidatas['X_IMAGE']
    Y = candidatas['Y_IMAGE']
    nrot = len(candidatas)
    X2 = X**2
    Y2 = Y**2
    XY = X * Y
    Sx = X.sum()
    Sy = Y.sum()
    Sxy = XY.sum()
    Sxx = X2.sum()
    Syy = Y2.sum()
    m = (nrot*Sxy - Sx*Sy) / (nrot*Sxx - Sx*Sx)
    n = (Sxx*Sy - Sx*Sxy) / (nrot*Sxx - Sx*Sx)
    muaj = n + 20.0

    varx = np.array([X.min(), X.max()])
    vary = varx * m + n

    # x y of the line to do the selection
    #x=range(int(MAGBEST.max()-MAGBEST.min()))+MAGBEST.min()
    zero = MAGBEST.min()
    x = range(15) + zero
    y = m * x + (muaj - 20.)

    anchomag = 0.4  #width in magnitudes
    ancho = anchomag

    #print 'ancho',anchomag, ancho
    mumin = 7.8#float(raw_input('Ingrese mu minimo '))
    mumax = 20.5
    stars = np.zeros((nobj,26),float)
    gx = np.zeros((nobj,26),float)
    g = 0
    i = 0

    mu = m * cat['MAG_BEST'] + (muaj - 20.)
    hi_mu = cat['MU_MAX'] > mumin
    hi_fwhm = cat['FWHM_IMAGE'] > fwhm - 0.5
    lo_flag = cat['FLAGS'] < 4.

    goodies = cat[hi_mu & hi_fwhm & lo_flag]

    upp_mu = goodies['MU_MAX'] < mu+ancho
    low_mu = goodies['MU_MAX'] > mu-ancho
    under_mu = goodies['MU_MAX'] < mumax
    small_fwhm = goodies['FWHM_IMAGE'] < fwhm+0.8

    stars = goodies[upp_mu & low_mu & under_mu & small_fwhm]

    pre_gx = goodies[!(upp_mu & low_mu & under_mu & small_fwhm)]
    gx = pre_gx[pre_gx['CLASS_STAR'] > 0.8 | pre_gx['MU_MAX'] > mu-ancho]

    for j in range(nobj):
        mu=m*MAGBEST[j]+(muaj-20.)
        if MUMAX[j] > mumin and FWHM[j] > fwhm-0.5 and FLAG[j] < 4.0:
            if MUMAX[j] < mu+ancho and MUMAX[j] < mumax and FWHM[j] < fwhm+0.8 \
                and MUMAX[j] > mu-ancho:
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

    starsmag = stars[:,11]
    starsmu=stars[:,23]
    starsfw=stars[:,20]

    gxmag = gx[:, 11]
    gxmu = gx[:, 23]
    gxfw = gx[:, 20]

    if plot in ('s', 'S', 'si', 'Si', 'SI'):
        print 'm,n',m, n
        plt.plot(MAGBEST,MUMAX, 'k.',
                 starsmag, starsmu, 'b.',
                 gxmag, gxmu, 'g.',
                 x, y, 'r',
                 varx, vary, 'b')

        plt.xlabel('MAG_BEST')
        plt.ylabel('MU_MAX')
        plt.show()

        plt.plot(FWHM,MAGBEST, 'k.',
                 starsfw, starsmag, 'b*',
                 gxfw, gxmag, 'g.')

        plt.ylabel('MAGBEST')
        plt.xlabel('FWHM')
        plt.show()


    return stars, gx
