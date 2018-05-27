#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  gaussian2d_test.py
#
#  Copyright 2017 Bruno S <bruno@oac.unc.edu.ar>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#
import numpy as np
from astropy.modeling import fitting, models
import random
import os

def GAUSS(e, x0, y0, th, noise):
    escx = 30
    escy = 30
    B = np.zeros((escx, escy), float)

    x = np.zeros(2,float)
    xg = np.array([x0, y0])
    C = np.zeros((2,2),float)

    ab = 10.0
    Amp = 20.0

    theta = th*np.pi/180.0
    a = (ab*((1.0 + e)/(1.0 - e)))**0.5
    b = (ab*((1.0-e)/(1.0+e)))**0.5
    C[0,0] = 2 * ( (np.cos(theta)**2)/a**2 + (np.sin(theta)**2)/b**2)
    C[1,1] = 2 * ( (np.cos(theta)**2)/b**2 + (np.sin(theta)**2)/a**2)
    C[0,1] = (1/(b**2) - 1/(a**2))*np.sin(2*theta)
    C[1,0] = (1/(b**2) - 1/(a**2))*np.sin(2*theta)
    A = Amp * 2 * np.pi * np.linalg.det(C)

    for j in range(escy):
        for i in range(escx):
            x[0] = i
            x[1] = j
            resta = x - xg
            B[j,i]=np.e**(-0.5*(np.dot(np.dot(resta,C),resta)))

    B = B + noise*np.random.random(B.shape) + 30.0
    return B * Amp


def FIT(B):

    fitter=fitting.LevMarLSQFitter()

    y2, x2 = np.mgrid[:30, :30]
    ampl = B.max()-B.min())+models.Const2D(amplitude=B.min())
    p = models.Gaussian2D(x_mean=15., y_mean=15.,
                          x_stddev=1., y_stddev=1.,
                          theta=np.pi/4.,
                          amplitude=ampl)
    out = fitter(p,x2,y2,B,maxiter=5000)
    return out

thetain=85.
print 'theta in',thetain


for ellip in np.arange(0.01,0.9,0.1):

    x=np.zeros(2,float)

    C=np.zeros((2,2),float)

    e=ellip
    ab=4.0
    Amp=2.00
    xg=np.array([x0,y0])
    theta=np.pi-(th*np.pi)/180.0
    a=(ab*((1.0+e)/(1.0-e)))**0.5
    b=(ab*((1.0-e)/(1.0+e)))**0.5
    C[0,0]=2*(((np.cos(theta))**2)/a**2+((np.sin(theta))**2)/b**2)
    C[1,1]=2*((((np.cos(theta))**2)/b**2)+(((np.sin(theta))**2)/a**2))
    C[0,1]=((1/(b**2))-(1/(a**2)))*np.sin(2*theta)
    C[1,0]=((1/(b**2))-(1/(a**2)))*np.sin(2*theta)
    A=Amp*(2*np.pi*np.linalg.det(C))

    for j in range(escy):
        for i in range(escx):
            x[0]=i
            x[1]=j
            resta=x-xg
            B[j,i]=(A/(2*np.pi*np.linalg.det(C)))*(np.e**(-0.5*(np.dot(np.dot(resta,C),resta))))+B[j,i]+noise*np.random.uniform(-1.,1.)+0.14
    return B


def iterar(eup,elow,thetain,noise):
    for ellip in np.arange(eup,elow,0.02):

        B=GAUSS(ellip,x0=15.,y0=15.,th=thetain,noise=noise)
        os.system('rm *fits')
        hdu = pyfits.PrimaryHDU(B)
        hdu.writeto('prumar.fits')
        out=FIT(B)[0]
        out2=FIT(B)[1]
        y2,x2=np.mgrid[:30,:30]

        elip = (out.x_stddev - out.y_stddev)/(out.x_stddev + out.y_stddev)
        th = np.rad2deg(out.theta)
        #~ print 'theta out ',th- int(th/180)*180.
        print 'elip out = ',ellip, abs(elip)
        if out.x_stddev > out.y_stddev:
            theta=(th- int(th/180)*180.)
        else:
            theta=((th- int(th/180)*180.)-90)
            #~ print 'CAMBIOOOOOOOOO-------------------'
        if theta > 0:
            theta=180.-theta
        else:
            theta=abs(theta)
        if theta >180.:
            theta=theta-180.
        theta=180.-theta
        print 'THETA',theta
        print out2
        print out.amplitude
    return B, out
        #~ print out.x_stddev,out.y_stddev

def read_im2(file_out):
    im2out=np.loadtxt(file_out,skiprows=1)
    e1=im2out[18]
    e2=im2out[19]
    e=(e1**2+e2**2)**0.5
    theta=np.rad2deg((np.arctan2(e2,e1)/2.0)+(np.pi/2.0))
    amp=im2out[11]
    back=im2out[5]
    noise=im2out[4]

    print 'e',e
    print 'theta',theta
    print 'amp',amp
    print 'back',back
    print 'noise',noise

    e = (out.x_stddev - out.y_stddev)/(out.x_stddev + out.y_stddev)
    th = np.rad2deg(out.theta)
    print 'Input parameters:'
    print 'theta={}, e={}'.format(thetain, ellip)
    print 'Output parameters:'
    print 'theta={}, e={}'.format(th, e)
