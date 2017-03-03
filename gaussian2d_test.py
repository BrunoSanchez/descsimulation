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

    fitter = fitting.LevMarLSQFitter()

    y2, x2 = np.mgrid[:30,:30]
    p = models.Gaussian2D(x_mean=15.,
                          y_mean=15.,
                          x_stddev=1.,
                          y_stddev=1.,
                          theta=np.pi/4.,
                          amplitude=B.max()-B.min())+models.Const2D(amplitude=B.min())

    out = fitter(p, x2, y2, B, maxiter=100000)
    return out


thetain=85.
print 'theta in',thetain


for ellip in np.arange(0.01,0.9,0.1):

    B=GAUSS(ellip, x0=15.,y0=15.,th=thetain,noise=0.)
    out=FIT(B)[0]
    y2,x2=np.mgrid[:30,:30]

    e = (out.x_stddev - out.y_stddev)/(out.x_stddev + out.y_stddev)
    th = np.rad2deg(out.theta)
    print 'Input parameters:'
    print 'theta={}, e={}'.format(thetain, ellip)
    print 'Output parameters:'
    print 'theta={}, e={}'.format(th, e)
