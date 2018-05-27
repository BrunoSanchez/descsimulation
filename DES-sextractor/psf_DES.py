#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  psf_DES.py
#
#  Copyright 2016 Bruno S <bruno@oac.unc.edu.ar>
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
from astropy.nddata.utils import extract_array
from properimage import single_image as si
from properimage import numpydb as npdb
import seeing_des
import star_gx


class SingleImageDES(si.SingleImage):
    """Class to analyze a single image using propercoadd.SingleImage class.
    It improves the star selection for psf determination.
    It makes use of sextractor and seeing_des, star_gx modules.

    Parameters
    ----------
    seeing_par: dictionary
        a dictionary containing the parameters for seeing_des.seeing

    Returns
    -------
    A SingleImage class instance with best sources determination modified for
    sextractor use.

    Example
    -------
    >>>seeing_par = {}
    >>>seeing_par['img'] = '../imagenes-LN-PSF/im_06_1.fits[1]'
    >>>seeing_par['pixsize'] = 0.27
    >>>seeing_par['run'] = 33
    >>>seeing_par['filter'] = 'R'
    >>>seeing_par['magmax'] = 23.
    >>>seeing_par['magmin'] = 18.
    >>>seeing_par['fwhmmax'] = 4.5
    >>>seeing_par['plot'] = 'SI'
    >>>
    >>> im = SingleImageDES(img, imagefile=True, seeing_par=seeing_par)

    """

    def __init__(self, seeing_par=None, *args, **kwargs):
        self.seeing_par = seeing_par
        super(SingleImageDES, self).__init__(*args, **kwargs)


    @property
    def best_sources(self):
        """Property, a table of best sources detected in the image.

        """
        if not hasattr(self, '_best_sources'):
            par = self.seeing_par
            self.seeing, self.salida_sex = seeing_des.seeing(par['img'],
                                                             par['pixsize'],
                                                             par['run'],
                                                             par['filter'],
                                                             par['magmax'],
                                                             par['magmin'],
                                                             par['fwhmmax'],
                                                             par['plot'])

            stars, gxs = star_gx.star_gx(self.salida_sex,
                                         self.seeing/par['pixsize'],
                                         par['plot'])

            stars.sort('MAG_BEST')

            if len(stars) > 1800:
                jj = np.random.choice(len(stars), 1800, replace=False)
                stars = stars[jj]

            print 'Sources good to calculate = {}'.format(len(stars))
            stars['x'] = stars['X_IMAGE']-1.
            stars['y'] = stars['Y_IMAGE']-1.

            self._best_sources = stars
            self._gxs = gxs
            print 'returning best sources'
        return self._best_sources

    @property
    def gxs(self):
        try:
            return self._gxs
        except AttributeError:
            s = self.best_sources
            return self._gxs

    @property
    def stamp_shape(self):
        return self.__stamp_shape

    @stamp_shape.setter
    def stamp_shape(self, shape):
        if not hasattr(self, '__stamp_shape'):
            if shape is None:
                #~ percent = np.percentile(self.best_sources['ISOAREA_IMAGE'], q=65)
                #~ p_sizes = 2.*np.sqrt(percent)
                percent = np.percentile(self.best_sources['FWHM_IMAGE'], q=65)
                p_sizes = 3.*percent

                if p_sizes >= 9:
                    dx = int(p_sizes)
                    if dx % 2 != 1: dx += 1
                    shape = (dx, dx)
                else:
                    shape = (9, 9)
                print(('stamps will be {} x {}'.format(*shape)))
        self.__stamp_shape = shape

