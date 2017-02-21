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
from properimage import propercoadd as pc
from properimage import numpydb as npdb
import seeing_des
import star_gx


class SingleImageDES(pc.SingleImage):
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
    def _best_srcs(self):
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
            p_sizes = np.sqrt(np.percentile(stars['ISOAREA_IMAGE'],
                                            q=[35, 55, 75]))
            if not p_sizes[1] < 11:
                dx = int(p_sizes[1])
                if dx % 2 == 0: dx += 1
                fitshape = (dx, dx)
            else:
                fitshape = (11, 11)

            print 'Sources good to calculate = {}'.format(len(stars))
            self._best_sources = {'sources': stars, 'fitshape': fitshape}

            self.db = npdb.NumPyDB_cPickle(self.dbname, mode='store')

            pos = []
            jj = 0
            for row in stars:
                position = [row['Y_IMAGE'], row['X_IMAGE']]
                sub_array_data = extract_array(self.bkg_sub_img,
                                               fitshape, position,
                                               fill_value=self.bkg.globalrms)
                sub_array_data = sub_array_data/np.sum(sub_array_data)

                # Patch.append(sub_array_data)
                self.db.dump(sub_array_data, jj)
                pos.append(position)
                jj += 1

            # self._best_sources['patches'] = np.array(Patch)
            self._best_sources['positions'] = np.array(pos)
            self._best_sources['n_sources'] = jj
            # self._best_sources['detected'] = srcs
            # self.db = npdb.NumPyDB_cPickle(self._dbname, mode='store')

            print 'returning best sources'

            self._gxs = gxs
        return self._best_sources

