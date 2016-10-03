#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  test_psf_DES.py
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
import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import psf_DES

seeing_par = {}
seeing_par['pixsize'] = 0.27
seeing_par['run'] = 33
seeing_par['filter'] = 'R'
seeing_par['magmax'] = 23.
seeing_par['magmin'] = 18.
seeing_par['fwhmmax'] = 4.5
seeing_par['plot'] = None


def main(image, path, plot=None):
    seeing_par['img'] = image
    if plot is not None:
        seeing_par['plot'] = plot

    im = psf_DES.SingleImageDES(seeing_par=seeing_par, img=image)
    test_dir = os.path.abspath(path)
    if not os.path.exists(test_dir):
        os.makedirs(test_dir)

    a_fields, psf_basis = im.get_variable_psf(delete_patches=False,
                                              pow_th=0.999)


    plt.figure()# figsize=(16, 16))
    plt.imshow(np.log10(fits.getdata(image)), interpolation='none')
    plt.plot(im._best_srcs['sources']['X_IMAGE'],
             im._best_srcs['sources']['Y_IMAGE'],
             'ro')
    plt.colorbar(orientation='horizontal')
    plt.savefig(os.path.join(test_dir, 'test_frame.png'))
    plt.close()

    #  Pathces are in im._best_srcs['patches']
    subplots = int(np.sqrt(len(im._best_srcs['patches'])) + 1)
    plt.figure(figsize=(20, 20))
    for i in range(len(im._best_srcs['patches'])):
        plt.subplot(subplots, subplots, i+1)
        plt.imshow(im._best_srcs['patches'][i], interpolation='none')
        # plt.colorbar(orientation='horizontal')
    plt.savefig(os.path.join(test_dir, 'psf_patches.png'))
    plt.close()

    subplots = int(np.sqrt(len(psf_basis)) + 1)
    plt.figure(figsize=(16, 16))
    for i in range(len(psf_basis)):
        plt.subplot(subplots, subplots, i+1)
        plt.imshow(psf_basis[i], interpolation='none')
        plt.colorbar(orientation='horizontal')
    plt.savefig(os.path.join(test_dir, 'psf_basis.png'))
    plt.close()

    x, y = np.mgrid[:im.imagedata.shape[0], :im.imagedata.shape[1]]
    plt.figure(figsize=(16, 16))
    for i in range(len(a_fields)):
        plt.subplot(subplots, subplots, i+1)
        plt.imshow(a_fields[i](x, y))
        plt.plot(im._best_srcs['sources']['X_IMAGE'],
                 im._best_srcs['sources']['Y_IMAGE'],
                 'ro')
        plt.colorbar(orientation='horizontal')
    plt.savefig(os.path.join(test_dir, 'a_fields.png'))
    plt.close()
    return a_fields, psf_basis, im

if __name__ == '__main__':
    import sys
    if len(sys.argv) == 4:
        sys.exit(main(sys.argv[1], sys.argv[2], sys.argv[3]))
    else:
        sys.exit(main(sys.argv[1], sys.argv[2]))

