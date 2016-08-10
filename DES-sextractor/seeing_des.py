#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import shlex
import subprocess

import jinja2
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii


def write_sexconf(dest_file, sexconf_dict):
    """Writes a sextractor configuration file using a template given in
    'sextractor-conf-template.j2', using jinja2

    Everything is parsed by a dictionary with the desired values settled.
    """
    loader = jinja2.FileSystemLoader(os.path.abspath('.'))
    jenv = jinja2.Environment(loader=loader,
                              trim_blocks=True,
                              lstrip_blocks=True)

    template = jenv.get_template('sextractor-conf-template.j2')

    if os.path.exists(sexconf_dict['salida_sex']):
        os.remove(sexconf_dict['salida_sex'])

    with open(dest_file, 'w') as f1:
        f1.write(template.render(sexconf_dict))


def seeing(imagen, pixsize, corrida, filtro, magmax,
           magmin, fwhmmax, plot=None):
    """Function to calculate the seeing level of an image, and separate the
    stars from the galaxies.

    This makes use of sextractor

    Parameters
    ----------

    Returns
    -------

    Example
    -------
    """

# =============================================================================
# WRITING SExtractor CONFIGURATION FILE
# =============================================================================
    # name of configuration file of SExtractor
    sexfile = 'first' + filtro + str(corrida) + '.sex'
    # exit catalogue from SExtractor
    salida_sex = 'run1-' + filtro + str(corrida) + '.cat'

    sexconf = {}
    sexconf['salida_sex'] = salida_sex
    sexconf['thresh'] = 5.0
    sexconf['param_name'] = 'first.param'
    sexconf['satur_lev'] = 50000.
    sexconf['zp_mag'] = 26.73
    sexconf['gain'] = 5.16
    sexconf['px_scale'] = pixsize
    sexconf['seeing_fwhm'] = 1.

    write_sexconf(sexfile, sexconf)

    callsex_cmd = shlex.split('sex '+imagen+' -c '+sexfile)
    subprocess.call(callsex_cmd)

    # print imagen
    # wait=raw_input('chequear y presionar ENTER')
    cat = ascii.read(salida_sex, format='sextractor')
    # reads the catalogue of the first run of sextractor

# =============================================================================
#    COMPUTE SATURATION LEVEL
# =============================================================================
    FLUXMAX = cat['FLUX_MAX']
    SATUR = 0.8 * FLUXMAX.max()
    # print 'SATUR_LEVEL =', SATUR

# =============================================================================
#  COMPUTE THE SEEING
# =============================================================================
    # import ipdb; ipdb.set_trace()
    FWHM = cat['FWHM_IMAGE']
    FLAG = cat['FLAGS']
    MAGBEST = cat['MAG_BEST']
    nobj = len(cat)
    fwhm = np.zeros((nobj), float)
    mag = np.zeros((nobj), float)

    if plot is not None:
        plt.plot(FWHM, MAGBEST, 'b.')
        plt.ylabel('MAG_BEST')
        plt.xlabel('FWHM')
        plt.show()

    # make the cuts in magnitude to get mostly the stars
    himag = cat['MAG_BEST'] < magmax
    lomag = cat['MAG_BEST'] > magmin
    lo_fw = cat['FWHM_IMAGE'] < fwhmmax

    fwhm = cat['FWHM_IMAGE'][himag & lomag & lo_fw]
    mag = cat['MAG_BEST'][himag & lomag & lo_fw]

    # and get the maximun value of the FWHM distribution
    step = 0.05
    intervalos = int((fwhm.max() - fwhm.min()) / step)
    interval1 = fwhm.min()
    interval2 = interval1+step
    maximo = 0
    moda = 90.
    for j in range(intervalos):
        contar = 0
        for x in fwhm:
            if interval1 < x < interval2:
                contar = contar+1
            if contar > maximo:
                maximo = contar
                moda = (interval1+interval2)/2.
        interval1 = interval1+step
        interval2 = interval2+step

    # print 'SEEING in pix',moda
    # print 'SEEING in arcsec',moda*pixsize
    seeing = moda * pixsize
    print ' '
    print ' ------------ seeing: ', seeing
    print ' ------------ imagen: ', imagen
    print ' '

    puntox = np.array([moda], float)
    puntoy = np.array([18.], float)

    if plot is not None:
        plt.plot(FWHM, MAGBEST, 'b.', puntox, puntoy, 'ro')
        plt.ylabel('MAG_BEST')
        plt.xlabel('FWHM')
        plt.axis([0, 30, 16, 26])
        plt.show()

# =============================================================================
#   WRITING SExtractor CONFIGURATION FILE
# =============================================================================
    # name of configuration file of SExtractor
    sexfile = 'second' + filtro + str(corrida) + '_pru.sex'

    sexconf = {}
    sexconf['salida_sex'] = 'run2-' + filtro + str(corrida) + '.cat'
    sexconf['thresh'] = 1.5
    sexconf['param_name'] = 'second.param'
    sexconf['satur_lev'] = int(SATUR)
    sexconf['zp_mag'] = 26.73
    sexconf['gain'] = 5.16
    sexconf['px_scale'] = pixsize
    sexconf['seeing_fwhm'] = seeing

    write_sexconf(sexfile, sexconf)
    callsex_cmd = shlex.split('sex '+imagen+' -c '+sexfile)
    subprocess.call(callsex_cmd)

    del(cat)
    del(FWHM)
    del(FLAG)
    del(MAGBEST)
    del(nobj)
    del(fwhm)
    del(mag)

    return seeing, sexconf['salida_sex']
