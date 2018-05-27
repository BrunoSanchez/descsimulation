import os
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
import sys
from astropy.visualization import (MinMaxInterval, ZScaleInterval, LogStretch,
                                   SqrtStretch, LinearStretch, ImageNormalize)

from properimage import single_image as si
from properimage import utils


frame = '/home/bruno/Devel/DESCSimulation/imagenes-LN-PSF/im_06_13.fits'

im = si.SingleImage(frame)

# =============================================================================
#    PSF spatially variant
# =============================================================================
a_fields, psf_basis = im.get_variable_psf(inf_loss=0.005)

test_dir = os.path.abspath('/home/bruno/Devel/DESCSimulation/psf_tests3/')
if not os.path.exists(test_dir):
    os.mkdir(test_dir)
dt = fits.getdata(frame)

plt.figure(figsize=(12,12))
norm = ImageNormalize(dt, interval=ZScaleInterval(), stretch=LogStretch())
plt.imshow(dt, #vmin=med-3.*s, vmax=med+3.*s,
           interpolation='none', origin='lower', norm=norm, cmap='Greys_r')
plt.colorbar(orientation='vertical')
plt.savefig(os.path.join(test_dir, 'test_frame.png'), dpi=720)
plt.close()

subplots = int(np.sqrt(len(psf_basis)) + 1)
plt.figure(figsize=(12, 12))
for i in range(len(psf_basis)):
    plt.subplot(subplots, subplots, i+1);
    plt.imshow(psf_basis[i], interpolation='none')
    plt.colorbar(orientation='vertical')
plt.tight_layout()
plt.savefig(os.path.join(test_dir, 'psf_basis.png'), dpi=720)
plt.close()

if len(psf_basis) >1:
    x, y = im.get_afield_domain()
    plt.figure(figsize=(12, 12))
    for i in range(len(a_fields)):
        plt.subplot(subplots, subplots, i+1);
        plt.imshow(a_fields[i](x, y), origin='lower', interpolation='none')
        plt.colorbar(orientation='vertical')
    plt.savefig(os.path.join(test_dir, 'a_fields.png'), dpi=720)
    plt.close()

im._clean()

