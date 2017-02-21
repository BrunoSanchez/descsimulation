import os
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import sys

from properimage import propercoadd as pc
from properimage import utils


frame = '/home/bruno/Devel/DESCSimulation/imagenes-LN-PSF/im_06_13.fits'

im = pc.SingleImage(frame, imagefile=True)

# =============================================================================
#    PSF spatially variant
# =============================================================================
a_fields, psf_basis = im.get_variable_psf(pow_th=0.1)

test_dir = os.path.abspath('/home/bruno/Devel/DESCSimulation/psf_tests3/')
if not os.path.exists(test_dir):
    os.mkdir(test_dir)
dt = fits.getdata(frame)

plt.figure(figsize=(16,16))
plt.imshow(np.log10(dt), interpolation='none')
plt.colorbar(orientation='horizontal')
plt.savefig(os.path.join(test_dir, 'test_frame.png'))
plt.close()

subplots = int(np.sqrt(len(psf_basis)) + 1)
plt.figure(figsize=(16, 16))
for i in range(len(psf_basis)):
    plt.subplot(subplots, subplots, i+1);
    plt.imshow(psf_basis[i], interpolation='none')
    plt.colorbar(orientation='horizontal')
plt.savefig(os.path.join(test_dir, 'psf_basis.png'))
plt.close()

if len(psf_basis) >1:
    x, y = np.mgrid[:4096, :2048]
    plt.figure(figsize=(16, 16))
    for i in range(len(a_fields)):
        plt.subplot(subplots, subplots, i+1);
        plt.imshow(a_fields[i](x, y))
        plt.colorbar(orientation='horizontal')
    plt.savefig(os.path.join(test_dir, 'a_fields.png'))
    plt.close()



