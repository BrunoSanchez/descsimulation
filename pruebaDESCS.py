import os
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import sys

sys.path.insert(0, os.path.abspath('..'))

import zackay_code.properimage.propercoadd.propercoadd as pc

frame = '/home/bruno/Devel/DESCSimulation/imagenes-LN-PSF/im_06_10.fits'

im = pc.SingleImage(frame, imagefile=True)

#fitted_models = sim.fit_psf_sep()

# =============================================================================
#    PSF spatially variant
# =============================================================================
a_fields, psf_basis = im.get_variable_psf()

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

x, y = np.mgrid[:4096, :2048]
plt.figure(figsize=(16, 16))
for i in range(len(a_fields)):
    plt.subplot(subplots, subplots, i+1);
    plt.imshow(a_fields[i](x, y))
    plt.colorbar(orientation='horizontal')
plt.savefig(os.path.join(test_dir, 'a_fields.png'))
plt.close()



