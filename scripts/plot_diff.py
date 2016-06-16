'''
Plot the difference image between two epochs of NGC 1404.
'''

from astropy.io import fits           # package for loading FITS files
from matplotlib import pyplot as plt  # plotting
from matplotlib import cm             # color maps
import numpy as np                    # math, stats, etc

f1 = fits.open('SN2011iv_B_SWO_DC_2011_12_11SN.fits')
f2 = fits.open('SN2011iv_B_template.fits')

image1 = f1[0].data      # Data  from the primary FITS HDU
image2 = f2[0].data      # (see http://fits.gsfc.nasa.gov/fits_primer.html)

# image data from the telescope often has saturated/bad pixels that make
# the real objects we're interested hard to see. Use percentile to figure
# out where most of the pixel values lie (1% -> 99%)
vmin1,vmax1 = np.percentile(np.ravel(image1), (1,99))
vmin2,vmax2 = np.percentile(np.ravel(image2), (1,99))

# Make a nice wide figure
fig = plt.figure(figsize=(18,6))
# 3x1 panels
ax1 = fig.add_subplot(131)
# Use the vmin,vmax calculated above. Also use a reverse greyscale color
# map and put the origin of the coordinates at lower left.
ax1.imshow(image1, cmap=cm.gray_r, vmin=vmin1, vmax=vmax1, origin='lower')
ax1.set_xlabel('X (pix)')
ax1.set_ylabel('Y (pix)')
ax1.set_title('NGC 1404  Dec. 11 2011')

ax2 = fig.add_subplot(132)
ax2.imshow(image2, cmap=cm.gray_r, vmin=vmin2, vmax=vmax2, origin='lower')
ax2.set_xlabel('X (pix)')
ax2.set_ylabel('Y (pix)')
ax2.set_title('NGC 1404  Stack 2007')

ax3 = fig.add_subplot(133)
ax3.imshow(image1-image2, cmap=cm.gray_r, vmin=vmin1, vmax=vmax1, 
   origin='lower')
ax3.set_xlabel('X (pix)')
ax3.set_ylabel('Y (pix)')
ax3.set_title('Difference')

fig.savefig('SN2011iv_diff.pdf')
