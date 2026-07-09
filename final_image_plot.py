#/usr/bash/python
# a python script
import matplotlib.pyplot as plt
import sys
import os
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
import colormaps as cmaps
path_sofia_templates='~/'



def crop_to_valid_data(fits_filepath):
    hduli=fits.open(fits_filepath)
    header = hduli[0].header
    data = hduli[0].data
    np.where
    wcs = WCS(header)

    valid_coords = np.argwhere(~np.isnan(data))

    # np.argwhere outputs coordinates in (row, col) i.e., (y, x) order
    ymin, xmin = valid_coords.min(axis=0)
    ymax, xmax = valid_coords.max(axis=0)

    height = (ymax - ymin) + 20
    width = (xmax - xmin) + 20

    center_x = xmin + width/2.0
    center_y = ymin + height/2.0

    cutout = Cutout2D(data, position=(center_x, center_y), size=(height, width), wcs=wcs)

    new_hdu = fits.PrimaryHDU(data=cutout.data)
    new_hdu.header.update(cutout.wcs.to_header())

    new_hdu.writeto(fits_filepath, overwrite=True)



with open('autoselfcal_par.txt') as f:
    lines = [line.rstrip('\n') for line in f if not line.startswith('#')]


img_name=lines[0]#'D100_7_12_autoselfcal'
phase_c=lines[1]#'J2000 12:59:02.5974 +027.38.25.374'
imsz=[int(lines[2]),int(lines[3])]
v_c= int(lines[4]) #5700 #km/s
deltaV=int(lines[5])

print(img_name,phase_c,imsz,v_c,deltaV)
suffix='_Robust2_dilation5'
img_name2=lines[0]+suffix
#### 1 RUN sofia per immagini pulite
if not img_name+'_plot_mom0.fits' in os.listdir() or not img_name2+'_plot_mom0.fits' in os.listdir():
    os.system('sofia '+path_sofia_templates+'plot_image_template.par input.data='+img_name+'_pbcorr.fits'+' input.primaryBeam='+img_name+'_pb.fits output.filename='+img_name+'_plot')
    os.system('sofia '+path_sofia_templates+'plot_image_template.par input.data='+img_name2+'_pbcorr.fits'+' input.primaryBeam='+img_name2+'_pb.fits output.filename='+img_name2+'_plot')


######## 2 correzione mom1 per v_c
def v_correction(img_name, v_c):
    hdul = fits.open(img_name)
    vel=hdul[0].data
    vel_cor=vel-v_c*np.ones_like(vel)*1e3
    fits.writeto(img_name.replace('.fits','_corr.fits'), vel_cor, hdul[0].header,overwrite=True)

v_correction(img_name+'_plot_mom1.fits',v_c)
v_correction(img_name2+'_plot_mom1.fits',v_c)
##### 3 triple panel mom0-mom1-mom2

def zero(img_name):
    hdul = fits.open(img_name)
    vel=hdul[0].data
    vel_cor=np.where(vel==0.0,np.nan,vel)
    fits.writeto(img_name.replace('.fits','_corr.fits'), vel_cor, hdul[0].header,overwrite=True)


zero(img_name+'_plot_mom0.fits')
zero(img_name2+'_plot_mom0.fits')

crop_to_valid_data(img_name+'_plot_mom0_corr.fits')
crop_to_valid_data(img_name+'_plot_mom1_corr.fits')
crop_to_valid_data(img_name+'_plot_mom2.fits')
#
#
crop_to_valid_data(img_name2+'_plot_mom0_corr.fits')
crop_to_valid_data(img_name2+'_plot_mom1_corr.fits')
crop_to_valid_data(img_name2+'_plot_mom2.fits')


def plotter(img_name):
    fig = plt.figure(figsize=(18, 6))



    mom0 = fits.open(img_name+'_mom0_corr.fits')
    mom0_data=mom0[0].data
    mom1 = fits.open(img_name+'_mom1_corr.fits')
    mom1_data=mom1[0].data
    mom2 = fits.open(img_name+'_mom2.fits')
    mom2_data=mom2[0].data



    ax1 = fig.add_subplot(1, 3, 1, projection=WCS(mom0[0].header))
    ax2 = fig.add_subplot(1, 3, 2, projection=WCS(mom1[0].header))
    ax3 = fig.add_subplot(1, 3, 3, projection=WCS(mom2[0].header))



    axes = [ax1, ax2, ax3]
    datasets = [mom0_data, mom1_data, mom2_data]
    titles = ['mom0', 'mom1', 'mom2']

    im1=ax1.imshow(np.log10(datasets[0]), cmap=cmaps.teal, origin='lower',vmin=-0.5,vmax=3.0)
    im2=ax2.imshow(datasets[1]/1e3, cmap=cmaps.rdbu.discrete(8), origin='lower',vmin=-deltaV/2.,vmax=deltaV/2.)
    im3=ax3.imshow(datasets[2]/1e3, cmap=cmaps.torch.discrete(8), origin='lower',vmin=1.,vmax=100.)

    plt.colorbar(im1,ax=ax1,location='bottom',label='log10I [Jy/beam]')
    plt.colorbar(im2,ax=ax2,location='bottom',label='V [km/s]')
    plt.colorbar(im3,ax=ax3,location='bottom',label='log10sigmaV [km/s]')


    # 3. Loop through panels to plot data and format coordinates
    for ax, data, title in zip(axes, datasets, titles):
        # Plot the image data
        ax.set_title(title, fontsize=14, pad=12)

        ax.tick_params(axis='both', direction='in')

        # Format labels
        ax.coords[0].set_axislabel('Right Ascension (J2000)', minpad=0.5)
        ax.coords[1].set_axislabel('Declination (J2000)', minpad=0.5)

        # Optional: Force specific tick formatting (e.g., hh:mm:ss and dd:mm:ss)
        ax.coords[0].set_major_formatter('hh:mm:ss')
        ax.coords[1].set_major_formatter('dd:mm:ss')

    # Adjust layout to prevent label clipping between panels
    plt.tight_layout()

    # Display or save the figure
    plt.savefig(img_name+'.png',dpi=200)
plotter(img_name+'_plot')
plotter(img_name2+'_plot')
##### 4 single panel mom0 figo without frame
fig = plt.figure()
mom0 = fits.open(img_name2+'_plot_mom0_corr.fits')
mom0_data=mom0[0].data
ax1 = fig.add_subplot(1, 1, 1, projection=WCS(mom0[0].header),aspect=1.0)
#ax1.imshow(np.log10(mom0_data), cmap=cmaps.guppy, origin='lower',vmin=0.,vmax=3.0)
ax1.contourf(np.log10(mom0_data), cmap=cmaps.gem, origin='lower',levels=np.linspace(0.1,4,15))
ax1.contour(np.log10(mom0_data), colors='white',lw=0.1,alpha=0.7, origin='lower',levels=np.linspace(0.1,4,15))

ax1.axis('off')
plt.savefig(img_name2+'_solo.png',dpi=200)
