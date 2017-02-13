""" Utility routines
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
import os
import glob
import pdb

from matplotlib import pyplot as plt

from astropy.table import Table
from astropy.io import fits

from xastropy.xutils import xdebug as xdb


def set_extraction_region(obj_tr, segm, coadd_corrtag_woPHA_file, apert=25., check=False):
    """
    Parameters
    ----------
    obj_tr
    segm
    apert
    corrtag_file : str
      For wavelength info
    ywidth : float, optional

    Returns
    -------
    ex_region : dict

    """
    # Load
    data = Table.read(coadd_corrtag_woPHA_file)
    wave = data['WAVELENGTH'].data
    # Set
    if segm == 'FUVA':
        x1=1200. #1315.   #more?
        x2=max(wave)      #2400.   #approx/
    elif segm == 'FUVB':
        x1=900. ##50.
        x2=max(wave) #1000.

    ex_region = {}
    ex_region['extraction'] = [x1,x2, obj_tr-apert, obj_tr+apert]

    # Write and Return
    #outfile = coadd_corrtag_woPHA_file.replace('.fits', '_exregion.json')

    if check:
        yfull = data['YFULL']
        plt.scatter(wave,yfull,s=1)
        # Region
        x1,x2,y1,y2 = ex_region['extraction']
        plt.plot([x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1],'b',linewidth=3.3)
        # Axes
        plt.xlim(x1-10,x2+10)
        if segm == 'FUVB':
            plt.xlim(50.,x2+10)
        plt.ylim(200., 800.)
        plt.show()

    # Return
    return ex_region

def coadd_exposures(x1d_files, segm, outfile, bin=None): #, checkdq0=False):
    """
    Parameters
    ----------
    x1d_files
    segm
    outfile

    Returns
    -------

    """

    from scipy.interpolate import interp1d
    if segm == 'FUVA':
        spec_row = 0
        subseg = 'a'
    elif segm == 'FUVB':
        spec_row = 1
        subseg = 'b'

    # Load
    xtbls = []
    dark_files = []
    for x1d_file in x1d_files:
        dark_files.append(x1d_file.replace('_x1d.fits','_{:s}_bkgd.fits'.format(subseg)))
        if not os.path.isfile(dark_files[-1]):
            print("No background file named {:s}".format(dark_files[-1]))
            raise IOError("Make it with dark_to_exposures()")
        #
        xtbl = Table.read(x1d_file)
        xtbls.append(xtbl[spec_row:spec_row+1])

    # Grab one wavelength array
    wave = xtbls[0]['WAVELENGTH'][0,:].data

    # Sum exposure time
    total_time = np.zeros_like(wave)
    for xtbl in xtbls:
        total_time += xtbl['DQ_WGT'][0,:]*xtbl['EXPTIME']

    # Find DQmin for all exposures -- Why are we doing this step??
    dqmin = np.ones_like(wave).astype(int) * 99999
    for xtbl in xtbls:
        # Reset DQ
        dq = xtbl['DQ'][0,:].data
        reset_1024 = dq == 1024
        dq[reset_1024] = 2
        dqmin = np.minimum(dq, dqmin)

    # Find DQ_WGT max for all exposures
    DQWmax = np.zeros_like(wave)
    for xtbl in xtbls:
        # Reset DQ
        dqw = xtbl['DQ_WGT'][0,:].data
        DQWmax = np.maximum(dqw, DQWmax)


    # ####################
    # CALIBRATION
    wave_calib, calib = [], []
    for xtbl in xtbls:
        gddq = (xtbl['DQ'] == 0) & (xtbl['FLUX'] > 0)  ## not dq > 0
        ### checkdq0
        #if checkdq0 == True:
        #    gddq = (xtbl['DQ'] == 0) & (xtbl['FLUX'] > 0)

        # Append
        wave_calib.append(xtbl['WAVELENGTH'][gddq].data.flatten())
        calib.append( (xtbl['NET'][gddq] / xtbl['FLUX'][gddq]).data)

    # arrays
    wave_calib = np.concatenate(wave_calib)
    calib = np.concatenate(calib)
    # sort
    srt = np.argsort(wave_calib)
    wave_calib = wave_calib[srt]
    calib = calib[srt]

    # Cut down
    gdwv = wave_calib < 2100.

    # Spline
    sens_func = interp1d(wave_calib[gdwv], calib[gdwv], bounds_error=False, fill_value=0.)  # cubic behaves badly


    # Total counts in science and background
    total_counts = np.zeros_like(wave)
    total_dark = np.zeros_like(wave)
    for ss, xtbl in enumerate(xtbls):
        # Science
        dqw = xtbl['DQ_WGT'][0,:].data
        total_counts += dqw * xtbl['GCOUNTS'][0,:]
        # Dark
        bkgd = Table.read(dark_files[ss])
        total_dark += dqw * bkgd['DARK'].data

    # Bin
    if bin is not None:
        # Check
        if bin not in [2,3]:
            raise IOError("Only ready for binning by 2 or 3 channels")
        # Ugly for loop
        nchannel = len(total_counts)
        new_tot_counts, new_tot_dark, new_tot_time, new_wave, new_dqmin, new_DQW = [], [], [], [], [], []
        for kk in np.arange(0, nchannel, bin):
            # Simple stuff sums
            new_tot_counts.append(np.sum(total_counts[kk:kk+bin]))
            new_tot_dark.append(np.sum(total_dark[kk:kk+bin]))
            new_tot_time.append(np.sum(total_time[kk:kk+bin]))
            new_dqmin.append(np.min(dqmin[kk:kk+bin]))
            new_DQW.append(np.max(DQWmax[kk:kk+bin]))
            # Wavelength
            new_wave.append(np.mean(wave[kk:kk+bin]))
        # Turn into arrays
        new_tot_counts = np.array(new_tot_counts)
        new_tot_dark = np.array(new_tot_dark)
        new_tot_time = np.array(new_tot_time)
        new_wave = np.array(new_wave)
    else:
        new_tot_counts, new_tot_dark, new_tot_time, new_wave = total_counts, total_dark, total_time, wave
        new_dqmin, new_DQW = dqmin, DQWmax


    # Flux array
    flux = np.zeros_like(new_tot_time)
    calib = sens_func(new_wave)
    gd_time_sens = (new_tot_time > 0.) & (calib > 0.)
    flux[np.where(gd_time_sens)[0]] = (new_tot_counts[gd_time_sens]-new_tot_dark[gd_time_sens]) / (calib[gd_time_sens] * new_tot_time[gd_time_sens])

    # Simple error estimate
    error = np.zeros_like(new_tot_time)
    gd_error = (new_tot_time > 0.) & (calib > 0.) & (new_tot_counts > 0)
    error[np.where(gd_error)[0]] = np.sqrt(new_tot_counts[gd_error]) / (calib[gd_error] * new_tot_time[gd_error])

    # Final spectral information
    coadd = Table()
    coadd['wave'] = new_wave
    coadd['flux'] = flux
    coadd['error'] = error
    coadd['counts'] = new_tot_counts
    coadd['bgkd'] = new_tot_dark
    coadd['eff_time'] = new_tot_time
    coadd['calib'] = calib
    coadd['DQ_MIN'] = new_dqmin
    coadd['DQW_max'] = new_DQW

    # Write
    coadd.write(outfile, overwrite=True)
    print("Wrote {:s}".format(outfile))
