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

from cosredux import utils as cr_utils


def set_background_region(obj_tr, segm, coadd_corrtag_woPHA_file, apert=25., ywidth=50., low_yoff=-10.,
                          check=False):
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
    bg_region : dict

    """
    # Load
    data = Table.read(coadd_corrtag_woPHA_file)
    wave = data['WAVELENGTH'].data
    # Set
    bg_region = {}
    if segm == 'FUVA':
        x1=1200. #1315.   #more?
        x2=max(wave)      #2400.   #approx/
    elif segm == 'FUVB':
        x1=50.
        x2=max(wave) #1000.

    bg_region['lower'] = (x1,x2, obj_tr-apert+low_yoff, obj_tr-apert-ywidth+low_yoff)
    if segm == 'FUVA':
        bg_region['upper'] = (x1,x2, obj_tr+apert, obj_tr+apert+ywidth)

    # Write and Return
    outfile = coadd_corrtag_woPHA_file.replace('.fits', '_bgregion.json')

    if check:
        yfull = data['YFULL']
        plt.scatter(wave,yfull,s=1)
        # Upper
        for key in ['lower', 'upper']:
            try:
                x1,x2,y1,y2 = bg_region[key]
            except KeyError:
                continue
            plt.plot([x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1],'r',linewidth=3.3)
        plt.xlim(x1-10,x2+10)
        plt.ylim(200., 800.)
        plt.show()

    # Return
    return bg_region



def get_pha_values_science(region, corrtagfile, background=True):
    """ Grab the PHA values in the input region
    Parameters
    ----------
    corrfile

    Returns
    -------
    phas : ndarray
    xdopp_min : float
    xdopp_max : float

    """
    if background:
        reg_keys = ['lower', 'upper']
    else:
        reg_keys = ['extraction']
    # Read
    data = Table.read(corrtagfile)
    wave = data['WAVELENGTH']
    xdopp = data['XDOPP'].data
    yfull = data['YFULL']
    dq = data['DQ']
    pha = data['PHA'].data

    # Loop
    all_phas = []
    all_xdopp = []
    for key in reg_keys:
        try:
            x1,x2,y1,y2 = region[key]
        except KeyError:
            print("No region={:s} provided.  Skipping it".format(key))
            continue
        iphas = (wave >= x1) & (wave <= x2) & (dq == 0) & (yfull > y1) & (yfull < y2)
        all_phas.append(pha[iphas])
        all_xdopp.append(xdopp[iphas])
    # Concatenate
    all_phas = np.concatenate(all_phas)
    all_xdopp = np.concatenate(all_xdopp)
    xdopp_min, xdopp_max = np.min(all_xdopp), np.max(all_xdopp)

    # Return
    return all_phas, xdopp_min, xdopp_max


def get_pha_values_dark(bg_region, dark_corrtag, xdopp_mnx):
    """ Grab the PHA values in the input region
    Parameters
    ----------
    corrfile

    Returns
    -------
    phas : ndarray
    xdopp_min : float
    xdopp_max : float

    """
    # Read
    data = Table.read(dark_corrtag)
    xdopp = data['XDOPP'].data
    yfull = data['YFULL']
    dq = data['DQ']
    pha = data['PHA'].data

    all_phas = []
    for key in ['lower', 'upper']:
        try:
            x1,x2,y1,y2 = bg_region[key]
        except KeyError:
            print("No region={:s} provided.  Skipping it".format(key))
            continue
        # Over-ride with xdopp
        x1, x2 = xdopp_mnx
        iphas = (xdopp >= x1) & (xdopp <= x2) & (dq == 0) & (yfull > y1) & (yfull < y2)
        all_phas.append(pha[iphas])
    # Concatenate
    all_phas = np.concatenate(all_phas)

    # Return
    return all_phas

def extract_dark_spectrum(coadd_dark_file, science_exp_file, obj_tr, segm, pha_mnx, apert=25., plot=False,
                          npix=16385):

    # Read coadded dark
    dark_data = Table.read(coadd_dark_file)
    xfull = dark_data['XFULL'].data
    dq=dark_data['DQ']
    yfull=dark_data['YFULL']
    pha=dark_data['PHA']

    # Get SHIFT1X from science header
    sci_header = fits.open(science_exp_file)[1].header
    if segm == 'FUVA':
        sci_shift = sci_header['SHIFT1A']
    elif segm == 'FUVB':
        sci_shift = sci_header['SHIFT1B']
    else:
        raise IOError("Bad segm input")

    # Extract
    in_extraction = (dq == 0) & (abs(yfull-obj_tr) <= apert) & (pha >= pha_mnx[0]) & (pha <= pha_mnx[1])
    pha_ex = pha[in_extraction]
    xfull_ex = xfull[in_extraction]

    # Histogram
    xshift = xfull_ex - sci_shift
    hbins=min(xshift)+np.arange(npix)*(max(xshift)-min(xshift))/(npix-1.)
    hist, edges = np.histogram(xshift, bins=hbins)

    # Plot
    if plot:
        plt.plot(hbins[:-1]+0.5, hist)
        plt.show()

    # Return
    return hist, pha_ex

def perform_kstest(sci_phas, dark_phas, criterion=0.1):
    """
    Parameters
    ----------
    sci_phas
    dark_phas
    criterion

    Returns
    -------

    """
    import scipy.stats as scc
    Dstat, PK_S = scc.ks_2samp(sci_phas, dark_phas)
    # Good?
    if PK_S > criterion:
        return True
    else:
        return False

def dark_to_exposures(exposures, bg_region, obj_tr, segm, defaults, min_ndark=4, show_spec=False,
                      N_smooth=500, verbose=True):
    """
    Parameters
    ----------
    exposures : list
      List of exposures without PHA filtering
    bg_region
    obj_tr
    segm
    defaults
    min_ndark
    show_spec
    N_smooth

    Returns
    -------

    """

    iend = exposures[0].rfind('/')
    sci_path = exposures[0][0:iend+1]
    dark_path = sci_path+'darks_'
    if segm == 'FUVA':
        dark_path += 'a/'
        sub_seg = 'a'
    elif segm == 'FUVB':
        dark_path += 'b/'
        sub_seg = 'b'
    # HVLEVELs
    hva, hvb = cr_utils.get_hvlevels(exposures)
    # Loop
    for ss, exposure in enumerate(exposures):
        if verbose:
            print("Working on exposure: {:s}".format(exposure))
        # Paths
        # PHA values in science region + xdopp values
        pha_values, xdopp_min, xdopp_max = get_pha_values_science(bg_region, exposure)
        # Find list of darks
        if segm == 'FUVA':
            sub_folder = 'a_{:d}/'.format(hva[ss])
        elif segm == 'FUVB':
            sub_folder = 'b_{:d}/'.format(hvb[ss])
        dark_list = glob.glob(dark_path+sub_folder+'*corrtag*')
        # Loop on Darks -- Keep list of good ones
        if verbose:
            print("Matching darks to the exposure")
        gd_darks= []
        for darkfile in dark_list:
            # PHAS
            drk_phas = get_pha_values_dark(bg_region, darkfile, (xdopp_min,xdopp_max))
            # KS
            if perform_kstest(pha_values, drk_phas):
                gd_darks.append(darkfile)
            ndark = len(gd_darks)
        if ndark < min_ndark:
            print("Need more darks.  Probably can just change criterion")
            return
        else:
            print("We have {:d} darks to use".format(ndark))

        # Coadd
        i0 = exposure.rfind('/')
        i1 = exposure.rfind('_corrtag')
        root_file = exposure[i0+1:i1]
        dark_coadd_file = dark_path+sub_folder+root_file+'_darks.fits'
        _ = cr_utils.coadd_bintables(gd_darks, outfile=dark_coadd_file)

        # Scaling
        coadd_drk_phas = get_pha_values_dark(bg_region, dark_coadd_file, (xdopp_min,xdopp_max))
        scale_sci_drk = float(pha_values.size) / coadd_drk_phas.size
        if verbose:
            print("Will use scale factor = {:g}".format(scale_sci_drk))

        # Extract
        spec, pha_ex = extract_dark_spectrum(dark_coadd_file, exposure, obj_tr, segm, defaults['pha_mnx'], plot=show_spec)
        if verbose:
            print("Extracting..")

        # Smooth
        smooth_spec = np.convolve(spec, np.ones((N_smooth,))/N_smooth, mode='same')

        # Generate a dark spectrum file
        outfile = sci_path+root_file+'_{:s}_bkgd.fits'.format(sub_seg)
        tbl = Table()
        tbl['DARK'] = smooth_spec
        tbl.write(outfile, overwrite=True)
        if verbose:
            print("Wrote background spectrum to {:s}".format(outfile))





