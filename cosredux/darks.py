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



def get_pha_values_science(bg_region, corrtagfile):
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
    data = Table.read(corrtagfile)
    wave = data['WAVELENGTH']
    xdopp = data['XDOPP'].data
    yfull = data['YFULL']
    dq = data['DQ']
    pha = data['PHA'].data

    # Loop
    all_phas = []
    all_xdopp = []
    for key in ['lower', 'upper']:
        try:
            x1,x2,y1,y2 = bg_region[key]
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

