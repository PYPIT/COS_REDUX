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
        x1=50.
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
        plt.ylim(200., 800.)
        plt.show()

    # Return
    return ex_region




