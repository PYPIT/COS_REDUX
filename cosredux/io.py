""" Handles I/O for COS products
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
from astropy.table import Table
from astropy.io import fits

from xastropy.xutils import xdebug as xdb

def modify_table_value(filename, column, row_dict, value, outfil=None, clobber=False):
    """ Open a file, grab rows of interest, update the value
    Parameters
    ----------
    filename
    column : str
    row_dict : dict
    value : float, str, int
    outfil : str, optional
      If provided, generates a new file with the Table (and header)

    Returns
    -------
    tbl : Table
      Modified as desired

    """
    # Read
    tbl = Table.read(filename)
    # Find the row(s)
    mask = np.array([True]*len(tbl))
    for key, item in row_dict.items():
        mask &= (tbl[key] == item)

    # Set value
    tbl[column][mask] = value

    # Write?
    if outfil is not None:
        hdu = fits.open(filename)
        phdu = fits.PrimaryHDU()
        phdu.header = hdu[0].header
        thdu = fits.table_to_hdu(tbl)
        thdu.header = hdu[1].header
        thdulist = fits.HDUList([phdu,thdu])
        thdulist.writeto(outfil, clobber=clobber)

    # Return
    return tbl
