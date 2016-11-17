""" Utility routines
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
import glob
import pdb

from astropy.table import Table
from astropy.io import fits

from xastropy.xutils import xdebug as xdb



def modify_rawtag_for_calcos(path):
    """ Open rawtag files, edit header
    WARNING: Over writes the header of the previous file

    Parameters
    ----------
    path : str

    Returns
    -------

    """
    # Find rawtag files
    rawfiles = glob.glob(path+'/*rawtag*')
    # Loop on the rawtag files
    for rawfile in rawfiles:
        print("Modifying header cards for rawtag file: {:s}".format(rawfile))
        with fits.open(rawfile, 'update') as f:
            hdu0 = f[0]
            hdu0.header['FLATCORR'] = 'OMIT' #[flatfielding of Poisson data creates fractional counts, which are hard to interpret, COS flatfielding is approximate anyhow]
            hdu0.header['PHACORR'] = 'OMIT'  #[initial setting, to be tuned later]
            hdu0.header['RANDSEED'] = 0      # [choose any non-negative value for reproducable results]
            hdu0.header['BACKCORR'] = 'OMIT' # [we will do this in post-processing]
            # - **Set your calibration files wisely:**
            hdu0.header['SPOTTAB'] = 'N/A'   # [This calibration step only works on select datasets where hotspots have been identified. How and by whom? Better get rid of this nonsense.]
            hdu0.header['GSAGTAB'] = 'lref$x6l1439el_gsag.fits'  #[This will flag all gain sag holes as of 2013Jun20, so mainly geocoronal gain sag holes at the overused LP1.
                                                           #    Keeps the LP1 and LP2 traces usable that partly overlap with the dark calibration regions. All this is safe and tested!
                                                           #    Download this file from https://hst-crds.stsci.edu.]
            hdu0.header['XTRACTAB'] = 'lref$x6q17586l_1dx.fits'   # [This will work for LP2 data with CALCOS v2.21. In newer versions STScI renamed some columns, so they might not work with CALCOS v2.21.
                                                           #    All entries are to be modified later (see below).
                                                           #    Download this file from https://hst-crds.stsci.edu.]


def modify_table_value(filename, column, row_dict, value, outfil=None, clobber=False):
    """ Open a file, grab rows of interest, update the values

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


def modify_LP2_1dx_calib(calib_path, OPT_ELEM='G140L', CENWAVE=1280):
    """  Modify WCA and PSA definitions in the 1dx calibration file
    of LP2 according to the values in the file for LP3
    Only necessary for CALCOS v2

    Step 3b of the procedure

    Parameters
    ----------

    Returns
    -------

    """
    lp3_dict = {}
    # Read LP3 file
    LP3_1dx_file = calib_path+'/z2d19237l_1dx.fits'
    lp3 = Table.read(LP3_1dx_file)
    for segment in ['FUVA', 'FUVB']:
        lp3_dict[segment] = {}
        for aperture in ['WCA', 'PSA']:
            # Grab values
            row = (lp3['OPT_ELEM'] == OPT_ELEM) & (lp3['CENWAVE']==CENWAVE) & (
                lp3['SEGMENT'] == segment) & (lp3['APERTURE'] == aperture)
            if np.sum(row) != 1:
                pdb.set_trace()
            # Fill
            lp3_dict[segment][aperture] = {}
            lp3_dict[segment][aperture]['SLOPE'] = lp3[row]['SLOPE'].data[0]
            lp3_dict[segment][aperture]['B_SPEC'] = lp3[row]['B_SPEC'].data[0]

    # Open LP2 file for editing; read in info
    LP2_1dx_file = calib_path+'/x6q17586l_1dx.fits'
    hdu = fits.open(LP2_1dx_file)
    lp2 = Table(hdu[1].data)
    hdu0 = hdu[0]

    # Update values in LP2 using modify_table_value
    for segment in lp3_dict.keys():
        for aperture in lp3_dict[segment].keys():
            # Find row
            row = (lp2['OPT_ELEM'] == OPT_ELEM) & (lp2['CENWAVE']==CENWAVE) & (
                lp2['SEGMENT'] == segment) & (lp2['APERTURE'] == aperture)
            # Set values
            for column in lp3_dict[segment][aperture].keys():
                if (aperture == 'PSA') and (column == 'SLOPE'):
                    lp2[column][row] = 0.
                elif (aperture == 'PSA') and (column == 'B_SPEC'):
                    if OPT_ELEM == 'G140L':
                        lp2[column][row] = int(lp3_dict[segment][aperture][column]) + 0.5  # Assumes odd height c
                    else:
                        pdb.set_trace()  # Not ready for this element
                else:
                    lp2[column][row] = lp3_dict[segment][aperture][column]
            if OPT_ELEM == 'G140L':
                lp2[row]['HEIGHT'] = 25
            else:
                pdb.set_trace()  # Not ready for this element

    # Write
    thdu = fits.table_to_hdu(lp2)
    thdulist = fits.HDUList([hdu0,thdu])
    thdulist.writeto(LP2_1dx_file, clobber=True)

    #


### COADD

from astropy.table import join


def coadd_fits(fa1, fa2, fa3, fa4, fa, clobber=False):
    """
    Parameters
    ----------
    fa1, fa2, fa3, fa4 - fits files to coadd
    clobber

    Returns
    -------
    fa - coadded fa1,fa2,fa3,fa4

    """

    # Read
    hdu = fits.open(fa1)
    tbl1 = Table(hdu[1].data)
    hdu.close()
    hdu = fits.open(fa2)
    tbl2 = Table(hdu[1].data)
    hdu.close()
    hdu = fits.open(fa3)
    tbl3 = Table(hdu[1].data)
    hdu.close()
    hdu = fits.open(fa4)
    tbl4 = Table(hdu[1].data)
    hdu.close()

    # Coadd tables
    tbltot = join(tbl1, tbl2, join_type='outer')
    tbltot = join(tbltot, tbl3, join_type='outer')
    tbltot = join(tbltot, tbl4, join_type='outer')

    # write new file (only hdu[1] for now, and header as in fa1 file)
    filename0 = fa1
    hdu = fits.open(filename0)
    phdu = fits.PrimaryHDU()
    phdu.header = hdu[0].header
    thdu = fits.table_to_hdu(tbltot)
    thdulist = fits.HDUList([phdu, thdu])
    thdulist.writeto(fa, clobber=clobber)
    hdu.close()

    # Return
    return tbltot



