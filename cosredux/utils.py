""" Utility routines
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
import os
import glob
import shutil
import pdb

from astropy.table import Table
from astropy.io import fits

from xastropy.xutils import xdebug as xdb





def modify_rawtag_for_calcos(path, verbose=False):
    """ Open rawtag files, edit header for running Calcos 1st time on a dataset
    WARNING: Over writes the header of the previous file

    Parameters
    ----------
    path : str
      path to files

    Returns
    -------

    """
    # Find rawtag files
    rawfiles = glob.glob(path+'/*rawtag*')
    # Loop on the rawtag files
    for rawfile in rawfiles:
        if verbose:
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
    filename : str
    column : str
    row_dict : dict
    value : float, str, int
      new value
    outfil : str, optional
      If provided, generates a new file with the Table (and header)
    clobber : bool, optional
      overwrite?

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
        thdulist.writeto(outfil, overwrite=clobber)

    # Return
    return tbl


def modify_LP2_1dx_calib(calib_path, LP3_1dx, OPT_ELEM='G140L', CENWAVE=1280, verbose=True):
    """  Modify WCA and PSA definitions in the 1dx calibration file
    of LP2 according to the values in the file for LP3
    Only necessary for CALCOS v2

    Step 3b of the procedure

    Parameters
    ----------
    calib_path : str
    LP3_1dx : str
    OPT_ELEM : str, optional (assumed: 'G140L')
    CENWAVE : int ; optional (assumed 1280)


    Returns
    -------

    """
    lp3_dict = {}
    # Read LP3 file
    LP3_1dx_file = calib_path+LP3_1dx  ## e.g. 'z2d19237l_1dx.fits'
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
    thdulist.writeto(LP2_1dx_file, overwrite=True)
    # Return
    return


def coadd_bintables(infiles, outfile=None, clobber=True):
    """ Coadd a set of input BINARY table FITS files

    Parameters
    ----------
    infiles : list
      list of files to coadd
    outfile : str, optional
      If given, generate a new FITS file
    clobber : bool, optional
    #add_sun_target_info : bool, optional
    #  Should be True for a corrtag combine

    Returns
    -------
    tbltot : Table
      combined table

    """
    from astropy.table import vstack

    tbllist = []
    head0 = None
    fileslist = infiles[0]
    for ifile in infiles:
        # Read
        hdu = fits.open(ifile)
        tbl1 = Table(hdu[1].data)
        if head0 is None:
            head0 = hdu[0].header
        else:
            pass
        if (ifile != infiles[0]):
            fileslist += ', '+ifile

        tbllist.append(tbl1)
        hdu.close()

    # Coadd tables
    tbltot = vstack(tbllist, join_type='exact')

    # write new file (header from first)
    if outfile is not None:
        phdu = fits.PrimaryHDU()
        phdu.header = head0
        phdu.header['history'] = 'Combined fits files '+fileslist   # Add  HISTORY lines to header
        thdu = fits.table_to_hdu(tbltot)
        thdulist = fits.HDUList([phdu, thdu])
        thdulist.writeto(outfile, overwrite=clobber)
        print("Wrote outfile {:s}".format(outfile))

    # Return
    return tbltot



'''

def add_sun_target_columns(corrtag_files_n, clobber=True):
    """  Add SUN_ALT and TARGET_ALT columns to coadded file
     WOULD NEED TO INTERPOLATE!
     SHOULD UPDATE TABLE1 in EACH CORRTAG FILE

    corrtag_files_n : list
      List of corrtag files
    coadded_file : str
      Filename of coadded corrtag file
    outfile : str
      Output filename
    :return:
    """
    sunalts = []
    limbangs = []
    for file in corrtag_files_n:
        hdu = fits.open(file)
        tbl = Table(hdu[3].data)
        sunaltf = tbl['SUN_ALT']
        limbangf = tbl['TARGET_ALT']
        hdu.close()
        sunalts = np.append(sunalts, sunaltf)
        limbangs = np.append(limbangs, limbangf)
    # print(max(sunalt),len(sunalt),len(sunalt1))

    data = fits.open(coadded_file)[1].data
    cols = []
    cols.append( fits.Column(name=str('SUN_ALT'), format='E', array=sunaltf) )
    cols.append( fits.Column(name=str('LimbAng'), format='E', array=limbangf) )
    orig_cols = data.columns
    new_cols = fits.ColDefs(cols)
    hdu = fits.BinTableHDU.from_columns(orig_cols + new_cols)
    hdu.writeto(outfile, overwrite=clobber)
    ##hdu.close()
'''


def get_hvlevels(files):
    """ Grab the HV Levels from a given file
    Parameters
    ----------
    files : list of files

    Returns
    -------
    hva : list
     list of HVLEVELA values
    hvb : list
     list of HVLEVELB values

    """
    hva, hvb = [], []
    for ifile in files:
        hdu = fits.open(ifile)
        ihva, ihvb = hdu[1].header['HVLEVELA'], hdu[1].header['HVLEVELB']
        hva.append(ihva)
        hvb.append(ihvb)
        print("HVA={:g} and HVB={:g} in file {:s}".format(ihva, ihvb, ifile))
        hdu.close()
    return hva, hvb





def change_pha(calibfld, lowa=2, upa=15, lowb=2, upb=15):
    """ Modify PHA values in calibration PHA files
    Worseck recommends doing this file by file..

    Parameters
    ----------
    calibfld : str
    lowa : int, optional
    upa : int, optional
    lowb : int, optional
    upb : int, optional
      lower and upper limits for PHA values in FUVA and FUVB

    Returns
    -------

    """
    phafiles = glob.glob(calibfld + '*pha.fits')
    for phafile in phafiles:
        print("Updating PHA values in {:s}".format(phafile))
        with fits.open(phafile, 'update') as f:
            head1 = f[1].header
            head1['PHALOWRA'] = lowa
            head1['PHALOWRB'] = lowb
            head1['PHAUPPRA'] = upa
            head1['PHAUPPRB'] = upb


def clean_for_calcos_phafiltering(redux_dir):
    """ Remove files before running calcos 2nd time (with PHACORR = PERFORM)
    Push corrtag files to woPHA

    Parameters
    ----------
    redux_dir : str

    Returns
    -------

    """
    # Rename corrtag files
    corrtag_files = glob.glob(redux_dir+'/*_corrtag_*')
    for cfile in corrtag_files:
        new_cfile = cfile.replace('corrtag_', 'corrtag_woPHA_')
        print("Renaming corrtag file to {:s}".format(new_cfile))
        os.rename(cfile, new_cfile)

    # Remove unwanted files
    for bad_exten in ['_flt', '_x1d', '_lampflash', '_counts', 'jnk']:
        bad_files = glob.glob(redux_dir+'/*'+bad_exten+'*')
        for bad_file in bad_files:
            print("Removing {:s}".format(bad_file))
            os.remove(bad_file)



def modify_phacorr(rawtag_path):
    """   Change PHACORR = PERFORM
    Parameters
    ----------
    rawtag_path : str
      Path to rawtag files

    Returns
    -------

    """
    # Find rawtag files
    rawfiles = glob.glob(rawtag_path+'/*rawtag*')
    # Loop on the rawtag files
    for rawfile in rawfiles:
        print("Modifying PHACORR header card to PERFORM for rawtag file: {:s}".format(rawfile))
        with fits.open(rawfile, 'update') as f:
            hdu0 = f[0]
            hdu0.header['PHACORR'] = 'PERFORM'


def change_dq_wgt(x1d_folder, clobber=True):
    """ Update DQ values in x1d frames

    Parameters
    ----------
    x1d_folder : str
       path to x1d files
    clobber : bool, optional
      overwrite

    Returns
    -------

    """
    x1dfiles = glob.glob(x1d_folder + '*_x1d*.fits')

    clobber = True
    for filename in x1dfiles:
        # read DQ
        hdu = fits.open(filename)
        tbl = Table(hdu[1].data)

        # without for:
        badDQ = (tbl['DQ'] > 2) & (tbl['DQ'] != 1024)
        tbl['DQ_WGT'][badDQ] = 0
        goodDQ = tbl['DQ'] == 1024
        tbl['DQ_WGT'][goodDQ] = 1

        # Write
        hdu = fits.open(filename)
        phdu = fits.PrimaryHDU()
        phdu.header = hdu[0].header
        thdu = fits.table_to_hdu(tbl)
        thdu.header = hdu[1].header
        thdulist = fits.HDUList([phdu, thdu])
        thdulist.writeto(filename, overwrite=clobber)

        hdu.close()



def copy_files(pths,ext):
    """ Make a copy of files in paths pths that are ending with ext+'.fits'
        Copies are in the same paths, and end with _copy1.fits

    Parameters
    ----------
    pths : list of str
       list of paths to the files
    ext : str
      files end with this string, e.g. '_1dx'

    Returns
    -------

    """
    for pth in pths:
        fls = glob.glob(pth+'*'+ext+'.fits')
        for fl in fls:
            flcp = fl[:-5]+'_copy1.fits'
            if len(glob.glob(flcp)) == 0:
                print('Copying file {:s}'.format(fl))
                shutil.copyfile(fl,flcp)
            else:
                print('File {:s} already exists'.format(flcp))
        print(' ')



def load_fnames(pth, lenname=7):
    """ Load paths to data and folder names

    Parameters
    ----------
    pth : str
       path to data folders, including * (for glob)
    lenname : str
      length of names of folders with data
      (should be the same for all data)

    Returns
    -------
    ff1 : list of str
      paths to data
    fldnames : list of str
      folder names
    """

    ff1 = []
    fldnames = []
    ff = np.sort(glob.glob(pth))
    for i in range(len(ff)):
        if len(ff[i]) == (len(pth)+lenname-2):
            ff1.append(ff[i])
            fldnames.append(ff[i][len(pth)-3:-1])
    return ff1, fldnames


def printallspoff(pths):
    """ Print SP_OFF values for all CORRTAG files (should be same for RAWTAG files)

    Parameters
    ----------
    pths : list of str
       list of paths to the files

    Returns
    -------

    """

    for ifldpth in pths:
        raws = glob.glob(ifldpth + '*corrtag_*.fits')
        for raw in raws:
            with fits.open(raw) as f:
                print(raw, f[1].header['SP_OFF_A'], f[1].header['SP_OFF_B'])



def summfile(ifldpth,idx):
    """ Find _x1dsum file ending with idx + '_x1dsum.fits'

    Parameters
    ----------
    ifldpth : str
       path to the files
    idx : str
      the file ends with idx + '_x1dsum.fits'
      e.g. idx = '10'

    Returns
    -------
      fname : str
        sum file name

    """

    fname = glob.glob(ifldpth+'*'+idx+'_x1dsum.fits')
    if len(fname) != 1:
        print('multiple or zero '+idx+'_x1dsum files: ', len(fname))
        print(ifldpth)
    else:
        fname = fname[0]
    return fname

