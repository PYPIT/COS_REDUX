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
    """ Defines background region (2 regions for FUVA, 1 region for FUVB)
    Parameters
    ----------
    obj_tr : float, int
      object trace
    segm : str
      segment
    coadd_corrtag_woPHA_file : str
      For wavelength info
    apert: float, optional
      aperture
    ywidth : float, optional
      width of the region
    low_yoff: float, int, optional
      offset from the aperture for the lower region
    check : bool, optional
      plot region

    Returns
    -------
    bg_region : dict
       dict with background region(s)

    """
    # Load
    data = Table.read(coadd_corrtag_woPHA_file)
    wave = data['WAVELENGTH'].data
    # Set
    bg_region = {}
    if segm == 'FUVA':
        x1=1200.
        x2=max(wave)
    elif segm == 'FUVB':
        x1=900.
        x2=max(wave)

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
        if segm == 'FUVB':
            plt.xlim(50.,x2+10)
        plt.ylim(min(yfull[wave > x1]),max(yfull[wave > x1]))
        plt.show()

    # Return
    return bg_region



def get_pha_values_science(region, corrtagfile, segm, background=True):
    """ Grab the PHA values in the input region
    Parameters
    ----------
    region : dict
      background or extraction region(s)
    corrtagfile: str
    segm: str
    background : bool, optional
      is it background region

    Returns
    -------
    all_phas : ndarray
       PHA values in the region
    xdopp_min : float
    xdopp_max : float
      min and max values for XDOPP in the region(s)

    """
    if background:
        reg_keys = ['lower', 'upper']
        if segm == 'FUVB':
            reg_keys = ['lower']
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
        if (y1 > y2):
            yn=y1
            y1=y2
            y2=yn
        if (x1 > x2):
            xn=x1
            x1=x2
            x2=xn
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
    bg_region : dict
       background region(s)
    dark_corrtag : str
       dark file
    xdopp_mnx : tuple
       min and max values for XDOPP in the region(s)

    Returns
    -------
    all_phas : ndarray
       PHA values in the region
    # xdopp_min : float
    # xdopp_max : float

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
        if (y1 > y2):
            yn=y1
            y1=y2
            y2=yn
        if (x1 > x2):
            xn=x1
            x1=x2
            x2=xn
        iphas = (xdopp >= x1) & (xdopp <= x2) & (dq == 0) & (yfull > y1) & (yfull < y2)
        all_phas.append(pha[iphas])
    # Concatenate
    all_phas = np.concatenate(all_phas)

    # Return
    return all_phas

def extract_dark_spectrum(coadd_dark_file, science_exp_file, obj_tr, segm, pha_mnx, apert=25., offs1=0., offs2=0., plot=False,
                          npix=16385):
    """ Extracts dark spectrum for a (coadded) dark file.

    Parameters
    ----------
    coadd_dark_file : str
      coadded dark file
    science_exp_file : str
      science exposure file
    obj_tr : float, int
      object trace
    segm: str
      segment
    pha_mnx: tuple
      min and max PHA values
    apert : float, int, optional
      aperture
    offs1 : float, int, optional
    offs2 : float, int, optional
      left and right offset from the aperture
    plot: bool, optional
    npix: int, optional
      number of pixels for xfull

    Returns
    -------
    hist : ndarray
      spectrum
    pha_ex : ndarray
      PHA values in the extraction region

    """


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
    in_extraction = (dq == 0) & (yfull <= (obj_tr + apert +offs2)) & (yfull >= (obj_tr - apert +offs1)) & (pha >= pha_mnx[0]) & (pha <= pha_mnx[1])
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
    """ Compare the cumulative pulse height histogram of the science exposure to the
    corresponding histograms of the dark exposures,
    and check if they are similar to a Kolmogorov-Smirnov test with some criterion

    Parameters
    ----------
    sci_phas : ndarray
    dark_phas : ndarray
       science and dark exposure PHAs
    criterion : float, optional
       for the KS test

    Returns
    -------
    True / False

    """
    import scipy.stats as scc
    Dstat, PK_S = scc.ks_2samp(sci_phas, dark_phas)
    # Good?
    if PK_S > criterion:
        return True
    else:
        return False


def dark_to_exposures(exposures, bg_region, obj_tr, segm, defaults, min_ndark=4, show_spec=False,
                      N_smooth=500, verbose=True, offs1=0.,offs2=0.):
    """ For each exposure coaddes dark spectra which are similar to science spectrum,
        and writes the spectrum in a file.
    Parameters
    ----------
    exposures : list
      List of exposures without PHA filtering
    bg_region : dict
      background region
    obj_tr : float, int
      object trace
    segm : str
      segment
    defaults : dict
      default values (here used: min and max PHAs)
    min_ndark : int, optional
      minimum number of darks to coadd
    show_spec : bool, optional
    N_smooth : int, optional
      smoothing
    verbose : bool, optional
    offs1 : float, int, optional
    offs2 : float, int, optional
      left and right offset from the aperture

    Returns
    -------

    """

    iend = exposures[0].rfind('/')
    sci_path = exposures[0][0:iend+1]
    hva, hvb = cr_utils.get_hvlevels(exposures)
    if segm == 'FUVA':
        sub_seg = 'a'
        hvl = hva
    elif segm == 'FUVB':
        sub_seg = 'b'
        hvl = hvb
    # HVLEVELs
    dark_path = sci_path+'darks_'+sub_seg
    # Loop
    for ss, exposure in enumerate(exposures):
        if verbose:
            print("Working on exposure: {:s}".format(exposure))
        # Paths
        # PHA values in science region + xdopp values
        pha_values, xdopp_min, xdopp_max = get_pha_values_science(bg_region, exposure, segm)
        # Find list of darks
        new_dark_path = dark_path+'_{:d}'.format(hvl[ss])+'/'
        dark_list = glob.glob(new_dark_path+'/*corrtag*')
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
            print("Changing criterion...")
            gd_darks = []
            pks_list = []
            import scipy.stats as scc
            for darkfile in dark_list:
                drk_phas = get_pha_values_dark(bg_region, darkfile, (xdopp_min, xdopp_max))
                Dstat, PK_S = scc.ks_2samp(pha_values, drk_phas)
                pks_list.append(PK_S)
            pks=np.asarray(pks_list)
            newcrit=np.partition(pks,len(pks)-min_ndark-1)[len(pks)-min_ndark-1]
            for i in np.arange(len(pks)):
                if pks[i] >= newcrit:
                    gd_darks.append(dark_list[i])
            print("New criterion: ", newcrit)
            ndark = len(gd_darks)
            #            return
        #else:
        #    print("We have {:d} darks to use".format(ndark))
        print("We have {:d} darks to use".format(ndark))

        # Coadd
        i0 = exposure.rfind('/')
        i1 = exposure.rfind('_corrtag')
        root_file = exposure[i0+1:i1]
        dark_coadd_file = new_dark_path+root_file+'_darks.fits'
        _ = cr_utils.coadd_bintables(gd_darks, outfile=dark_coadd_file)

        # Scaling
        coadd_drk_phas = get_pha_values_dark(bg_region, dark_coadd_file, (xdopp_min,xdopp_max))
        scale_sci_drk = float(pha_values.size) / coadd_drk_phas.size
        if verbose:
            print("Will use scale factor = {:g}".format(scale_sci_drk))

        # Extract
        spec, pha_ex = extract_dark_spectrum(dark_coadd_file, exposure, obj_tr, segm, defaults['pha_mnx'], plot=show_spec, offs1 = offs1, offs2 = offs2)
        if verbose:
            print("Extracting..")

        # Smooth
        smooth_spec = np.convolve(spec, np.ones((N_smooth,))/N_smooth, mode='same')

        # Generate a dark spectrum file
        outfile = sci_path+root_file+'_{:s}_bkgd.fits'.format(sub_seg)
        tbl = Table()
        tbl['DARK'] = smooth_spec * scale_sci_drk
        tbl.write(outfile, overwrite=True)
        if verbose:
            print("Wrote background spectrum to {:s}".format(outfile))



def dark_calcos_script(dark_files, segm, science_folder):
    """ Generate a simple script for running calcos on the series of dark frames
    Also copy in .spt files from science folder

    Parameters
    ----------
    dark_files : list
    segm : str
    science_folder : str
      path to where science reduction is performed
    """
    # Create dark folder
    try:
        os.mkdir(science_folder+'/darks_{:s}/'.format(segm))
    except OSError: # likely already exists
        pass
    #
    clfile = science_folder+'/darks/calcos_darkscript_{:s}.cl'.format(segm)
    with open(clfile, 'w') as f:
        for dark_file in dark_files:
            f.write('calcos {:s}'.format(dark_file))


'''
def separate_darks(darksfiles, path):
    """

    Parameters
    ----------
    darksfiles - list
       List with darks files - corrtag files?
    path - str
       path

    Returns
    -------

    """
    # hvlevela, hvlevelb values
    import shutil
    # loop through darks
    hva, hvb = [], []
    seg = []
    for ifile in darksfiles:
        hdu = fits.open(ifile)
        ihva, ihvb = hdu[1].header['HVLEVELA'], hdu[1].header['HVLEVELB']
        hva.append(ihva)
        hvb.append(ihvb)
        iseg = hdu[0].header['SEGMENT']
        seg.append(iseg)
        hdu.close()

    # find unique hvlevela, hvlevelb
    uniq_hva=np.unique(hva)
    uniq_hvb=np.unique(hvb)

    # generate folder for each
    for ihva in uniq_hva:
        dirname=path+'a_'+str(ihva)
        try:
            os.stat(dirname)
        except:
            os.mkdir(dirname)

    for ihvb in uniq_hvb:
        dirname = path +'b_' + str(ihvb)
        try:
            os.stat(dirname)
        except:
            os.mkdir(dirname)

    # files names without path
    filesonly = []
    for dstr in darksfiles:
        strnew = dstr.split("/")[-1]
        filesonly.append(strnew)

    # move corrtags into folders
    for i in np.arange(len(darksfiles)):
        if seg[i] == 'FUVA':
            dirname = path + str('a_') + str(hva[i]) + str('/')
        if seg[i] == 'FUVB':
            dirname = path + str('b_') + str(hvb[i]) + str('/')
        # move file in the folder
        shutil.move(darksfiles[i], dirname + filesonly[i])
'''

def find_darks(darksfld, scifile, segm, hvlvl, ndays=90):
    """  Finds darks taken within ndays around the date
         of your science observations,
         with the HVlevel = hvlvl
    Parameters
    ----------
    darksfld : str
      Path to the darks folder
    scifile : str
      Science frame (could be raw, processed, etc.)
      Requires date information in the header
    segm : str
      COS segment -- 'a' or 'b'
    hvlvl : int or str
       HVlevel
    ndays : int, optional
      Number of days to search within to match to darks


    Returns
    -------
    dlist : list
      List of darks within +/- ndays

    """
    # data: date
    hdu = fits.open(scifile)
    head0 = hdu[0].header
    ftime = head0['DATE']

    npdays = np.timedelta64(ndays, 'D')

    # find darks
    darkfiles = glob.glob(darksfld + '*' + segm + '.fits')

    dlist = []
    for ifile in darkfiles:
        hdu = fits.open(ifile)
        head1 = hdu[1].header
        # Query on HVL
        ihva, ihvb = head1['HVLEVELA'], head1['HVLEVELB']
        if segm == 'a':
            if ihva != hvlvl:
                continue
        elif segm == 'b':
            if ihvb != hvlvl:
                continue
        # Query on DATE
        itime = head1['DATE-OBS']
        dd1 = np.datetime64(itime) - np.datetime64(ftime)
        if (np.abs(dd1) < npdays):
            dlist.append(ifile)

    # new array with darks
    return dlist

def setup_for_calcos(darksfld, scifile, segm, **kwargs):
    """  Identify darks to process, generate sub-folder,
    copy over, modify headers, generate calcos script

    Parameters
    ----------
    darksfld : str
      path to darks
    scifile : str
      one corrtag file from a given visit
      includes path to the scifile
    segm : str
    #ndays

    Returns
    -------
    d_subf : str

    """
    from shutil import copyfile as shutilcp
    # Path
    rdxsci_path = scifile[0:scifile.rfind('/')+1]

    # Get HVLVL
    hva, hvb = cr_utils.get_hvlevels([scifile])

    if segm == 'FUVA':
        hvl = hva[0]
        subseg = 'a'
    elif segm == 'FUVB':
        hvl = hvb[0]
        subseg = 'b'

    # Folders
    d_subf = 'darks_'+subseg+'_{:d}'.format(hvl)
    try:
        os.stat(rdxsci_path+d_subf)
    except:
        os.mkdir(rdxsci_path+d_subf)

    # Find darks
    darks = find_darks(darksfld, scifile, subseg, hvl, **kwargs)

    # Copy to folder (create if needed)
    print("Copying dark frames into sub folder: {:s}".format(d_subf))
    new_darks = []
    for darkfile in darks:
        darkname = darkfile.split("/")[-1]
        darkfile2 = rdxsci_path+d_subf+'/'+darkname
        # rawtag file
        shutilcp(darkfile,darkfile2)
        new_darks.append(darkname)
        # spt file
        spt_file = darkfile.replace('rawtag_{:s}.fits'.format(subseg), 'spt.fits')
        sptname = spt_file.split("/")[-1]
        sptfile2 = rdxsci_path+d_subf+'/'+sptname
        shutilcp(spt_file,sptfile2)

    # Edit headers
    print("Editing dark frame headers")
    cr_utils.modify_rawtag_for_calcos(rdxsci_path+d_subf)

    # Generate calcos script
    clfile = rdxsci_path+d_subf+'/'+d_subf+'.cl'
    f = open(clfile, 'w')
    for ifile in new_darks:
        f.write('calcos ' + ifile + '\n') # calcos [dark1]_rawtag_a.fits ...
    print("Wrote calcos script: {:s}".format(clfile))

    # Return sub-folder path
    return d_subf


def clean_after_calcos(path):
    """ Remove unwanted (large) files from the dark sub-folder
    after running calcos.
    Saves the corrtag files and the .cl script

    Parameters
    ----------
    path : str
      Full path to sub-folder for processed dark frames
    """
    all_files=glob.glob(path+'/*.fits')
    corr_files=glob.glob(path+'/*corrtag*.fits')
    cl_files=glob.glob(path+'/*.cl') # Save script
    keep_files = corr_files+cl_files
    for dfile in all_files:
        if dfile not in keep_files:
            os.remove(dfile)
