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
from pypeit.core import coadd as coadd

from linetools.spectra import utils as spltu
import linetools.spectra.xspectrum1d as xspec
from linetools.spectra import utils as ltsu

import elg_cgm.sn as ecsn

from xastropy.xutils import xdebug as xdb


def set_extraction_region(obj_tr, segm, coadd_corrtag_woPHA_file, apert=25., offs1=0., offs2=0., check=False):
    """ Defines extraction region
    Parameters
    ----------
    obj_tr : float, int
      object trace
    segm : str
      segment
    apert : float, int, optional
    coadd_corrtag_woPHA_file : str
      For wavelength info
    offs1 : float, int, optional
    offs2 : float, int, optional
      left and right offsets from the aperture
       could be used for FUVB
    check : bool, optional
      show extraction region
    #ywidth : float, optional

    Returns
    -------
    ex_region : dict

    """
    # Load
    data = Table.read(coadd_corrtag_woPHA_file)
    wave = data['WAVELENGTH'].data
    # Set
    if segm == 'FUVA':
        x1=1200.
        x2=max(wave)
    elif segm == 'FUVB':
        x1=900.
        x2=max(wave)

    ex_region = {}
    ex_region['extraction'] = [x1, x2, obj_tr - apert+offs1, obj_tr + apert+offs2]

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
        plt.ylim(min(yfull[wave > x1]),max(yfull[wave > x1]))
        plt.show()

    # Return
    return ex_region

def coadd_exposures(x1d_files, segm, outfile, bin=None):
    """ Coadd exposures (step 9 of the procedure)
    Parameters
    ----------
    x1d_files : list of str
    segm : str
    outfile : str
      file with coadded exposures
    bin : 2, 3, None
      bin the output spectrum

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
        gddq = (xtbl['DQ'] == 0) & (xtbl['FLUX'] > 0)

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


def combinespectfiles(spfile_a, spfile_b, file_ab):
    """ Coadd two spectra, and write output in a file

    Parameters
    ----------
    spfile_a : str
    spfile_b : str
       .fits files with spectra
    file_ab : str
       output .fits file with combined spectra

    Returns
    -------

    """
    from linetools.spectra import io as tio
    from linetools.spectra import utils as spltu
    file_a = tio.readspec(spfile_a)
    file_b = tio.readspec(spfile_b)
    spliced_sp=spltu.splice_two(file_b, file_a, chk_units=False)
    spliced_sp.write(file_ab)
    #print("Wrote {:s}".format(file_ab))







def flxwave(fname1, fname2s, xlim=(1900, 1950), ylim=None, norm=True, nseg=0):
    """ Plot flux vs wavelength for default and customized data reduction

    Parameters
    ----------
    fname1 : str
      a summ file that is output of the default data reduction
    fname2s : list of str
      summ files that are output ot customized data reductions.
      It is needed that these files have spectra taken with the same optical element and central wavelenght as fname1.
    xlim : float
      displayed wavelength range
    ylim : float
      dislayed flux range
    norm : bool
      normalize fluxes to the same median flux
    nseg : 0,1,2 (for segments A, B, C)
      plot flux for segment nseg

    Returns
    -------

    """

    # default data reduction
    fnm = fname1
    tbl = Table.read(fnm)
    wave1 = tbl['WAVELENGTH'][nseg]
    flx1 = tbl['FLUX'][nseg]
    err1 = tbl['ERROR'][nseg]

    plt.figure(figsize=(6, 5))
    plt.plot(wave1, flx1, color='blue')
    plt.plot(wave1, err1, '--', color='blue', label='previous')

    # customized data reduction
    colors = ['red', 'orange', 'limegreen']
    for i in range(len(fname2s)):
        fname2 = fname2s[i]
        fnm = fname2
        tbl = Table.read(fnm)
        wave2 = tbl['WAVELENGTH'][nseg]
        flx2 = tbl['FLUX'][nseg]
        err2 = tbl['ERROR'][nseg]

        # normalize to median flux
        if norm:
            sc = np.median(flx1) / np.median(flx2)
        else:
            sc = 1
        print(sc)
        flx2 = flx2 * sc
        err2 = err2 * sc

        # plot
        plt.plot(wave2, flx2, color=colors[i])
        plt.plot(wave2, err2, '--', color=colors[i], label='new ' + str(i+1))

    plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)
    plt.xlabel('Wavelength')
    plt.ylabel('Flux')
    plt.legend()
    plt.show()






def findsp(fldpth, verbose=True):
    """ Find all NUV and FUV spectra

    Parameters
    ----------
    fldpth : str
       path to the files

    Returns
    -------
    fuvsp : list of XSpectrum1D objects
      FUV spectra
    nuvsp : list of XSpectrum1D objects
      NUV spectra

    """

    # find all spectra
    sumfiles = glob.glob(fldpth + '*sum.fits')

    fuvsp = []
    nuvsp = []
    for i in range(len(sumfiles)):
        spf = sumfiles[i]
        hdu = fits.open(spf)
        # FUV or NUV
        det = hdu[0].header['DETECTOR']
        hdu.close()

        # read data
        tbl = Table.read(spf)
        for j in range(len(tbl)):
            wave = tbl['WAVELENGTH'][j]
            flx = tbl['FLUX'][j]
            err = tbl['ERROR'][j]

            # append xspectrum
            xsp = xspec.XSpectrum1D.from_tuple((wave, flx, err))
            if det == 'NUV':
                nuvsp.append(xsp)
            elif det == 'FUV':
                fuvsp.append(xsp)
            else:
                print('Not defined opt_elem?')

    if verbose:
        print('Number of NUV and FUV spectra: ', len(nuvsp), len(fuvsp))
    nuvsp = np.asarray(nuvsp)
    fuvsp = np.asarray(fuvsp)

    return nuvsp, fuvsp




def findspsn(spectra,det,minsn=1,verbose=True):
    """ Find all spectra, with median S/N per resolution element greater than minsn.
        Spectra are taken with detector det.

    Parameters
    ----------
    fldpth : str
       path to the files
    det : 'FUV' or 'NUV'
      detector
    minsn : float
      minimum S/N per resolution element, default = 1

    Returns
    -------
    spectrasn : list of XSpectrum1D objects
      spectra with median S/N per resolution element greater than minsn

    """

    spectrasn = []
    for ispec in spectra:
        # S/N
        if det == 'NUV':
            sn = ecsn.snf1(ispec.wavelength,ispec.flux,ispec.sig)
        if det == 'FUV':
            sn = ecsn.snf2(ispec.wavelength,ispec.flux,ispec.sig)
        medsn = np.median(sn)

        # append
        if medsn > minsn:
            spectrasn.append(ispec)
            if verbose:
                print('S/N:',medsn)
        else:
            if verbose:
                print('S/N low ',medsn)

    if verbose:
        print('Number of spectra with S/N > {:f}: {:f}'.format(minsn,len(spectrasn)))

    spectrasn = np.asarray(spectrasn)

    return spectrasn











def coaddsp(spects,plotsp=False,outf=None,overwrite=False):
    """ Coadd spectra

    Parameters
    ----------
    spects : list of XSpectrum1D objects
       list of spectra
    plotsp : bool
      if true, plot the coadded spectrum
    outf : str
      output file

    Returns
    -------
    sptot : XSpectrum1D
      coadded spectrum

    """

    # sort spectra by minimum wavelength
    wvmins = []
    for i in range(len(spects)):
        wvmins.append(spects[i].wvmin.value)
    ii = np.argsort(wvmins)
    spects = spects[ii]

    # Create a list of spectra (outpspects) in which spectra do not overlap in wavelenght with each other
    # Coadd spectra whose wavelengths overlap (spects1).
    outpspects = []
    spects1 = [spects[0]] # a temporary list in which every next spectrum overlaps in wavelenght with the previous one
    wvmax1 = spects[0].wvmax.value # max wavelenght in the previous spectrum
    for i in range(len(spects) - 1):
        sp = spects[i + 1]
        if sp.wvmin.value < wvmax1:
            spects1.append(sp)
        else:
            mspec = ltsu.collate(spects1)
            coaddsp = coadd.coadd_spectra(mspec)
            outpspects.append(coaddsp)
            #
            spects1 = [sp]
        wvmax1 = sp.wvmax.value

    mspec = ltsu.collate(spects1)
    coaddsp = coadd.coadd_spectra(mspec)
    outpspects.append(coaddsp)

    # coadd outpspects spectra (which do not overlap in wavelenght) by using spltu.splice_two
    sptot = outpspects[0]
    for i in range(len(outpspects) - 1):
        sptot = spltu.splice_two(sptot, outpspects[i + 1], chk_units=False)

    # plot sptot
    if plotsp:
        plt.figure(figsize=(17,4))
        plt.plot(sptot.wavelength.value, sptot.flux.value, color='black')
        plt.plot(sptot.wavelength.value, sptot.sig.value, color='blue')
        plt.plot([sptot.wvmin.value, sptot.wvmax.value], [0, 0], '--',color='lightgray')
        plt.xlabel('Wavelength')
        plt.ylabel('Flux')

    # write to file
    if outf is not None:
        sptot.write(outf,clobber=overwrite)

    return sptot


def medsn(ispec, det):
    """ Find S/N per resolution element for a given spectrum

    Parameters
    ----------
    ispec : XSpectrum1D object
      spectrum
    det : 'NUV' or 'FUV'

    Returns
    -------
    medsn : float
      median S/N per resolution element

    """

    # S/N
    if det == 'NUV':
        sn = ecsn.snf1(ispec.wavelength, ispec.flux, ispec.sig)
    elif det == 'FUV':
        sn = ecsn.snf2(ispec.wavelength, ispec.flux, ispec.sig)
    else:
        print('Det options are only NUV and FUV')
        pdb.set_trace()

    medsn = np.median(sn)

    return medsn




def smoothsp(spects, det, snmin):
    """ Smooth noisy spectra, and return smoothed spectra with S/N per pixel higher than snmin

    Parameters
    ----------
    spects : list of XSpectrum1D objects
      spectra
    det : 'NUV' or 'FUV'
      'NUV' or 'FUV' detector
    snmin : float
      median S/N per pixel needs to be > snmin

    Returns
    -------
    smspects : list of XSpectrum1D objects
      smoothed spectra with S/N > snmin

    """

    outpspec = []
    print('Total number of spectra: ',len(spects))

    # For each spectrum find median S/N, and boxcar smooth over ism pixels
    # Return only spectra that have S/N > snmin
    for ispec in spects:
        print(ispec)

        # S/N
        imedsn = medsn(ispec, det)
        print(imedsn)

        # sm - lists of the number of pixels to smooth over
        if det == 'NUV':
            sm = [2,3]
        if det == 'FUV':
            sm = [2,4,8,12,16,20]
            # or define any other list with sm, e.g.:
            # sm = [2,4,6,8,10,12,14,16,18,20]

        # smooth over ism pixels, until the spectrum has S/N > snmin
        if imedsn > snmin:
            outpspec.append(ispec)
        else:
            for ism in sm:
                ispecsm = ispec.box_smooth(ism)
                imedsn = medsn(ispecsm, det)
                print(ism,imedsn)
                if imedsn > snmin:
                    outpspec.append(ispecsm)
                    print('Smoothing with {:f} pixels'.format(ism))
                    break

    return outpspec

 