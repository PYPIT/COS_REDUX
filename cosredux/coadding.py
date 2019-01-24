""" Module for coadding spectra, and comparing customized and default data reduction
"""

from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
import glob
import pdb

from matplotlib import pyplot as plt

from astropy.table import Table
from astropy.io import fits
from pypeit.core import coadd as coadd

from linetools.spectra import utils as spltu
import linetools.spectra.xspectrum1d as xspec
from linetools.spectra import utils as ltsu

from xastropy.xutils import xdebug as xdb


def flxwave(fname1, fname2s, xlim=(1900, 1950), ylim=None, norm=True, nseg=0, figsize=(6, 5)):
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

    plt.figure(figsize=figsize)
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






def findsp(fldpth, verbose=True, emax=None):
    """ Find all NUV and FUV spectra

    Parameters
    ----------
    fldpth : str
       path to the files
    emax : int
      read extensions 0, 1, ..., emax-1 of sum.fits files (corresponding to segments A,B,C)

    Returns
    -------
    fuvsp : array of XSpectrum1D objects
      FUV spectra
    nuvsp : array of XSpectrum1D objects
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
        jmax = len(tbl)
        if (emax is not None):
            if (emax < len(tbl)):
                jmax = emax
        for j in range(jmax):
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
                raise IOError("Detector not well defined in ", spf)

    if verbose:
        print('Number of NUV and FUV spectra: ', len(nuvsp), len(fuvsp))
    nuvsp = np.asarray(nuvsp)
    fuvsp = np.asarray(fuvsp)

    return nuvsp, fuvsp




def snf1(wave,flux,error):
    """ Calculate S/N per NUV resolution element = 3 pixels

    Parameters
    ----------
    wave : list of floats
      wavelength
    flux : list of floats
      flux
    error : list of floats
      error in flux
    Returns
    -------
    sn : list of floats
      S/N per resolution element, for each pixel

    """

    sn = []
    for i in np.arange(len(wave)):
        if flux[i] == 0:
            sn.append(0)
        else:
            s = 0
            n = 0
            for j in [-1,0,1]: # resolution element = 3 pixels
                ij = i+j
                if ij < 0:
                    ij = 0
                if ij >= len(flux):
                    ij = len(flux)-1
                s = s + flux[ij]
                n = n + error[ij]**2
            isn = s/(n**0.5)
            sn.append(isn)
    return sn


def snf2(wave,flux,error, R = 19000.):
    """ Calculate S/N per resolution R (for FUV)

    Parameters
    ----------
    wave : list of floats
      wavelength
    flux : list of floats
      flux
    error : list of floats
      error in flux
    Returns
    -------
    sn : list of floats
      S/N per  R = 19000, for each pixel

    """

    sn = []
    for i in np.arange(len(wave)):
        if flux[i] == 0:
            sn.append(0)
        else:
            s = 0
            n = 0
            dlam = wave[i]/R
            lam1 = wave[i] - dlam/2
            lam2 = wave[i] + dlam/2
            pix1 = np.argmin(abs(wave - lam1))
            pix2 = np.argmin(abs(wave - lam2))
            for j in (np.arange(pix2-pix1+1)+pix1):
                s = s + flux[j]
                n = n + error[j]**2
            isn = s/(n**0.5)
            sn.append(isn)
    return sn





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
    spectrasn : array of XSpectrum1D objects
      spectra with median S/N per resolution element greater than minsn

    """

    spectrasn = []
    for ispec in spectra:
        # S/N
        if det == 'NUV':
            sn = snf1(ispec.wavelength,ispec.flux,ispec.sig)
        elif det == 'FUV':
            sn = snf2(ispec.wavelength,ispec.flux,ispec.sig)
        else:
            raise IOError("det (Detector) could be 'NUV' or 'FUV' only ")
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
    for sp in spects[1:]:
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
    for outsp in outpspects[1:]:
        sptot = spltu.splice_two(sptot, outsp, chk_units=False)

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
        #pdb.set_trace()
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
        sn = snf1(ispec.wavelength, ispec.flux, ispec.sig)
    elif det == 'FUV':
        sn = snf2(ispec.wavelength, ispec.flux, ispec.sig)
    else:
        print('Det options are only NUV and FUV')
        pdb.set_trace()

    medsn = np.median(sn)

    return medsn




def smoothsp(spects, det, snmin, outf=None):
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
        elif det == 'FUV':
            sm = [2,4,8,12,16,20]
            # or define any other list with sm, e.g.:
            # sm = [2,4,6,8,10,12,14,16,18,20]
        else:
            raise IOError("Detector could be 'NUV' or 'FUV' only ")

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
            if imedsn <= snmin:
                outpspec.append(ispecsm)
                print('Smoothing with {:f} pixels. S/N still low.'.format(ism))

    if outf is not None:
        outpspec[0].write(outf)

    return outpspec