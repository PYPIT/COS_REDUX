""" Trace the source in Data file
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
import pdb
import glob
from astropy.table import Table
from astropy.io import fits
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

from xastropy.xutils import xdebug as xdb

import pypeit.utils as pyutils

from cosredux import utils
from cosredux import io as cr_io


def crude_histogram(yfull, ymin=300, ymax=700, ytl=550, pk_window=4.5, verbose=False):
    """ Derive a trace given an input COS image

    Parameters
    ----------
    yfull : ndarray
      y pixel values of the counts
        (typically taken from YFULL column of data file)
    ymin : int, optional
      Minimum search window for the trace
    ymax : int, optional
      Minimum search window for the trace
    ytl : int, optinal
      Estimated y-position to separate arc from object

    Returns
    -------
    obj_y : float
      Estimated y position of the object
    arc_y : float
      Estimated y position of the arc

    """

    nh=ymax-ymin    # bins of width = 1

    # Generate histogram
    ycen=np.arange(nh)+ymin+0.5
    yedges=np.arange(nh+1)+ymin
    yhist, edges = np.histogram(yfull, bins=yedges)

    # Find Object
    iobj = np.where(ycen < ytl)[0]
    ihist_pk = np.argmax(yhist[iobj])
    obj_y_guess = ycen[iobj[ihist_pk]]


    # Find arc
    iarc = np.where(ycen > ytl)[0]
    ihist_pk = np.argmax(yhist[iarc])
    arc_y_guess = ycen[iarc[ihist_pk]]
    if verbose:
        print("Crude: Obj={:g}, Arc={:g}".format(obj_y_guess,arc_y_guess))

    # Refine
    obj_y = refine_peak(yfull, obj_y_guess)
    arc_y = refine_peak(yfull, arc_y_guess)
    if verbose:
        print("Refined: Obj={:g}, Arc={:g}".format(obj_y,arc_y))

    # Return
    return obj_y, arc_y



def crude_histogram_multi(yfull, ypeaks, verbose=False, pk_window=5.):
    """ Derive multiple traces given an input COS image and expected peaks positions

    Parameters
    ----------
    yfull : ndarray
      y pixel values of the counts
        (typically taken from YFULL column of data file)
    ypeaks : int array
      expected approximate positions of peaks

    Returns
    -------
    obj_y : float array
      Estimated y positions of the object(s)
    """

    ymax = np.max(yfull)
    ymin = np.min(yfull)
    nh=ymax-ymin    # bins of width = 1

    # Generate histogram
    ycen=np.arange(nh)+ymin+0.5
    yedges=np.arange(nh+1)+ymin
    yhist, edges = np.histogram(yfull, bins=yedges)

    # Find Objects
    ypeaks = np.asarray(ypeaks)
    idx = np.argsort(ypeaks)
    ypeaks = ypeaks[idx]
    crudeobj = []
    refinedobj = []

    for i in range(len(ypeaks)):
        if i == 0:
            yb1 = ymin
        else:
            yb1 = (ypeaks[i-1]+ypeaks[i])/2
        if i == (len(ypeaks)-1):
            yb2 = ymax
        else:
            yb2 = (ypeaks[i+1]+ypeaks[i])/2

        ###

        yb1 = ypeaks[i] - pk_window
        yb2 = ypeaks[i] + pk_window


        iobj = np.where((ycen > yb1) & (ycen < yb2))[0]
        ihist_pk = np.argmax(yhist[iobj])
        obj_y_guess = ycen[iobj[ihist_pk]]
        crudeobj.append(obj_y_guess)
        # Refine
        obj_y = refine_peak(yfull, obj_y_guess)
        refinedobj.append(obj_y)

    if verbose:
        print("Crude Obj: ", crudeobj[idx])
        print("Refined: Obj: ", refinedobj[idx])

    # Return
    refinedobj = np.asarray(refinedobj)

    return refinedobj[idx]




def refine_peak(yfull, y_guess, pk_window=5., per_lim=(0.25,0.75)):
    """  Refines peak: calculates average value of yfull between percentiles per_lim
    Parameters
    ----------
    yfull : ndarray
    y_guess : float, int
      peak guess
    pk_window : float, int, optional
    per_lim : tuple, optional
      percentiles limits

    Returns
    -------
    pk : float
      refined peak value

    """
    # Cut to y values near the peak
    cut = np.abs(yfull-y_guess) < pk_window
    ycut = yfull[cut]

    # Sort
    ycut.sort()

    # Cumulative sum
    ysum = np.cumsum(ycut)

    left = np.argmin(np.abs(per_lim[0]-ysum/ysum[-1]))
    right = np.argmin(np.abs(per_lim[1]-ysum/ysum[-1]))

    # Average
    pk = (ycut[left]+ycut[right])/2.

    # Return
    return pk



def show_traces(wave, yfull, obj_y, arc_y, segm, plottype='all', outfilpl=None):
    """ Shows traces of object and lamp, or only of object
    Parameters
    ----------
    wave : ndarray
      WAVELENGTH
    yfull : ndarray
      YFULL
    obj_y : float, int
    arc_y : float, int
      traces of object and lamp
    segm : str
      segment
    plottype : str
      'all'  - display both object and lamp trace
      'obj'  - y axis: display only object
      'objb' - x,y axes: display only object
    outfilpl : str
      if not None: save figure in this file

    Returns
    -------

    """
    from matplotlib import pyplot as plt

    plt.clf()
    # figsize
    if plottype == 'all':
        plt.figure(figsize=(15,15))
    elif plottype == 'obj':
        plt.figure(figsize=(15, 4))
    elif plottype == 'objb':
        plt.figure(figsize=(8,4))
    ax = plt.gca()
    # plot
    ax.scatter(wave,yfull,s=1)
    wvmin = np.maximum(np.min(wave), 1200)
    if segm == 'FUVB':
        wvmin = np.min(wave)
    wvmax = np.max(wave)

    ax.set_xlabel('Wavelength')
    ax.set_ylabel('YFULL')
    ax.set_xlim(wvmin, wvmax)
    # Y range
    gdp = (wave > wvmin) & (wave < wvmax)
    if segm == 'FUVB':
        gdp = (wave > 200) & (wave < wvmax)
    ymin = np.min(yfull[gdp])
    ymax = np.max(yfull[gdp])
    # plot traces? display only object?
    if plottype == 'all':
        ax.set_ylim(ymin, ymax)
        ax.plot([wvmin, wvmax], [obj_y, obj_y], color='blue')
        ax.plot([wvmin, wvmax], [arc_y, arc_y], color='red')
    elif plottype == 'obj':
        ax.set_ylim(obj_y - 40, obj_y + 40)
        if segm == 'FUVB':
            ax.set_xlim(0, np.max(wave))
            ax.set_ylim(obj_y - 40 -10, obj_y + 40 -10)
    elif plottype == 'objb':
        ax.set_ylim(obj_y - 40 - 10, obj_y + 40 - 10)
        ax.set_xlim(900, np.max(wave))

    # save figure
    if outfilpl is not None:
        plt.savefig(outfilpl)
    plt.show()



def traces(filename, calib_path, segment, row_dict=None, LP='LP3',
           outfil=None, clobber=False, show=False, calcos_version='v2', plottype='all', outfilpl=None):
    """ Finds traces of the object and the lamp, and updates these values in the calibration file
    Parameters
    ----------
    filename : str
      File for which we want to find the trace. E. g. it could be corrtag file.
    calib_path : str
      Path to calibration files
    segment : str
    row_dict : dict, optional
      Dict that describes which row(s) we want to modify
      Defaulted to G140L, CENWAVE=1280
    LP
    #ymin : float, optional
    #  define search window for trace. Default values for FUVA on LP3
    #ymax : float, optional
    #  define search window for trace. Default values for FUVA on LP3
    #ytl : float, optional
    #outfil : str, optional
    clobber : bool, optional
      overwrite object trace value in the calibration file
    show : bool, optional
      show traces?
    calcos_version : str, optional
    plottype : str, optional
      display both object and lamp, or only object
    outfilpl : str, optional
      output file for the plot

    Returns
    -------
    obj_y : float
      Estimated y position of the object
    arc_y : float
      Estimated y position of the arc
    """
    assert LP == 'LP3'  # Only one coded
    if segment == 'FUVA':
        ymin, ymax, ytl = 300, 700, 550
    elif segment == 'FUVB':
        # ymin, ymax, ytl = 360, 760, 550
        ymin, ymax, ytl = 500, 650, 550
    else:
        raise IOError("Not ready for this segment")
    # Prepare to modify table
    if row_dict is None:
        row_dict = dict(OPT_ELEM='G140L', CENWAVE=1280, APERTURE='PSA')
    # Add segment
    row_dict['SEGMENT'] = segment
    # FITS
    data = Table.read(filename)
    wave = data['WAVELENGTH']
    yfull = data['YFULL']
    obj_y, arc_y = crude_histogram(yfull, ymin=ymin, ymax=ymax, ytl=ytl)
    if show:
        show_traces(wave, yfull, obj_y, arc_y, segment, plottype=plottype, outfilpl=outfilpl)
    # Update trace value
    if calcos_version == 'v2':
        filecal = calib_path+'x6q17586l_1dx.fits'  # WHEN RUNNING calcos 2
    else:
        raise IOError("Not ready for another calcos version")
    # Modify
    #pdb.set_trace()
    utils.modify_table_value(filecal, 'B_SPEC', row_dict, obj_y, outfil=filecal, clobber=clobber)
    print('Updated trace for segment={:s} in {:s}'.format(segment, filecal))

    # Write to hard-drive
    outfile = filename.replace('.fits', '_traces.json')
    cr_io.write_traces(obj_y, arc_y, outfile)

    return obj_y, arc_y


def traceshist(file_tr,traces_n,plottype='all',ymin=300,ymax=700, offs1ob = 0., offs2ob = 0., apert=25):
    """ Displays histogram of yfull

    Parameters
    ----------
    file_tr : str
      fits file
    traces_n : float, int
      traces in the fits file
    #segm : str
    #  segment
    plottype : str
      'all' - histogram for all wave
      'objb' - histogram only for wave > 1000 (approx. for object trace location on segment B)
    ymin : float, int, optional
    ymax : float, int, optional
      histogram for yfull between ymin and ymax
    offs1ob : float, int, optional
    offs2ob : float, int, optional
      offset from aperture
    apert : float, int, optional
      aperture

    Returns
    -------

    """
    from matplotlib import pyplot as plt

    data = Table.read(file_tr)
    yfull = data['YFULL']
    wave = data['WAVELENGTH']
    nh=ymax-ymin    # bins of width = 1

    yfullhist = yfull
    #if segm == 'FUVB':
    if plottype == 'objb':
        yfullhist = yfull[wave > 1000.]

    # Generate histogram
    ycen=np.arange(nh)+ymin+0.5
    yedges=np.arange(nh+1)+ymin
    yhist, edges = np.histogram(yfullhist, bins=yedges)

    plt.plot(ycen,yhist)
    plt.plot([traces_n[0],traces_n[0]],[0,max(yhist)],'r')
    plt.plot([traces_n[1],traces_n[1]],[0,max(yhist)],'b')
    plt.plot([traces_n[0]+offs1ob-apert,traces_n[0]+offs1ob-apert],[0,max(yhist)],'r--')
    plt.plot([traces_n[0]+offs2ob+apert,traces_n[0]+offs2ob+apert],[0,max(yhist)],'r--')
    plt.show()








def findpeaks(fldpth, fnamecorrs, verbose=True, printorig=True, xmin=0, xmax=None, height = None,pk_window=5,apertures=None):
    """ Find new traces (bspecs and slopes; YFULL = bspec + slope * XFULL)
        XFULL and YFULL are combined from all fnamecorrs.
        First, we fit a line YFULL = bspec + slope * XFULL, and find new slopes. For FUV we ignore this step,
        and slopes are zero.
        Next, we correct YFULL for the slopes, and find new bspec using cosredux.trace.refine_peak.
        Using the code to change WCA traces might not be reliable in some cases (for default input parameters).
        It is recomended using plottraces and plothist to check the results.

    Parameters
    ----------
    fldpth : str
       path to the files
    fnamecorrs : list of str
      corrtag files. These files need to have the same settings (optical element, central wavelength, FUVA/FUVB).
    printorig :bool
      print original peaks and slopes
    xmin : float
      minimum XFULL that is taken into account when finding traces
    xmax : float
      maximum XFULL that is taken into account when finding traces
    height : float
      height (YFULL range) used when finding traces
    apertures : list of str
      apertures

    Returns
    -------
    allnewypeaks : list of lists of floats
      new yspecs
    allnewslopes : list of lists of floats
      new slopes
    ypeaks : list of floats
      previous yspecs
    slopes : list of floats
      previous slopes

    """

    # read fnamecorrs header: cenwave, opt_elem, detector, LP2_1dx_file
    hdu = fits.open(fnamecorrs[0])
    hd = hdu[0].header
    hdu.close()
    cenwave = hd['CENWAVE']
    opt_elem = hd['OPT_ELEM']

    if hd['DETECTOR'] == 'NUV':
        segments = ['NUVA', 'NUVB', 'NUVC']
        LP2_1dx_file = fldpth + hd['XTRACTAB'][5:-5] + '_copy1.fits'
        if apertures is None:
            apertures = ['PSA', 'WCA']
    elif hd['DETECTOR'] == 'FUV':
        segments = [hd['SEGMENT']]
        LP2_1dx_file = fldpth + hd['TWOZXTAB'][5:-5] + '_copy1.fits'
        if apertures is None:
            apertures = ['PSA']
    else:
        raise IOError("Detector not well defined in ",fnamecorrs[0])

    # LP2_1dx_file: ypeaks, slopes, heights  (should be LP3, etc. This is just notation.)
    hdu = fits.open(LP2_1dx_file)
    lp2 = Table(hdu[1].data)
    #
    ypeaks = []
    slopes = []
    heights = []
    for aperture in apertures:
        for segment in segments:
            irow = np.where((lp2['CENWAVE'] == cenwave) &
                            (lp2['APERTURE'] == aperture) &
                            (lp2['SEGMENT'] == segment) &
                            (lp2['OPT_ELEM'] == opt_elem))[0]
            if len(irow) != 1:
                print('irow error')
                pdb.set_trace()
            else:
                irow = irow[0]
            if verbose:
                if hd['DETECTOR'] == 'NUV':
                    print(lp2[irow]['B_SPEC'], lp2[irow]['SLOPE'])
                if hd['DETECTOR'] == 'FUV':
                    print(lp2[irow]['B_SPEC'])
            ypeaks.append(lp2[irow]['B_SPEC'])
            if hd['DETECTOR'] == 'NUV':
                slopes.append(lp2[irow]['SLOPE'])
            else:
                slopes.append(0)
            if height is not None:
                heights.append(height)
            else:
                heights.append(lp2[irow]['HEIGHT'])
    if printorig:
        print('Original values:')
        print('slopes: ', slopes)
        print('bspecs: ', ypeaks)

    ####################################################
    # read and combine xfull, yfull for corrtag files
    xfull = []
    yfull = []
    for fnamecorr in fnamecorrs:
        hdu = fits.open(fnamecorr)
        data1 = hdu[1].data
        if len(data1['XFULL']) < 1:
            print(' No data for ', fnamecorr)
        else:
            xfull = xfull + list(data1['XFULL'])
            yfull = yfull + list(data1['YFULL'])
    xfull = np.asarray(xfull)
    yfull = np.asarray(yfull)

    if len(yfull) < 1:
        print(' No data available for', fldpth)
        return

    ##################################################

    # find new slopes
    newslopes1 = []
    newypeaks1 = []
    for i in range(len(slopes)):
        if xmax is None:
            xmax = np.max(xfull)
        k = np.where((abs(yfull - slopes[i] * xfull - ypeaks[i]) < heights[i] / 2) &
                     (xfull >= xmin) & (xfull <= xmax))[0]
        yfull1 = yfull[k]
        xfull1 = xfull[k]
        # fit
        if hd['DETECTOR'] == 'NUV':
            x0, x1 = pyutils.robust_polyfit(xfull1, yfull1, 1)[1]
            newslopes1.append(x1)
            # assuming that central point in the trace did not change, i.e.:
            #     newslopes1 * (np.max(xfull)+np.min(xfull))*0.5 + newypeaks1 = const
            newypeaks1.append( (slopes[i] - x1) * (np.max(xfull)+np.min(xfull))*0.5 + ypeaks[i] )

        else:
            newypeaks1.append(ypeaks[i])
            newslopes1.append(slopes[i])
    newypeaks1 = np.asarray(newypeaks1)
    newslopes1 = np.asarray(newslopes1)

    # find new bspecs
    dyfulls = []
    for i in range(len(newslopes1)):
        dyfulls.append(abs(yfull - newslopes1[i] * xfull - newypeaks1[i]))
    jj = np.argmin(dyfulls, axis=0)
    #
    inewslopes = newslopes1[jj]
    yfull2 = yfull - inewslopes * xfull # correct YFULL for the slope
    newypeaks = crude_histogram_multi(yfull2, newypeaks1, pk_window=pk_window)
    newslopes = newslopes1

    #
    if verbose:
        print(newypeaks)
        print(ypeaks)
        print(newslopes)
        print(slopes)

    return newypeaks, newslopes, ypeaks, slopes




def modifyxtractab(fldpth, fnamecorrs1, new_ebh=None, new_slopes=None, new_bspecs=None, new_loout=None, new_upout=None,
                   new_loinn=None, new_upinn=None, verbose=True,overwrite=True,apertures=None):
    """ Modify XTRACTAB traces.

    Parameters
    ----------
    fldpth : str
       path to the files
    fnamecorrs1 : list of str
      corrtag files. These files need to have the same settings (optical element, central wavelength, FUVA/FUVB).
    new_ebh : float
      new height
    new_slopes : list of floats
      new slopes
    new_bspecs : list of floats
      new bspecs
    overwrite : bool
      True - save changes, False - only print what would change in the XTRACTAB file
    apertures : list of str
      apertures for which traces will be modified
    new_loout : float
      new LOWER_OUTER in TWOZXTAB file
    new_upout : float
      new UPPER_OUTER in TWOZXTAB file
    new_loinn : float
      new LOWER_INNER in TWOZXTAB file
    new_upinn : float
      new UPPER_INNER in TWOZXTAB file

    Returns
    -------

    """

    # read fnamecorrs header: cenwave, opt_elem, detector, LP2_1dx_file
    fnamecorr = fnamecorrs1[0]
    hdu = fits.open(fnamecorr)
    hd = hdu[0].header
    cenwave = hd['CENWAVE']
    opt_elem = hd['OPT_ELEM']
    hdu.close()

    if hd['DETECTOR'] == 'NUV':
        segments = ['NUVA', 'NUVB', 'NUVC']
        LP2_1dx_file = fldpth + hd['XTRACTAB'][5:]
        if apertures is None:
            apertures = ['PSA', 'WCA']
    elif hd['DETECTOR'] == 'FUV':
        segments = [hd['SEGMENT']]
        LP2_1dx_file = fldpth + hd['TWOZXTAB'][5:]
        if apertures is None:
            apertures = ['PSA']
    else:
        raise IOError("Detector not well defined in ", fnamecorrs1[0])
    if verbose:
        print(LP2_1dx_file)

    # open and modify LP2_1dx_file
    with fits.open(LP2_1dx_file) as f:
        lp2 = f[1].data

        i = 0
        for aperture in apertures:
            for segment in segments:
            # the order for apertures and segments is the same as the output of findpeaks()
                irow = np.where((lp2['CENWAVE'] == cenwave) &
                                (lp2['APERTURE'] == aperture) &
                                (lp2['SEGMENT'] == segment) &
                                (lp2['OPT_ELEM'] == opt_elem))[0]
                if len(irow) != 1:
                    print('irow error')
                    pdb.set_trace()
                else:
                    irow = irow[0]
                if verbose:
                    if new_bspecs is not None:
                        print(i, 'bspec', lp2[irow]['B_SPEC'], new_bspecs[i])
                    if (hd['DETECTOR'] == 'NUV') & (new_slopes is not None):
                        print(i, 'slope', lp2[irow]['SLOPE'], new_slopes[i])

                # modify bspec, slope and height
                if new_bspecs is not None:
                    lp2[irow]['B_SPEC'] = new_bspecs[i]
                if (hd['DETECTOR'] == 'NUV') & (new_slopes is not None):
                    lp2[irow]['SLOPE'] = new_slopes[i]

                # height
                if (aperture == 'PSA') & (new_ebh is not None):
                    if verbose:
                        print('ebh', lp2[irow]['HEIGHT'], new_ebh)
                    lp2[irow]['HEIGHT'] = new_ebh

                # twozone
                if new_loout is not None:
                    if verbose:
                        print('lo_out', lp2[irow]['LOWER_OUTER'], new_loout)
                    lp2[irow]['LOWER_OUTER'] = new_loout
                if new_upout is not None:
                    if verbose:
                        print('up_out', lp2[irow]['UPPER_OUTER'], new_upout)
                    lp2[irow]['UPPER_OUTER'] = new_upout
                if new_loinn is not None:
                    if verbose:
                        print('lo_inn', lp2[irow]['LOWER_INNER'], new_loinn)
                    lp2[irow]['LOWER_INNER'] = new_loinn
                if new_upinn is not None:
                    if verbose:
                        print('up_inn', lp2[irow]['UPPER_INNER'], new_upinn)
                    lp2[irow]['UPPER_INNER'] = new_upinn

                i = i + 1

        if overwrite:
            f.writeto(LP2_1dx_file, overwrite=overwrite)
            print('Wrote new traces in ',LP2_1dx_file)

        if verbose:
            print(' ')



def modifyspoff(fldpth, seg, new_bspec, previous_bspec, verbose=True, overwrite=True):
    """ Modify SP_SET values in RAWTAG files

    Parameters
    ----------
    fldpth : str
       path to the files
    seg : str
      segment; options: 'a', 'b'
    new_bspec : float
      new bspec that we want to use
    previous_bspec: float
      default bspec

    Returns
    -------

    """

    # segment
    if seg == 'a':
        Seg = 'A'
    if seg == 'b':
        Seg = 'B'

    # modify rawtag files
    raws = glob.glob(fldpth + '*rawtag_' + seg + '.fits')
    for raw in raws:
        with fits.open(raw) as f:
            if verbose:
                print('Replacing', f[1].header['SP_OFF_' + Seg], ' with ', previous_bspec - new_bspec, ' in ', raw)
            f[1].header['SP_SET_' + Seg] = previous_bspec - new_bspec
            if overwrite:
                f.writeto(raw, overwrite=overwrite)






def plottraces(fldpth, corrtag, newypeaks=None, newslopes=None, verbose=False, dylim=None, showorigtraces=True,
               apertures=None):
    """ Plot YFULL vs XFULL, and show original and new traces.
     The traces are shown for object and lamp, and for all segments. For each case, a separate figure is shown.

    Parameters
    ----------
    fldpth : str
       path to the files
    corrtag : list of str
      list of corrtag files. (XFULL, YFULL) data are combined for all of these files.
      These files need to have the same settings (optical element, central wavelength, FUVA/FUVB).
    newypeaks : list of floats
      bspecs for the new traces
    newslopes : list of floats
      slopes for the new traces
    dylim : float
      YFULL range to show around traces. ylim = (bspec - dylim, bspec + dylim)
    showorigtraces : bool
      If true, show original traces
    apertures : list of str
      apertures

    Returns
    -------

    """

    # read and combine xfull, yfull for corrtag files
    xfull = []
    yfull = []
    for icorrtag in corrtag:
        hdu = fits.open(icorrtag)
        data1 = hdu[1].data
        xfull = xfull + list(data1['XFULL'])
        yfull = yfull + list(data1['YFULL'])
    xfull = np.asarray(xfull)
    yfull = np.asarray(yfull)

    # read corrtag file header data
    hdu = fits.open(corrtag[0])
    hd = hdu[0].header
    hdu.close()
    cenwave = hd['CENWAVE']
    opt_elem = hd['OPT_ELEM']
    if hd['DETECTOR'] == 'NUV':
        segments = ['NUVA', 'NUVB', 'NUVC']
        LP2_1dx_file = fldpth + hd['XTRACTAB'][5:-5] + '_copy1.fits'
        if apertures is None:
            apertures = ['PSA', 'WCA']
    elif hd['DETECTOR'] == 'FUV':
        segments = [hd['SEGMENT']]
        LP2_1dx_file = fldpth + hd['TWOZXTAB'][5:-5] + '_copy1.fits'
        if apertures is None:
            apertures = ['PSA']
    else:
        raise IOError("Detector not well defined in ", corrtag[0])

    # read ypeaks and slopes from a copy of the LP2_1dx_file
    hdu = fits.open(LP2_1dx_file)
    lp2 = Table(hdu[1].data)
    #
    ypeaks = []
    slopes = []
    for aperture in apertures:
        for segment in segments:
            irow = np.where((lp2['CENWAVE'] == cenwave) &
                            (lp2['APERTURE'] == aperture) &
                            (lp2['SEGMENT'] == segment) &
                            (lp2['OPT_ELEM'] == opt_elem))[0]
            if len(irow) != 1:
                print('irow error')
                pdb.set_trace()
            else:
                irow = irow[0]

            # print ypeaks and slopes (optional)
            if verbose:
                if hd['DETECTOR'] == 'NUV':
                    print(lp2[irow]['B_SPEC'], lp2[irow]['SLOPE'])
                if hd['DETECTOR'] == 'FUV':
                    print(lp2[irow]['B_SPEC'])

            # append ypeaks and slopes
            ypeaks.append(lp2[irow]['B_SPEC'])
            if hd['DETECTOR'] == 'NUV':
                slopes.append(lp2[irow]['SLOPE'])
            else:
                slopes.append(0)

    # plot yfull vs xfull, and show traces
    xmnmx = np.asarray([np.min(xfull), np.max(xfull)])
    slopes = np.asarray(slopes)
    for i in range(len(ypeaks)):
        k = np.where(abs(yfull - ypeaks[i]) <= dylim*2)[0]
        plt.figure(figsize=(14, 4))
        plt.scatter(xfull[k], yfull[k], s=0.005)
        # traces
        if showorigtraces:
            plt.plot(xmnmx, ypeaks[i] + slopes[i] * xmnmx, color='red')
        if newypeaks is not None:
            plt.plot(xmnmx, newypeaks[i] + newslopes[i] * xmnmx, color='orange')
        #
        plt.xlabel('XFULL')
        plt.ylabel('YFULL')
        if dylim is not None:
            plt.ylim(newypeaks[i]-dylim  +slopes[i] * np.average(xmnmx) ,newypeaks[i]+dylim +slopes[i] * np.average(xmnmx))
        plt.show()


def plothist(fldpth, corrtag, newypeaks=None, slopes=None, newheights=None, iheight=57,
             verbose=False, nhres=10, fitgauss=False, percs=None, showorigperc=False,apertures=None):
    """ Displays histogram of YFULL corrected for the slope, and shows original and new traces.
    For each case (object or lamp; different segments), a separate figure is shown.

    Parameters
    ----------
    fldpth : str
      path to the files
    corrtag : list of str
      list of corrtag files. (XFULL, YFULL) data are combined for all of these files.
      These files need to have the same settings (optical element, central wavelength, FUVA/FUVB).
    newypeaks : list of floats
      bspecs for the new traces
    newslopes : list of floats
      slopes for the new traces
    newheights : list of floats
      list of new heights to show (assumed to be the same for all cases). Could contain any number of elements.
    iheight : float
      height used to calculate the histogram
    nhres : float
      histogram resolution. Bin width is 1/nhres
    fitgauss : bool
      if True, fit a Gaussian to the data
    percs : list of floats
      list of percentiles to show; used for FUV twozone.
      It is needed first to define iheight to be equal to the height from the twozxtab
    showorigperc : bool
      Show original outer region boundaries (from twozxtab file)
    apertures : list of str
      apertures

    Returns
    -------

    """

    # read corrtag file header data
    hdu = fits.open(corrtag[0])
    hd = hdu[0].header
    hdu.close()
    #
    cenwave = hd['CENWAVE']
    opt_elem = hd['OPT_ELEM']
    print(corrtag[0], opt_elem, cenwave)
    #
    if hd['DETECTOR'] == 'NUV':
        segments = ['NUVA', 'NUVB', 'NUVC']
        LP2_1dx_file = fldpth + hd['XTRACTAB'][5:-5] + '_copy1.fits'
        if apertures is None:
            apertures = ['PSA', 'WCA']
    elif hd['DETECTOR'] == 'FUV':
        segments = [hd['SEGMENT']]
        LP2_1dx_file = fldpth + hd['TWOZXTAB'][5:-5] + '_copy1.fits'  # '_copy1.fits'
        if apertures is None:
            apertures = ['PSA']
    else:
        raise IOError("Detector not well defined in ", corrtag[0])

    if fitgauss:
        apertures = ['PSA']

    # read ypeaks and heights from a copy of the LP2_1dx_file
    hdu = fits.open(LP2_1dx_file)
    lp2 = Table(hdu[1].data)

    ypeaks = []
    heights = []
    #
    apertures1 = []
    segments1 = []
    # loop apertures and segments
    for aperture in apertures:
        for segment in segments:
            irow = np.where((lp2['CENWAVE'] == cenwave) &
                            (lp2['APERTURE'] == aperture) &
                            (lp2['SEGMENT'] == segment) &
                            (lp2['OPT_ELEM'] == opt_elem))[0]
            if len(irow) != 1:
                print('irow error')
                pdb.set_trace()
            else:
                irow = irow[0]
            if verbose:
                if hd['DETECTOR'] == 'NUV':
                    print(lp2[irow]['B_SPEC'], lp2[irow]['SLOPE'])
                if hd['DETECTOR'] == 'FUV':
                    print(lp2[irow]['B_SPEC'])
            ypeaks.append(lp2[irow]['B_SPEC'])
            heights.append(lp2[irow]['HEIGHT'])
            apertures1.append(aperture)
            segments1.append(segment)

    # read and combine xfull, yfull for corrtag files
    xfull = []
    yfull = []
    for icorrtag in corrtag:
        hdu = fits.open(icorrtag)
        data1 = hdu[1].data
        xfull = xfull + list(data1['XFULL'])
        yfull = yfull + list(data1['YFULL'])
    xfull = np.asarray(xfull)
    yfull = np.asarray(yfull)



    for i in range(len(ypeaks)):
        # info
        print(apertures1[i], segments1[i])
        # data
        yfullcorr = yfull - slopes[i] * xfull
        k = np.where(abs(yfullcorr - ypeaks[i]) < iheight / 2)[0]
        yfull1 = yfullcorr[k]

        if fitgauss:
            yfull1 = yfull1 - ypeaks[i]

        # Generate histogram
        ymax = np.max(yfull1)
        ymin = np.min(yfull1)
        nh = (ymax - ymin) * nhres  # number of bins
        wbin = 1 / nhres # bin width
        ycen = np.linspace(ymin + wbin * 0.5, ymax - wbin * 0.5, nh)
        yedges = np.linspace(ymin, ymax, nh + 1)
        yhist, edges = np.histogram(yfull1, bins=yedges)

        # plot
        plt.figure(figsize=(14, 4))
        plt.plot(ycen, yhist, color='black')
        hmax = np.max(yhist)

        # ypeaks, heights
        if not fitgauss:
            # original ypeaks and heights
            plt.plot([ypeaks[i], ypeaks[i]], [0, hmax], color='red')
            plt.plot([ypeaks[i] - heights[i]/2, ypeaks[i] - heights[i]/2], [0, hmax], '--', color='red')
            plt.plot([ypeaks[i] + heights[i]/2, ypeaks[i] + heights[i]/2], [0, hmax], '--', color='red')

            # new ypeaks and heights
            col2 = 'blue'
            if newypeaks is not None:
                plt.plot([newypeaks[i], newypeaks[i]], [0, hmax], color=col2)
                if newheights is not None:
                    for newheight in newheights:
                        plt.plot([newypeaks[i] - newheight/2, newypeaks[i] - newheight/2], [0, hmax], '--',
                                 color=col2)
                        plt.plot([newypeaks[i] + newheight/2, newypeaks[i] + newheight/2], [0, hmax], '--',
                                 color=col2)

            # show percentiles
            col3 = 'limegreen'
            col3b = 'orange'
            if (percs is not None) | showorigperc:
                # data
                k = np.where(abs(yfull - newypeaks[0]) < iheight / 2)[0]
                yfull2 = yfull[k]
            if percs is not None:
                for iperc in percs:
                    yperc = np.percentile(yfull2,iperc,overwrite_input=False)
                    plt.plot([yperc,yperc],[0,hmax],'--',color=col3)
            if showorigperc:
                origperc = [lp2['LOWER_OUTER'][irow]*100, lp2['UPPER_OUTER'][irow]*100]
                for iperc in origperc:
                    yperc = np.percentile(yfull2, iperc, overwrite_input=False)
                    plt.plot([yperc, yperc], [0, hmax], '--', color=col3b)


        # fit a Gaussian
        if fitgauss:
            mean = 0
            sigma = 1
            cons = 1
            popt, pcov = curve_fit(gaus, ycen, yhist, p0=[1, mean, sigma, cons])
            print(mean, sigma, popt)
            plt.plot(ycen, gaus(ycen, *popt), color='red')

        plt.plot([ymin + wbin * 0.5, ymax - wbin * 0.5],[0,0],color='lightgray')
        plt.xlabel('YFULL_corr')
        plt.ylabel('N(YFULL_corr)')
        plt.show()



def gaus(x,a,x0,sigma,cons):
    """ For input array (or list) x, calculate a Gaussian, and with cons added.

    Parameters
    ----------
    x : array or list of floats
    a : float
    x0 : float
    sigma : float
    cons : float

    Returns
    -------
    Gaussian, multiplied with a constant cons, and with norm added

    """

    return (a*np.exp(-((x-x0)**2)/(2*sigma**2))) + cons


