""" Trace the source in Data file
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
import pdb
from astropy.table import Table
from astropy.io import fits

from xastropy.xutils import xdebug as xdb

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


