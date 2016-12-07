# Module to run tests on spectra.io
from __future__ import print_function, absolute_import, \
     division, unicode_literals
import os
import pytest
from astropy.table import Table
import numpy as np

from linetools.spectra import io
from linetools.spectra.xspectrum1d import XSpectrum1D
from linetools.spectra import utils as ltsu

from cosredux import utils

from cosredux import trace


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_crude():
    # FITS
    filename = data_path('l01corrtagsapp_a.fits')
    data = Table.read(filename)
    wave = data['WAVELENGTH']
    yfull = data['YFULL']
    # Crude
    obj_y, arc_y = trace.crude_histogram(yfull)
    # Refine object
    trace.refine_peak(yfull, obj_y)
    # Plot
    trace.show_traces(wave, yfull, obj_y, arc_y)


def test_find_dark():
    dstr = 'lcya'
    datastr = '01'
    datafld0 = 'new_redux_new/'
    cosfile='/home/marijana/Marijana/COS/LCYA01010/'
    darksfld=cosfile+'darksall2/'
    fa1r=cosfile+datafld0+'lcya'+datastr+'fyq_rawtag_a.fits'
    darks_a = utils.find_darks(darksfld, fa1r, 'a')


