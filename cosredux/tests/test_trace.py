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

tst_path = os.getenv('DROPBOX_DIR')+'/COS-LRG/tests/'

def test_crude():
    # FITS
    filename = tst_path+'l01corrtagsapp_a.fits'
    data = Table.read(filename)
    wave = data['WAVELENGTH']
    yfull = data['YFULL']
    # Crude
    obj_y, arc_y = trace.crude_histogram(yfull)
    # Refine object
    pk = trace.refine_peak(yfull, obj_y)
    np.testing.assert_allclose(pk, 466.02716064453125)
    # Plot
    show = False
    if show:
        trace.show_traces(wave, yfull, obj_y, arc_y)


def test_find_dark():
    # Setup
    darksfld = tst_path+'darks/'
    scifilea = tst_path+'raw/lcya01fyq_rawtag_a.fits'
    scifileb = tst_path+'raw/lcya01fyq_rawtag_b.fits'
    # Run
    darks_a = utils.find_darks(darksfld, scifilea, 'a')
    darks_b = utils.find_darks(darksfld, scifileb, 'b')
    # Test
    assert len(darks_a) == 1
    assert len(darks_b) == 1


