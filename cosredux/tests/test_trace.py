# Module to run tests on spectra.io
from __future__ import print_function, absolute_import, \
     division, unicode_literals
import os
import pytest

import numpy as np
from shutil import copyfile
import glob

from astropy.io import fits
from astropy.table import Table

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


def test_modify_raw():
    # Copy raw into files/raw
    try:
        os.mkdir(data_path('raw'))
    except OSError: # likely already exists
        pass
    datafld0=tst_path+'raw/'
    raw_files = glob.glob(datafld0+'*rawtag*')
    new_files = []
    for kk,raw_file in enumerate(raw_files[0:2]):
        root = raw_file[raw_file.rfind('/')+1:]
        new_file = data_path('raw/'+root)
        copyfile(raw_file,new_file)
        new_files.append(new_file)
        # Grab for testing
        if kk == 0:
            with fits.open(raw_file) as f:
                head = f[0].header
                fcorr = head['FLATCORR']
    #
    assert fcorr == 'PERFORM'
    utils.modify_rawtag_for_calcos(data_path('raw/'))
    # Test
    with fits.open(new_files[0]) as f:
        newhead = f[0].header
        newfcorr = newhead['FLATCORR']
    assert newfcorr == 'OMIT'


def test_modify_calibs():
    try:
        os.mkdir(data_path('calibs'))
    except OSError: # likely already exists
        pass
    calibfld=tst_path+'calibs/'
    calib_files = glob.glob(calibfld+'*1dx*')
    new_files = []
    for kk,calib_file in enumerate(calib_files[0:2]):
        root = calib_file[calib_file.rfind('/')+1:]
        new_file = data_path('calibs/'+root)
        copyfile(calib_file,new_file)
        new_files.append(new_file)
        # Grab for testing
        if kk == 0:
            with fits.open(calib_file) as f:
                head = f[0].header
    utils.modify_LP2_1dx_calib(data_path('calibs/'))
    # Test
    OPT_ELEM='G140L'
    CENWAVE=1280
    segment = 'FUVA'
    aperture = 'PSA'
    LP3_1dx_file = calibfld+'/z2d19237l_1dx.fits'
    LP2_1dx_file = calibfld+'/x6q17586l_1dx.fits'
    new_LP2_file = data_path('calibs/x6q17586l_1dx.fits')
    B_SPEC_vals = []
    for ifile in [LP3_1dx_file, LP2_1dx_file, new_LP2_file]:
        tbl = Table.read(ifile)
        row = (tbl['OPT_ELEM'] == OPT_ELEM) & (tbl['CENWAVE']==CENWAVE) & (
            tbl['SEGMENT'] == segment) & (tbl['APERTURE'] == aperture)
        B_SPEC_vals.append(tbl[row]['B_SPEC'].data[0])
    assert B_SPEC_vals[0] != B_SPEC_vals[1]
    np.testing.assert_allclose(B_SPEC_vals[2], 459.5)


def test_coadd():
    datafld=tst_path+'corrtag/'
    datastr='01'
    fa = data_path('l' + datastr + 'corrtagsapp_a.fits')
    fb = data_path('l' + datastr + 'corrtagsapp_b.fits')
    corrtag_files_a = glob.glob(datafld + '*_corrtag_a.fits')
    corrtag_files_b = glob.glob(datafld + '*_corrtag_b.fits')
    # Test on a
    nrows = []
    for afile in corrtag_files_a:
        with fits.open(afile) as f:
            nrows.append(len(f[1].data))
    utils.coadd_bintables(corrtag_files_a, outfile=fa)
    utils.coadd_bintables(corrtag_files_b, outfile=fb)
    # Test on A
    with fits.open(fa) as f:
        new_rows = len(f[1].data)
    assert np.sum(np.array(nrows)) == new_rows


def test_traces_2():
    datastr='01'
    fa = data_path('l' + datastr + 'corrtagsapp_a.fits')
    fb = data_path('l' + datastr + 'corrtagsapp_b.fits')
    fcd=data_path('calibs/x6q17586l_1dx.fits')   ###calibfld+'x6q17586l_1dx.fits'
    row_dict_a = {'OPT_ELEM': 'G140L', 'CENWAVE': '1280', 'SEGMENT': 'FUVA', 'APERTURE': 'PSA'}
    row_dict_b = {'OPT_ELEM': 'G140L', 'CENWAVE': '1280', 'SEGMENT': 'FUVB', 'APERTURE': 'PSA'}
    traces_a = trace.traces(fa, fcd, row_dict_a, outfil=fcd, clobber=True)
    traces_b = trace.traces(fb, fcd, row_dict_b, ymin=360, ymax=760, ytl=610, outfil=fcd, clobber=True)
    np.testing.assert_allclose(traces_a[0], 466.14, rtol=1e-4)
    np.testing.assert_allclose(traces_b[0], 525.9885, rtol=1e-4) #525.9885
   ## pytest.set_trace()

