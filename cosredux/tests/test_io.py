# Module to run tests on io
from __future__ import print_function, absolute_import, \
     division, unicode_literals
import os
import pytest
from astropy.table import Table
import numpy as np

from cosredux import io


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_modify_table():
    # FITS
    row_dict = {'SEGMENT': 'FUVA', 'OPT_ELEM': 'G140L', 'APERTURE': 'PSA', 'CENWAVE': 1280}
    filename = data_path('z2d19237l_1dx.fits')
    io.modify_table_value(filename, 'B_SPEC', row_dict, 1.234)

