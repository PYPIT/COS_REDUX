""" Handles I/O for COS products
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
from astropy.table import Table
from astropy.io import fits

from linetools import utils as ltu

from xastropy.xutils import xdebug as xdb

def write_bg_regions(bg_region, outfile):
    """ Write background regions to a simple JSON file

    Parameters
    ----------
    bg_region : dict
    outfile : str

    Returns
    -------

    """
    jdict = ltu.jsonify(bg_region)
    # Write
    ltu.savejson(outfile, jdict, easy_to_read=True, overwrite=True)
    print("Wrote Background Regions to {:s}",outfile)


def read_traces(coadd_file):
    """ Read traces from a .json file with the same name as a .fits file (coadd_file)
    Parameters
    ----------
    coadd_file : str

    Returns
    -------

    """
    trc_file = coadd_file.replace('.fits', '_traces.json')
    tdict = ltu.loadjson(trc_file)
    # Return
    return tdict['obj'], tdict['arc']


def write_traces(obj, arc, outfile):
    """ Write a simple JSON file
    Parameters
    ----------
    obj : float
    arc : float
    outfile : str

    Returns
    -------

    """
    tdict = dict(obj=obj, arc=arc)
    jdict = ltu.jsonify(tdict)
    # Write
    ltu.savejson(outfile, jdict, easy_to_read=True, overwrite=True)
    print("Wrote Traces to {:s}",outfile)

