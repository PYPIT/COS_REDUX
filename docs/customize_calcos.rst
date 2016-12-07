.. highlight:: rest

****************
Customize CALCOS
****************


THIS IS STILL A DRAFT VERSION!
==============================




This document will describe data reduction of COS spectra, by using CALCOS v2.21 and Python v3.5.2.

This is a brief recipe that explains how to reduce the data. More detailed explanations will be added later.

Primary reference to be cited: Worseck et al. 2016, ApJ, 825, 144.



1) Data needed
==============

- MAST raw and calibration data:
------------------------------

  Download from **https://archive.stsci.edu/** all data except calibrated.

  It might be needed to previously open an account at STScI.
  When available, the data could be retrieved by typing in the terminal:
    **ftp stdatu.stsci.edu**        (Then enter your username and password)
    cd /stage/.../.../          (Enter the location of your data)
    ls                          (List your files)
    prompt
    binary
    mget *
    quit

Create a separate directory for calibration files (the files downloaded).
Organize the raw data by object as you wish (see next step).


- Other uncalibrated data
-----------------------

  In order to (properly) run CALCOS, other needed files are:
  Uncalibrated Support Data (**asn, jit, jif, spt, trl**) and Uncalibrated Science Data (for time tag spectra: **rawtag** files).
  Names of needed files could be found in the COS Data Handbook (Table 2.1).
  No needed to copy jit and jif files, since they are in calibrated files. ?

  Copy these files in the directory with raw data (in a separate directory).
  Do not have Intermediate and Final data products in the same directory.
  There are also intermediate trl files.

- Dark frames
-----------

Download from MAST:

  a) At http://archive.stsci.edu/hst/search.php **search for darks** by entering:
     Target Name = DARK
     Resolver = Don't resolve
     Imagers = COS
     Observations = Calibration
     User-specified field 1 = Instrument Config = COS/FUV
     Sort By = Start Time

     This gives a list of all COS FUV darks sorted by observation date.

  b) From this list **download all exposures taken within ~3 months** around the
      date of your science observations. Regular dark monitoring exposures are
      1330s long.

  (tbd: obs --> darks; later)


2a) Get CALCOS Running (v2.21)
==============================

- **Download UREKA 1.4.1.2 from http://ssb.stsci.edu/ureka/1.4.1.2**

  Backup your login scripts: e.g. .bashrc and IRAF's login.cl (you may see the name of this file with ls -a).
  Install UREKA as explained in the webpage above.
  In the folder where UREKA is downloaded you can copy IRAF's login.cl file. (if it is not working)

- **Start UREKA:**

  In the terminal from the folder where UREKA is downloaded type:
  xgterm (or xterm)
  ur_setup (in the xgterm window)
  mkiraf (if you copied login.cl file)
  export lref=/data/cos/calib/ [this sets the directory with the calibration files]
                               [or do ...]
  pyraf

  Within pyraf type:
     stsdas
     hst_calib
     hstcos

  After ur_setup, you can also test UREKA installation: ur_test


2b) Python codes
================

Github...



Follow the steps below for each segment (a, b) separately.
==========================================================

3) Customize and run CALCOS (1st Pass)
======================================

(a) Original Notes
------------------

- **cd folder_with_raw_data**

- **For all FUV data set the following CALCOS calibration switches within pyraf:**

  thedit *rawtag*.fits[0] FLATCORR 'OMIT' [flatfielding of Poisson data creates fractional counts, which are hard to interpret, COS flatfielding is approximate anyhow]

  thedit *rawtag*.fits[0] PHACORR 'OMIT'  [initial setting, to be tuned later]

  thedit *rawtag*.fits[0] RANDSEED '0'    [choose any non-negative value for reproducable results]

  thedit *rawtag*.fits[0] BACKCORR 'OMIT' [we will do this in post-processing]

- **Set your calibration files wisely:**

  thedit *rawtag*.fits[0] SPOTTAB 'N/A'   [This calibration step only works on select datasets where hotspots have been identified. How and by whom? Better get rid of this nonsense.]

  thedit *rawtag*.fits[0] GSAGTAB 'lref$x6l1439el_gsag.fits'  [This will flag all gain sag holes as of 2013Jun20, so mainly geocoronal gain sag holes at the overused LP1.
                                                               Keeps the LP1 and LP2 traces usable that partly overlap with the dark calibration regions. All this is safe and tested!
                                                               Download this file from https://hst-crds.stsci.edu.]

  thedit *rawtag*.fits[0] XTRACTAB 'lref$x6q17586l_1dx.fits'  [This will work for LP2 data with CALCOS v2.21. In newer versions STScI renamed some columns, so they might not work with CALCOS v2.21.
                                                               All entries are to be modified later (see below).
                                                               Download this file from https://hst-crds.stsci.edu.]

  Downloading files from https://hst-crds.stsci.edu :
  Do not just use right click to save files.

(a) COS_REDUX Script
--------------------

The above thedit commands can be executed on a set of rawtag files
in a given folder with the modify_rawtag_for_calcos() method.  Here is an example::

    from cosredux import utils as crdxu
    crdxu.modify_rawtag_for_calcos(path_to_files)

This over-writes various header cards in the primary HDU.

(b) Original Notes
------------------

- **Define your WCA and PSA traces in CALCOS v2.21:** In the calibration directory type:

  tedit [your LP3 1dx calibration file]

          - Note: the WCA (the listed values -- Slope) and PSA (B_SPEC parameter) trace definitions of your setup.
          - (Note that the ycorr (B_SPEC) coordinates of the two detector segments differ by ~60 pixels.)
          - How to use tedit:
                      - exit: ctrl+D, type quit, write changes? yes/no
                      - (help: ctrl+D, type help)
                      - change a single value: just type, press space, and go to another cell (this could be also done in Python)

  tedit lref$x6q17586l_1dx.fits

          - Change the WCA and PSA trace definitions for your setup to the LP3 values from above.
          - For the PSA set SLOPE=0 (to avoid fractional counts)
                        and set B_SPEC to the next to the next .5 of the pixel given in the LP3 trace table.
          - For G140L data set the PSA HEIGHT=25.
          - (Change values only for G140L ?)

(b) COS_REDUX Script
--------------------

The LP3 trace solution may be applied to the LP2 with the modify_LP2_1dx_calib() script.
Here is an example::

    crdxu.modify_LP2_1dx_calib(cpath)

where cpath it the path to the calibration directory.

(c) Both
--------

- **Run CALCOS (1st Pass):**

  cd science_directory and type in pyraf:
  calcos [dataset prefix]_asn.fits

  (to reduce the data obtained in a particular visit at a particular setting.)

          - Look at the terminal output to see whether stim pulses and WCA traces are found.
            (Running CALCOS on _asn files issues a warning that the wavelength calibration
             lamp is not on. This will be automatically reset to Lamp 1, which is used for
             science, so this warning can be ignored. Also, warnings about stim pulses not
             being found for short time intervals (few seconds) should not cause trouble.)






4) Coadd and Inspect CORRTAG Files, and rerun CALCOS
====================================================

(a) Coadd files: in pyraf or in COS_Redux
-----------------------------------------

Pyraf
+++++

**a) Pyraf:** http://stsdas.stsci.edu/cgi-bin/gethelp.cgi?tmerge.hlp

tmerge lcya01fyq_corrtag_a.fits,lcya01ggq_corrtag_a.fits,lcya01g0q_corrtag_a.fits,lcya01gjq_corrtag_a.fits l01corrtagsapp_a.fits append

tmerge lcya01fyq_corrtag_b.fits,lcya01ggq_corrtag_b.fits,lcya01g0q_corrtag_b.fits,lcya01gjq_corrtag_b.fits l01corrtagsapp_b.fits append


COS_Redux
+++++++++

Input a list of corrtag files to combine with astropy.table.vstack.
Here is an example::

    corrtag_files = glob.glob(os.getenv('DROPBOX_DIR')+'/COS-LRG/LCYA01010/*corrtag_a.fits')
    coadd_bintables(corrtag_files, outfile='tst.fits')




(b) Find trace of object and arc
--------------------------------

hdu = fits.open(fa)
tbl = Table(hdu[1].data)
wave = tbl['WAVELENGTH']
yfull = tbl['YFULL']
hdu.close()

A: ymin=300, ymax=700, ytl=550
B: ymin, ymax , ytl

peaks=crude_histogram(yfull, ymin=300, ymax=700, ytl=550, pk_window=4.5, verbose=False)
obj=peaks[0]
arc=peaks[1]

obj_y=refine_peak(yfull, obj, pk_window=5., per_lim=(0.25,0.75))

arc_y=refine_peak(yfull, arc, pk_window=5., per_lim=(0.25,0.75))

(show_traces(wave, yfull, obj_y, arc_y))


**change value of trace:**

modify_table_value(filename, column, row_dict, value, outfil=None, clobber=False)



**- other:**

  slope

  pha


**- pyraf:**

  thedit *rawtag*.fits[0] PHACORR 'PERFORM'

  delete

  calcos [dataset prefix]_asn.fits

  (stim pulses)


**- files**

**- Change the DQ_WGT column in the _x1d tables according to DQ**





5) Darks
========






6)







This document will describe how to install PYPIT.

Installing Dependencies
=======================
Though we have tried to keep the number of dependencies low,
there are a few packages that need to be installed (various python packages,
GSL, and linetools).

In general, we recommend that you use Anaconda for the majority
of these installations.

Detailed installation instructions are presented below:

Python Dependencies
-------------------

PYPIT depends on the following list of Python packages. 

We recommend that you use `Anaconda <https://www.continuum.io/downloads/>`_ to install and/or update these packages.

* `python <http://www.python.org/>`_ versions 2.7, or 3.3 or later
* `numpy <http://www.numpy.org/>`_ version 1.10 or later
* `astropy <http://www.astropy.org/>`_ version 1.1 or later
* `scipy <http://www.scipy.org/>`_ version 0.17 or later
* `matplotlib <http://matplotlib.org/>`_  version 1.4 or later
* `PyQT4 <https://wiki.python.org/moin/PyQt/>`_ version 4 (needed for linetools)
* `Ginga <https://ginga.readthedocs.io/en/latest/>`_ latest version (highly recommended; essentially required)
* `h5py <https://www.h5py.org/>`_ version 2.6 (for data I/O)
*  yaml -- On Python 3 (at least), you may need to install pyyaml

If you are using Anaconda, you can check the presence of these packages with::

	conda list "^python$|numpy|astropy$|scipy$|matplotlib|PyQT|ginga|yaml|h5py"

If the packages have been installed, this command should print out all the packages and their version numbers.  

If any of these packages are missing you can install them with a command like::

	conda install PyQT

If any of the packages are out of date, they can be updated with a command like::

	conda update scipy

Installing Linetools
--------------------
The latest version of `linetools <https://github.com/linetools/linetools/>`_ is
also required for PYPIT.
Linetools is a package designed for the analysis of 1-D spectra.
The installation steps for linetools are provided
`here <http://linetools.readthedocs.io/en/latest/install.html/>`_.

According to the linetools documentation page, "If you wish to have
full functionality of the GUIs and are using MacOSX, then you probably
need to change your backend from macosx to TkAgg in the matplotlibrc file."


GSL
---

GSL installation
++++++++++++++++

The package complies Cython code that links to gsl routines.
These must be installed on your system prior to PYPIT installation.
We recommend that if you need to install GSL that you use Anaconda,
e.g.::

    conda install -c https://conda.anaconda.org/asmeurer gsl

You are also required to point the ENVIRONMENTAL variable
GSL_PATH to the path above the lib/ and include/ directories
You may determine this path with::

    gsl-config --prefix

It is possible you will also need to set the
LD_LIBRARY_PATH environmental variable to the gsl lib directory,
e.g.::

    export LD_LIBRARY_PATH=/u/xavier/anaconda/lib

.. _GSLELCAPITAN:

GSL on Mac OSX El Capitan
+++++++++++++++++++++++++
.. warning::

	**The above method for installing GSL with Anaconda will not work
	if you are using Mac OSX El Capitan!**

The Mac OSX El Capitan operating system introduced
"Sytem Integrity Protection" (SIP), which restricts root access to as well
as the creation of symlinks in SIP-protected folders (ex: /usr, /bin etc).
The /Users folder, where Anaconda generally installs packages,
is also SIP-protected. This means that the relative paths produced by
some of our Cython code are interfered with by SIP and will cause PYPIT to crash.

Here are some hacks to make the anaconda installation work as
well as some alternate installation instructions:

**1) Replace relative paths in compiled Cython files with full path** 
::

	 #in this example, GSL is installed in '/Users/USERNAME/anaconda/lib/'
	 cd PYPIT/pypit/
	 install_name_tool -change "@rpath/./libgsl.0.dylib" "/Users/USERNAME/anaconda/lib/libgsl.0.dylib" arcyextract.so
	 install_name_tool -change "@rpath/./libgslcblas.0.dylib" "/Users/USERNAME/anaconda/lib/libgslcblas.0.dylib" arcyextract.so
	 install_name_tool -change "@rpath/./libgsl.0.dylib" "/Users/USERNAME/anaconda/lib/libgsl.0.dylib" arcytrace.so
	 install_name_tool -change "@rpath/./libgslcblas.0.dylib" "/Users/USERNAME/anaconda/lib/libgslcblas.0.dylib" arcytrace.so
	 install_name_tool -change "@rpath/./libgsl.0.dylib" "/Users/USERNAME/anaconda/lib/libgsl.0.dylib" arcycomb.so
	 install_name_tool -change "@rpath/./libgslcblas.0.dylib" "/Users/USERNAME/anaconda/lib/libgslcblas.0.dylib" arcycomb.so
	 install_name_tool -change "@rpath/./libgsl.0.dylib" "/Users/USERNAME/anaconda/lib/libgsl.0.dylib" arcyproc.so
	 install_name_tool -change "@rpath/./libgslcblas.0.dylib" "/Users/USERNAME/anaconda/lib/libgslcblas.0.dylib" arcyproc.so
	 install_name_tool -change "@rpath/./libgsl.0.dylib" "/Users/USERNAME/anaconda/lib/libgsl.0.dylib" arcyutils.so
	 install_name_tool -change "@rpath/./libgslcblas.0.dylib" "/Users/USERNAME/anaconda/lib/libgslcblas.0.dylib" arcyutils.so
	 install_name_tool -change "@rpath/./libgsl.0.dylib" "/Users/USERNAME/anaconda/lib/libgsl.0.dylib" arcyarc.so
	 install_name_tool -change "@rpath/./libgslcblas.0.dylib" "/Users/USERNAME/anaconda/lib/libgslcblas.0.dylib" arcyarc.so
	 

**2) Disable System Integrity Protection**

This is a last resort solution and we do not
recommend it due to security concerns. Instructions for how
to do this can be
found `here <https://www.quora.com/How-do-I-turn-off-the-rootless-in-OS-X-El-Capitan-10-11/>`_.


**3) Install GSL with Homebrew instead of Anaconda**

Since `Homebrew <http://brew.sh/>`_ installs programs in /usr/local , which is not SIP protected, this should work without additional hacks.::

  brew install gsl

in which case the ``GSL_PATH`` variable should be set to ``/usr/local/Cellar/gsl/1.16/``, where ``1.16`` might have to
be replaced with whatever version number you have installed.

Since Homebrew installs programs in /usr/local , which is not
SIP protected, this should work without additional hacks.


Installing PYPIT
================

We recommend that you grab the code from github::

	#go to the directory where you would like to install PYPIT.
	git clone https://github.com/PYPIT/PYPIT.git

From there, you can build and install either with install or develop, e.g.::

	cd PYPIT
	python setup.py develop

or::

	cd PYPIT
	python setup.py install

This should compile all the necessary Cython files, etc.

Tests
=====
In order to assess whether PYPIT has been properly installed,
we suggest you run the following tests:

1. Ensure run_pypit works
-------------------------
Go to a directory outside of the PYPIT directory (e.g. your home directory),
then type run_pypit.::

	cd
	run_pypit


2. Run the PYPIT unit tests
---------------------------

Enter the PYPIT directory and do::

    python setup.py test


3. Try the test suite
---------------------
We have provided a suite of tests that you can download and run via this Repo:
`TestSuite <https://github.com/PYPIT/PYPIT-development-suite>`_

It can be installed as follows::

	# we suggest installing this in the directory above PYPIT
	git clone https://github.com/PYPIT/PYPIT-development-suite.git

To run the test::

	cd PYPIT-development-suite
	./pypit_test all

.. note::

	pypit_test can also take the argument kast instead of all. 


The test takes a while to run but should run without issue if all the packages have been properly installed. 


**If you installed GSL with anaconda, a common error from running ./pypit_test all is:**

|[BUG]     :: There appears to be a bug on Line 7 of arproc.py with error:

| dlopen(/Users/USERNAME/software/PYPIT/pypit/arcyextract.so, 2): Library not loaded: @rpath/./libgsl.0.dylib

| Referenced from: /Users/USERNAME/software/PYPIT/pypit/arcyextract.so


**To fix this bug:**

a) Make sure GSL_PATH and LD_LIBRARY_PATH are defined in your .bashrc or .tcshrc file and that the appropriate rc file has been sourced

b) If that does not work, check out :ref:`GSLELCAPITAN`.
