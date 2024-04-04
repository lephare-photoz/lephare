List of keywords available for configuration of the LePhare executables 
=======================================================================

All the input ascii files can have comment lines starting with #. A single configuration file can be used for all the executables.

filter
------

.. list-table:: Keywords
   :widths: 25 25 50
   :header-rows: 1

   * - Keyword name
     - values
     - Comment
   * - FILTER_REP
     - 
     - Repository in which the filters are stored.
   * - FILTER_LIST
     - 
     - Comma-separated list of filter file paths relative to ``FILTER_REP``. These files are expected to have two columns: wavelengths and corresponding transmission.
   * - TRANS_TYPE
     - 0[def] or 1 
     - Define the transmission as being of type energy (0), or photon count (1). See ``flt::trans``
   * - FILTER_CALIB
     - 0[def],1,2,3,4,5 
     - comma-separated list of integer corresponding in order exactly to ``FILTER_LIST``, and defining which kind of calibration to apply to the filter transmission curves. 
   * - FILTER_FILE
     -  
     - output filename written in ``$LEPHAREWORK/filt``

sedtolib
--------

.. list-table:: Keywords
   :widths: 25 25 50
   :header-rows: 1

   * - Keyword name
     - values
     - Comment
   * - GAL_SED
     -  
     - name of file containing the list of Galaxy SED files to include
   * - GAL_FSCALE
     -  
     - Arbitrary Flux Scale
   * - GAL_LIB
     -  
     - basename of the Galaxy SED library output bin and doc files under ``$LEPHAREWORK/lib_bin``
   * - SEL_AGE
     -  
     - filename of the file containing the ages to consider|
   * - AGE_RANGE
     -  min,max 
     - age range to be considered
   * - QSO_SED
     -  
     - name of file containing the list of AGN/QSO SED files to include|
   * - QSO_FSCALE
     -  
     - Arbitrary Flux Scale
   * - QSO_LIB
     -  
     - basename of the AGN/QSO SED library output bin and doc files under ``$LEPHAREWORK/lib_bin``
   * - STAR_SED
     -  
     - name of file containing the list of Star SED files to include|
   * - STAR_LIB
     -  
     - basename of the Star SED library output bin and doc files under ``$LEPHAREWORK/lib_bin``
   * - STAR_FSCALE
     -  
     - Arbitrary Flux Scale

mag_gal
-------

.. list-table:: Keywords
   :widths: 25 25 50
   :header-rows: 1

   * - Keyword name
     - values
     - Comment
   * - COSMOLOGY
     -  H0,OmegaM,OmegaL
     -  LCDM relevant cosmological keywords
   * - FILTER_FILE
     -  
     -  same as above
   * - MAGTYPE
     -  AB[def] or VEGA 
     -  Magnitude type (AB or VEGA) of the library
   * - EXTINC_LAW
     -  calzetti.dat[def] 
     -   Extinction laws stored in ext/. Could include several files separated by comma.
   * - EB_V
     -  0[def] 
     -  Reddening color excess E(B-V) values to be applied. Values separated by comma.
   * - MOD_EXTINC
     -  0,0[def] 
     -  Range of models for which extinction will be applied. Two values per attenutation curves, separated by comma.
   * - ZGRID_TYPE
     -  0[def] or 1
     -  constant step in redshift if 0; evolving step in redshift as (1+z) if 1.
   * - Z_STEP
     -  0.04,0,6 [def] 
     -  dz,zmin,zmax: redshift step (dz), the minimum (zmin) and the maximum redshift (zmax)
   * - GAL_LIB_IN
     -  
     -  Input  galaxy library (in $ZPHOTWORK/lib_bin/)
   * - QSO_LIB_IN
     -  
     -   Input  AGN library (in $ZPHOTWORK/lib_bin/)
   * - STAR_LIB_IN
     -  
     -  Input  stellar library (in $ZPHOTWORK/lib_bin/)
   * - GAL_LIB_OUT
     -  
     -  Output galaxy library (in $ZPHOTWORK/lib_mag/)
   * - QSO_LIB_OUT
     -  
     -  Output AGN library (in $ZPHOTWORK/lib_mag/)
   * - STAR_LIB_OUT
     -  
     -  Output stellar library (in $ZPHOTWORK/lib_mag/)
   * - LIB_ASCII
     -  NO[def] or YES 
     -  Writing the magnitude library in ASCII file
   * - EM_LINES
     -  NO[def] or EMP_UV or EMP_SFR or PHYS 
     -  choice between emission line prescription
   * - EM_DISPERSION
     -  1[def] 
     -  Dispersion allowed in the emission line flux factor. Example: 0.5,0.75,1,1.5,2

zphota
------

.. list-table:: Keywords
   :widths: 25 25 50
   :header-rows: 1

   * - Keyword name
     - values
     - Comment
   * - CAT_IN
     -  
     -  photometric catalog input
   * - INP_TYPE
     -  F or M 
     -  Input  values:  Flux (F) or Magnitude (M)
   * - CAT_TYPE
     -  LONG[def] or SHORT 
     -   Input catalog format (long requires at minimum spec-z and context)
   * - CAT_MAG
     -  AB[def] or VEGA 
     -  Input magnitude type
   * - CAT_FMT
     -  MEME[def] or MMEE  
     -  Input format for photometry (MEME alternate mag-flux with error)
   * - CAT_LINES
     -  0,10000000000[def] 
     -   Min and max rows read in input catalog (to run only a subsample)
   * - PARA_OUT
     -  
     -   Name of the file with selected output parameters|
   * - CAT_OUT
     -  
     -  Name of the output file
   * - ZPHOTLIB
     -  
     -   Library names   (with no extension). Could have several separated by comma. Should be in LEPHAREWORK/lib_mag.
   * - ADD_EMLINES
     -  
     -   Range of galaxy models in which considering emission lines contribution.
   * - ERR_SCALE
     -  0[def] 
     -   Systematic errors (in mag) add in quadrature to the observations. One per filter, separated by comma.
   * - ERR_FACTOR
     -  1.0[def]  
     -  Scaling factor to the errors (in flux)
   * - BD_SCALE
     -  0[def] 
     -  Band used for scaling  the models to the observations (sum of $2^i$, as context). 0 means all.
   * - GLB_CONTEXT
     -  0[def] 
     -  Forces the context of all objects (sum of 2^i, as context). 0 means all.
   * - FORB_CONTEXT
     -  0[def] 
     -  context for removing some bands from the fit (sum of 2^i, as context). 0 means inactive.
   * - MASS_SCALE
     - 0, 0[def]
     -  Prior: allowed range in mass
   * - MAG_ABS
     -  0, 0[def]
     -  Prior: Absolute magnitude range allowed for the model)
   * - MAG_REF
     - 
     - Reference filter for the prior in abs. mag. (start at 1)
   * - NZ_PRIOR
     -  -1[def] 
     -   N(z) prior as function of i-band. Give the rank of i-band filter in input (start at 1). Negative value means no prior.
   * - ZFIX
     -  NO[def] or YES 
     -  Fixed redshift with the spec-z value|
   * - EXTERNALZ_FILE
     -  NONE[def] 
     -  Use the spec-z from an extrenal file (format Id,zs)
   * - Z_INTERP
     -  NO[def] 
     -   Parabolic interpolation between original step (dz)
   * - DZ_WIN
     -  0.25[def] 
     -  smoothing  window function for 2nd peak search in L(z)
   * - MIN_THRES
     -  0.1[def]
     -  threshold for the detection of 2nd peak in normalised L(z) (between 0 and 1)
   * - SPEC_OUT
     -  NO[def]  
     -  Output files with Gal/Star/QSO spectra (one file per object)
   * - CHI2_OUT
     -  NO[def]  
     -  Output files with the chi2 for the full library (one file per object)
   * - PDZ_OUT
     -  NONE[.pdz] 
     -  Output file name in which PDZ will be stored (full path). The code will add automatically the extension[.pdz]
   * - PDZ_TYPE
     -  BAY[def]  or MIN
     -   Define how the PDF is given in output of OUT\_PDZ.
   * - FIR_LIB
     -  NONE[def] 
     -   Far-IR libraries separated by comma
   * - FIR_LMIN
     -  7[def] 
     -   $\lambda$ min for FIR analysis (in $\mu m$)
   * - FIR_CONT
     -  -1[def] 
     -   Context for bands to be used in Far-IR
   * - FIR_SCALE
     -  -1[def]  
     -  Context for bands to be used for scaling
   * - FIR_FREESCALE
     -  NO[def] 
     -  Allows for free scaling
   * - FIR_SUBSTELLAR
     -  NO[def] 
     -  Removing stellar component from best optical fit
   * - MABS_METHOD
     -  0[def] or 1 or 2 or 3 or 4 
     -  Method used for absolute magnitudes in each filter
   * - Z_METHOD
     -  BEST[def] or ML 
     -   compute the absolute magnitude at a given redshift solution ML or BEST
   * - MABS_CONTEXT
     -  0[def]  
     -  Context for the bands used to derive  Mabs.
   * - MABS_REF
     -   1[def]  
     -  Filter used to derive the Mabs if method=2
   * - MABS_FILT
     -  1[def]   
     -  For method 4: list of  fixed filters chosen to derive Mabs in all bands according to the redshift bins
   * - MABS_ZBIN
     -  0,6[def] 
     -  For method 4: list of Redshift bins associated with  MABS_FILT. Even number of values.
   * - ADDITIONAL_MAG
     -  
     -   Name of file with filters, to derive Mabs in additional filters
   * - RF_COLORS
     -  -1,-1,-1,-1[def]  
     -  When computing uncertainties on abs. mag., do it in two colors (4 filters)
   * - M_REF
     -  
     -   Filter in which to compute the absolute magnitudes and associated errorbars
   * - APPLY_SYSSHIFT
     -  0[def] 
     -  Apply systematic shifts in each bands. Number of values must correspond to the number of filters.
   * - AUTO_ADAPT
     -  NO[def] 
     -   Optimize zero-points with spec-z 
   * - ADAPT_BAND
     -  1[def] 
     -  Reference band for the selection in magnitude
   * - ADAPT_LIM
     -  15,35[def] 
     -  Mag range for spectro in reference band
   * - ADAPT_ZBIN
     -  0.01,6[def] 
     -   Redshift's interval used for training
   * - LIMITS_ZBIN
     -  0.0,90.[def] 
     -  Redshift limits used to split in N bins, separated by a coma.
   * - LIMITS_MAPP_REF
     -  1[def] 
     -  Compute z-max. Band in which the absolute magnitude is computed.
   * - LIMITS_MAPP_SEL
     -  1[def] 
     -  Compute z-max. Give the selection band in each redshift bin.  Need 1 or N values.
   * - LIMITS_MAPP_CUT
     -  90[def] 
     -  Compute z-max. Magnitude cut used in each redshift bin. Need 1 or N values.
   * - RM_DISCREPENT_BD
     -  200,3[def] 
     -  define the threshold in chi2 to consider and the minimum improvement  in chi2 expected
   * - Z_RANGE
     -  
     -  Not implemented in c++
   * - EBV_RANGE
     -  
     -  Not implemented in c++
   * - ADAPT_CONTEXT
     -  
     -  Not implemented in c++
   * - ADAPT_MODBIN
     -  
     -    Not implemented in c++
   * - PROB_INTZ
     -  
     -  Not implemented in c++
   * - MAG_ABS_AGN
     -  
     -  Not implemented in c++
   * - CHI2_OUT
     -  
     -  Not implemented in c++
   * - PDZ_MABS_FILT
     -  
     -  Not implemented in c++
   * - LIR_PRIOR
     -  
     -  Not implemented in c++
   * - LIR_PRIOR_NEW
     -  
     -  Not implemented in c++
   * - ZFORM_MIN
     -  
     -  Not implemented in c++
   * - PHYS_LIB
     -  
     -  Not implemented in c++
   * - PHYS_CONT
     -  
     -  Not implemented in c++
   * - PHYS_SCALE
     -  
     -  Not implemented in c++
   * - PHYS_NMAX
     -  
     -  Not implemented in c++
   * - FAST_MODE
     -  
     -  Not implemented in c++
   * - COL_NUM
     -  
     -  Not implemented in c++
   * - COL_SIGMA
     -  
     -  Not implemented in c++
   * - COL_SEL
     -  
     -  Not implemented in c++
   * - ERROR_ADAPT
     -  
     -  Not implemented in c++ 