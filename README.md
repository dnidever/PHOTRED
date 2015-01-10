# PHOTRED
Automated DAOPHOT-based PSF photometry pipeline

PHOTRED is an automated DAOPHOT-based PSF photometry pipeline.  It's very generic but does required reduced flat images.
PHOTRED does WCS fitting, aperture photometry, single-image PSF photometry (ALLSTAR), source matching across multiple images,
forced PSF photometry across multiple exposures using a master source list created from a deep stack of all exposures
(ALLFRAME), aperture correction, calibration using photometric transformation equations, and dereddening.
STDRED is also included and can be used to derive the transformation equations from exposures of standard star fields.

PHOTRED does require the following that need to be installed separately:
-IRAF
-IDL
-The IDL Astronomy Uers's Library
-The stand-alone fortran version of DAOPHOT/ALLFRAME
-SExtractor

PDF manuals for PHOTRED and STDRED are included in the doc/ directory.

IDL "sav" files of all the main PHOTRED programs are in the sav/ directory.  These can be used with the IDL virtual
machine for people who don't have an IDL license.
