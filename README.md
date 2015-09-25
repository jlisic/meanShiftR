meanShiftR
==========

R package for mean shift classification.

This R package requires that the flann library and development libraries are installed.  On OS X flann is available throught brew or ports, under Linux you should refer to your package management system repository.

If you are having issues with hdf5.h not being found on OS X, and are using brew.  You might want to just do 

  R CMD INSTALL --configure-vars='INCLUDE_DIR=/usr/local/include' meanShiftR 

This avoids using pkg-config, which at least on my laptop causes some odities.


