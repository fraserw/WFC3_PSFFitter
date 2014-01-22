
The WFC3PSFFits package is intended to provide all the necessary tools to distort, 
plant, and fit TinyTim PSFs to WFC3 images (both ir and UV1/2 detectors). WFC3_phot
is intended to provide the necessary tools for aperture photometry of those same
data.

Note: The documentation is **woefully** incomplete. This project is fully functional
but not yet ready for novice consumption. Give it a month to get the documentation
in order.

For basic examples of how to use the code, take a look at:

aperPhot.py:
-in this example you will see how to call functions in WFC3_phot.py to get a good
psf centroid using a refinement grid. One should ignore all the data copying and
such, as this is local to my own particular data.

doMagFit.py:
-this code uses the centroids (provided by aperPhot.py) to perform both psf and
aperture photometry. The psf fitting is done with ROOT's minuit in the WFC3_PSFFits
module. The aperture photometry is done using pyraf in the WFC3_phot module.

For further questions, please contact the author, Wesley Fraser (westhefras<SPAMSUCKS>@<SPAMSUCKS>gmail.com.
  