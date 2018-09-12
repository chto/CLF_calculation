#!/usr/bin/python

import numpy as np

#Function for converting 10-sigma limiting magnitudes to 5-sigma
def get_limmag_5sigma(lim_limmag,lim_exptime):
    n2 = 25.
    zp = 22.5
    flim10 = 10.**((lim_limmag-zp)/(-2.5))
    fn = ((flim10**2.*lim_exptime)/(10.0**2.) - flim10)
    flim5 = (n2 + np.sqrt(n2**2. + 4*lim_exptime*n2*fn))/(2.*lim_exptime)

    limmag_new = zp - 2.5*np.log10(flim5)

    return limmag_new
