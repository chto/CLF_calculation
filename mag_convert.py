#!/u/ki/dapple/bin/python

#!/usr/bin/env python

import sys

import numpy as np

#Routines for converting back and forth between magnitudes and
#solar luminosities
#Default -- set for z=0.1 kcorrect zero points, i-band

def mag_to_Lsolar(mag,use_des=0,abs_solar=4.57595):
    '''
    mag -- input magnitudes
    use_des -- switch to DES i-band instead of SDSS (default) if ==1
    abs_solar -- Solar magnitude in input band.  Default to sdss i-band at 
                 z=0.1
    Note:  Other band corrections at z=0.1 are:
           SDSS: u 6.77927
                 g 5.43365
                 r 4.75539
                 i 4.57495
                 z 4.51551
           DES:  g 5.36500
                 r 4.72384
                 i 4.55972
                 z 4.51189
                 Y 4.50623
    All values listed drawn from k_solar_magnitudes from kcorrect
    '''
    
    if (use_des == 1) & (abs_solar==4.57595):
        abs_solar = 4.55972

    return 10.**(0.4*(abs_solar-mag))

def Lsolar_to_mag(lum,use_des=0,abs_solar=4.57595):
    '''
    lum -- input solar luminosities
    use_des -- switch to DES i-band instead of SDSS (default) if ==1
    abs_solar -- Solar magnitude in input band.  Default to sdss i-band at 
                 z=0.1
    Note:  Other band corrections at z=0.1 are:
           SDSS: u 6.77927
                 g 5.43365
                 r 4.75539
                 i 4.57495
                 z 4.51551
           DES:  g 5.36500
                 r 4.72384
                 i 4.55972
                 z 4.51189
                 Y 4.50623
    All values listed drawn from k_solar_magnitudes from kcorrect
    '''
    
    if (use_des == 1) & (abs_solar==4.57595):
        abs_solar = 4.55972

    return abs_solar - 2.5*np.log10(lum)
