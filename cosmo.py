#!/usr/bin/env python

from math import *
import numpy as np
from chtolib import cosmoCon as cn

#H(z)
def Hubble(z, H0=cn.H0, omegaM=cn.omegaM, omegaL=cn.omegaL):
    return H0*sqrt(omegaM*(1.+z)**3+omegaL)

def rco(z, H0=cn.H0, omegaM=cn.omegaM, omegaL=cn.omegaL):
    '''
    #comoving distance as a function of z
    #Defaults set to Planck cosmology
    ##In MPC unit
    '''
    c = cn.c

    #perform numerical integration
    dz = 0.001
    if len(z) > 1:
        zmax = max(z)
    else:
        zmax = z[0]
    nbins = int(ceil(zmax/dz))
    if zmax == 0:
        return z
    dz = zmax/(nbins-1)
    if (nbins == 0):
        return np.zeros(len(z))

    myr = np.zeros(nbins)
    myz = np.zeros(nbins)
    for i in range(nbins):
        myz[i] = dz*i

    for i in range(nbins):
        if myz[i] == 0:
            continue
        myr[i] = myr[i-1] + c*(1./Hubble(dz*(i+1), H0=H0, omegaM=omegaM, omegaL=omegaL)+
                    1./Hubble(dz*i, H0=H0, omegaM=omegaM, omegaL=omegaL) )*dz/2.

    #now interpolate to get all desired values
    #print z, myz, myr, nbins
    r = np.interp(z,myz,myr)
    return r

#luminosity distance modulus given distance and (true, non-peculiar) redshift
#Assumes that r is in Mpc
def lm_distmod(r,z):
    dm = 5.*np.log10((1.+np.array(z))*np.array(r)*1e6/10.)
    return dm

#Gives volume in a spherical shell
def comoving_volume(zmin,zmax,H0=cn.H0, omegaM=cn.omegaM, omegaL=cn.omegaL):
    '''
      in MPC^3

    '''
    [rmin] = rco(np.array([zmin]),H0=H0,omegaM=omegaM,omegaL=omegaL)
    [rmax] = rco(np.array([zmax]),H0=H0,omegaM=omegaM,omegaL=omegaL)
    
    vol = (rmax**3-rmin**3)*4.*np.pi/3.

    return vol
 

#Calculate virial overdensity at given redshift -- uses Bryan + Norman 1998
def virial_dens(z, H0=cn.H0, omegaM=cn.omegaM, omegaL=cn.omegaL):
    Esq = omegaM*(1+z)**3 + omegaL
    x = omegaM*(1+z)**3/Esq-1

    delta = 18*np.pi**2 + 82*x - 39*x**2

    return delta

#Critical density in Msolar/Mpc^3
def rho_crit(z,H0=cn.H0, omegaM=cn.omegaM, omegaL=cn.omegaL):
    rho = 3*Hubble(z, H0=H0, omegaM=omegaM, omegaL=omegaL)**2/8./np.pi/G
    return rho

#Calculation of sigma(R)^2 from input power spectrum
#Note that this requires input of the power spectrum, and includes conversion of the 
#power spectrum values to remove h's
def sigma_sq_matter_perturb(R, k, Pk, h=cn.h):
    #Make the Gaussian filter
    #fk2 = np.exp( -(k/h)*(k/h)*R*R )
    #Make top-hat filter
    y = k*R
    fk = 3/y**3 * (np.sin(y) - y*np.cos(y))
    fk2 = fk*fk

    #Get dlnk, assuming things have log-bins
    dlnk = np.log(k[1]/k[0])

    #And integrate
    sigma2 = np.sum(fk2*Pk*h**3*(k/h)*(k/h)*(k/h*dlnk))*4.*np.pi/(2.*np.pi)**3

    return sigma2
