#/u/ki/dapple/bin/python

##!/usr/bin/env python

import sys

import pyfits
import numpy as np

import cosmo
import pz_utils

#Note this function is very basic, and has not been updated
def redm_nlambda_calc(cat,area,outdir,zmin,zmax):
    '''
    Calculates n_clusters(lambda) for several preset bins in redshift,
    then outputs the results to  files
    Incorporates P(z), although the current redshift bins are fairly broad
    Output units are in [Mpc/h]^-3
    '''
    cat_zmax = np.max(cat['z_lambda'])
    
    nlbins = 91
    lmbins = 10+np.array(range(nlbins))

    zlist = np.where(zmin < cat_zmax)[0]
    if len(zlist) == 0:
        print >> sys.stderr, "ERROR:  Catalog zmax < 0.15 ?!"
        print >> sys.stderr, "        Unable to calculate n(lambda)"
        return

    zmin = zmin[zlist]
    zmax = zmax[zlist]
    
    n_of_lambda = np.zeros([len(zmin),len(lmbins)])

    for i in range(len(zmin)):
        p_in_bin = pz_utils.p_in_zbin(cat['pz'],cat['pzbins'],zmin[i],zmax[i])
        for j in range(nlbins):
            mylist = np.where( (cat['lambda_chisq'] >= lmbins[j]) & (cat['lambda_chisq'] < lmbins[j] + 1) )[0]
            n_of_lambda[i,j] = np.sum(p_in_bin[mylist])

        #Calculate volume for this redshift range, with default cosmology
        #assumes are given in square degrees
        vol = cosmo.comoving_volume(zmin[i],zmax[i])*area/41252.96

        n_of_lambda[i,:] = n_of_lambda[i,:]/vol

        #print output
        outfile = outdir+"nlambda_z_"+str(zmin[i])+"_"+str(zmax[i])+".dat"
            
        f = open(outfile,'w')
        for j in range(nlbins):
            print >> f, lmbins[j], n_of_lambda[i,j]
        f.close()

    return

def redm_nlambda(lambda_chisq,lambda_err,z,zmin,zmax,area,p_zbin,use_pz=False):
    '''
    Basic n(lambda) calculation; incorporates P(z) and P(lambda) if requested
    Note that this takes the p_in_zbin vector as an input to speed matters up
    '''
    cat_zmax = np.max(z)
    
    nlbins = 91
    lmbins = 10+np.array(range(nlbins))

    zlist = np.where(zmin < cat_zmax)[0]
    if len(zlist) == 0:
        print >> sys.stderr, "ERROR:  Catalog zmax < ",zmin[0]," ?!"
        print >> sys.stderr, "        Unable to calculate n(lambda)"
        return -1

    #print >> sys.stderr, "TEST2: ",len(zmin),len(zlist),cat_zmax
    zmin = zmin[zlist]
    zmax = zmax[zlist]
    
    n_of_lambda = np.zeros([len(zmin),len(lmbins)])

    for i in range(len(zmin)):
        if use_pz:
            for j in range(nlbins):
                p_in_lbin = pz_utils.p_in_lmbin(lambda_chisq,lambda_err,lmbins[j],lmbins[j]+1)
                n_of_lambda[i,j] = np.sum(p_zbin[i]*p_in_lbin)
        else:
            for j in range(nlbins):
                mylist = np.where( (lambda_chisq >= lmbins[j]) & (lambda_chisq < lmbins[j] + 1) & (z >= zmin[i]) & (z < zmax[i]) )[0]
                n_of_lambda[i,j] = len(mylist)

        #Calculate volume for this redshift range, with default cosmology
        #assumes are given in square degrees
        vol = cosmo.comoving_volume(zmin[i],zmax[i])*area((zmin[i]+zmax[i])/2.)/41252.96

        n_of_lambda[i,:] = n_of_lambda[i,:]/vol

    return n_of_lambda

#Updated to include widths in lambda
def redm_nlambda_err(lambda_chisq,lambda_err,z,pz,pzbins,
                     bootlist,outdir,zmin,zmax,area):
    '''
    Calculates n_clusters(lambda) for several preset bins in redshift,
    then outputs the results to  files
    Incorporates P(z), although the current redshift bins are fairly broad
    Output units are in [Mpc/h]^-3
    Includes error calculation
    '''

    nlbins = 91
    lmbins = 10+np.array(range(nlbins))
    nz = len(zmin)
    nclusters = len(lambda_chisq)
    
    #Set up the p_zbin array to save time
    p_zbin = np.zeros([nz,nclusters])
    for i in range(nz):
        p_zbin[i] = pz_utils.p_in_zbin(pz,pzbins,zmin[i],zmax[i])

    #Get the initial results first
    n_of_lambda = redm_nlambda(lambda_chisq,lambda_err,z,zmin,zmax,area,p_zbin,use_pz=True)

    nboot = len(bootlist)
    nl_boot = np.zeros([nboot,len(zmin),nlbins])
    for i in range(nboot):
        nl_boot[i,:,:] = redm_nlambda(lambda_chisq[bootlist[i]],lambda_err[bootlist[i]],z[bootlist[i]],zmin,zmax,area,p_zbin[:,bootlist[i]],use_pz=True)
        
    #Calculate the estimated error and covariance
    nl_err = np.zeros([len(zmin),nlbins])

    for i in range(len(zmin)):
        for j in range(nlbins):
            nl_err[i,j] = np.sum( (n_of_lambda[i,j] - nl_boot[:,i,j])**2 )/(nboot-1.) 
    nl_err = np.sqrt(nl_err)

    #Make the covariance matrix
    nl_covar = np.zeros([len(zmin),nlbins,nlbins])
    for i in range(len(zmin)):
        for j in range(nlbins):
            for k in range(nlbins):
                nl_covar[i,j,k] = np.sum( (n_of_lambda[i,j] - nl_boot[:,i,j])*
                                          (n_of_lambda[i,k] - nl_boot[:,i,k]) )/(nboot-1.) 

    #And print out the results
    for i in range(len(zmin)):
        outfile = outdir+"nlambda_z_"+str(zmin[i])+"_"+str(zmax[i])+".dat"
        f = open(outfile,'w')
        for j in range(nlbins):
            print >> f, lmbins[j],n_of_lambda[i,j],nl_err[i,j]
        f.close()

    #And print out the full covariance matrix
    for i in range(len(zmin)):
        outfile = outdir+"nlambda_covar_z_"+str(zmin[i])+"_"+str(zmax[i])+".dat"
        f = open(outfile,'w')
        for j in range(nlbins):
            for k in range(nlbins):
                f.write(str(nl_covar[i,j,k])+" ")
            f.write("\n")
        f.close()

    return
