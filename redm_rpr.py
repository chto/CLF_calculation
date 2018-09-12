#!/u/ki/dapple/bin/python

##!/usr/bin/env python

import sys

import pyfits
import numpy as np

from collections import Counter
import time

#Routines for calculating the radial profiles

#Weights galaxies by p
#-- histogram with essentially CLFs of individual clusters
def count_galaxies_p_r(c_mem_id,scaleval,g_mem_id,r,p,mag,rbins,minlum,maxlum):
    nclusters = len(cat)
    nlum = len(minlum)

    nr = len(rbins)
    dr = np.log10(rbins[1]/rbins[0])
    minr = rbins[0]/10**(dr/2.)
    maxr = rbins[-1]*10**(dr/2.)

    #Make the array that will contain the full count of galaxies
    #in each cluster
    count_arr = np.zeros([nclusters,nlum,nr])

    #Make the index that takes mem_match_id to cluster
    index = np.zeros(np.max(g_mem_id)+1) - 100
    index[c_mem_id] = range(len(c_mem_id))

    #Select only galaxies that are in the given r, luminosity range
    min_minlum = min(minlum)
    max_maxlum = max(maxlum)
    mylist = np.where( (r <= maxr) & (r >= minr) & (mag >= min_minlum) & (mag < max_maxlum) )[0]
    if len(mylist)==0:
        print >> sys.stderr, "WARNING:  No galaxies found in range ",minr,maxr,min_minlum,max_maxlum
        return count_arr

    #Loop over galaxies, counting each in turn
    for i in range(len(mylist)):
        mycluster = index[g_mem_id[mylist[i]]]
        mybin = np.floor(np.log10(r[mylist[i]]/minr)/dr)
        mybin = mybin.astype(int)
        if (mybin < 0) | (mybin >= nr):
            continue
        mylbins = np.where( (mag[mylist[i]] > minlum) & (mag[mylist[i]] < maxlum) )[0]
        if len(mylbins) == 0:
            continue
        count_arr[mycluster,mylbins,mybin] += p[mylist[i]]*scaleval[mycluster]
    
    return count_arr

#Alternative random galaxy counting setup, which runs on all bootstrap samples at once
def count_galaxies_rand_all_r(bootlist, gboot, c_mem_id, scaleval, g_mem_id, minlum, maxlum, rbins, mag, p, r):
    nclusters = len(c_mem_id)
    ngals = len(g_mem_id)

    nlum = len(minlum)

    nr = len(rbins)
    dr = np.log10(rbins[1]/rbins[0])
    minr = rbins[0]/10**(dr/2.)
    maxr = rbins[-1]*10**(dr/2.)

    nboot = len(bootlist)

    count_arr = np.zeros([nclusters, nlum, nr])
    
    #Pick out which bin is needed for each galaxy
    mylist = np.where( (mag < max(maxlum) ) & (mag >= min(minlum)) & (r < maxr) & (r >= minr) )[0]
    if len(mylist)==0:
        print >> sys.stderr, "WARNING:  No galaxies found in range ",minr,maxr,minlum,maxlum
        return count_arr
    mybin = np.floor(np.log10(r[mylist]/minr)/dr)
    mybin = mybin.astype(int)
    start = time.clock()

    #All galaxiest listed for each cluster in gboot are included with weight 1
    #So, loop over clusters first
    for i in range(nclusters):
        glist = gboot[i]
        if len(glist)==0:
            continue
        binlist = np.where( (mybin[glist] > 0) & (mybin[glist] < nr) )[0]
        if len(binlist)==0:
            continue
        glist = glist[binlist]
        for j in range(len(glist)):
            lbins = np.where( (mag[glist[j]] >= minlum) & (mag[glist[j]] < maxlum) )[0]
            if len(lbins)==0:
                continue
            count_arr[i,lbins,j] = count_arr[i,lbins,j] + scaleval[bootlist[i]]
        
    print >> sys.stderr, "Time to complete loop: ",time.clock()-start,(time.clock()-start)/len(mylist)

    return count_arr


#Function for summing up a count_arr to make a set of radial profiles
#For a set of limits on luminosity, and fixed limits on redshift, richness
def make_single_rpr(lm,z,minlum,maxlum,count_arr,lm_min,lm_max,zmin,zmax,rbins):
    dr = np.log10(rbins[1]/rbins[0])
    rpr = np.zeros([len(minlum),len(rbins)])

    clist = np.where( (z >= zmin) & (z < zmax) & (lm >= lm_min) & (lm < lm_max) )[0]
    
    if len(clist) == 0:
        print >> sys.stderr, "WARNING: no clusters found for limits of: "
        print >> sys.stderr, lm_min,lm_max,zmin,zmax
        return rpr
    
    #Add up radial profiles for all relevant clusters, for each luminosity range
    for i in range(len(minlum)):
        if len(clist) == 0:
            continue
        rpr[i,:] = np.sum(count_arr[clist,i,:],0)/len(clist)/( (rbins*10.**(dr/2.))**2 - (rbins/10.**(dr/2.))**2 )/(4*np.pi)
    
    return rpr

#Prints CLFs and Covar to file for given redshift and lambda bin
def print_rpr_covar(rbins,rpr,covar,outfile,covarfile):
    nlum = len(rbins)
    f = open(outfile,'w')
    for i in range(nlum):
        f.write(str(rbins[i])+" "+str(rpr[i])+" "+str(np.sqrt(covar[i,i]))+"\n")
    f.close()

    f = open(covarfile,'w')
    for i in range(nlum):
        for j in range(nlum):
            f.write(str(covar[i,j])+" ")
        f.write("\n")
    f.close()
    return

#Main radial profiles function that does all requested operations and also
#outputs results to files
def redm_rpr(cat,mem,mag,lm_min,lm_max,zmin,zmax,minlum,maxlum,
             bootlist,gboot,outdir):

    #First, define the limits in radius
    #Note that this uses logarithmic binning
    dr = 0.08
    nr = 30
    rbins = 10.**(np.array(range(nr))*dr-2+dr/2.)

    nlambda = len(lm_min[0])
    nz = len(zmin)
    nlum = len(minlum)
    nboot = len(bootlist)

    #Make the empty radial profile array
    rpr = np.zeros([nz,nlambda,nlum,nr])

    #Set up all necessary counting arrays for the main radial profiles
    #Counting all galaxies
    count_arr = count_galaxies_p_r(cat['mem_match_id'],cat['scaleval'],mem['mem_match_id'],mem['r'],mem['p'],mag,rbins,minlum,maxlum)
    
    #Total up the radial profile for each bin in radius and redshift
    for i in range(nz):
        for j in range(nlambda):
            #print np.max(cat['lambda_chisq']),np.min(cat['lambda_chisq'])
            #print >> sys.stderr, i, j
            rpr[i,j,:,:] = make_single_rpr(cat['lambda_chisq'],cat['z_lambda'],
                                           minlum,maxlum,count_arr,lm_min[i,j],lm_max[i,j],
                                           zmin[i],zmax[i],rbins)

    print >> sys.stderr, "Beginning radial profile covariance calculations..."
    #Now, run through all of the bootstrap samples to make errors and covariance matrices    
    rpr_boot = np.zeros([nboot,nz,nlambda,nlum,nr])
    
    for i in range(nboot):
        #Use the all-counting function to count galaxies
        count_arr_b = count_galaxies_rand_all_r(bootlist[i],gboot[i], cat['mem_match_id'], cat['scaleval'],mem['mem_match_id'], minlum, maxlum, rbins, mag, mem['p'], mem['r'])
        for j in range(nz):
            for k in range(nlambda):
                rpr_boot[i,j,k] = make_single_rpr(cat['lambda_chisq'][bootlist[i]], cat['z_lambda'][bootlist[i]],minlum,maxlum, count_arr_b,lm_min[j,k],lm_max[j,k], zmin[j],zmax[j],rbins) 
    
    #Covariance matrices -- currenly only within a single z, lambda, luminosity bin
    covar_rpr = np.zeros([nz,nlambda,nlum,nr,nr])
    for i in range(nz):
        for j in range(nlambda):
            for k in range(nlum):
                for l in range(nr):
                    for m in range(nr):
                        covar_rpr[i,j,k,l,m] = np.sum( ( rpr_boot[:,i,j,k,l] - rpr[i,j,k,l] )*
                                                     ( rpr_boot[:,i,j,k,m] - rpr[i,j,k,m] ) )/(nboot-1.)

                #And print the resulting covariance matrices, and profiles with errors
                outfile = outdir+"rpr_z_"+str(zmin[i])+"_"+str(zmax[i])+"_lm_"+str(lm_min[i,j])[0:5]+"_"+str(lm_max[i,j])[0:5]+"_m_"+str(minlum[k])+"_"+str(maxlum[k])+".dat"
                covarfile = outdir+"rpr_covar_z_"+str(zmin[i])+"_"+str(zmax[i])+"_lm_"+str(lm_min[i,j])[0:5]+"_"+str(lm_max[i,j])[0:5]+"_m_"+str(minlum[k])+"_"+str(maxlum[k])+".dat"
                print_rpr_covar(rbins,rpr[i,j,k],covar_rpr[i,j,k],outfile,covarfile)

    #Test printing of bootstrap measurements
    for i in range(nboot):
        for j in range(nlum):
            outfile = outdir+"rpr_boot_"+str(i)+"_z_"+str(zmin[3])+"_"+str(zmax[3])+"_lm_"+str(lm_min[3,1])[0:5]+"_"+str(lm_max[3,1])[0:5]+"_m_"+str(minlum[j])+"_"+str(maxlum[j])+".dat"
            print_rpr_covar(rbins,rpr_boot[i,3,1,j],covar_rpr[3,1,j],outfile,"test.dat")
    
    return
