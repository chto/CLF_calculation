#!/usr/bin/python

import numpy as np
import sys

def get_p_bcg_not_cen(cat,cenmag,cindex,pcen,g_mem_id,mag,p,
                      zmin,zmax,outdir,use_lum=1,weight_cen=0):
    '''
    Note that this produces P(BCG!=central) as a function of lambda
    in each redshift bin

    Note that this also assumes that the clusters and galaxies are ordered by mem_match_id
    '''
    ncl = len(cat)
    nz = len(zmin)
    #set up the arrays for holding our data
    lm_bins = np.array(range(100))+5.
    count_clusters = np.zeros([nz,len(lm_bins)])
    count_brighter = np.copy(count_clusters)

    zlist = np.where( (cat['z_lambda'] >= np.min(zmin)) & (cat['z_lambda'] < np.max(zmax)))[0]

    #Loop over all clusters in the given redshift range
    count = 0L
    count_last = 0L
    ncen = 1
    #print len(zlist),np.max(zlist),len(cat)
    for i in range(len(zlist)):
        zbin = np.where( (cat['z_lambda'][zlist[i]] >= zmin) & (cat['z_lambda'][zlist[i]] < zmax) )[0]
        if len(zbin)==0:
            print >> sys.stderr, "SKIP: ",i,cat['z_lambda'][zlist[i]]
            continue
        zbin = zbin[0]
        lmbin = np.where( (cat['lambda_chisq'][zlist[i]] >= lm_bins) & (cat['lambda_chisq'][zlist[i]] < lm_bins+1.) )[0]
        if len(lmbin)==0:
            continue
        lmbin = lmbin[0]
    
        #Add the cluster
        count_clusters[zbin,lmbin] = count_clusters[zbin,lmbin] + 1.

        if weight_cen==1:
            ncen = cat['ncent_good'][zlist[i]]
        if ncen == 0:
            #Oops, no centrals!  Set BCG>central to 1 and move on
            count_brighter[zbin,lmbin] = count_brighter[zbin,lmbin] + 1
            continue
        brighter = 0.
        
        #Now, find the first galaxy that matches our list
        while g_mem_id[count] < cat['mem_match_id'][zlist[i]]:
            count = count+1
        #print >> sys.stderr, "TEST: ",i, count, count_last, len(zlist),zbin, lmbin, g_mem_id[0], cat['mem_match_id'][zlist[i]]

        #Find the last galaxy listed in the cluster
        count_last = count
        while g_mem_id[count_last] == cat['mem_match_id'][zlist[i]]:
            count_last = count_last+1
            if count_last==len(g_mem_id):
                break
        #Get the list of all galaxies, sorted on brightness
        glist = np.argsort(mag[count:count_last])+count
        #print >> sys.stderr, "TEST2: ",glist,count,count_last
        if use_lum:
            glist = glist[::-1]

        #Run through the four brightest galaxies, and add up those that are
        #Brighter than the central, not the central, with appropriate probabilities
        ngals_check = np.min([5,len(glist)])
        brighter = 0.

        p_cen = cat['p_cen']
        if len(p_cen.shape)==1:
            p_cen = p_cen.reshape(-1,1)
        for j in range(ncen):
            if use_lum == 1:
                clist = np.where( (glist!=cindex[zlist[i]][j]) & (mag[glist] > cenmag[zlist[i]][j]) )[0]
            else:
                clist = np.where( (glist!=cindex[zlist[i]][j]) & (mag[glist] < cenmag[zlist[i]][j]) )[0]
            if len(clist) == 0:
                brighter = brighter + p_cen[zlist[i]][j]
                continue
            brighter = brighter + p_cen[zlist[i]][j]*np.prod(1-p[glist[clist]])

        #Add it all up
        count_brighter[zbin,lmbin] = count_brighter[zbin,lmbin] + (1-brighter)

    #Print out the results
    print >> sys.stderr, "    Printing P(BCG!=central)..."
    for i in range(nz):
        outfile = outdir + "pbcg_cen_z_"+str(zmin[i])+"_"+str(zmax[i])+".dat"
        #print >> sys.stderr, outfile
        f = open(outfile,'w')
        for j in range(len(lm_bins)):
            print >> f, lm_bins[j], count_clusters[i,j], count_brighter[i,j]#, count_brighter[i,j]/count_clusters[i,j]
        f.close()
        
    print >> sys.stderr, "    Done printing"
        
    return


#A variant of the above function, which, rather than tabulating, just
#returns the luminosities of the top four satellites in each cluster, their
#probabilities, and which central they match up with (if any)
#Used for initial check of lambda-Lcen-Lsat(brightest) relation
def get_brightest_sat(cat,cindex,g_mem_id,mag,p,use_lum=1):
    '''
    Testing function for lambda-Lcen-Lsat(bright) relation
    Returns luminosities, probabilities, and central id (if any) for 
    each cluster
    '''
    ncl = len(cat)
    nsat = 4

    sat_mag = np.zeros([ncl,nsat])
    sat_p = np.zeros([ncl,nsat])
    cen_id = np.zeros([ncl,nsat]).astype(int)-1

    #Loop over all clusters
    count = 0L
    count_last = 0L
    ncen = 1
    for i in range(ncl):
        ncen = cat['ncent_good'][i]
        
        #Now, find the first galaxy that matches our list
        while g_mem_id[count] < cat['mem_match_id'][i]:
            count = count+1
        count_last = count
        #And now the last galaxies that matches
        while g_mem_id[count_last] == cat['mem_match_id'][i]:
            count_last = count_last+1
            if count_last==len(g_mem_id):
                break
        glist = np.argsort(mag[count:count_last])+count
        if use_lum:
            glist = glist[::-1]

        for k in range(np.min([len(glist),nsat])):
            cmatch = np.where(glist[k] == cindex[i])[0]
            if len(cmatch) > 0:
                cen_id[i][k] = cmatch[0]
            sat_mag[i,k] = mag[glist[k]]
            sat_p[i,k] = p[glist[k]]
        

    return sat_mag, sat_p, cen_id
