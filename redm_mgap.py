#/u/ki/dapple/bin/python

import os
import sys

import pyfits
import numpy as np

#Various routines necessary for calculating magnitude gaps

#Count weighted gap measure,including central weighting
def gaps_weight(c_mem_id,cenmag,cengalindex,pcen,g_mem_id,mag,p,use_lum=False):
    clist = np.argsort(c_mem_id)
    glist = np.argsort(g_mem_id)

    gaps_count = np.zeros([len(clist),20])
    gaps_val = np.zeros([len(clist),20])
    count_lo = 0L
    count_hi = 0L
    for i in range(len(clist)):
        count_lo = count_hi
        #Getting the range of matched satellites
        while c_mem_id[clist[i]] > g_mem_id[glist[count_lo]]:
            count_lo = count_lo+1
            count_hi = count_lo
        while c_mem_id[clist[i]] == g_mem_id[glist[count_hi]]:
            count_hi = count_hi+1
            if count_hi >= len(glist):
                break

        #Set up the temporary array to hold probabilities, etc., while we collect them
        ntestsat = 10
        p_temp = np.zeros([5,ntestsat])
        gap_temp = 0.*p_temp

        #Get the sorted list of luminosities for each cluster
        slist = np.argsort(mag[glist[count_lo:count_hi]])
        if use_lum:
            slist = slist[::-1]
        slist = glist[slist+count_lo]
        for j in range(len(pcen[clist[i]])):
            if pcen[clist[i]][j] == 0:
                continue
            mbin = 0
            mcount = 0
            while mbin < ntestsat:
                if mcount >= len(slist):
                    break
                #Skip the magnitude that corresponds to the central
                if cengalindex[clist[i],j] == slist[mcount]:
                    mcount = mcount+1
                    continue
                
                #Probabilistic counting
                if mbin == 0:
                    p_temp[j,mbin] = p[slist[mcount]]
                else:
                    fac = 1.
                    for k in range(mcount):
                        if cengalindex[clist[i],j] != slist[k]:
                            fac = fac*(1-p[slist[k]])
                    #fac = 1-fac
                    p_temp[j,mbin] = p[slist[mcount]]*fac
                p_temp[j,mbin] = p_temp[j,mbin]*pcen[clist[i],j]
                gap_temp[j,mbin] = cenmag[clist[i],j] - mag[slist[mcount]]

                mcount = mcount+1
                mbin = mbin+1
        #Now leaving the loop over all centrals for this cluster
        p_temp = p_temp.flatten()
        gap_temp = gap_temp.flatten()
        #Sort, high to low probability
        plist = np.argsort(p_temp)[::-1]
        #Save the most probable values
        for j in range(len(gaps_count[clist[i]])):
            gaps_count[clist[i],j] = p_temp[plist[j]]
            gaps_val[clist[i],j] = gap_temp[plist[j]]

    return [ gaps_count, gaps_val ]

#Count weighted gap measure, using most likely central only
def gaps_noweight(c_mem_id,cenmag,cengalindex,g_mem_id,mag,p,use_lum=False):
    clist = np.argsort(c_mem_id)
    glist = np.argsort(g_mem_id)
    
    count_lo = 0L
    count_hi = 0L
    gaps_count = np.zeros([len(clist),10])
    gaps_val = np.zeros([len(clist),10])
    for i in range(len(clist)):
        mbin = 0
        mcount = 0
        count_lo = count_hi
        while c_mem_id[clist[i]] > g_mem_id[glist[count_lo]]:
            count_lo = count_lo+1
            count_hi = count_lo
        while c_mem_id[clist[i]] == g_mem_id[glist[count_hi]]:
            count_hi = count_hi+1
            if count_hi >= len(glist):
                break
        
        slist = np.argsort(mag[glist[count_lo:count_hi]])
        if use_lum:
            slist = slist[::-1]
        slist = glist[slist+count_lo]
        while mbin < 10:
            if mcount >= len(slist):
                break
            #Skip the magnitude that corresponds to the central
            if cengalindex[clist[i]] == slist[mcount]:
                mcount = mcount+1
                continue
            
            #Correct probabilistic counting
            if mbin == 0:
                gaps_count[clist[i]][mbin] = p[slist[mcount]]
            else:
                fac = 0.
                for j in range(mbin):
                    fac = fac+gaps_count[clist[i]][j]
                fac = 1. - fac
                gaps_count[clist[i]][mbin] = p[slist[mcount]]*fac
            gaps_val[clist[i]][mbin] = cenmag[clist[i]] - mag[slist[mcount]]
            #if i < 3 and mbin>0:
            #    print i,mbin,p[slist[mcount]],gaps_count[clist[i]][mcount],np.sum(gaps_count[clist[i]][0:mbin]),fac
            mbin = mbin+1
            mcount=mcount+1

    return [gaps_count, gaps_val]

#Count weighted gap measure, using most likely central only
#Also returns an array of galaxy indices for alternative uses
def gaps_noweight_index(c_mem_id,cenmag,cengalindex,g_mem_id,mag,p,use_lum=False):
    clist = np.argsort(c_mem_id)
    glist = np.argsort(g_mem_id)
    
    count_lo = 0L
    count_hi = 0L
    gaps_count = np.zeros([len(clist),10])
    gaps_val = np.zeros([len(clist),10])
    g_index = np.zeros([len(clist),10]).astype(long)-1
    for i in range(len(clist)):
        mbin = 0
        mcount = 0
        count_lo = count_hi
        while c_mem_id[clist[i]] > g_mem_id[glist[count_lo]]:
            count_lo = count_lo+1
            count_hi = count_lo
        while c_mem_id[clist[i]] == g_mem_id[glist[count_hi]]:
            count_hi = count_hi+1
            if count_hi >= len(glist):
                break
        
        slist = np.argsort(mag[glist[count_lo:count_hi]])
        if use_lum:
            slist = slist[::-1]
        slist = glist[slist+count_lo]
        while mbin < 10:
            if mcount >= len(slist):
                break
            #Skip the magnitude that corresponds to the central
            if cengalindex[clist[i]] == slist[mcount]:
                mcount = mcount+1
                continue
            
            #Correct probabilistic counting
            if mbin == 0:
                gaps_count[clist[i]][mbin] = p[slist[mcount]]
            else:
                fac = 0.
                for j in range(mbin):
                    fac = fac+gaps_count[clist[i]][j]
                fac = 1. - fac
                gaps_count[clist[i]][mbin] = p[slist[mcount]]*fac
            gaps_val[clist[i]][mbin] = cenmag[clist[i]] - mag[slist[mcount]]
            g_index[clist[i]][mbin] = slist[mcount]

            mbin = mbin+1
            mcount=mcount+1

    return [gaps_count, gaps_val, g_index]

#Printing function (note lack of error bars here)
def mgap_print(lbins,gaps,outfile):
    f = open(outfile,'w')
    for i in range(len(lbins)):
        print >> f, lbins[i], gaps[i]
    f.close()
    return

#Main magnitude gap routine
#Calculates weighted magnitude gap distribution
def redm_mgap(cat,mem,cenmag,cengalindex,mag,zmin,zmax,
              lm_min,lm_max,
              bootlist,gboot,outdir,
              use_lum=False,use_obs=False,weight_cen=False):

    #Make binning ranges -- initial test ranges
    if use_lum:
        dlum = 0.05
        minlum = np.array(range(80))*dlum-2
        maxlum = minlum+dlum
    else:
        dlum = 0.125
        minlum = np.array(range(80))*dlum-5
        maxlum = minlum+dlum

    nz = len(zmin)
    nlm = len(lm_min[0])
    nlbins = len(minlum)

    #Make array for counting
    gaps = np.zeros([nz, nlm, nlbins])

    #Split up depending on whether using most likely central only
    if weight_cen:
        #Use p_cen
        [gaps_count, gaps_val] = gaps_weight(cat['mem_match_id'],cenmag,cengalindex,
                                             cat['p_cen'],mem['mem_match_id'],mag,
                                             mem['p'],use_lum=use_lum)
    else:
        #Don't use p_cen -- most likely central only
        [gaps_count, gaps_val] = gaps_noweight(cat['mem_match_id'],cenmag,cengalindex,
                                               mem['mem_match_id'],mag,mem['p'],
                                               use_lum=use_lum)

    #Total them up
    for i in range(len(cat)):
        if (cat[i]['z_lambda'] < np.min(zmin)) | (cat[i]['z_lambda'] > np.max(zmax)):
            continue
        if (cat[i]['lambda_chisq'] < np.min(lm_min)) | (cat[i]['lambda_chisq'] > np.max(lm_max)):
            continue
        my_zbin = np.where( (cat[i]['z_lambda'] >= zmin) & (cat[i]['z_lambda'] < zmax))[0][0]
        my_lmbin = np.where( (cat[i]['lambda_chisq'] >= lm_min[my_zbin]) & (cat[i]['lambda_chisq'] < lm_max[my_zbin]))[0]
        my_lbin = np.floor( (gaps_val[i]-np.min(minlum))/dlum ).astype(int)
        mlist = np.where( (my_lbin >= 0) & (my_lbin < nlbins) )[0]
        
        for j in range(len(my_lmbin)):
            for k in range(len(mlist)):
                gaps[my_zbin,my_lmbin[j],my_lbin[mlist[k]]] = gaps[my_zbin,my_lmbin[j],my_lbin[mlist[k]]] + gaps_count[i][mlist[k]]

    #Divide out by the number of clusters
    for i in range(nz):
        for j in range(nlm):
            ncl_bin = len( np.where( (cat['z_lambda'] >= zmin[i]) & (cat['z_lambda'] < zmax[i]) & 
                                     (cat['lambda_chisq'] >= lm_min[i][j]) & 
                                     (cat['lambda_chisq'] < lm_max[i][j] ) )[0] )
            gaps[i,j,:] = gaps[i,j,:]/ncl_bin/dlum

    #Print the outputs
    for i in range(nz):
        for j in range(nlm):
            mgap_print(minlum+dlum/2.,gaps[i,j,:],
                       outdir+"mgap_z_"+str(zmin[i])+"_"+str(zmax[i])+"_lm_"+
                       str(lm_min[i,j])+"_"+str(lm_max[i,j])+".dat")
            
    return
