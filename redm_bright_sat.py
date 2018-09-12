#!/usr/bin/python

import numpy as np
import sys

def get_brightest_satellite(pcen,cindex,mag,p,start_index,
                            lumbins,use_lum=1):
    '''
    Function that finds the brightest satellite, probabilistically, for a single
    input cluster and list of associated galaxies
    '''
    nlum = len(lumbins)
    dlum = lumbins[1]-lumbins[0]
    minlum = lumbins[0]-dlum/2.
    maxlum = lumbins[-1]+dlum/2.

    clf_bsat = 0*lumbins

    slist = np.argsort(mag)
    if use_lum == 1:
        slist = slist[::-1]

    ngals = np.min([len(slist),7])

    #Match each of the ngals galaxy to a central, if possible
    pcen_sat = np.zeros(ngals)
    for i in range(ngals):
        cplace = np.where(start_index + slist[i] == cindex)[0].astype(int)
        if len(cplace) > 0:
            pcen_sat[i] = pcen[cplace]


    pvals = np.zeros(ngals)
    for i in range(ngals):
        mybin = np.floor((mag[slist[i]] - minlum)/dlum)
        mybin = mybin.astype(int)
        
        myp = p[slist[i]]*(1 - pcen_sat[i])

        if i > 0:
            myp = myp*np.prod(1-p[slist[:i]]+p[slist[:i]]*pcen_sat[:i] )
            
        pvals[i] = myp
        #print i, myp, mybin, pcen_sat[i],np.prod(1-p[slist[:i]]+p[slist[:i]]*pcen_sat[:i] )
        if (mybin < 0) | (mybin > nlum-1):
            continue
        clf_bsat[mybin] = clf_bsat[mybin]+myp

    return clf_bsat

def get_brightest_satellite_all(cat,mem,mag,cindex,lm_min,lm_max,zmin,zmax,
                                bootlist,gboot,match_index,
                                outdir,weight_cen=0,obs_clf=0,use_lum=1):
    '''
    Function that loops over all clusters to get the distribution for the
    brightest satellite galaxy; outputs to files

    Does calculate covariance matrices
    '''

    #First, define the limits in magnitude/luminosity
    #Default magnitudes
    lumbins = np.array(range(35))*0.2-25.
    if (obs_clf==1) & (use_lum==0):
        #Observed magnitudes -- not recommended
        lumbins = np.array(range(60))*0.2+12
    if use_lum==1:
        lumbins = np.array(range(35))*0.08+9
    nlum = len(lumbins)
    dlum = lumbins[1]-lumbins[0]

    nlambda = len(lm_min[0])
    nz = len(zmin)
    nboot = len(bootlist)
    nclusters = len(cat)
    ngals = len(mem)

    #Make the empty brightest sat CLF array
    bsatclf = np.zeros([nz,nlambda,nlum])
    nboot = len(bootlist)

    #Set up the initial counting array -- note that this loops over clusters
    count_arr = np.zeros([nclusters,nlum])
    count = 0L
    count_hi = 0L
    for i in range(nclusters):
        #if i % 100 == 0:
            #print i, " Running counting..."
        if weight_cen==1:
            if len(cat['p_cen'].shape)==1:
               pcen=np.array([cat[i]['p_cen']])
            else:
               pcen = cat[i]['p_cen']
        else:
            pcen = np.array([1])
        
        #Getting the matching list of satellites
        while cat['mem_match_id'][i] > mem['mem_match_id'][count]:
            count = count+1
        count_hi = count
        while cat['mem_match_id'][i] == mem['mem_match_id'][count_hi]:
            count_hi = count_hi + 1
            if count_hi >= ngals:
                break

        count_arr[i] = get_brightest_satellite(pcen,cindex[i],mag[count:count_hi],
                                               mem['p'][count:count_hi],count,
                                               lumbins,use_lum=use_lum)

    #Now, total up the count
    #print lm_min[0][0]
    for i in range(nz):
        for j in range(nlambda):
            clist = np.where( (cat['lambda_chisq'] >= lm_min[i][j]) & (cat['lambda_chisq'] < lm_max[i][j]) & (cat['z_lambda'] >= zmin[i]) & (cat['z_lambda'] < zmax[i]) )[0]
            bsatclf[i,j] = np.sum( count_arr[clist], axis=0 )/float(len(clist))/dlum
            #print i, j, len(clist)

    #Now get the results for each bootstrap sample
    bsatcov = np.zeros([nz,nlambda,nlum,nlum])
    bsat_temp = np.zeros([nboot,nz,nlambda,nlum])
    for i in range(nboot):
        for j in range(nz):
            for k in range(nlambda):
                clist = bootlist[i][np.where( (cat['lambda_chisq'][bootlist[i]] >= lm_min[j][k]) & (cat['lambda_chisq'][bootlist[i]] < lm_max[j][k]) & (cat['z_lambda'][bootlist[i]] >= zmin[j]) & (cat['z_lambda'][bootlist[i]] < zmax[j]) )[0]]
                bsat_temp[i,j,k] = np.sum(count_arr[clist], axis=0)/float(len(clist))/dlum

    #And now calculate the covariance matrix
    for i in range(nz):
        for j in range(nlambda):
            for k in range(nlum):
                for l in range(nlum):
                    bsatcov[i,j,k,l] = np.sum( (bsat_temp[:,i,j,k] - bsatclf[i,j,k])*(bsat_temp[:,i,j,l]-bsatclf[i,j,l]) )/(nboot-1.)


    #Printout section
    for i in range(nz):
        for j in range(nlambda):
            f = open(outdir+"clf_sat_bright_z_"+str(zmin[i])+"_"+str(zmax[i])+"_lm_"+
                     str(lm_min[i][j])[0:5]+"_"+str(lm_max[i][j])[0:5]+".dat",'w')
            for k in range(nlum):
                f.write(str(lumbins[k])+" "+str(bsatclf[i,j,k])+" "+str(np.sqrt(np.diag(bsatcov[i,j])[k]))+"\n")
            f.close()

    #And output the full covariance matrices
    for i in range(nz):
        for j in range(nlambda):
            f = open(outdir+"clf_sat_bright_covar_z_"+str(zmin[i])+"_"+str(zmax[i])+"_lm_"+
                     str(lm_min[i][j])[0:5]+"_"+str(lm_max[i][j])[0:5]+".dat",'w')
            for k in range(nlum):
                for l in range(nlum):
                    f.write(str(bsatcov[i,j,k,l])+" ")
                f.write("\n")
            f.close()

    return #lumbins, bsatclf, count_arr


#This version, rather than getting the brightest sat distribution as a function
#of luminosity only, gets the distribution of probability in Lcen, Lsat(brightest) space
def get_bright_sat_cen(pcen,cenmag,cindex,mag,p,start_index,
                       lumbins,use_lum=1):
    '''
    Inputs are: pcen, cindex, mag, p, start_index, lumbins

    Runs calculation for a single cluster at a time
    Note that lumbins is used as the binning in both the Lcen and Lsat dimensions
    '''
    nlum = len(lumbins)
    dist = np.zeros([nlum,nlum])
    dlum = lumbins[1]-lumbins[0]
    minlum = lumbins[0]-dlum/2.
    maxlum = lumbins[-1]+dlum/2.

    #Sorting on luminosity
    slist = np.argsort(mag)
    if use_lum==1:
        slist = slist[::-1]

    ngals = np.min([len(slist),7])

    #First, loop over centrals
    ncen = len(pcen)
    for i in range(ncen):
        #Get appropriate luminosity bin for the central
        cbin = int(np.floor((cenmag[i]-minlum)/dlum))
        if (cbin < 0) | (cbin > nlum-1):
            continue

        #check to see if any of the galaxies match this central
        cplace = np.where(start_index + slist == cindex[i])[0]
        if len(cplace) == 0:
            cplace = -1
        else:
            cplace = cplace[0]

        #Now loop over the satellites, giving a probability for each
        sbin = np.floor((mag[slist]-minlum)/dlum).astype(int)
        for j in range(ngals):
            #Skip if the particular satellite is the central
            if j == cplace:
                continue
            #Skip if the satellite is too bright/dim to count
            if (sbin[j] < 0) | (sbin[j] > nlum-1):
                continue
            if j == 0:
                myp = pcen[i]*p[slist[0]]
            else:
                glist = np.where( ( np.array(range(ngals)) != cplace) & (np.array(range(ngals)) < j))[0]
                myp = pcen[i]*np.prod(1-p[slist[glist]])*p[slist[j]]
            #print cbin, sbin, nlum
            dist[cbin,sbin[j]] = dist[cbin,sbin[j]] + myp
            #print i, j, pcen[i], p[slist[j]], myp, cbin, sbin[j], cplace, ngals
    
    #print "TEST: ",np.sum(dist)

    return dist

#Run the distribution for central/brightest satellite luminosity for all clusters
def get_bright_sat_cen_all(cat,mem,mag,cindex,cenmag,lm_min,lm_max,zmin,zmax,
                           bootlist,gboot,match_index,
                           outdir,weight_cen=0,obs_clf=0,use_lum=1):
    '''
    Function that loops over all clusters to get the joint distribution of the
    central and brightest satellite galaxies; outputs to files

    Does not yet calculate covariance matrices
    '''

    #First, define the limits in magnitude/luminosity
    lumbins = np.array(range(35))*0.2-25.
    if (obs_clf==1) & (use_lum==0):
        #Observed magnitudes -- not recommended
        lumbins = np.array(range(60))*0.2+12
    if use_lum==1:
        lumbins = np.array(range(35))*0.08+9
    nlum = len(lumbins)
    dlum = lumbins[1]-lumbins[0]

    nlambda = len(lm_min[0])
    nz = len(zmin)
    nboot = len(bootlist)
    nclusters = len(cat)
    ngals = len(mem)

    #Make the emtpy array
    dist = np.zeros([nz,nlambda,nlum,nlum])
    nboot = len(bootlist)

    #Set up the initial counting array -- note that this loops over clusters
    count_arr = np.zeros([nclusters,nlum,nlum])
    count = 0L
    count_hi = 0L
    for i in range(nclusters):
        if weight_cen==1:
            if len(cat['p_cen'].shape)==1:
               pcen=np.array([cat[i]['p_cen']])
            else:
               pcen = cat[i]['p_cen']
        else:
            pcen = np.array([1])
        
        #Getting the matching list of satellites
        while cat['mem_match_id'][i] > mem['mem_match_id'][count]:
            count = count+1
        count_hi = count
        while cat['mem_match_id'][i] == mem['mem_match_id'][count_hi]:
            count_hi = count_hi + 1
            if count_hi >= ngals:
                break
        count_arr[i] = get_bright_sat_cen(pcen,cenmag[i],cindex[i],mag[count:count_hi],
                                          mem['p'][count:count_hi],count,
                                          lumbins,use_lum=use_lum)
        
    #Now, total everything up as desired
    for i in range(nz):
        for j in range(nlambda):
            clist = np.where( (cat['lambda_chisq'] >= lm_min[i][j]) & (cat['lambda_chisq'] < lm_max[i][j]) & (cat['z_lambda'] >= zmin[i]) & (cat['z_lambda'] < zmax[i]) )[0]
            dist[i,j] = np.sum( count_arr[clist], axis=0 )/float(len(clist))/dlum/dlum

    #Print the outputs
    #Note that this prints lumbins and the 2D distribution matrix
    for i in range(nz):
        for j in range(nlambda):
            f = open(outdir+"dist_cen_sat_bright_z_"+str(zmin[i])+"_"+str(zmax[i])+"_lm_"+
                     str(lm_min[i][j])[0:5]+"_"+str(lm_max[i][j])[0:5]+".dat",'w')
            for k in range(nlum):
                f.write(str(lumbins[k])+" ")
                for l in range(nlum):
                    f.write(str(dist[i,j,k,l])+" ")
                f.write("\n")
            f.close()

    #Make the matrix of errors (not bother with full covariance matrix right now -- it's a mess)
    #First, get the bootstrap sample results
    dist_err = np.zeros([nz,nlambda,nlum,nlum])
    dist_temp = np.zeros([nboot,nz,nlambda,nlum,nlum])
    for i in range(nboot):
        for j in range(nz):
            for k in range(nlambda):
                clist = bootlist[i][np.where( (cat['lambda_chisq'][bootlist[i]] >= lm_min[j][k]) & (cat['lambda_chisq'][bootlist[i]] < lm_max[j][k]) & (cat['z_lambda'][bootlist[i]] >= zmin[j]) & (cat['z_lambda'][bootlist[i]] < zmax[j]) )[0]]
                dist_temp[i,j,k] = np.sum( count_arr[clist], axis=0 )/float(len(clist))/dlum/dlum

    for i in range(nz):
        for j in range(nlambda):
            for k in range(nlum):
                for l in range(nlum):
                    dist_err[i,j,k,l] = np.sum( (dist[i,j,k,l] - dist_temp[:,i,j,k,l] )*(dist[i,j,k,l] - dist_temp[:,i,j,k,l] ) )/(nboot-1.)

    #And print the error estimates
    for i in range(nz):
        for j in range(nlambda):
            f = open(outdir+"dist_err_cen_sat_bright_z_"+str(zmin[i])+"_"+str(zmax[i])+"_lm_"+
                     str(lm_min[i][j])[0:5]+"_"+str(lm_max[i][j])[0:5]+".dat",'w')
            for k in range(nlum):
                f.write(str(lumbins[k])+" ")
                for l in range(nlum):
                    f.write(str(dist_err[i,j,k,l])+" ")
                f.write("\n")
            f.close()

    return count_arr
