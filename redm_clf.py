#!/u/ki/dapple/bin/python

##!/usr/bin/env python

import sys

import pyfits
import numpy as np

from collections import Counter
import time
import os
import pz_utils
import matplotlib.pyplot as plt
#Routines for calculating the CLF

#Weights galaxies by p
#-- histogram with essentially CLFs of individual clusters
def count_galaxies_p(c_mem_id,scaleval,g_mem_id,p,mag,lumbins):
    nclusters = len(c_mem_id)
    
    nlum = len(lumbins)
    dlum = lumbins[1]-lumbins[0]
    minlum = lumbins[0]-dlum/2.
    maxlum = lumbins[-1]+dlum/2.

    #Make the array that will contain the full count of galaxies
    #in each cluster
    count_arr = np.zeros([nclusters,nlum])

    #Make the index that takes mem_match_id to cluster
    max_id = np.max(g_mem_id)
    index = np.zeros(max_id+1) - 100
    index[c_mem_id] = range(len(c_mem_id))

    #Select only galaxies that are in the giving luminosity range
    mylist = np.where( (mag <= maxlum) & (mag >= minlum) )[0]
    if len(mylist)==0:
        print >> sys.stderr, "WARNING:  No galaxies found in range ",minlum,maxlum
        return count_arr

    #Loop over galaxies, counting each in turn
    for i in range(len(mylist)):
        #Catch for case where cluster mem_match_id too high to be in list
        if g_mem_id[mylist[i]] > max_id:
            continue
        mycluster = int(index[g_mem_id[mylist[i]]])
       # print(mycluster)
        mybin = np.floor((mag[mylist[i]]-minlum)/dlum)
        mybin = mybin.astype(int)
        if (mybin < 0) | (mybin >= nlum) | (mycluster == -100):
            continue
        count_arr[mycluster,mybin] += p[mylist[i]]*scaleval[mycluster]
    
    return count_arr

#Weights galaxies by p and p_cen -- counts centrals only
'''
def count_galaxies_p_cen(cat,cenmag,lumbins,weight_cen=0):
    nclusters = len(cat)
    
    nlum = len(lumbins)
    dlum = lumbins[1]-lumbins[0]
    minlum = lumbins[0]-dlum/2.

    #Make the array that will contain the full count of galaxies
    #in each cluster
    count_arr = np.zeros([nclusters,nlum])

    #Loop over central galaxies, counting each in turn
    ncen=1
    for i in range(nclusters):
        #If weighting, should check all good centrals
        if weight_cen==1:
            ncen = cat['ncent_good'][i]
        mybin = np.floor((cenmag[i,0:ncen]-minlum)/dlum)
        mybin = mybin.astype(int)
        myskip=False
        if ncen==0:
            myskip=True
        if ncen==1:
            if (mybin < 0) | (mybin >= nlum):
                myskip=True
        if ncen>1:
            count = 0
            for j in range(ncen):
                if (mybin[j] < 0) | (mybin[j] >= nlum):
                    count = count+1
            if count == ncen:
                myskip=True
            else:
                sublist = np.where( (mybin >= 0 ) & (mybin < nlum))[0]

        if myskip:
            continue
        
        if ncen==1:
            count_arr[i,mybin] += cat['p_cen'][i][0:ncen]
        else:
            for j in range(len(sublist)):
                count_arr[i,mybin[sublist[j]]] += cat['p_cen'][i][sublist[j]]
        #If not weighting, give best central all probability
        if weight_cen==0:
            count_arr[i,mybin] += cat['p_sat'][i][0:ncen]

    return count_arr
'''
###Chto's implementattion on count_galaxies_rand_cen
##I would say it is much much faster
def count_galaxies_p_cen(cat,cenmag,lumbins,weight_cen=0):
    nlum = len(lumbins)
    nclusters = len(cat)
    dlum = lumbins[1]-lumbins[0]
    minlum = lumbins[0]-dlum/2.
    chto_countArray=np.zeros([len(cat),nlum])
    mybin = np.floor((cenmag[:,:]-minlum)/dlum).astype(np.int)
    p_cen = cat['p_cen']
    if len(p_cen.shape)==1:
       p_cen = p_cen.reshape(-1,1)
    p_sat = cat['p_sat']
    if len(p_sat.shape)==1:
       p_sat = p_sat.reshape(-1,1)
    if weight_cen:
        ncen = np.repeat(cat['ncent_good'].reshape(-1,1),p_cen.shape[1],axis=1)
        weight = p_cen
    else:
        ncen = np.zeros(p_cen.shape)
        ncen = np.hstack(((np.ones(p_cen.shape[0])).reshape(-1,1),ncen[:,:-1])).astype(np.int) 
        weight = p_sat+p_cen
    
    ys = np.outer(np.ones(p_cen.shape[0]),np.arange(p_cen.shape[1]))
    mask = ((mybin>=0) & (mybin<nlum)& (ys<ncen)).astype(np.float64)
    newbin=mybin*mask
    #index=np.arange(nclusters)
    indexArray=np.arange(nlum)
    for i in indexArray:
        masknew=(mybin==i).astype(np.float64)*mask
        newNewbin=np.sum(masknew*weight,axis=1)
        chto_countArray[:,i]=newNewbin[:]
    return chto_countArray

#Includes/excludes galaxies at random, with probability p 
#-- histogram with essentially CLFs of individual clusters
#NOTE THIS FUNCTION NOT IN USE -- do not use without careful forethought
def count_galaxies_rand(c_mem_id,g_mem_id,scaleval,lumbins,mag,select):
    nclusters = len(c_mem_id)
    ngals = len(g_mem_id)

    nlum = len(lumbins)
    dlum = lumbins[1]-lumbins[0]
    minlum = lumbins[0]-dlum/2.
    maxlum = lumbins[-1]+dlum/2.

    #Make the array that will contain the full count of galaxies
    #in each cluster
    count_arr = np.zeros([nclusters,nlum])

    #Make the index that takes mem_match_id to cluster
    #index = np.zeros(np.max(mem['mem_match_id'])+1) - 100
    #index[cat['mem_match_id']] = range(len(cat))
    #Check to get the max number of times mem_match_id repeats
    mem_id  = Counter(c_mem_id)
    rep_max = mem_id.most_common(1)[0][1]

    #First value in each row is # of valid clusters with that mem_match_id
    index  = np.zeros([np.max(g_mem_id)+1,rep_max+1],long)-1
    index[:,0] = np.zeros(len(index))
    #Loop over clusters to assign to index
    #print len(index)
    for i in range(nclusters):
        index[c_mem_id[i],0] = index[c_mem_id[i],0] + 1
        #if i < 10:
        #    print >> sys.stderr, i, cat['mem_match_id'][i], index[cat['mem_match_id'][i],0]
        #    print >> sys.stderr, "    ",index[cat['mem_match_id'][i]]
        #if index[cat['mem_match_id'][i],0] > rep_max:
        #    print >> sys.stderr, "Oops",i,rep_max,cat['mem_match_id'][i],index[cat['mem_match_id'][i],0]
        #    print >> sys.stderr, index[cat['mem_match_id'][i]]
        index[c_mem_id[i],index[c_mem_id[i],0]] = i
    
    #Select only galaxies that are in the given luminosity range
    #And which are randomly selected
    mylist = np.where( (mag <= maxlum) & (mag >= minlum) & (select < mem['p']))[0]
    if len(mylist)==0:
        print >> sys.stderr, "WARNING:  No galaxies found in range ",minlum,maxlum
        return count_arr

    print >> sys.stderr, "Done making index",len(mylist)
    #Loop over galaxies, counting each in turn
    start = time.clock()
    mybin = np.floor((mag[mylist]-minlum)/dlum)
    mybin = mybin.astype(int)
    for i in range(len(mylist)):
        #if i==1000:
        #    print >> sys.stderr, i, time.clock()-start, (time.clock()-start)/1000, len(mylist)
        #if i == 10000:
        #    print >> sys.stderr, i, time.clock()-start, (time.clock()-start)/10000,len(mylist)
        #if i==100000:
        #    print >> sys.stderr, i, time.clock()-start, (time.clock()-start)/100000, len(mylist)
        if (mybin[i] < 0) | (mybin[i] >= nlum):
            continue
        mycluster = index[g_mem_id[mylist[i]]]
        if mycluster[0] == 0:
            continue
        count_arr[mycluster[1:mycluster[0]],mybin[i]] = count_arr[mycluster[1:mycluster[0]],mybin[i]] + scaleval[mycluster[1:mycluster[0]]]
                                                                                                      
    print >> sys.stderr, "Time for counting: ", time.clock()-start, len(mylist),(time.clock()-start)/len(mylist)
    return count_arr

#Alternative random galaxy counting setup, which runs on only one bootstrap sample at a time
def count_galaxies_rand_all(bootlist, gboot, c_mem_id, scaleval, g_mem_id, lumbins, mag, p):
    nclusters = len(c_mem_id)
    ngals = len(g_mem_id)

    nlum = len(lumbins)
    dlum = lumbins[1]-lumbins[0]
    minlum = lumbins[0]-dlum/2.
    maxlum = lumbins[-1]+dlum/2.

    count_arr = np.zeros([nclusters, nlum])

    #Pick out what bin is needed for each galaxy
    mylist = np.where( (mag < maxlum) & (mag >= minlum))[0]
    if len(mylist)==0:
        print >> sys.stderr, "WARNING:  No galaxies found in range ",minlum,maxlum
        return count_arr
    mybin = np.floor((mag-minlum)/dlum)
    mybin = mybin.astype(int)

    #All galaxies listed for each cluster in gboot are included with weight 1
    #So, loop over clusters first
    start = time.clock()
    for i in range(nclusters):
        glist = gboot[i]
        if len(glist)==0:
            continue
        binlist = np.where( (mybin[glist] > 0) & (mybin[glist] < nlum) )[0]
        if len(binlist)==0:
            continue
        glist = glist[binlist]
        for j in range(len(glist)):
            count_arr[i,mybin[glist[j]]] = count_arr[i,mybin[glist[j]]] + scaleval[bootlist[i]]

    print >> sys.stderr, "Time to complete loop: ",time.clock()-start,(time.clock()-start)/len(mylist)

    return count_arr


#Includes/excludes galaxies at random, with probability p and p_cen 
#-- histogram with essentially CLFs of individual clusters
#centrals only; takes input of which galaxies are "selected" in sample
#extremely slow need to rewrite
##
'''
def count_galaxies_rand_cen(cat,cenmag,cengalindex,lumbins,weight_cen=0):
    nclusters = len(cat)
    
    nlum = len(lumbins)
    dlum = lumbins[1]-lumbins[0]
    minlum = lumbins[0]-dlum/2.

    #Make the array that will contain the full count of galaxies
    #in each cluster
    count_arr = np.zeros([nclusters,nlum])

    #randoms for choosing which central
    cen_select = np.random.uniform(size=nclusters)
    p_cen = cat['p_cen']
    if len(p_cen.shape)==1:
       p_cen = p_cen.reshape(-1,1)

    #Loop over central galaxies, counting each in turn
    ncen=1
    for i in range(nclusters):
        #If weighting, should check all good centrals
        if weight_cen==1:
            ncen = cat['ncent_good'][i]
        mybin = np.floor((cenmag[i,0:ncen]-minlum)/dlum)
        mybin = mybin.astype(int)
        myskip=False
        if ncen==0:
            myskip=True
        if ncen==1:
            if (mybin < 0) | (mybin >= nlum):
                myskip=True
        if ncen>1:
            count = 0
            for j in range(ncen):
                if (mybin[j] < 0) | (mybin[j] >= nlum):
                    count = count+1
            if count == ncen:
                myskip=True
            else:
                sublist = np.where( (mybin >= 0 ) & (mybin < nlum))[0]

        if myskip:
            continue

        #If weighting, do a random selection for the central
        if cen_select[i] > np.sum(p_cen[i][0:ncen]):
            continue
        for j in range(ncen):
            if (cen_select[i] < np.sum(p_cen[i][0:j+1])) & (mybin[j] >= 0) & (mybin[j] < nlum):
                count_arr[i,mybin[j]] += 1.
                break
    
    return count_arr
'''
###Chto's implementattion on count_galaxies_rand_cen
##I would say it is much much faster
def count_galaxies_rand_cen(cat,cenmag,cengalindex,lumbins,weight_cen=0):
    nlum = len(lumbins)
    nclusters = len(cat)
    cen_select=np.random.uniform(size=len(cat))
    dlum = lumbins[1]-lumbins[0]
    minlum = lumbins[0]-dlum/2.
    chto_countArray=np.zeros([len(cat),nlum])
    p_cen = cat['p_cen']
    if len(p_cen.shape)==1:
       p_cen=p_cen.reshape(-1,1)
    cumCat_p=np.cumsum(p_cen,axis=1)
    mybin = np.floor((cenmag[:,:]-minlum)/dlum).astype(np.int)
    if weight_cen:
        ncen = np.repeat(cat['ncent_good'].reshape(-1,1),p_cen.shape[1],axis=1)
    else:
        ncen = np.repeat(np.ones(p_cen.shape[0]).reshape(-1,1),p_cen.shape[1],axis=1)
    
    cen_select2 = np.repeat(cen_select.reshape(-1,1),p_cen.shape[1],axis=1)
    ys = np.outer(np.ones(p_cen.shape[0]),np.arange(p_cen.shape[1]))

    mask = ((cumCat_p >= cen_select2)&(mybin>0) & (mybin<nlum)& (ys<ncen)).astype(np.float64)
    mask2 = np.cumsum(mask, axis=1)
    mask3 = np.cumsum(mask2, axis=1)
    mask3[np.where(mask3>1)]=0.
    newbin=np.sum((mybin*mask3), axis=1)
    mask_new=np.sum(mask3,axis=1)[:, np.newaxis]
    indexArray=np.hstack((np.arange(0,p_cen.shape[0]).reshape(-1,1),newbin.reshape(-1,1))).astype(np.int)
    chto_countArray[indexArray[:,0],indexArray[:,1]]=1.
    chto_countArray =  chto_countArray*mask_new
    return chto_countArray


########

#Function for calculating the number of clusters above a given limiting magnitude
def cluster_Lcount(lumbins,limmag,use_lum=0,p=[]):
    nclusters_lum = np.zeros_like(lumbins)
    binlist = np.zeros_like(limmag).astype(long)
    dlum = lumbins[1]-lumbins[0]
    minlum = lumbins[0]-dlum/2.
    nlumbins = len(lumbins)

    if len(p) == 0:
        p = np.zeros_like(limmag)+1
    #print(binlist)
    for i in range(len(limmag)):
        if use_lum == 1:
            mybin = np.ceil( (limmag[i] - minlum)/dlum )
            mybin = mybin.astype(int)
            binlist[i] = nlumbins-1
            if mybin <= 0:
                nclusters_lum = nclusters_lum + p[i]
                binlist[i] = 0
            if (mybin > 0) & (mybin < nlumbins):
                nclusters_lum[mybin:] = nclusters_lum[mybin:] + p[i]
                binlist[i] = mybin
        else:
            mybin = np.floor( (limmag[i] - minlum)/dlum )
            mybin = mybin.astype(int)
            binlist[i] = 0
            if mybin > nlumbins:
                nclusters_lum = nclusters_lum + p[i]
                binlist[i] = nlumbins
            if (mybin > 0) & (mybin <= nlumbins):
                nclusters_lum[:mybin] = nclusters_lum[:mybin] + p[i]
                binlist[i] = mybin

    return nclusters_lum, binlist

#Function for summing up a count_arr to make a single CLF for specified ranges
#In cluster properties
def make_single_clf(lm,z,lumbins,count_arr,lm_min,lm_max,zmin,zmax,
                    limmag=[],use_lum=0,p=[]):
    dlum = lumbins[1]-lumbins[0]
    clf = np.zeros_like(lumbins)

    if len(p)==0:
    #Version using peak, single-value of P(z) distribution
        clist = np.where( (z >= zmin) & (z < zmax) & (lm >= lm_min) & (lm < lm_max) )[0]
    
        if len(clist) == 0:
            print >> sys.stderr, "WARNING: no clusters found for limits of: "
            print >> sys.stderr, lm_min,lm_max,zmin,zmax
            return clf
    
        #Add CLFs for all relevant clusters
        if len(limmag)==0:
            clf = np.sum(count_arr[clist,:],0)/dlum
            clf = clf/len(clist)
        else:
            #Use limiting magnitudes to determine effective number of clusters in each LF bin
            [nclusters_lum, binlist] = cluster_Lcount(lumbins,np.array(limmag[clist]),use_lum=use_lum)
            #Total the CLF only where the galaxy count for the cluster is complete
            for i in range(len(clist)):
                if use_lum==1:
                    clf[binlist[i]:] = clf[binlist[i]:] + count_arr[clist[i],binlist[i]:]
                else:
                    clf[:binlist[i]] = clf[:binlist[i]] + count_arr[clist[i],:binlist[i]]
                
            #Now divide by number of clusters that are okay in each bin
            clf = clf/nclusters_lum/dlum

    else:
    #Version where p(in z bin) and p(in lambda bin) is given using p
        clist = np.where( p>0.01 )[0]
        
        if len(clist) == 0:
            print >> sys.stderr, "WARNING: no clusters found for limits of: "
            print >> sys.stderr, lm_min,lm_max,zmin,zmax
            return clf
        #Add CLFs for all relevant clusters
        if len(limmag)==0:
            for i in range(len(clist)):
                clf = clf + count_arr[clist[i],:]*p[clist[i]]
            clf = clf/np.sum(p[clist])/dlum
        else:
            #Use limiting magnitudes to determine effective number of clusters in each LF bin
            [nclusters_lum, binlist] = cluster_Lcount(lumbins,limmag[clist],use_lum=use_lum,p=p[clist])
            #Total the CLF only where the galaxy count for the cluster is complete
            for i in range(len(clist)):
                if use_lum==1:
                    clf[binlist[i]:] = clf[binlist[i]:] + count_arr[clist[i],binlist[i]:]*p[clist[i]]
                else:
                    clf[:binlist[i]] = clf[:binlist[i]] + count_arr[clist[i],:binlist[i]]*p[clist[i]]
                
            #Now divide by number of clusters that are okay in each bin
            clf = clf/nclusters_lum/dlum
        

    return clf
    
#PlotCLF
def plot_clf(lm_min,lm_max,zmin,zmax,indir,outdir="plots/"):
    os.system("mkdir -p "+indir+outdir)
    nz = len(zmin)
    nlambda = len(lm_min[0])
    #print(nz,nlambda,lm_min.shape,zmin.shape,lm_max.shape,zmax.shape)
    for i in range(nz):
        for j in range(nlambda):
            #cengal
            clf_cen_name = "clf_cen_z_"+str(zmin[i])+"_"+str(zmax[i])+"_lm_"+str(lm_min[i,j])[0:5]+"_"+str(lm_max[i,j])[0:5]+".dat"
            covar_cen_name = "clf_cen_covar_z_"+str(zmin[i])+"_"+str(zmax[i])+"_lm_"+str(lm_min[i,j])[0:5]+"_"+str(lm_max[i,j])[0:5]+".dat"
            clf_cen=np.loadtxt(indir+clf_cen_name)
            clf_cov_cen=np.loadtxt(indir+covar_cen_name)
            #satgal
            clf_sat_name = "clf_sat_z_"+str(zmin[i])+"_"+str(zmax[i])+"_lm_"+str(lm_min[i,j])[0:5]+"_"+str(lm_max[i,j])[0:5]+".dat"
            covar_sat_name = "clf_sat_covar_z_"+str(zmin[i])+"_"+str(zmax[i])+"_lm_"+str(lm_min[i,j])[0:5]+"_"+str(lm_max[i,j])[0:5]+".dat"
            clf_sat=np.loadtxt(indir + clf_sat_name)
            clf_cov_sat=np.loadtxt(indir + covar_sat_name)
            def determinecolor(matrix):
                mask = np.diag(matrix)< np.median(np.diag(matrix))+5*np.std(np.diag(matrix))
                colormax = np.median(np.diag(matrix)[mask])+3*np.std(np.diag(matrix)[mask])
                return colormax
            def setfigure():
                fig=plt.figure(figsize=(5,5))
                ax = plt.gca();
                ax.set_xticks(np.arange(9, 12, 0.5));
                ax.set_yticks(np.arange(9, 12, 0.5));
                plt.xlabel(r'$L[L_\odot]$')
                plt.ylabel(r'$L[L_\odot]$')
                extent = (np.min(clf_cen[:,0]),np.max(clf_cen[:,0]),np.max(clf_cen[:,0]),np.min(clf_cen[:,0]))
                return extent
            colormax = determinecolor(clf_cov_cen)
            extent=setfigure()
            plt.imshow(clf_cov_cen,vmax=colormax, extent=extent)
            plt.colorbar(label=r"$\Phi^2[(dlogL)^{-2}]$")
            plt.savefig(indir+outdir+covar_cen_name[:-4]+".png")
            
            colormax = determinecolor(clf_cov_sat)
            extent=setfigure()
            plt.imshow(clf_cov_sat,vmax=colormax, extent=extent)
            plt.colorbar(label=r"$\Phi^2[(dlogL)^{-2}]$")
            plt.savefig(indir+outdir+covar_sat_name[:-4]+".png")
            #plt.clf()
            

#Prints CLFs and Covar to file for given redshift and lambda bin
def print_clf_covar(lumbins,clf,covar,outfile,covarfile):
    nlum = len(lumbins)
    f = open(outfile,'w')
    for i in range(nlum):
        f.write(str(lumbins[i])+" "+str(clf[i])+" "+str(np.sqrt(covar[i,i]))+"\n")
    f.close()

    f = open(covarfile,'w')
    for i in range(nlum):
        for j in range(nlum):
            f.write(str(covar[i,j])+" ")
        f.write("\n")
    f.close()
    return

#Testing output -- prints all bootstrap samples
def print_boot_test(lumbins,clf,outfile):
    nboot = len(clf)
    nlum = len(lumbins)
    f = open(outfile,'w')
    for i in range(nlum):
        f.write(str(lumbins[i])+" ")
        for j in range(nboot):
            f.write(str(clf[j,i])+" ")
        f.write("\n")
    f.close()
    return

#Main CLF function that does all requested operations and also
#outputs CLF results to files
def redm_clf(cat,mem,mag,cenmag,cengalindex,lm_min,lm_max,zmin,zmax,
             pcen_all,bootlist,gboot,match_index,
             outdir,weight_cen=0,obs_clf=0,use_lum=0,
             limmag=[], trough=False, stellarMass=False,  jackknifeerr=0):

    #First, define the limits in magnitude/luminosity
    #Default magnitudes
    print mag
    lumbins = np.array(range(35))*0.2-25.
    if (obs_clf==1) & (use_lum==0):
        #Observed magnitudes -- not recommended
        lumbins = np.array(range(60))*0.2+12
    if use_lum==1:
        lumbins = np.array(range(35))*0.08+9
    if stellarMass:
        lumbins = np.array(range(35))*0.08+9.8
    nlum = len(lumbins)

    nlambda = len(lm_min[0])
    nz = len(zmin)
    nboot = len(bootlist)
    nclusters = len(cat)

    #Make the empty CLF arrays
    cenclf = np.zeros([nz,nlambda,nlum])
    satclf = np.zeros([nz,nlambda,nlum])

    #Set up all necessary counting arrays for the main CLF
    #Counting all galaxies
    sat_count_arr = count_galaxies_p(cat['mem_match_id'],cat['scaleval'],mem['mem_match_id'],
                                 mem['p']*(1-pcen_all),mag,lumbins)
#    np.save("/u/ki/chto100/code/desy3workshop/clf/makeclf/chtolog/sat_count.py", sat_count_arr)
#    print "done satcount arr\n\n\n\n\n\n"
    #Counting centrals
    cen_count_arr = count_galaxies_p_cen(cat,cenmag,lumbins,weight_cen=weight_cen)

    #Set up arrays to collect the probability of being in the given redshift, lambda bin
    #These will be used later with the bootstraps, to avoid repeating computations
#    p_zbin = np.zeros([nz,nclusters])
#    p_lbin = np.zeros([nz,nlambda,nclusters])
    
    #Total up the CLF for each bin in lambda and redshift
    for i in range(nz):
        #Get probability that each cluster is in this bin
#        p_zbin[i] = pz_utils.p_in_zbin(cat['pz'],cat['pzbins'],zmin[i],zmax[i])

        for j in range(nlambda):
            #Get probability that the cluster is in this richness bin
#            p_lbin[i,j] = pz_utils.p_in_lmbin(cat['lambda_chisq'],cat['lambda_chisq_e'],
#                                                lm_min[i,j],lm_max[i,j])
#            start=time.clock()
            cenclf[i,j] = make_single_clf(cat['lambda_chisq'],cat['z_lambda'],
                                          lumbins,cen_count_arr,lm_min[i,j],lm_max[i,j],
                                          zmin[i],zmax[i],limmag=limmag,
                                          use_lum=use_lum)
            satclf[i,j] = make_single_clf(cat['lambda_chisq'],cat['z_lambda'],
                                          lumbins,sat_count_arr,lm_min[i,j],lm_max[i,j],
                                          zmin[i],zmax[i],limmag=limmag,
                                          use_lum=use_lum)
#            if (i>1) and (j==1):
#               print satclf[i,j]
#            end=time.clock() 
#            print("total time use for single_clf: {0} s".format(end-start))
    #And print out the CLFs
    for i in range(nz):
        for j in range(nlambda):
            f = open(outdir+"clf_cen_z_"+str(zmin[i])+"_"+str(zmax[i])+"_lm_"+
                     str(lm_min[i,j])[0:5]+"_"+str(lm_max[i,j])[0:5]+".dat",'w')
            for k in range(nlum):
                f.write(str(lumbins[k])+" "+str(cenclf[i,j,k])+"\n")
            f.close()
            f = open(outdir+"clf_sat_z_"+str(zmin[i])+"_"+str(zmax[i])+"_lm_"+
                     str(lm_min[i,j])[0:5]+"_"+str(lm_max[i,j])[0:5]+".dat",'w')
            for k in range(nlum):
                f.write(str(lumbins[k])+" "+str(satclf[i,j,k])+"\n")
            f.close()
    
    #If nboot == 0, call it done
    if nboot == 0:
        return
            
    print >> sys.stderr, "Beginning CLF covariance calculations..."
    #Now, run through all of the bootstrap samples to make errors and covariance matrices    
    cenclf_boot = np.zeros([nboot,nz,nlambda,nlum])
    satclf_boot = np.zeros([nboot,nz,nlambda,nlum])
    #Use the all-counting function to count galaxies
    if jackknifeerr:
       for i in range(nboot):

           start=time.clock()
           print >> sys.stderr, "running on jacknife: {0}".format(i)
           gboot = pz_utils.getjackgal(bootlist[i], cat['mem_match_id'], mem['mem_match_id'], match_index)
           print(np.where(cat['mem_match_id'][bootlist[i]]>np.max(mem['mem_match_id'][gboot])))
           sat_count_arr_b = count_galaxies_p(cat['mem_match_id'][bootlist[i]],
                                 cat['scaleval'][bootlist[i]],mem['mem_match_id'][gboot],
                                 mem['p'][gboot]*(1-pcen_all[gboot]),mag[gboot],lumbins)
           #Counting centrals
           cen_count_arr_b = count_galaxies_p_cen(cat[bootlist[i]],cenmag[bootlist[i]],lumbins,weight_cen=weight_cen) 

           if len(limmag)==0:
               my_limmag = []
           else:
               my_limmag = limmag[bootlist[i]]

           for j in range(nz):
               for k in range(nlambda):                
                   cenclf_boot[i,j,k] = make_single_clf(cat['lambda_chisq'][bootlist[i]],
                                                        cat['z_lambda'][bootlist[i]],
                                                        lumbins,cen_count_arr_b,lm_min[j,k],lm_max[j,k],
                                                        zmin[j],zmax[j],limmag=my_limmag,
                                                        use_lum=use_lum)
                   satclf_boot[i,j,k] = make_single_clf(cat['lambda_chisq'][bootlist[i]],
                                                        cat['z_lambda'][bootlist[i]],
                                                        lumbins,sat_count_arr_b,lm_min[j,k],lm_max[j,k],
                                                        zmin[j],zmax[j],limmag=my_limmag,
                                                        use_lum=use_lum)

           end = time.clock() 
           print("total time used: {0} s".format(end-start))

    else:
       for i in range(nboot):
           print >> sys.stderr, "running on bootlist: {0}".format(i)
#        start=time.clock()
           sat_count_arr_b = count_galaxies_rand_all(bootlist[i],gboot[i], cat['mem_match_id'],
                                                     cat['scaleval'],
                                                     mem['mem_match_id'], lumbins, mag, mem['p']*(1-pcen_all))
#        end=time.clock()
#        print("total time use sat_count: {0} s".format(end-start))
           start=time.clock()
           cen_count_arr_b = count_galaxies_rand_cen(cat[bootlist[i]],cenmag[bootlist[i]],
                                                     cengalindex[bootlist[i]],
                                                     lumbins,weight_cen=weight_cen)
           end = time.clock() 
           print("total time use cen_count: {0} s".format(end-start))
           if len(limmag)==0:
               my_limmag = []
           else:
               my_limmag = limmag[bootlist[i]]

#        start=time.clock()
           for j in range(nz):
               for k in range(nlambda):                
                   cenclf_boot[i,j,k] = make_single_clf(cat['lambda_chisq'][bootlist[i]],
                                                        cat['z_lambda'][bootlist[i]],
                                                        lumbins,cen_count_arr_b,lm_min[j,k],lm_max[j,k],
                                                        zmin[j],zmax[j],limmag=my_limmag,
                                                        use_lum=use_lum)
                   satclf_boot[i,j,k] = make_single_clf(cat['lambda_chisq'][bootlist[i]],
                                                        cat['z_lambda'][bootlist[i]],
                                                        lumbins,sat_count_arr_b,lm_min[j,k],lm_max[j,k],
                                                        zmin[j],zmax[j],limmag=my_limmag,
                                                        use_lum=use_lum)

                   if trough:
                     cenclf_boot[i,j,k]/=np.max(cenclf_boot[i,j,k])
                     satclf_boot[i,j,k]/=np.max(satclf_boot[i,j,k])

#        end = time.clock() 
#        print("total time use single_clf: {0} s".format(end-start))
       #print satclf[0,0],satclf_boot[0,0,0],satclf[0,0]-satclf_boot[0,0,0]
       
       #Covariance matrices -- currenly only within a single z, lambda bin
    covar_cen = np.zeros([nz,nlambda,nlum,nlum])
    covar_sat = np.zeros([nz,nlambda,nlum,nlum])
    meanboot_cen = np.zeros([nz,nlambda,nlum])
    meanboot_sat = np.zeros([nz,nlambda,nlum])
    for i in range(nz):
        for j in range(nlambda):
            for k in range(nlum):
                for l in range(nlum):
                    if jackknifeerr:
                        mean_boot_cen_k = np.mean(cenclf_boot[:,i,j,k])
                        mean_boot_cen_l = np.mean(cenclf_boot[:,i,j,l])
                        mean_boot_sat_k = np.mean(satclf_boot[:,i,j,k])
                        mean_boot_sat_l = np.mean(satclf_boot[:,i,j,l])
                        meanboot_cen[i,j,k] = mean_boot_cen_k
                        meanboot_sat[i,j,k] = mean_boot_sat_k
                        covar_cen[i,j,k,l] = np.sum( ( cenclf_boot[:,i,j,k] - mean_boot_cen_k )*
                                                 ( cenclf_boot[:,i,j,l] - mean_boot_cen_l) )/(nboot)*(nboot-1)
                        covar_sat[i,j,k,l] = np.sum( ( satclf_boot[:,i,j,k] - mean_boot_sat_k )*
                                                 ( satclf_boot[:,i,j,l] - mean_boot_sat_l ) )/(nboot)*(nboot-1)
                    else:

                        covar_cen[i,j,k,l] = np.sum( ( cenclf_boot[:,i,j,k] - cenclf[i,j,k] )*
                                                 ( cenclf_boot[:,i,j,l] - cenclf[i,j,l] ) )/(nboot-1.)
                        covar_sat[i,j,k,l] = np.sum( ( satclf_boot[:,i,j,k] - satclf[i,j,k] )*
                                                 ( satclf_boot[:,i,j,l] - satclf[i,j,l] ) )/(nboot-1.)
            if jackknifeerr:
               f = open(outdir+"clf_cen_jack_mean_z_"+str(zmin[i])+"_"+str(zmax[i])+"_lm_"+
                        str(lm_min[i,j])[0:5]+"_"+str(lm_max[i,j])[0:5]+".dat",'w')
               for k in range(nlum):
                  f.write(str(lumbins[k])+" "+str(meanboot_cen[i,j,k])+"\n")
               f.close()
               f = open(outdir+"clf_sat_jack_mean_z_"+str(zmin[i])+"_"+str(zmax[i])+"_lm_"+
                     str(lm_min[i,j])[0:5]+"_"+str(lm_max[i,j])[0:5]+".dat",'w')
               for k in range(nlum):
                   f.write(str(lumbins[k])+" "+str(meanboot_sat[i,j,k])+"\n")
               f.close() 
         #And print the resulting covariance matrices, and CLFs with errors
            outfile = outdir+"clf_cen_z_"+str(zmin[i])+"_"+str(zmax[i])+"_lm_"+str(lm_min[i,j])[0:5]+"_"+str(lm_max[i,j])[0:5]+".dat"
            covarfile = outdir+"clf_cen_covar_z_"+str(zmin[i])+"_"+str(zmax[i])+"_lm_"+str(lm_min[i,j])[0:5]+"_"+str(lm_max[i,j])[0:5]+".dat"
            print_clf_covar(lumbins,cenclf[i,j],covar_cen[i,j],outfile,covarfile)

            outfile = outdir+"clf_sat_z_"+str(zmin[i])+"_"+str(zmax[i])+"_lm_"+str(lm_min[i,j])[0:5]+"_"+str(lm_max[i,j])[0:5]+".dat"
            covarfile = outdir+"clf_sat_covar_z_"+str(zmin[i])+"_"+str(zmax[i])+"_lm_"+str(lm_min[i,j])[0:5]+"_"+str(lm_max[i,j])[0:5]+".dat"
            print_clf_covar(lumbins,satclf[i,j],covar_sat[i,j],outfile,covarfile)

         #PRINT TESTING RESULTS
         #outfile = outdir+"test_sat_z_"+str(zmin[i])+"_"+str(zmax[i])+"_lm_"+str(lm_min[i,j])[0:5]+"_"+str(lm_max[i,j])[0:5]+".dat"
         #print_boot_test(lumbins,satclf_boot[:,i,j,:],outfile)

    return
