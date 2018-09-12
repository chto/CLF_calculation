#!/u/ki/dapple/bin/python

import os
import sys

import numpy as np

#Functions to output a list of the N most massive clusters in
#each redshift bin (default N=50)
def redm_bigcount(lambda_chisq,z,zmin,zmax,outdir,ncl=50):
    nz = len(zmin)

    #Sorted lambda values
    slist = np.argsort(lambda_chisq)
    slist = slist[::-1]

    #Run for each redshift range
    for i in range(nz):
        clist = slist[np.where( (z[slist]>zmin[i]) & (z[slist] < zmax[i]) )[0]]

        #Print the resulting list to a file
        outfile = outdir+"lambda_max_list_z_"+str(zmin[i])+"_"+str(zmax[i])+".dat"
        f = open(outfile,'w')
        for j in range(ncl):
            print >> f, lambda_chisq[clist[j]]
        f.close()

    return
