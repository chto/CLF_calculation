#!/u/ki/dapple/bin/python

##!/usr/bin/env python

import sys

import numpy as np

#Routines for doing conversion to abundance-matched lambda limits

def abm_limit_single(lambda1,z1,area1,lambda2,z2,area2,zmin,zmax,lcut):
    '''
    Routine for abundance matching sample 1 to one lambda limit 
    applied to sample 2
    Assumes that lambda1 has already been sorted on lambda
    Inputs:  lambda1 z1 area1 lambda2 z2 area2 zmin zmax lcut
    '''
    
    nclusters = len(np.where( (z2>=zmin) & (z2<zmax) & (lambda2>=lcut) )[0])
    nval = int(np.floor(nclusters/area2*area1))
    list1 = np.where( (z1 >= zmin) & (z1<zmax) )[0]
    lcut_abm = lambda1[list1[nval]]

    return lcut_abm

def abm_limits(cat,cat_abm,area,abm_area,zmin,zmax,lcut):
    lambda1 = np.copy(cat['lambda_chisq'])
    z1 = np.copy(cat['z_lambda'])

    lambda2 = np.copy(cat_abm['lambda_chisq'])
    z2 = np.copy(cat_abm['z_lambda'])
    
    ncuts = len(lcut)
    nz = len(zmin)
    lcut_abm = np.zeros([nz,ncuts])

    sortlist = np.argsort(lambda1)
    sortlist = sortlist[::-1]
    lambda1 = lambda1[sortlist]
    z1 = z1[sortlist]

    for i in range(nz):
        for j in range(ncuts):
            zmid=(zmin[i]+zmax[i])/2.
            lcut_abm[i,j] = abm_limit_single(lambda1,z1,area(zmid),lambda2,z2,
                                             abm_area(zmid),zmin[i],zmax[i],lcut[j])

    return lcut_abm
