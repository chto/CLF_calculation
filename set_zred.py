#!/u/ki/dapple/bin/python

##!/usr/bin/env python

import os
import sys

import pyfits
import numpy as np

#Various functions for convering cluster catalogs for using zred and
#lambda_zred instead of z_lambda, lambda_chisq
def set_zred(cat,mem):
    #Set catalog values
    #Should always use z_lambda value -- better cluster photoz
    #cat['z_lambda'][:] = cat['zred'][:]
    cat['lambda_chisq'][:] = cat['lambda_zred'][:]

    #Set cluster member values if included
    index = np.zeros(max(cat['mem_match_id'])+1,dtype=int)-100
    index[cat['mem_match_id']] = range(len(cat))
    #mem['z'][:] = cat['z_lambda'][index[mem['mem_match_id']]][:]

    return cat, mem
