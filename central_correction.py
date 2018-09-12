#!/usr/bin/python

import numpy as np
import pyfits

import esutil

def get_cindex_and_cenmag(id_cent,g_id,g_imag):
    cindex = np.zeros_like(id_cent)
    cenmag = np.zeros_like(cindex).astype(float)
    #Hash table taking galaxy ID to index
    g_index = np.zeros(np.max(g_id)+1).astype(long)-1
    g_index[g_id] = np.array(range(len(g_imag)))
    for i in range(len(cindex[0])):
        cindex[:,i] = g_index[id_cent[:,i]]
        cenmag[:,i] = g_imag[g_index[id_cent[:,i]]]

    return cindex, cenmag

def get_sindex(cindex,ncent_good,g_imag):
    sindex = np.array(range(len(g_imag)))
    
    #Exclude first two centrals
    sindex[cindex[:,0]] = 0*cindex-1
    slist = np.where(ncent_good >=2)[0]
    sindex[cindex[:,1][slist]] = 0*slist-1

    sindex = sindex[np.where(sindex > 0)[0]]
    satmag = g_imag[sindex]

    return sindex, satmag

def correct_dr8_cen(param,z):
    '''
    Empirical correction of DR8 most likely centrals magnitude
    '''
    x = param[0]*(z-0.15)+param[1]
    xlist = np.where(x > 0)[0]
    if len(xlist) > 0:
        x[xlist] = 0*xlist
    return x

if __name__ == "__main__":
    #Main run sequence for estimating corrections to DR8 luminosities
    
    c_dr8 = pyfits.open("/nfs/slac/g/ki/ki19/des/erykoff/clusters/process/dr8_redmapper_v5.10/run_ubermem/dr8_run_redmapper_v5.10_lgt5_catalog.fit")
    g_dr8 = pyfits.open("/nfs/slac/g/ki/ki19/des/erykoff/clusters/process/dr8_redmapper_v5.10/run_ubermem/dr8_run_redmapper_v5.10_lgt5_catalog_members_mod.fit")
    c_dr8 = c_dr8[1].data
    g_dr8 = g_dr8[1].data

    c_s82 = pyfits.open("/nfs/slac/g/ki/ki19/des/erykoff/clusters/process/stripe82_redmapper_v5.10/run_ubermem/stripe82_run_redmapper_v5.10_lgt5_catalog.fit")
    g_s82 = pyfits.open("/nfs/slac/g/ki/ki19/des/erykoff/clusters/process/stripe82_redmapper_v5.10/run_ubermem/stripe82_run_redmapper_v5.10_lgt5_catalog_members_mod.fit")
    c_s82 = c_s82[1].data
    g_s82 = g_s82[1].data

    print "Done with read-ins"
    
    #Pick out lists of all galaxies
    #DR8
    cindex_dr8,cenmag_dr8 = get_cindex_and_cenmag(c_dr8['id_cent'],
                                                  g_dr8['id'],g_dr8['imag'])

    #S82
    cindex_s82,cenmag_s82 = get_cindex_and_cenmag(c_s82['id_cent'],
                                                  g_s82['id'],g_s82['imag'])

    print "Done getting cindex, cenmag data"
    
    #Setup matching between most likely centrals
    h = esutil.htm.HTM(11)
    m1, m2, d12 = h.match(g_dr8[cindex_dr8[:,0]]['ra'],
                          g_dr8[cindex_dr8[:,0]]['dec'],
                          g_s82[cindex_s82[:,0]]['ra'],
                          g_s82[cindex_s82[:,00]]['dec'],
                          2./3600)

    #Match the second most likely centrals
    clist_dr8 = np.where(c_dr8['ncent_good'] >= 2)[0]
    clist_s82 = np.where(c_s82['ncent_good'] >= 2)[0]    
    m1_c2, m2_c2, d12_c2 = h.match(g_dr8[cindex_dr8[clist_dr8,1]]['ra'],
                                   g_dr8[cindex_dr8[clist_dr8,1]]['dec'],
                                   g_s82[cindex_s82[clist_s82,1]]['ra'],
                                   g_s82[cindex_s82[clist_s82,1]]['dec'],
                                   2./3600)
    m1_c2 = clist_dr8[m1_c2]
    m2_c3 = clist_s82[m2_c2]

    print "Centrals matching done, now working with satellites"

    #Make the sindex, which pulls only satellite galaxies (not 1st or 3nd central)
    sindex_dr8, satmag_dr8 = get_sindex(cindex_dr8,c_dr8['ncent_good'],
                                        g_dr8['imag'])

    sindex_s82, satmag_s82 = get_sindex(cindex_s82,c_s82['ncent_good'],
                                        g_s82['imag'])
    
    #Trim down -- can't handle so many galaxies at once
    mcut = 19. #Apparent magnitude cut
    slist = np.where(satmag_dr8 < mcut)[0]
    sindex_dr8 = sindex_dr8[slist]
    satmag_dr8 = satmag_dr8[slist]

    slist = np.where(satmag_s82 < mcut)[0]
    sindex_s82 = sindex_s82[slist]
    satmag_s82 = satmag_s82[slist]

    print "Doing satellite matching..."

    #And now run the matching
    m1_s, m2_s, d12_s = h.match(g_dr8[sindex_dr8]['ra'],g_dr8[sindex_dr8]['dec'],
                                g_s82[sindex_s82]['ra'],g_s82[sindex_s82]['dec'],
                                2./3600)
