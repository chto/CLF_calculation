#!/u/ki/dapple/bin/python

import os
import sys
import numpy as np
import scipy.integrate
import scipy.special
import pyfits

import pz_utils

#Functions for correcting the p_ext probabilities; 
#Correction mostly to account for presence of blue galaxies
#in reality
#See Eduardo's work for more details

#Note that this may be updated in the future
def pext_correct(p_old, chi2):
    '''
    Currently, p_old (initial p_ext) is the only input

    Calculates delta_p based on fixed parameter set

    May later add L dependency?

    Note that this requires input of chi^2 for each galaxy
    Note that this requires read-in of epsilon_functions.fit file

    Includes correction for p_new > 1 or p_new < 0
    '''

    #Old version
    #Parameters -- drawn from GAMA
    #p0 = 0.00494789
    #p1 = 0.655880
    #p2 = -2.86675
    #p3 = 2.19239

    #Change in probabilities
    #delta = p0 + p1*p_old + p2*p_old*p_old + p3*p_old**3

    #p_new = p_old + delta

    #Updated, current version
    #Read in fit data
    if (os.path.isfile("epsilon_functions.fit")==False):
        print >> sys.stderr, "WARNING: Unable to perform correction to probabilities."
        print >> sys.stderr, "Running without correction."
        return p_old

    dat = pyfits.open("epsilon_functions.fit")
    dat = dat[1].data

    #Interpolation to obtain the epsilon values
    lnchi2 = np.log(chi2)
    e_chi2 = np.interp(lnchi2,dat['lnchisq'],dat['epschi'])
    e_blue = np.interp(lnchi2,dat['lnchisq'],dat['epsblue'])

    #Constants
    #c = 0.061
    #delta = -0.0168
    #Updated constant
    #c = 0.066
    #delta = -0.02123
    #Constants updated again -- and again
    c = 0.070
    delta = -0.075

    #And putting it all together -- see Rozo et al for details
    p_new = p_old * (1+e_chi2)/(1+p_old*e_chi2) / (1+p_old*c) / (1+p_old*e_blue) * (1+delta) / (1+p_old*delta)

    #Checks for impossibly high/low p_new
    alist = np.where(p_new > 1)[0]
    if len(alist) > 0:
        p_new[alist] = 0*alist + 1
    blist = np.where(p_new < 0)[0]
    if len(blist) > 0:
        p_new[blist] = 0*blist

    return p_new

def pext_correct_partial(p_old):
    c = 0.070

    #First correction for correlated structure
    p_next = p_old * 1./(1+p_old*c)

    #Estimate the delta correction
    delta = np.sum(p_next)/np.sum(p_old)-1.

    #And correct
    p_new = p_next/(1+p_next*delta)
    return p_new

#Latest correction version -- run with v5.10 redmapper
#Function for conversion so this can run on S82
def remap_three_colors_to_four(chisq):

    nbins = 1000
    x = np.array(range(nbins))*30./999-0.
    
    dof = 4.0
    y=x**(0.5*dof-1.0)*np.exp(-0.5*x)/np.math.gamma(0.5*dof)/2.0**(0.5*dof)	
    dof=3.0
    y3=x**(0.5*dof-1.0)*np.exp(-0.5*x)/np.math.gamma(0.5*dof)/2.0**(0.5*dof)
    cdf4=0.0*y
    cdf3=cdf4
    
    for i in np.array(range(nbins-1))+1:
        cdf4[i] = scipy.integrate.simps(y[0:i+1],x[0:i+1])
        cdf3[i] = scipy.integrate.simps(y3[0:i+1],x[0:i+1])
    xeq = np.interp(cdf3,cdf4,x)

    chisq_new = np.interp(chisq,x,xeq)

    return chisq_new

#Making correction section
def pext_correct_blue(p_old,cat_chisq,cat_chisq_mod,r,ncolors=4.):
    if ncolors != 4. and ncolors != 3.:
        print >> sys.stderr, "WARNING: Only able to correct for n=3 or n=4 colors, not ",ncolors
        return p_old
    
    struct_blue = pyfits.open("epsilon_blue.fit")
    struct_blue = struct_blue[1].data

    xblue = struct_blue['lnchisq']
    epsblue = struct_blue['epsblue']
    epsc = 0.045
    delta = -0.093

    s = np.log(cat_chisq/cat_chisq_mod)
    bad = np.where(np.isfinite(s)==False)[0]
    if len(bad) > 0:
        s[bad] = 0*bad

    x = np.array(range(1000))*14/999.-6.
    dof = ncolors
    chisq = np.exp(x)
    y=1.0/2.0**(dof/2.0)/np.math.gamma(dof/2.0)*chisq**(dof/2.0)*np.exp(-chisq/2.0)
    
    xvec = np.log(cat_chisq_mod)
    bad = np.where(cat_chisq_mod <= np.exp(-6))[0]
    if len(bad) > 0:
        xvec[bad] = -6
    xvec1 = xvec+s

    rho = np.exp(-s)*np.interp(xvec,x,y)
    rho0 = np.interp(xvec1,x,y)

    bad = np.where(rho < 0)[0]
    if len(bad) > 0:
        rho[bad] = np.min(y)
    bad = np.where(rho0 < 0)[0]
    if len(bad) > 0:
        rho0[bad] = np.min(y)

    gal_epschi=rho/rho0-1.0

    bad = np.where(gal_epschi >= 1e4)[0]
    if len(bad) > 0:
        gal_epschi[bad] = 1e4
    #Remap if ncolors=3
    if ncolors==3:
        cat_chisq = remap_three_colors_to_four(cat_chisq)
        
    gal_epsblue = np.interp(np.log(cat_chisq+1e-5),xblue,epsblue)
    
    b = delta + gal_epschi+delta*gal_epschi
    a = gal_epsblue + epsc

    p_new = p_old*(1.+b)/(1.+p_old*(a+b+a*b))
    #Checking for centrals -- this part likely to change
    bad = np.where(r < 1e-5)[0]
    if len(bad) > 0:
        p_new[bad] = 1.0

    return p_new

def predfunc_fiducial(x):
    mu = 2.44 + 1.00
    sigma = 0.28+0.00
    
    y = (x-mu)/sigma
    y = 0.5*(1.0-scipy.special.erf(y/np.sqrt(2)))

    return y

#Most updated version of probability correction; requires predfunc_fiducial
def pext_correct_analytic(p_old,cat_chisq,cat_chisq_mod,r,ncolors=4.):
    epsc = 0.045
    delta = -0.093

    s = np.log(cat_chisq/cat_chisq_mod)
    bad = np.where(np.isfinite(s)==False)[0]
    if len(bad) > 0:
        s[bad] = 0.

    x = np.array(range(1000))*14/999.-6.
    dof = ncolors
    chisq = np.exp(x)
    y=1.0/2.0**(dof/2.0)/np.math.gamma(dof/2.0)*chisq**(dof/2.0)*np.exp(-chisq/2.0)

    xvec = np.log(cat_chisq_mod)
    bad = np.where(cat_chisq_mod <= np.exp(-6))[0]
    if len(bad) > 0:
        xvec[bad] = -6
    xvec1 = xvec+s

    rho = np.exp(-s)*np.interp(xvec,x,y)
    rho0 = np.interp(xvec1,x,y)

    bad = np.where(rho < 0)[0]
    if len(bad) > 0:
        rho[bad] = np.min(y)
    bad = np.where(rho0 < 0)[0]

    gal_epschi=rho/rho0-1.0

    bad = np.where(gal_epschi >= 1e4)[0]
    if len(bad) > 0:
        gal_epschi[bad] = 1e4
    #Remap if ncolors=3
    if ncolors==3:
        cat_chisq_mod = remap_three_colors_to_four(cat_chisq_mod)

    fblue=1.0-predfunc_fiducial(np.log(cat_chisq_mod))
    gal_epsblue=fblue/(1.0-fblue)
    bad=np.where(1.0-fblue < 1e-5)[0]
    if len(bad) > 0:
        gal_epsblue[bad]=1e5	

    b=delta+gal_epschi+delta*gal_epschi
    a=gal_epsblue+epsc

    p_new = p_old*(1.0+b)/(1.0+p_old*(a+b+a*b))
    bad = np.where(r < 1e-5)[0]
    if len(bad) > 0:
        p_new[bad] = 1.0

    return p_new


#Overall pext_correct routine -- handles routine switching
#This is the part that should be called by redm_full
def pext_correct_full(cat,mem,pcen_all,use_p_ext,ncolors):
    if use_p_ext == 0:
        return mem['p']
    if use_p_ext < 4:
        #Include correction to p_ext values, based on GAMA
        p_ext = pext_correct(mem['p_ext'],mem['chisq'])
        
        #excluding those with r > r_lambda if requested (use_p_ext==2)
        p = pz_utils.dimmer_rlambda_p(use_p_ext,
                                      cat['mem_match_id'],
                                      cat['lambda_chisq'],cat['r_lambda'],
                                      mem['mem_match_id'],
                                      mem['p'],p_ext,mem['r'])
        return p
    if use_p_ext == 4:
        #New redmapper format
        p = mem['p'][:]*mem['pfree'][:]
        p = pext_correct(p,mem['chisq'])
        #Cut to only keep galaxies with r<r_lambda
        p = pz_utils.dimmer_rlambda_p_new(use_p_ext,
                                          cat['mem_match_id'],
                                          cat['lambda_chisq'],
                                          cat['r_lambda'],
                                          mem['mem_match_id'],
                                          p,mem['r'])
        return p
    if use_p_ext == 5:
        #New redmapper format, partial correction for S82
        p = mem['p'][:]*mem['pfree'][:]
        #Run the partial correction -- correlated structure
        p = pext_correct_partial(p)
        #Cut to keep only galaxies with r<r_lambda
        p = pz_utils.dimmer_rlambda_p_new(use_p_ext,
                                          cat['mem_match_id'],
                                          cat['lambda_chisq'],cat['r_lambda'],
                                          mem['mem_match_id'],
                                          p,mem['r'])
        return p
    if use_p_ext == 6:
        #Includes color corrections
        #Note this WILL fail if run on catalog without chisq_mod tag
        m_names = mem.columns.names
        fail = True
        for name in m_names:
            if name == 'CHISQ_MOD':
                fail = False
                break
        if fail:
            print >> sys.stderr, "WARNING: Unable to locate necessary tag chisq_mod.  Exitting..."
            print >> sys.stderr, "         No correction made"
            return mem['p']
        #Found the tag okay, so let's go
        p = pext_correct_blue(mem['p'],mem['chisq'],
                              mem['chisq_mod'],mem['r'],
                              ncolors=ncolors)
        #Cut to keep only galaxies with r<r_lambda
        p = pz_utils.dimmer_rlambda_p_new(use_p_ext,
                                          cat['mem_match_id'],
                                          cat['lambda_chisq'],cat['r_lambda'],
                                          mem['mem_match_id'],
                                          p,mem['r'])
        p = p*mem['pfree']

        #Get the central galaxy satellite membership
        clist = np.where(mem['r'] < 1e-5)[0]
        if len(clist) > 0:
            p[clist] = mem['p'][clist]*(1-pcen_all[clist])
            
        return p

    if use_p_ext == 7:
        #More recent version; is purely analytic
        #Requires a catalog with the chisq_mod tag
        m_names = mem.columns.names
        fail = True
        for name in m_names:
            if name == 'CHISQ_MOD':
                fail = False
                break
        if fail:
            print >> sys.stderr, "WARNING: Unable to locate necessary tag chisq_mod.  Exitting..."
            print >> sys.stderr, "         No correction made"
            return mem['p']
        #Found the tag okay, so let's go
        p = pext_correct_analytic(mem['p'],mem['chisq'],
                                  mem['chisq_mod'],mem['r'],
                                  ncolors=ncolors)
        #And do the radius trimming
        p = pz_utils.dimmer_rlambda_p_new(use_p_ext,
                                          cat['mem_match_id'],
                                          cat['lambda_chisq'],cat['r_lambda'],
                                          mem['mem_match_id'],
                                          p,mem['r'])

        p = p*mem['pfree']
        #Get central galaxy satellite membership
        clist = np.where(mem['r'] < 1e-5)[0]
        if len(clist) > 0:
            p[clist] = mem['p'][clist]*(1-pcen_all[clist])

        return p

    if use_p_ext == 8:
        #Simple version, that does NOT perform more than a basic correction;
        #Includes extension of galaxies to dim values, but no chisq correction
        p = mem['p']*mem['pfree']
        p = pz_utils.dimmer_rlambda_p_new(0,cat['mem_match_id'],
                                 cat['lambda_chisq'],cat['r_lambda'],
                                 mem['mem_match_id'],
                                 p,mem['r'])
        #Get central galaxy satellite membership
        clist = np.where(mem['r'] < 1e-5)[0]
        if len(clist) > 0:
            p[clist] = mem['p'][clist]*(1-pcen_all[clist])

    print >> sys.stderr, "Failed to match use_p_ext -- skipping correction"
    return mem['p']
