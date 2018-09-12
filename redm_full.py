#!/usr/bin/env python
import dodotune
import os
import sys
import gc
import pyfits
import numpy as np

#various local packages
import cosmo
from get_central_mag import *
from mag_convert import *
import set_zred
from redm_nz_calc import *
from redm_nlambda_calc import *
from redm_bigcount import *
from abm_limits import *
import pext_correct
from central_correction import correct_dr8_cen
from get_limmag_5sigma import get_limmag_5sigma
import healpy as hp
from glob import glob
import redm_clf
import redm_rpr
import redm_mgap
import redm_pbcg
import redm_bright_sat
import copy
from scipy.sparse import dok_matrix
from scipy import interpolate
import kmeans_radec 
from kmeans_radec import KMeans, kmeans_sample
####For Trough by chto####
def IndexToDeclRa(ipix,nside):
    theta, phi = hp.pix2ang(nside, ipix)
    ra = np.rad2deg(phi)
    dec = np.rad2deg(0.5 * np.pi - theta)
    return [ra, dec]
def DeclRaToIndex(ra,dec,nside):
    theta = 0.5 * np.pi - np.deg2rad(dec)
    phi = np.deg2rad(ra)
    return hp.ang2pix(nside, theta, phi)



#####

#reads parameter file, outputs dict of values
#Also checks that has all necessary parameters
def read_redm_param(filename):
    keys = []
    vals = []

    f = open(filename)
    for line in f:
        if line.strip() == '':
            continue
        if line[0] == '#':
            continue
        entries = line.split()
        keys.append(entries[0])
        if len(entries[1:]) == 1:
            vals.append(entries[1])
        else:
            vals.append(entries[1:])
    f.close()


    #print keys, vals
    #flag variable declaration
    skip_rpr = False
    
    #First check: if repeated keys, then fail w/ error
    setup = []
    for i in range(len(keys)):
        if keys.count(keys[i]) > 1:
            print >> sys.stderr, "ERROR:  Input value "+keys[i]+" is repeated."
            sys.exit(1)
        setup.append((keys[i],vals[i]))
    #make the preliminary dictionary
    params = dict(setup)

    #Format check on all keys
    if keys.count('cluster_file') == 0:
        print >> sys.stderr, "ERROR:  No input cluster_file"
        sys.exit(1)
    if keys.count('member_file') == 0:
        print >> sys.stderr, "ERROR:  No input member_file"
        sys.exit(1)
    if keys.count('kcorr_file') == 0:
        print >> sys.stderr, "ERROR:  No input kcorr_file"
        sys.exit(1)
    if keys.count('cindex_file') == 0:
        print >> sys.stderr, "ERROR:  No input cindex_file"
        sys.exit(1)
    if keys.count('area') == 0:
        print >> sys.stderr, "ERROR:  No input area"
        sys.exit(1)
    if keys.count('outdir') == 0:
        print >> sys.stderr, "ERROR:  No output directory outdir"
        sys.exit(1)
    
    if keys.count('use_des') == 0:
        print >> sys.stderr, "WARNING:  Default use_des to 0"
        setup.append(('use_des','0'))

    if keys.count('stellarmass') == 0:
        print >> sys.stderr, "WARNING:  Default stellarmass to 0"
        setup.append(('stellarmass',0))



    if keys.count('no_uband') == 0:
        print >> sys.stderr, "WARNING:  Default no_uband to 1"
        setup.append(('no_uband','1'))
    if keys.count('use_zred') == 0:
        print >> sys.stderr, "WARNING:  Default use_zred to 0"
        setup.append(('use_zred','0'))
    if keys.count('zcut_max') == 0:
        setup.append(('zcut_max','3'))
    if keys.count('LRG') == 0:
        print >> sys.stderr, "WARNING:  Default LRG to 1"
        setup.append(('LRG','1'))
    if keys.count('bandshift') == 0:
        print sys.stderr, "WARNING:  Default bandshift is z=0.1"
        setup.append(('bandshift','0.1'))
    if keys.count('abs_solar') == 0:
        print >> sys.stderr, "WARNING:  Default abs_solar is for SDSS i-band at z=0.3"
        #setup.append(('abs_solar','4.71493'))
        setup.append(('abs_solar','4.71718'))
    if keys.count('use_scaleval') == 0:
        print >> sys.stderr, "WARNING:  Default use_scaleval to 1"
        setup.append(('use_scaleval','1'))
    if keys.count('use_p_ext') == 0:
        print >> sys.stderr, "WARNING:  Default is to not use p_ext (use_p_ext=0)"
        setup.append(('use_p_ext','0'))
    if keys.count('ncolors') == 0:
        print >> sys.stderr, "WARNING:  Default number of colors is 4"
        print >> sys.stderr, "          Tag does nothing unless use_p_ext>=6"
        setup.append(('ncolors','4'))
    if keys.count('use_dr8_cen_corr') == 0:
        print >> sys.stderr, "WARNING:  Default is NOT to use central mag correction"
        print >> sys.stderr, "          use_dr8_cen_corr = 0"
        setup.append(('use_dr8_cen_corr','0'))
    if keys.count('p_offset') == 0:
        print >> sys.stderr, "WARNING:  Default probability offset is zero"
        setup.append(('p_offset','0.'))
    if keys.count('use_lambda_err') == 0:
        print >> sys.stderr, "WARNING: Default is to not use lambda errors"
        setup.append(('use_lambda_err','0'))
    if keys.count('do_nz') == 0:
        print >> sys.stderr, "WARNING:  Default skip do_nz (0)"
        setup.append(('do_nz','0'))
    if keys.count('nz_descale') == 0:
        print >> sys.stderr, "WARNING:  Default nz_descale (0)"
        setup.append(('nz_descale','0'))
    if keys.count('do_nlambda') == 0:
        print >> sys.stderr, "WARNING:  Default skip do_nlambda (0)"
        setup.append(('do_nlambda','0'))
    if keys.count('zmin') == 0:
        print >> sys.stderr, "ERROR: No min z specified (zmin)"
        sys.exit(1)
    if keys.count('zmax') == 0:
        print >> sys.stderr, "ERROR: No max z specified (zmax)"
        sys.exit(1)
    if len(params['zmin']) != len(params['zmax']):
        print >> sys.stderr, "ERROR: Number of min and max z limits unequal"
        sys.exit(1)
    if keys.count('lm_min') == 0:
        print >> sys.stderr, "ERROR: No min lambda specified (lm_min)"
        sys.exit(1)
    if keys.count('lm_max') == 0:
        print >> sys.stderr, "ERROR: No max lambda specified (lm_max)"
        sys.exit(1)
    if len(params['lm_min']) != len(params['lm_max']):
        print >> sys.stderr, "ERROR: Number of min and max lambda limits unequal"
        sys.exit(1)
    if len(params['use_lum']) == 0:
        print >> sys.stderr, "WARNING: Default use magnitudes (use_lum==0)"
        setup.append(('use_lum','0'))
    if keys.count('ABM') == 0:
        print >> sys.stderr, "WARNING:  Default to no abundance matching (ABM=0)"
        setup.append(('ABM','0'))
    if keys.count('Analband') == 0:
        print >> sys.stderr, "WARNING:  Default is to use i band"
        setup.append(('Analband','i'))

    if keys.count('jackknife_error') == 0:
        print >> sys.stderr, "WARNING:  Default is to do bootstrap error"
        setup.append(('jackknife_error','0'))

    if keys.count('plot_clf') == 0:
        print >> sys.stderr, "WARNING:  Default is not to plot CLF"
        setup.append(('plot_clf','0'))

    if keys.count('use_REFMAG') == 0:
        print >> sys.stderr, "WARNING:  Default use_REFMAG to {0}".format(params['use_des'])
        setup.append(('use_REFMAG',params['use_des']))
    #print params['ABM']
    if (params['ABM'] == '1') & (keys.count('abm_file') == 0):
        print >> sys.stderr, "ERROR:  To do abundance matching, require abm_file"
        sys.exit(1)
    if (params['ABM'] == '1') & (keys.count('abm_area') == 0):
        print >> sys.stderr, "ERROR:  To do abundance matching, require abm_area"
        sys.exit(1)
    if keys.count('do_clf')==0:
        print >> sys.stderr, "WARNING:  Default skip do_clf (0)"
        setup.append(('do_clf','0'))
    if keys.count('weight_cen')==0:
        print >> sys.stderr, "WARNING:  Default not use weight_cen (0)"
        setup.append(('weight_cen','0'))
    params = dict(setup)
    if keys.count('do_env')==0:
        print >> sys.stderr, "WARNING:  Default skip do_env (0)"
        setup.append(('do_env','0'))
    if keys.count('env_proj')==0:
        setup.append(('env_proj','1'))
    if keys.count('env_rmax')==0:
        setup.append(('env_rmax','100'))
    if keys.count('do_rpr')==0:
        print >> sys.stderr, "WARNING:  Default skip do_rpr (0)"
        setup.append(('do_rpr','0'))
    #Make adjustments
    if keys.count('rpr_minlum')==0 or keys.count('rpr_maxlum')==0:
        print >> sys.stderr, "WARNING:  Require rpr_minlum, rpr_maxlum, skipping radial profiles"
        skip_rpr = True
    if keys.count('do_mgap')==0:
        print >> sys.stderr, "WARNING:  Default skip do_mgap (0)"
        setup.append(('do_mgap','0'))
    if keys.count('do_pbcg') == 0:
        print >> sys.stderr, "WARNING:  Default skip do_pbcg (0)"
        setup.append(('do_pbcg','0'))
    if keys.count('do_bsat') == 0:
        print >> sys.stderr, "WARNING:  Default skip do_bsat (0)"
        setup.append(('do_bsat','0'))
    if keys.count('troughNames')==0:
        print >> sys.stderr, "Default is not to do trough"
        setup.append(('troughNames',None))
    if keys.count('dohaloCat')==0:
        print >> sys.stderr, "Default is not to do halo catalog"
        setup.append(('dohaloCat',None))

    if keys.count('get_limmag_5sigma')==0:
        print >> sys.stderr, "Default is do limmag 5 sigma correction"
        setup.append(('get_limmag_5sigma',1))
 
 
    #Remake the dictionary as needed
    params = dict(setup)
    #print params


###############
   #For area issue:
    if os.path.isfile(params['area']):
      area=pyfits.open(params['area'])[1]
      areaf=interpolate.interp1d(area.data['z'], area.data['area'],bounds_error=False,fill_value="extrapolate")
      params['area']=areaf
    else:
      area_exclusive=float(params['area'])
      def f_area(x):
         try:
            return np.array([area]*len(x))
         except:
            return area_exclusive
      params['area']=f_area
   
    if (params['ABM'] == '1'):

      if os.path.isfile(params['abm_area']):
         area=pyfits.open(params['abm_area'])[1]
         areaf=interpolate.interp1d(area.data['z'], area.data['area'],bounds_error=False,fill_value="extrapolate")
         params['abm_area']=areaf
      else:
         area_abm_exclusive=float(params['abm_area'])
         def f_abm_area(x):
            try:
               return np.array([area]*len(x))
            except:
               return area_abm_exclusive
         params['abm_area']=f_abm_area



############

    #If minlum/maxlum not supplied, do_rpr is overridden
    if skip_rpr:
        params['do_rpr'] = '0'

    return params

def main(filename, zmin=None,zmax=None):

    #Read the parameter file
    params = read_redm_param(filename)
    if zmin != None and zmax != None:
       params['zmin'] = zmin
       params['zmax'] = zmax
    print >> sys.stderr, "Parameter file read successful"

    #Read in the data
    print >> sys.stderr, "Reading clusters file..."
    if os.path.isfile(params['cluster_file'])==False:
        print >> sys.stderr, "ERROR:  Cluster catalog file "+params['cluster_file']+" not found"
        sys.exit(1)
    cat = pyfits.open(params['cluster_file'])
    cat = cat[1].data
    print >> sys.stderr, len(cat), " clusters read in"
    print >> sys.stderr, np.max(cat['z'])
    print >> sys.stderr, "Reading members file..."
    if os.path.isfile(params['member_file'])==False:
        print >> sys.stderr, "ERROR:  Member catalog file "+params['member_file']+" not found"
        sys.exit(1)
    mem = pyfits.open(params['member_file'])
    mem = mem[1].data
    print >> sys.stderr, len(mem)," galaxies read in"

    #If asked to use zred instead of z_lambda, set up for that
    if int(params['use_zred']) == 1:
        [cat, mem] = set_zred.set_zred(cat,mem)
    
    #If not asked to use lambda_errors in calculations, set them to zero
    if int(params['use_lambda_err']) == 0:
        cat['lambda_chisq_e'][:] = cat['lambda_chisq_e'][:]*0.

    #Also read in data to abundance match to, if requested
    if int(params['ABM']) == 1:
        print >> sys.stderr, "Reading in clusters to use for abundance matching..."
        if os.path.isfile(params['abm_file'])==False:
            print >> sys.stderr, "ERROR:  ABM catalog file "+params['abm_file']+" not found"
            sys.exit(1)
        abm = pyfits.open(params['abm_file'])
        abm = abm[1].data
        print >> sys.stderr, len(abm)," clusters for ABM"
    if params['stellarmass']:
       mag = mem['MSTAR_50']
    else:
    #Read in amag data if it exists, otherwise run kcorrect
       if os.path.isfile(params['kcorr_file'])==False:
           #Running kcorrect -- note this requires running IDL
           #os.system('setenv CATFILE '+params['cluster_file'])
           os.environ['CATFILE'] = params['cluster_file']
           os.environ['MEMFILE'] = params['member_file']
           os.environ['KCORRFILE'] = params['kcorr_file']
           os.environ['USE_DES'] = params['use_des']
           os.environ['NO_UBAND'] = params['no_uband']
           os.environ['LRG'] = params['LRG']
           os.environ['BANDSHIFT'] = params['bandshift']
           os.environ['USE_ZRED'] = params['use_zred']
           os.system('/afs/slac/g/ki/software/idl/idl70/bin/idl < /afs/slac.stanford.edu/u/ki/chto100/code/redmapper_clf/redmapper/rm_kcorr_wrapper.pro')
           if os.path.isfile(params['kcorr_file'])==False:
               print >> sys.stderr, "ERROR:  Kcorrect failed, no file found"
               sys.exit(1)
       kcorr = pyfits.open(params['kcorr_file'])
       kcorr = kcorr[0].data
       print >> sys.stderr, len(kcorr), " galaxies kcorrected"
       if len(kcorr) != len(mem):
           print >> sys.stderr, "ERROR!   Number of galaxies != number of kcorrect values"
           sys.exit(1)

       #Pick up the mags we want for all galaxies
       if (int(params['obs_clf'])==1) & (int(params['use_lum'])==0):
           mag = mem['imag']
       else:
           if int(params['use_REFMAG'])==1:
               if params['Analband'] == "r": 
                  print >> sys.stderr,"using r band"
                  mag = mem['REFMAG'] + (kcorr[:,1] - mem['model_mag'][:,3])
               elif params['Analband'] == "i":
                  mag = mem['REFMAG'] + (kcorr[:,2] - mem['model_mag'][:,3])
               else:
                  print >> sys.stderr, "ERRPR: OOPS I don't know the band you specified (because I am too lazy), which is {0}".format(params['Analband'])
           else:
               print(params['use_des'])
               if len(mem['model_mag'][0])==4:
                   mag = mem['imag'] + (kcorr[:,2] - mem['model_mag'][:,2])
               else:
                   mag = mem['imag'] + (kcorr[:,2] - mem['model_mag'][:,3])
    
    #Pulling limits from the parameters structure
    lm_min = np.array(params['lm_min'])
    lm_max = np.array(params['lm_max'])
    zmin = np.array(params['zmin'])
    zmax = np.array(params['zmax'])
    lm_min = lm_min.astype(float)
    lm_max = lm_max.astype(float)
    zmin = zmin.astype(float)
    zmax = zmax.astype(float)
    #If requested, switch to ABM lambda limits
    if int(params['ABM']) == 1:
        print >> sys.stderr, "Running abundance matching..."
        lm_min = abm_limits(cat,abm,params['area'],
                            params['abm_area'],zmin,zmax,lm_min)
        lm_max = abm_limits(cat,abm,params['area'],
                            params['abm_area'],zmin,zmax,lm_max)
        print "abundance matching: lm_min",lm_min
        print "abundance matching: lm_max",lm_max
    else:
        my_nz = len(zmin)
        lm_min = np.repeat([lm_min],my_nz,axis=0)
        lm_max = np.repeat([lm_max],my_nz,axis=0)
    
    #Get central absolute mags
    #Will only recalcuate index if it doesn't exist; otherwise, read in the index
    #and get the magnitudes from the mag values already available
    #Updated version -- for redmapper v5.10 and greater -- uses indices supplied
    #in the catalog files for finding the necessary index
    c_names = cat.columns.names
    use_id_cent = False
    #Check to see if the central IDs are available
    for name in c_names:
        if name == 'ID_CENT':
            use_id_cent = True
            break
    ''' 
    if use_id_cent:
        cengalindex = np.zeros_like(cat['id_cent'])
        cenmag = np.zeros_like(cengalindex).astype(float)
        #Hash table taking galaxy ID to index
        offset = np.min(mem['id'])
        g_index = np.zeros(np.max(mem['id'])-offset+1)-1
        g_index[mem['id']-offset] = np.array(range(len(mem)))
#modified by chto to handle that cat['id_cent']=0
        for i in range(len(cengalindex[0])):  
            index = np.where(cat['id_cent'][:,i]-offset>=0)
            cengalindex[index,i] = g_index[cat['id_cent'][index,i]-offset]
            cenmag[index,i] = mag[cengalindex[index,i]]
        if int(params['weight_cen'])==0:
            cenmag = cenmag[0]
            cengalindex = cengalindex[0]
        del g_index
    '''
   #############
   #Algorithm changed by Chto
   #The reason is that is np.zeros of a large number will cause mem error in Python
   #############
    if use_id_cent:
        cengalindex = np.zeros_like(cat['id_cent'])
        cenmag = np.zeros_like(cengalindex).astype(float)
        #Hash table taking galaxy ID to index
        offset = np.min(mem['id'])
        g_index = dok_matrix((np.max(mem['id'])-offset+1,1), dtype=np.int)
        g_index[mem['id']-offset] = np.array(range(len(mem)))[:,np.newaxis]
        if len(cenmag.shape)==1:
            index = np.where(cat['id_cent'][:]-offset>=0)[0]
            cengalindex[index] = g_index[cat['id_cent'][index]-offset].toarray().flatten()
            cenmag[index] = mag[cengalindex[index].astype(int)]
            cengalindex=cengalindex.astype(int)
            cengalindex=cengalindex.reshape(-1,1)
            cenmag=cenmag.reshape(-1,1)
        else:
            for i in range(len(cengalindex[0])):
               index = np.where(cat['id_cent'][:,i]-offset>=0)[0]
               cengalindex[index,i] = g_index[cat['id_cent'][index,i]-offset].toarray().flatten()
               cenmag[index,i] = mag[cengalindex[index,i]]
        if int(params['weight_cen'])==0:
            cenmag = cenmag[:,0]
            cengalindex = cengalindex[:,0]
            cengalindex=cengalindex.reshape(-1,1)
            cenmag=cenmag.reshape(-1,1)
        del g_index

        gc.collect()
 #########################
    else:
        #No central IDs found -- doing it the hard way
        if os.path.isfile(params['cindex_file'])==False:
            print >> sys.stderr, "Getting central magnitudes..."
            [cenmag, cengalindex] = get_central_mag(cat,mem,mag,weight_cen=int(params['weight_cen']))
            cengalindex = cengalindex.astype(long)
            hdu = pyfits.PrimaryHDU(cengalindex)
            hdu.writeto(params['cindex_file'])
        else:
            print >> sys.stderr, "Reading in central magnitudes..."
            cengalindex = pyfits.open(params['cindex_file'])
            cengalindex = cengalindex[0].data
            cengalindex = cengalindex.astype(long)
            cenmag = np.zeros_like(cengalindex).astype(float)
            print len(cenmag)
            for i in range(len(cenmag)):
                cenmag[i] = mag[cengalindex[i]]

    #If available, get the limiting magnitude information
    #Currently set for the more generous cut
    use_limmag = False
    for name in c_names:
        if name == 'LIM_LIMMAG':
            use_limmag = True
            limmag = cat['LIM_LIMMAG']
            break
    #Convert limiting magnitude to absolute if needed; based on central
    #galaxy's k-correction
    if use_limmag:
        if params['dohaloCat']:
            limmag = cat['LIM_LIMMAG_DERED']
        else:
            if int(params['get_limmag_5sigma'])==1:
               limmag = get_limmag_5sigma(limmag,cat['lim_exptime'])
            if int(params['use_des'])==1:
               if params['Analband'] == "r": 
                  limmag = limmag - ( cat['MODEL_MAG'][:,1]- cenmag[:,0] )
               elif params['Analband'] == "i":
                  limmag = limmag - ( cat['MODEL_MAG'][:,2]- cenmag[:,0] )
               else:
                  print >> sys.stderr, "ERRPR: OOPS I don't know the band you specified (because I am too lazy), which is {0}".format(params['Analband'])
            else:
               limmag = limmag - ( cat['imag']- cenmag[:,0] )
    else:
        limmag = []
    if params['stellarmass']:
       limmag=[] 

    #For the rest of the calculations, since these are predicated on the
    #given lambda/z thresholds, remove all the clusters we don't care about.
    #This should speed things up significantly/avoid memory issues.
    print >> sys.stderr, "Ultimate max z cut is at: ", float(params['zcut_max'])
    clist = np.where( cat['z_lambda'] < float(params['zcut_max']) )[0]
    cmlist = np.where( mem['z'] < float(params['zcut_max']))[0]
    cat = cat[clist]
    cenmag = cenmag[clist]
    cengalindex = cengalindex[clist]
    if len(limmag)!=0:
         limmag = limmag[clist]
    #Note that we also must trim galaxies
    mem = mem[cmlist]
    if not params['stellarmass']:
      kcorr = kcorr[cmlist]
    idn_list = np.zeros(len(mag))
    idn_list[cmlist] = np.array(range(len(cmlist)))
    mag = mag[cmlist]
    cengalindex = idn_list[cengalindex]
    cengalindex = cengalindex.astype(long)
    del cmlist, clist
    print "cmlist, clist: ", gc.collect()
####
    ##Try to know what is the magnitude. 
    #WARNING: Cuts should be VERY generous, otherwise, may have issues with P(z) tails...

    #Make the main output directory
    os.system("mkdir -p "+params['outdir'])

    #Now that we have our magnitudes, add central corrections if necessary
    if int(params['use_dr8_cen_corr']) == 1:
        if int(params['weight_cen'])==1:
            cenmag[:,0] = cenmag[:,0] + correct_dr8_cen([0.213,-0.08],cat['z_lambda'])
            cenmag[:,1] = cenmag[:,1] + correct_dr8_cen([0.104,-0.036],cat['z_lambda'])
            mag[cengalindex[:,0]] = mag[cengalindex[:,0]] + correct_dr8_cen([0.213,-0.08],mem['z'][cengalindex[:,0]])
            corrlist = np.where(cat['ncent_good'] >=2)[0]
            if len(corrlist) > 0:
                mag[cengalindex[corrlist,1]] = mag[cengalindex[corrlist,1]] + correct_dr8_cen([0.104,-0.036],mem['z'][cengalindex[corrlist,1]])
        else:
            cenmag = cenmag + correct_dr8_cen([0.213,-0.08],cat['z_lambda'])
            mag[cengalindex] = mag[cengalindex] + correct_dr8_cen([0.213,-0.08],cat['z_lambda'])

    #Convert everything to log(L) if requested
    if int(params['use_lum'])==1:
        #Value currently hard-coded to offset for z=0.3 bandshift in SDSS i-band
        abs_solar = float(params['abs_solar'])
        mag = np.log10(mag_to_Lsolar(mag,use_des=int(params['use_des']),abs_solar=abs_solar))
        cenmag = np.log10(mag_to_Lsolar(cenmag,use_des=int(params['use_des']),
                                        abs_solar=abs_solar))
        if use_limmag:
            limmag = np.log10(mag_to_Lsolar(limmag,use_des=int(params['use_des']),
                                            abs_solar=abs_solar))
    np.save(params['outdir']+"limmag.npy", limmag)
    np.save(params['outdir']+"cengalindex.npy",cengalindex)
    print "I save the limiting magnitude.... HAAH"
    print >> sys.stderr, "Finished converting mags to log(Lsolar)"
    np.save(params['outdir']+"cen_mag.npy",cenmag)
#    limmag=[] #so hacky :( this is not my style

    print "magnitude distribution: "
    print np.histogram(mag, bins=np.array(range(35))*0.08+9.3)
    np.save(params['outdir']+"mag.npy", mag)
    #Fix the normaliztion of the p(z) so that triangular integration works okay in pz_utils
    ##Modified by chto@@
    dz = cat['pzbins'][:,1] - cat['pzbins'][:,0]
    weight = np.sum(cat['pz'],axis=1)*dz
    cat['pz'] = cat['pz']/weight[:,None]
    del dz,weight
    print "dz, weight: ", gc.collect()
    print >> sys.stderr, "Done renormalizing P(z)"

    print >> sys.stderr, "Max lambda: ",np.max(cat['lambda_chisq'])

    #Now produce the galaxy samples and a hash table that takes cluster ID to the first
    #of its listed galaxies
    #Create an array that gives p_cen for all member galaxies
    pcen_all = np.zeros(len(mem))
    if int(params['weight_cen'])==0:
        pcen_all[cengalindex] = 0*cengalindex + 1.
    else:
        for i in range(len(cengalindex[0])):
            clist = np.where(cengalindex[:,i] != -1)[0]
            #print len(cengalindex),cengalindex[i][0],clist[0]
            pcen_all[cengalindex[clist,i]] = pcen_all[cengalindex[clist,i]] + cat['p_cen'][clist,i]

    np.save(params['outdir']+"pcen_all.py",pcen_all)

    #Reassign p with p_ext if requested
    if int(params['use_p_ext'])>0:
        mem['p'][:] = pext_correct.pext_correct_full(cat,mem,pcen_all,int(params['use_p_ext']),
                                                     int(params['ncolors']))

    #Add a systematic probability offset if requested
    if float(params['p_offset']) != 0:
        mem['p'][:] = mem['p'][:] + float(params['p_offset'])
        #Fix any objects with p>1 or p<0
        plist = np.where(mem['p'] > 1.)[0]
        if len(plist) > 0:
            mem['p'][plist] = 0.*plist+1.
        plist = np.where(mem['p'] < 0.)[0]
        if len(plist) > 0:
            mem['p'][plist] = 0.*plist

    #Now make the satellite lists    
    print >> sys.stderr, "PCEN: ",np.max(pcen_all),len(cat),len(np.where(pcen_all > 0)[0]), len(np.where(pcen_all > 0.9)[0]), np.sum(pcen_all)

    #Make necessary bootstrap samples
    #Includes redshifts taken from P(z)
    #Only want to do this once
    #Also allows covariance estimates between measuresments, which are
    #Not currently implemented
    #Number of bootstrap samples -- current hard-coded
    nboot = dodotune.njack
    bootlist = pz_utils.make_boot_samples_simple(nboot,cat)
    [match_index, gboot] = pz_utils.make_boot_samples_gal_full(bootlist,cat['mem_match_id'],
                                                               mem['mem_match_id'],mem['p']*(1-pcen_all))

    print >> sys.stderr, "Done setting up bootstrap samples"
    assert params['jackknife_error']
    if float(params['jackknife_error']) !=0:
       print "doing jackknife CLF" 
       njack = dodotune.njack
       jacklist = pz_utils.make_jack_samples_simple(njack,cat)

       match_index_jack = pz_utils.getjackgal(jacklist[0], cat['mem_match_id'], mem['mem_match_id'])
       assert (mem['mem_match_id'][match_index_jack[cat['mem_match_id']][:,0]] == cat['mem_match_id']).all()

#       [match_index_jack, gjack] = pz_utils.make_jack_samples_gal_full(jacklist,cat['mem_match_id'],
#                                                               mem['mem_match_id'],mem['p']*(1-pcen_all))
       print >> sys.stderr, "Done setting up jackknife samples"


    #If requested, calculate n(z)
    if int(params['do_nz'])==1:
        print >> sys.stderr, "Calculating n(z)..."
        redm_nz_calc(cat,params['area'],params['outdir'],bootlist,
                     descale=bool(int(params['nz_descale'])))
        print >> sys.stderr, "Done calculating n(z)"

    #If requested, calculate n(lambda)
    if int(params['do_nlambda'])==1:
        print >> sys.stderr, "Calculating n(lambda)..."
        redm_nlambda_err(cat['lambda_chisq'],cat['lambda_chisq_e'],cat['z_lambda'],cat['pz'],cat['pzbins'],
                         bootlist,params['outdir'],zmin,zmax,params['area'] )
        print >> sys.stderr, "Done calculating n(lambda)"
        
        redm_bigcount(cat['lambda_chisq'],cat['z_lambda'],zmin,zmax,params['outdir'])

        print >> sys.stderr, "Done calculating bonus listing of massive clusters"
        
    print >> sys.stderr,"lm_min: ", lm_min, "lm_max: ", lm_max,"zmin: ", zmin,"zmax: ", zmax
    
    #Calculate the CLF
    if int(params['do_clf'])==1:
        #np.save("mag.npy",mag)
        #np.save("cenmag.npy",cenmag)
        #np.save("cat.npy",cat)
        #np.save("mem.npy", mem)
        print >> sys.stderr, "Calculating CLF..."
        if params["troughNames"] is not None:
            trough_data=np.array([hp.read_map(params['troughNames'].format(i))for i in range(5)])
            nside=hp.npix2nside(trough_data.shape[1])
            def getTroughProb(RaDec):
                ra,dec=RaDec
                return trough_data[:,DeclRaToIndex(ra,dec,nside)]
            positionArray=np.array([cat['RA'],cat['DEC']]).T  
            positionArray_mem=np.array([mem['RA'],mem['DEC']]).T
            troughProb=np.array(map(getTroughProb, positionArray))
            troughProb_mem=np.array(map(getTroughProb, positionArray_mem))
            oldmemP=copy.deepcopy(mem['p'])
            oldcatpcen=copy.deepcopy(cat['p_cen'])
            oldcatpsat=copy.deepcopy(cat['p_sat'])
            for i in range(5):
                mem['p']*=troughProb_mem[:,i]
                cat['p_cen']*=troughProb[:,i].reshape(-1,1)
                cat['p_sat']*=troughProb[:,i].reshape(-1,1)
                outdir=params['outdir']+"trough_{0}/".format(i)
                os.system("mkdir -p "+outdir)
                redm_clf.redm_clf(cat,mem,mag,cenmag,cengalindex,lm_min,lm_max,zmin,zmax,
                           pcen_all,bootlist,gboot,match_index,outdir,
                           weight_cen=int(params['weight_cen']),
                           obs_clf=int(params['obs_clf']),use_lum=int(params['use_lum']),
                           limmag=limmag, trough=True,stellarMass=params['stellarmass'])
                mem['p']=oldmemP
                cat['p_cen']=oldcatpcen
                cat['p_sat']=oldcatpsat
        else:     
            #print limmag
            if float(params['jackknife_error']) !=0:
                print "doing jackknife CLF" 

                redm_clf.redm_clf(cat,mem,mag,cenmag,cengalindex,lm_min,lm_max,zmin,zmax,
                           pcen_all,jacklist,gboot=None,match_index=match_index_jack, 
                           outdir=params['outdir'],
                           weight_cen=int(params['weight_cen']),
                           obs_clf=int(params['obs_clf']),use_lum=int(params['use_lum']),
                           limmag=limmag, stellarMass=params['stellarmass'], jackknifeerr=int(params['jackknife_error']))
            else:
                print "doing bootstrap CLF"
                redm_clf.redm_clf(cat,mem,mag,cenmag,cengalindex,lm_min,lm_max,zmin,zmax,
                           pcen_all,bootlist,gboot,match_index,params['outdir'],
                           weight_cen=int(params['weight_cen']),
                           obs_clf=int(params['obs_clf']),use_lum=int(params['use_lum']),
                           limmag=limmag, stellarMass=params['stellarmass'])
        print >> sys.stderr, "Done calculating CLF"

    #PLot CLF:
    if int(params['plot_clf'])==1:
       redm_clf.plot_clf(lm_min,lm_max,zmin,zmax, indir=params['outdir'])
#    assert False

    #Calculate the radial profiles
    if int(params['do_rpr']) == 1:
        print >> sys.stderr, "Calculating radial profiles..."
        
        #First, get necessary input parameter limits on magnitude
        rpr_minlum = np.array(params['rpr_minlum'])
        rpr_maxlum = np.array(params['rpr_maxlum'])
        rpr_minlum = rpr_minlum.astype(float)
        rpr_maxlum = rpr_maxlum.astype(float)
        
        #Basic weirdness checking for luminosity limits
        error_check = 0
        elist = np.where(rpr_minlum > rpr_maxlum)[0]
        if len(elist) > 0:
            print >> sys.stderr, "ERROR:  Require rpr_minlum < rpr_maxlum"
            error_check = 1
        if int(params['obs_clf']) == 1 or int(params['use_lum']) == 1 and min(rpr_minlum) < 0:
            print >> sys.stderr, "ERROR:  Radial profiles are using app mags or solar lum,"
            print >> sys.stderr, "        but rpr_minlum/rpr_maxlum < 0"
            error_check = 1
        if int(params['obs_clf']) == 0 and int(params['use_lum']) == 0 and max(rpr_maxlum) > 0:
            print >> sys.stderr, "ERROR:  Radial profiles are using abs mags,"
            print >> sys.stderr, "        but rpr_minlum/rpr_maxlum > 0"
            error_check = 1

        if error_check == 0:
            redm_rpr.redm_rpr( cat,mem,mag,lm_min,lm_max,zmin,zmax,rpr_minlum,rpr_maxlum,
                               bootlist,gboot,params['outdir'] )
        else:
            print >> sys.stderr, "SKIPPING RADIAL PROFILES"
            

    #Calculate magnitude gaps
    if int(params['do_mgap']) == 1:
        print >> sys.stderr, "Calculating magnitude gaps..."

        redm_mgap.redm_mgap(cat,mem,cenmag,cengalindex,mag,zmin,zmax,
                            lm_min,lm_max,
                            bootlist,gboot,params['outdir'],
                            use_lum=bool(int(params['use_lum'])),
                            use_obs=bool(int(params['obs_clf'])),
                            weight_cen=bool(int(params['weight_cen'])) )

    #calculate the probability that the brightest galaxy is not the central galaxy
    if int(params['do_pbcg'])==1:
        print >> sys.stderr, "Calculating P(BCG!=central)..."
        
        redm_pbcg.get_p_bcg_not_cen(cat,cenmag,cengalindex,cat['p_cen'],mem['mem_match_id'],mag,mem['p'],
                                    zmin,zmax,params['outdir'],use_lum=int(params['use_lum']),weight_cen=int(params['weight_cen']))

    #Calculate the distribution of the brightest satellite galaxy
    if int(params['do_bsat'])==1:
        print >> sys.stderr, "Calculating brightest satellite clf..."
        redm_bright_sat.get_brightest_satellite_all(cat,mem,mag,cengalindex,lm_min,lm_max,zmin,zmax,
                                                    bootlist,gboot,match_index,
                                                    params['outdir'],weight_cen=int(params['weight_cen']),
                                                    use_lum=int(params['use_lum']),obs_clf=int(params['obs_clf']))
        
        print >> sys.stderr, "Calculating joint brightest sat-central distribution..."
        count_arr = redm_bright_sat.get_bright_sat_cen_all(cat,mem,mag,cengalindex,
                                                           cenmag,lm_min,lm_max,
                                                           zmin,zmax,
                                                           bootlist,gboot,match_index,
                                                           params['outdir'],
                                                           weight_cen=int(params['weight_cen']),
                                                           use_lum=int(params['use_lum']),
                                                           obs_clf=int(params['obs_clf']))
    #output the params
    np.save(params['outdir']+"params.npy", params)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print >> sys.stderr, "ERROR:  Required input format is"
        print >> sys.stderr, "        paramfile"
        sys.exit(1)

    filename = sys.argv[1]
    if len(sys.argv) == 4 :
        zmin = eval(sys.argv[2])
        zmax = eval(sys.argv[3])
    else:
        zmin=None
        zmax=None
    main(filename, zmin=zmin, zmax=zmax)

