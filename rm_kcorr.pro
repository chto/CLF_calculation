pro rm_kcorr,dat,dat_mem,bandshift=bandshift,$
                 amag,kcorr=kcorr,des=des,no_uband=no_uband,$
                 LRG=LRG,use_zred=use_zred ;set to use LRG parameters
if n_elements(bandshift) eq 0 then bandshift = 0.5
if n_elements(des) eq 0 then des = 0
if n_elements(no_uband) eq 0 then no_uband = 0
if n_elements(LRG) eq 0 then LRG=0
if n_elements(use_zred) eq 0 then use_zred=0

;Set to use zred instead of z
if use_zred then begin
   dat.z_lambda = dat.zred
   dat.lambda_chisq = dat.lambda_zred
   index = lonarr(max(dat.mem_match_id)+1)
   index[dat.mem_match_id] = lindgen(n_elements(dat))
   dat_mem.z = dat[index[dat_mem.mem_match_id]].z_lambda
endif

;comment by chto: filter path = /afs/slac.stanford.edu/u/ki/rmredd/code/idl/kcorrect/data/filters
filterlist = ['sdss_u0.par','sdss_g0.par','sdss_r0.par','sdss_i0.par','sdss_z0.par']
if n_elements(dat_mem[0].model_mag) eq 4 then filterlist = filterlist[1:*]
nfilter = n_elements(filterlist)
if des then filterlist = ['DES_g.par','DES_r.par','DES_i.par','DES_z.par']

if nfilter eq 5 and (no_uband or des) then begin
    mag = dat_mem.model_mag[1:4]
    magerr = dat_mem.model_magerr[1:4]
    if not des then filterlist = ['sdss_g0.par','sdss_r0.par','sdss_i0.par','sdss_z0.par']
    nfilter=4
endif else begin
    mag = dat_mem.model_mag
    magerr = dat_mem.model_magerr
endelse

nclusters = n_elements(dat)
ngals = n_elements(dat_mem)

if ngals lt 1e6 then begin
   if LRG then begin
        ;print,size(mag),size(magerr),size(dat_mem.z),size(filterlist),bandshift
      print,filterlist,min(dat_mem.z),max(dat_mem.z)
      kcorrect,mag,magerr,dat_mem.z,kcorr,/magnitude,/stddev,$
               band_shift=bandshift,filterlist=filterlist,vname='lrg1'
    endif else begin
        kcorrect,mag,magerr,dat_mem.z,kcorr,/magnitude,/stddev,$
          band_shift=bandshift,filterlist=filterlist
    endelse
endif else begin
    setsize = 100000
    nsets = ceil(ngals/setsize)
    kcorr = fltarr(nfilter,ngals)
    for i=0,nsets-1 do begin
        print,"Running set ",i," of ",nsets
        list = lindgen(setsize)+i*setsize
        if i eq nsets-1 then list = lindgen(ngals-i*setsize)+i*setsize
        if LRG then begin
            kcorrect,mag[*,list],magerr[*,list],dat_mem[list].z,$
              kcorr_set,/magnitude,/stddev,$
              band_shift=bandshift,filterlist=filterlist,vname='lrg1'
        endif else begin
            kcorrect,mag[*,list],magerr[*,list],dat_mem[list].z,$
              kcorr_set,/magnitude,/stddev,$
              band_shift=bandshift,filterlist=filterlist
        endelse
        kcorr[*,list] = kcorr_set
    endfor
endelse

amag = mag

for i=0,nfilter-1 do amag[i,*] = amag[i,*] - lf_distmod(dat_mem.z) - kcorr[i,*]

end

pro run_mb_kcorr,g,amag,kcorr=kcorr,coeffs=coeffs
if n_elements(coeffs) eq 0 then coeffs = 0

filterlist = ['sdss_u0.par','sdss_g0.par','sdss_r0.par','sdss_i0.par','sdss_z0.par']

ngals = n_elements(g)

;k_load_vmatrix, vmatrix, lambda

print,"Starting kcorrection..."
if coeffs then begin
    kcorrect,omag,omagerr,g.z,kcorr,coeffs=g.coeffs,absmag=amag,$
      band_shift=0.5,filterlist=filterlist,$
      vmatrix=vmatrix,lambda=lambda,/magnitude,/stddev
endif else begin
    kcorrect,g.tmag,g.omagerr,g.z,kcorr,absmag=amag,$
      band_shift=0.5,filterlist=filterlist,$
      /magnitude,/stddev
endelse

end

;gets abs mag of the sun in the given bands
function get_solar_masses,bandshift=bandshift
if n_elements(bandshift) eq 0 then bandshift = 0.5

filterlist = ['sdss_u0.par','sdss_g0.par','sdss_r0.par','sdss_i0.par','sdss_z0.par']


return,k_solar_magnitudes(band_shift=bandshift,filterlist=filterlist)
end
