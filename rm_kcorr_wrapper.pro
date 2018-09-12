;batch file for running rm_kcorr
;to be called from python by redm_full

.comp rm_kcorr

spawn,"echo $CATFILE",catfile
spawn,"echo $MEMFILE",memfile
spawn,"echo $KCORRFILE",kcorrfile
spawn,"echo $USE_DES",use_des
spawn,"echo $NO_UBAND",no_uband
spawn,"echo $LRG",LRG
spawn,"echo $BANDSHIFT",bandshift
spawn,"echo $USE_ZRED",use_zred

use_des = long(use_des[0])
no_uband = long(no_uband[0])
LRG = long(LRG[0])
bandshift = double(bandshift[0])
use_zred = long(use_zred[0])

print,catfile
print,memfile

print,use_des,no_uband,LRG,bandshift,use_zred

;IDL readin
dat = mrdfits(catfile[0],1)
dat_mem = mrdfits(memfile[0],1)

rm_kcorr,dat,dat_mem,bandshift=bandshift,$
         amag,kcorr=kcorr,des=use_des,$
         no_uband=no_uband,$
         LRG=LRG,use_zred=use_zred

mwrfits,amag,kcorrfile[0]
