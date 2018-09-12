#!/bin/bash
QUEUE=long
Time=120:00
#Main runs
bsub -q $QUEUE -W $Time -W $Time  -oo dr8_finaltest_eli_new_test_0.1.log "python ../redm_full.py param_DES_y1a1_real_test.dat [0.1] [0.2]"
bsub -q $QUEUE -W $Time -W $Time  -oo dr8_finaltest_eli_new_test_0.2.log "python ../redm_full.py param_DES_y1a1_real_test.dat [0.2] [0.3]"
bsub -q $QUEUE -W $Time -W $Time  -oo dr8_finaltest_eli_new_test_0.3.log "python ../redm_full.py param_DES_y1a1_real_test.dat [0.3] [0.4]"
bsub -q $QUEUE -W $Time -W $Time  -oo dr8_finaltest_eli_new_test_0.4.log "python ../redm_full.py param_DES_y1a1_real_test.dat [0.4] [0.5]"
bsub -q $QUEUE -W $Time -W $Time  -oo dr8_finaltest_eli_new_test_0.5.log "python ../redm_full.py param_DES_y1a1_real_test.dat [0.5] [0.6]"
bsub -q $QUEUE -W $Time -W $Time  -oo dr8_finaltest_eli_new_test_0.6.log "python ../redm_full.py param_DES_y1a1_real_test.dat [0.6] [0.7]"
bsub -q $QUEUE -W $Time -W $Time  -oo dr8_finaltest_eli_new_test_0.7.log "python ../redm_full.py param_DES_y1a1_real_test.dat [0.7] [0.8]"
bsub -q $QUEUE -W $Time -W $Time  -oo dr8_finaltest_eli_new_test_0.8.log "python ../redm_full.py param_DES_y1a1_real_test.dat [0.8] [0.9]"

