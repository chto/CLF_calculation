#!/bin/bash
QUEUE=kipac-ibq
Time=20:00
#Main runs
#bsub -x -q $QUEUE -W $Time -W $Time  -oo chtolog/param_dr8_v5.10_paper_Ber_v1.log "python redm_full.py /u/ki/chto100/code/redmapper_clf/redmapper/chtoParam/SDSSPaper/Bernardi/chto_correct_v1.dat"
#bsub -x -q $QUEUE -W $Time -W $Time  -oo chtolog/param_dr8_v5.10_paper_Ber_v2.log "python redm_full.py /u/ki/chto100/code/redmapper_clf/redmapper/chtoParam/SDSSPaper/Bernardi/chto_correct_v2.dat"
#bsub -x -q $QUEUE -W $Time -W $Time  -oo chtolog/param_dr8_v5.10_paper_Ber_v3.log "python redm_full.py /u/ki/chto100/code/redmapper_clf/redmapper/chtoParam/SDSSPaper/Bernardi/chto_correct_v3.dat"
#bsub -x -q $QUEUE -W $Time -W $Time  -oo chtolog/param_dr8_v5.10_paper_Ber_v4.log "python redm_full.py /u/ki/chto100/code/redmapper_clf/redmapper/chtoParam/SDSSPaper/Bernardi/chto_correct_v4.dat"
#bsub -x -q $QUEUE -W $Time -W $Time  -oo chtolog/param_dr8_v5.10_paper_Ber_v5.log "python redm_full.py /u/ki/chto100/code/redmapper_clf/redmapper/chtoParam/SDSSPaper/Bernardi/chto_correct_v5.dat"
#bsub -x -q $QUEUE -W $Time -W $Time  -oo chtolog/param_dr8_v5.10_paper_Ber_v6.log "python redm_full.py /u/ki/chto100/code/redmapper_clf/redmapper/chtoParam/SDSSPaper/Bernardi/chto_correct_v6.dat"
bsub -x -q $QUEUE -W $Time -W $Time  -oo chtolog/param_dr8_v5.10_paper_Ber_v7.log "python redm_full.py /u/ki/chto100/code/redmapper_clf/redmapper/chtoParam/SDSSPaper/Bernardi/chto_correct_v7_nocencorr.dat"
