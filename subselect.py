#!/u/ki/dapple/bin/python

import numpy as np
import pyfits

#Various routins for subset selection
def match_g_to_c(g,csub):
    glist = []

    count = 0L
    for i in range(len(g)):
        if i % 100000==0:
            print "Now at ",i
        while g['mem_match_id'][i] > csub['mem_match_id'][count]:
            count = count+1
            if count >= len(csub):
                break
        if count >= len(csub):
            break
        if g['mem_match_id'][i] == csub['mem_match_id'][count]:
            glist.append(i)
            
    return glist
