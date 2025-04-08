import sys
import bz2
import numpy as np

def remora_call_freqs_bz2(path_methcalls):
    llrs = None
    print(path_methcalls)
    region_start, region_end = [int(x) for x in (path_methcalls.split('/')[-4].split('_')[1].split('-'))]
    with bz2.open(path_methcalls,mode='rb') as methcalls:
        header = methcalls.readline().decode().split()
        # print(header)
        idx_log_lik_mod = header.index('mod_log_prob')
        idx_log_lik_can = header.index('can_log_prob')

        idx_start = header.index('pos')

        llrs = np.array([float(line.split()[idx_log_lik_mod])-float(line.split()[idx_log_lik_can]) for line in methcalls if int(line.split()[idx_start])>region_end or int(line.split()[idx_start])+1<region_start])

    print('Cutoff\tNegatives\tUnsure\tPositives\tPos%Kept\tPos%All\tNeg%All\tUns%All')
    # Determine frequency of positive calls for various cutoff values
    for i in np.arange(0,4,1):
        neg = np.sum(llrs<=-i)
        pos = np.sum(llrs>=i)
        print(i,neg,llrs.shape[0]-neg-pos,pos,"{:.2f}".format(pos/(neg+pos)*100)+'%',"{:.2f}".format(pos/llrs.shape[0]*100)+'%',"{:.2f}".format(neg/llrs.shape[0]*100)+'%',"{:.2f}".format((llrs.shape[0]-neg-pos)/llrs.shape[0]*100)+'%')

    # return llrs

remora_call_freqs_bz2(sys.argv[1])
