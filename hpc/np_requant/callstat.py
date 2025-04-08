import sys
import numpy as np
import bz2

def get_log_contexts(path_kmer_log):
    contexts = []
    kmers = []
    with open(path_kmer_log) as kmer_log:
        # Skip anything before
        for line in kmer_log:
            if line.startswith('>contexts'):
                break
        # Read the actual contexts
        for line in kmer_log:
            # Quit if we find something else
            if line.startswith('>'):
                break
            contexts.append([int(x) for x in line.split()[1:3]])
        for line in kmer_log:
            # Quit if we find something else
            if line.startswith('>'):
                break
            kmers.append(line.strip())
    return contexts,kmers


def np_call_counts(path_methcalls,motif_sites):
    llrs = None
    # print(path_methcalls)
    with bz2.open(path_methcalls,mode='rb') as methcalls:
        header = methcalls.readline().decode().split()
        # print(header)
        idx_log_lik_ratio = header.index('log_lik_ratio')
        idx_start = header.index('start')

        # Slow-ish implementation to check for overlaps 
        # Perhaps a bisect on the sorted side is faster?
        # Alternatively sorting the calls first might be faster with premade implementations
        keep_calls = []
        drop_calls = []
        for line in methcalls:
            # print(line.split())
            call_start = int(line.split()[idx_start])
            call_end   = int(line.split()[idx_start+1])
            call_llr   = float(line.split()[idx_log_lik_ratio])

            for site in motif_sites:
                # Next site is beyond this call, must be safe to add
                if site[0] > call_end:
                    # print(site,call_start,call_end)
                    keep_calls.append(call_llr)
                    break
                # Overlap between site and call, involved in training
                if site[1] > call_start and site[0] < call_end:
                    # print(site,call_start,call_end,'drop')
                    drop_calls.append(call_llr)
                    break

    llrs = np.array(keep_calls)
    llrs_all = np.array(keep_calls+drop_calls)

    threshold = 2
    pos = np.sum(llrs>=threshold)
    neg = np.sum(llrs<=-threshold)

    pos_all = np.sum(llrs_all>=threshold)
    neg_all = np.sum(llrs_all<=-threshold)
    return [pos,neg,llrs.shape[0],pos_all,neg_all,llrs_all.shape[0]]


path_calls = sys.argv[1]
log_path = sys.argv[2]

file_calls = path_calls.split('/')[-1].split('.')[:3]

test_case,impute,modstate = file_calls
motif_sites, kmers_trained = get_log_contexts(log_path+'/'+test_case+'.tsv')

mod_pos, mod_neg, mod_len, mod_pos_all, mod_neg_all, mod_len_all = np_call_counts(path_calls,motif_sites)
print(test_case,modstate,impute,len(motif_sites),len(kmers_trained),
    ' ',"{:.2f}".format(mod_pos/(mod_neg+mod_pos)*100)+'%',mod_pos,mod_neg,mod_len,
    ' ',"{:.2f}".format(mod_pos_all/(mod_neg_all+mod_pos_all)*100)+'%',mod_pos_all,mod_neg_all,mod_len_all)