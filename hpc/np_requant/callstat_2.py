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
    return contexts,set(kmers)


def np_call_counts(path_methcalls,motif_sites,kmers_trained={},kmer_size=6):
    llrs = None
    # print(path_methcalls)
    with bz2.open(path_methcalls,mode='rb') as methcalls:
        header = methcalls.readline().decode().split()
        # print(header)
        idx_log_lik_ratio = header.index('log_lik_ratio')
        idx_start = header.index('start')
        idx_sequence = header.index('sequence')

        # Slow-ish implementation to check for overlaps 
        # Perhaps a bisect on the sorted side is faster?
        # Alternatively sorting the calls first might be faster with premade implementations
        keep_calls = [] # not involved in train
        drop_calls = [] # involved in train
        putr_calls = [] # pure train
        for line in methcalls:
            call_start = int(line.split()[idx_start])
            call_end   = int(line.split()[idx_start+1])
            call_llr   = float(line.split()[idx_log_lik_ratio])
            call_seq   = line.split()[idx_sequence].decode()

            # Last motif site is before this call, must be safe to add
            if motif_sites == [] or call_start > motif_sites[-1][1]:
                kmer_overlap = len(kmers_trained & set([call_seq[idx:idx + kmer_size] for idx in range(len(call_seq) - kmer_size + 1)]))
                keep_calls.append((call_llr,kmer_overlap/len(call_seq)))
            else:
                for site in motif_sites:
                    # Next site is beyond this call, must be safe to add
                    if site[0] > call_end:
                        kmer_overlap = len(kmers_trained & set([call_seq[idx:idx + kmer_size] for idx in range(len(call_seq) - kmer_size + 1)]))
                        keep_calls.append((call_llr,kmer_overlap/len(call_seq)))
                        break
                    # Overlap between site and call, involved in training
                    if site[1] > call_start and site[0] < call_end:
                        drop_calls.append(call_llr)
                        break

            # if len(keep_calls)>100000:
            #     # print(keep_calls)
            #     break

    return np.array(keep_calls), np.array(drop_calls)

def stat_call_counts(desc,llrs,llr_threshold=2):
    pos = np.sum(llrs>=llr_threshold)
    neg = np.sum(llrs<=-llr_threshold)

    print(test_case,len(motif_sites),len(kmers_trained),modstate,impute,desc,llr_threshold,
    "{:.2f}".format(pos/(neg+pos)*100)+'%',pos,neg,llrs.shape[0])


path_calls = sys.argv[1]
log_path = sys.argv[2]

file_calls = path_calls.split('/')[-1].split('.')[:3]

test_case,impute,modstate = file_calls
motif_sites, kmers_trained = get_log_contexts(log_path+'/'+test_case+'.tsv')

keep_calls, train_calls = np_call_counts(path_calls,motif_sites,kmers_trained)
for t in np.arange(0,3.5,0.5):
    stat_call_counts('train',train_calls,t) # Only sites used for training

    test_calls = keep_calls[:,0]
    stat_call_counts('test',test_calls,t) # Only sites not used for training

    all_calls = np.concatenate((test_calls,np.array(train_calls)))
    stat_call_counts('all',all_calls,t) # All sites

    stat_call_counts('zero',keep_calls[keep_calls[:,1]==0],t) # Test sites without any k-mers that were used for training
    stat_call_counts('max',keep_calls[keep_calls[:,1]==1],t) # Test sites with any k-mers that were used for training
