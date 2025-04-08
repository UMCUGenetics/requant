import sys
import numpy as np
# import numpy.ma as ma
import random
import re
import pysam
import os

'''
This script marks the reference genome with motifs for testing.
Sites are randomly sampled from the reference genome and stored,
and a copy of the reference genome is made with the motifs replaced
by modified motifs.
'''

path_fasta = '../tamatoa/data/lambda_phage_unix.fa'
fasta = pysam.FastaFile(path_fasta)


motif, motif_mod = 'CG', 'MG'
motif, motif_mod = 'GC', 'GM'
mer_size = 6
padding = mer_size - len(motif) # Padding differs a bit between fwd and rev sides later
max_picks = int(np.power(4,mer_size - len(motif))*5*0.10)
contexts = set()
last_end = 0

regions = []

revcomp_lookup = {'A':'T','C':'G','G':'C','T':'A'}
def rev_comp(x):
    return ''.join([revcomp_lookup[bp] for bp in x][::-1])


def get_motif_matches():
    all_contexts = []
    for refseq in fasta.references:
        faseq = fasta.fetch(refseq)
        mod_iter = (re.finditer(motif,faseq))

        for mpos in mod_iter:
            # Start of refseq match
            if mpos.start() < mer_size:
                continue

            env_wide = faseq[mpos.start()-padding:mpos.end()+padding]

            # End of refseq match
            if len(env_wide) != 10:
                continue

            # Ignore multiple motifs too close together
            if env_wide.count(motif) > 1:
                print(env_wide)
                continue

            all_contexts.append([refseq,mpos.start(),mpos.end(),env_wide])
    return all_contexts


def get_random_contexts(all_contexts, kmer_count):
    random.shuffle(all_contexts)
    kmers = set()
    context_sites = []
    contexts = set()
    for x in all_contexts:
        if len(kmers) >= kmer_count:
            break
        # print(x)

        # Only add context if not already completely added before
        add_mers = [x[3][i:i+mer_size] for i in range(mer_size-1)]
        if set(add_mers) - kmers:
            kmers.update(add_mers)
            add_mers_rev = [rev_comp(x[3])[i:i+mer_size] for i in range(mer_size-1)]
            kmers.update(add_mers_rev)
            context_sites.append(x)
            contexts.add(x[3])
            contexts.add(rev_comp(x[3]))
    return kmers,context_sites


def get_minimal_contexts(all_contexts, coverage):
    # Effectively this function picks a random motif context and looks for a
    # another motif context with 0 matching bases other than the motif.
    # This is repeated until enough such cases were matched to get the intended
    # coverage for each base at each position in the k-mer.

    def does_overlap(A,B):
        # Check if any character between these two strings overlaps (same char @ same index).
        for i,c in enumerate(A):
            if B[i] == c:
                return True
        return False

    def get_mers(y):
        # Get all k-mers from context sequence, including reverse complement options.
        return set([y[3][i:i+mer_size] for i in range(mer_size-1)]+[rev_comp(y[3])[i:i+mer_size] for i in range(mer_size-1)])

    kmers = set()
    context_sites = []
    contexts = set()

    picks = []

    while len(picks) < 2 * coverage:
        random.shuffle(all_contexts)
        # Pick a context
        first_pick = all_contexts[0]
        # Already picked this context before? Skip this first pick and try again
        if first_pick[3] in contexts or len(get_mers(first_pick))<10 or len(get_mers(first_pick)-kmers)<10:
            continue
        # Test all other contexts for any non-motif overlap
        for x in all_contexts:
            if not does_overlap(first_pick[3][:4]+first_pick[3][6:],x[3][:4]+x[3][6:]) and \
               not does_overlap(first_pick[3][:4]+first_pick[3][6:],rev_comp(x[3][:4]+x[3][6:])):
                # Already picked this context before? Skip this match and look for something better
                if x[3] in contexts or len(get_mers(x))<10 or len(get_mers(x)-kmers)<10:
                    continue

                for y in [first_pick,x]:
                    kmers.update(get_mers(y))
                    context_sites.append(y)
                    contexts.add(y[3])
                    contexts.add(rev_comp(y[3]))
                    picks.append(y[3])

                break

    return kmers,context_sites


def write_fasta(path_fasta, context_sites):
    context_sites.sort()
    refseq = context_sites[0][0] # one ref chromosome in lambda anyway

    seq_all = fasta.fetch(refseq)
    last = 0
    seq_comb = ''
    for context in context_sites:
        # print(last,context[1],seq_all[context[1]:context[2]])
        seq_comb = seq_comb+seq_all[last:context[1]]+motif_mod
        last = context[2]
    seq_comb = seq_comb+seq_all[last:]
    # print(last,len(seq_all))
    # print(len(seq_all),len(seq_comb))
    # print(seq_comb)

    with open(path_fasta, 'w') as f:
        f.write('>'+refseq+'\n')

        n = 80
        chunks = [seq_comb[i:i+n] for i in range(0, len(seq_comb), n)]
        f.writelines('\n'.join(chunks))


def write_contexts(path_contexts,kmers,context_sites):
    with open(path_contexts, 'w') as f:
        f.write('>contexts\n')
        for context in context_sites:
            f.write('\t'.join([str(x) for x in context])+'\n')

        f.write('>kmers\n')
        for kmer in kmers:
            f.write(kmer+'\n')


all_contexts = get_motif_matches()
print('Detected motifs in reference:',len(all_contexts))


# for cov in range(1,17):
# # for cov in [3]:
#     for i in range(5):
#         # Seems there is some low chance edge case we didn't catch yet, just check and try again if that happens
#         kmers=[]
#         while len(kmers)!=20*cov:
#             kmers,context_sites = get_minimal_contexts(all_contexts,cov)
#         print(cov,i,len(kmers),len(context_sites))
#         write_fasta('data/refcuts_'+motif+'_minimal_2/p'+str(cov)+'_'+str(i)+'.fa', context_sites)
#         write_contexts('data/refcuts_'+motif+'_minimal_2/p'+str(cov)+'_'+str(i)+'.tsv',kmers,context_sites)
# exit()

for perc in [90]:#[5,10,25,50,75,90]:
    max_picks = int(np.power(4,mer_size - len(motif))*5*perc/100)
    for i in range(5):
        kmers,context_sites = get_random_contexts(all_contexts,max_picks)
        print(max_picks,len(kmers),len(context_sites),len(kmers)/int(np.power(4,mer_size - len(motif))*5))#,kmers,context_sites)
        # print(context_sites)
        write_fasta('data/refcuts_'+motif+'_test_3/p'+str(perc)+'_'+str(i)+'.fa', context_sites)
        write_contexts('data/refcuts_'+motif+'_test_3/p'+str(perc)+'_'+str(i)+'.tsv',kmers,context_sites)
exit()

# Make output dir to store regions to
dir_out = './data/refcuts_'+motif+'_'+str(mer_size)+'_'+str(max_picks)+'_m'
os.mkdir(dir_out)

# Write region to separate fasta files
for reg in regions:
    seq_name = '>'+refseq+':'+str(reg[0])+'-'+str(reg[1])+' '+str(reg[2])
    print(seq_name)

    with open(dir_out+'/'+refseq+'_'+str(reg[0])+'-'+str(reg[1])+'.fa', 'w') as f:
        f.write(seq_name+'\n')

        seq = fasta.fetch(refseq,reg[0],reg[1])

        seq_mod = seq.replace('CG','MG')
        seq_all = fasta.fetch(refseq)
        seq_comb = seq_all[:reg[0]]+seq_mod+seq_all[reg[1]:]
        seq = seq_comb

        n = 80
        chunks = [seq[i:i+n] for i in range(0, len(seq), n)]
        f.writelines('\n'.join(chunks))

    with open(dir_out+'/'+refseq+'_'+str(reg[0])+'-'+str(reg[1])+'.log', 'w') as f:
        f.write('>contexts\n')
        f.writelines('\n'.join(reg[3]))
        f.write('\n>kmers\n')
        f.writelines('\n'.join(reg[4]))
        f.write('\n')