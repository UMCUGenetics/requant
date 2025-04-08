import sys
import numpy as np
# import numpy.ma as ma
import random
import re
import pysam
import os

'''
Similar to mark_fasta.py, but as Remora doesn't want to work with
marking motifs as modified or not, we had to just cut out regions
from the reference genome with the intended k-mer coverage.
'''

path_fasta = '../tamatoa/data/lambda_phage_unix.fa'
fasta = pysam.FastaFile(path_fasta)


motif = 'CG'
mer_size = 6
padding = mer_size - len(motif) # Padding differs a bit between fwd and rev sides later
max_picks = int(np.power(4,mer_size - len(motif))*5*0.90)
contexts = set()
kmers = set()
last_end = 0

regions = []

revcomp_lookup = {'A':'T','C':'G','G':'C','T':'A'}
def rev_comp(x):
    return ''.join([revcomp_lookup[bp] for bp in x][::-1])

for refseq in fasta.references:
    faseq = fasta.fetch(refseq)
    mod_iter = (re.finditer('CG',faseq))

    for mpos in mod_iter:
        if mpos.start() < mer_size:
            continue
        print('---',mpos.start(),mpos.end())

        if not kmers:
            cur_start = mpos.start()

        env_wide = faseq[mpos.start()-padding-1:mpos.end()+padding+1]
        env_seq = env_wide.replace('CG','MG')[:-1]
        contexts.add(env_seq)
        add_mers = [env_seq[i:i+mer_size] for i in range(mer_size)]
        kmers.update(add_mers)

        # Add reverse contexts as well
        env_seq_rev = rev_comp(env_wide).replace('CG','MG')[:-1]
        contexts.add(env_seq_rev)
        add_mers_rev = [env_seq_rev[i:i+mer_size] for i in range(mer_size)]
        kmers.update(add_mers_rev)

        # Break if we have more than specified
        # Technically we may break after a few more kmers more than intended
        if len(kmers) >= max_picks:
            # TODO? Adding the padding may mean we hit another CG motif close by,
            # leaving it out means we cannot always complete the context if we really make a cut
            regions.append([cur_start,mpos.end(),len(kmers),contexts,kmers])
            contexts = set()
            kmers = set()
            print('-')

# Make output dir to store regions to
dir_out = './data/refcuts_'+motif+'_'+str(mer_size)+'_'+str(max_picks)#+'_m'
os.mkdir(dir_out)

# Write region to separate fasta files
for reg in regions:
    seq_name = '>'+refseq+':'+str(reg[0])+'-'+str(reg[1])+' '+str(reg[2])
    print(seq_name)

    with open(dir_out+'/'+refseq+'_'+str(reg[0])+'-'+str(reg[1])+'.fa', 'w') as f:
        f.write(seq_name+'\n')

        seq = fasta.fetch(refseq,reg[0],reg[1])

        # seq_mod = seq.replace('CG','MG')
        # seq_all = fasta.fetch(refseq)
        # seq_comb = seq_all[:reg[0]]+seq_mod+seq_all[reg[1]:]
        # seq = seq_comb

        n = 80
        chunks = [seq[i:i+n] for i in range(0, len(seq), n)]
        f.writelines('\n'.join(chunks))

    with open(dir_out+'/'+refseq+'_'+str(reg[0])+'-'+str(reg[1])+'.log', 'w') as f:
        f.write('>contexts\n')
        f.writelines('\n'.join(reg[3]))
        f.write('\n>kmers\n')
        f.writelines('\n'.join(reg[4]))
        f.write('\n')