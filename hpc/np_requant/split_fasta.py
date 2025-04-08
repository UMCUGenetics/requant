import sys
import numpy as np
# import numpy.ma as ma
import random
import re
import pysam
import os

path_fasta = '../tamatoa/data/lambda_phage_unix.fa'

fasta = pysam.FastaFile(path_fasta)

motif = 'CG'
mer_size = 6
padding = mer_size - len(motif)
max_picks = int(np.power(4,padding)*5*0.25)
contexts = set()
kmers = set()
last_end = 0

regions = []

for refseq in fasta.references:
    faseq = fasta.fetch(refseq)
    mod_iter = (re.finditer('CG',faseq))

    for mpos in mod_iter:
        if mpos.start() < mer_size:
            continue
        print('---',mpos.start(),mpos.end())

        if not kmers:
            cur_start = mpos.start()

        env_seq = faseq[mpos.start()-padding:mpos.end()+padding]
        contexts.add(env_seq)
        add_mers = [env_seq[i:i+mer_size] for i in range(mer_size-1)]
        kmers.update(add_mers)

        # Break if we have more than specified
        # Technically we may break after a few more kmers more than intended
        if len(kmers) >= max_picks:
            # TODO? Adding the padding may mean we hit another CG motif close by
            regions.append([cur_start-padding,mpos.end()+padding,len(kmers),contexts,kmers])
            contexts = set()
            kmers = set()
            print('-')

# Make output dir to store regions to
dir_out = './data/refcuts_'+motif+'_'+str(mer_size)+'_'+str(max_picks)
os.mkdir(dir_out)

# Write region to separate fasta files
for reg in regions:
    seq_name = '>'+refseq+':'+str(reg[0])+'-'+str(reg[1])+' '+str(reg[2])
    print(seq_name)

    with open(dir_out+'/'+refseq+'_'+str(reg[0])+'-'+str(reg[1])+'.fa', 'w') as f:
        f.write(seq_name+'\n')

        seq = fasta.fetch(refseq,reg[0],reg[1])
        n = 80
        chunks = [seq[i:i+n] for i in range(0, len(seq), n)]
        f.writelines('\n'.join(chunks))

    with open(dir_out+'/'+refseq+'_'+str(reg[0])+'-'+str(reg[1])+'.log', 'w') as f:
        f.write('>contexts\n')
        f.writelines('\n'.join(reg[3]))
        f.write('\n>kmers\n')
        f.writelines('\n'.join(reg[4]))
        f.write('\n')