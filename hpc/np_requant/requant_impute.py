import sys
import requant as rq
import matplotlib.pyplot as plt
import numpy as np

f_np_model = sys.argv[1]
f_kmer_log = sys.argv[2]
f_imputed = sys.argv[3]

dict_canon, dict_modif, dict_main_sd = rq.load_polishmodel(f_np_model)

# for kmer in dict_modif:
#     print(kmer,dict_modif[kmer]==dict_canon[kmer.replace('MG','CG')])

motif, motif_mod = 'CG', 'MG'
# motif, motif_mod = 'GC', 'GM'

print(np.sum([dict_modif[k]!=dict_canon[k.replace(motif_mod,motif)] for k in dict_modif]))

if f_kmer_log != '-':
    kmers_kept = set()
    with open(f_kmer_log) as kmer_log:
        # Skip the context sequences
        for line in kmer_log:
            if line.startswith('>kmers'):
                break
        # Read the actual kmers
        for line in kmer_log:
            kmers_kept.add(line.strip())

    # We stored CG kmer motifs so need to replace with modified motifs (MG) now
    kmers_kept = list(kmers_kept)
    # print(len(kmers_kept))
    kmers_kept = [k.replace(motif, motif_mod) for k in kmers_kept]
    values = [dict_modif[k.replace(motif, motif_mod)] for k in kmers_kept]
    sub_dict_modif = dict(zip(kmers_kept, values))
else:
    sub_dict_modif = dict_modif

# Not all kmers are really trained, identify those by current values
untrained_kmers = set()
for key,value in sub_dict_modif.items():
    if value == dict_canon[key.replace(motif_mod, motif)]:
        untrained_kmers.add(key)
        #print(key,value,dict_canon[key.replace(motif_mod, motif)])

# Remove kmers with no difference between canon and mod
for kmer in untrained_kmers:
    sub_dict_modif.pop(kmer)

# Actual number of kmers we have left now:
print('Trained kmers available:',len(sub_dict_modif))


# Determine the effect of on each base by position in the kmer
eff_dict  = rq.my_shifter(dict_canon,sub_dict_modif)
rep_table = rq.generate_table(dict_canon,eff_dict,dict_main_sd)

# Write replaced model to a table
rq.write_model(rep_table,f_imputed+'.replaced.model')

# Impute the full model but retain trained kmers
add_table = rq.generate_table(dict_canon,eff_dict,dict_main_sd,sub_dict_modif)

# Write imputed model to a table
rq.write_model(add_table,f_imputed+'.added.model')


sorted_eff_dict = sorted(eff_dict.items(), key=lambda x: abs(x[1]), reverse=True)

for x in sorted_eff_dict: print(x[0]+'\t'+str(x[1]))