import sys
import numpy as np
import random

import itertools

import matplotlib.pyplot as plt
# import matplotlib.patches as patches
# from matplotlib import cm

mer_size = 6
mod_motif, mod_motif_M, mod_alphabet = 'CG', 'MG', 'cpg'
# mod_motif, mod_motif_M, mod_alphabet = 'GC', 'GM', 'gpc'
base_colors = {
    'A':'C2',
    'C':'C0',
    'G':'C1',
    'T':'C3',
    'M':'C7'}

# This should be replaced with a proper method that isn't CG specific
# Disabled/deprecated, was to catch half motif at end of kmer e.g. xxxxxM(G)
def replace_edge(kmer):
    # if kmer[-1] == 'M':
    #     return kmer[:-1]+'C'
    return kmer


def load_polishmodel(model_path):
    dict_main = {}
    dict_main_sd = {}
    dict_alt = {}
    with open(model_path) as fp:
        for line in fp:
            if line[0] == '#':
                continue
            splitline = line.split()
            if splitline[0] == 'kmer':
                continue
            kmer = splitline[0]

            # Does this kmer have a modified base
            if 'M' in kmer:
                # That kmer wasn't trained as it should have both a modified and an unmodified base
                if kmer.count(mod_motif) >= 1:
                    continue

                # Replace all mod motif occurrences with canonical motifs
                canon = replace_edge(kmer.replace(mod_motif_M,mod_motif))
                # kmer still has an M somewhere so it wasn't really a proper modification motif
                if 'M' in canon:
                    continue

                # Looks like a properly trained modified kmer, keep it
                # dict_alt[canon] = float(splitline[1])
                dict_alt[kmer] = float(splitline[1])

            # This is a canonical kmer, other dict
            else:
                dict_main[kmer] = float(splitline[1])
            dict_main_sd[kmer] = float(splitline[2])

    return dict_main, dict_alt, dict_main_sd


def find_strongest_eff(sub_diff):
    effects = []
    for base in 'ACGT':
        for index in range(mer_size):
            subsub = [val for kmer,val in sub_diff.items() if kmer[index]==base]
            if subsub:
                # print(base,index,np.median(subsub))
                effects.append([np.abs(np.median(subsub)),index,base,np.median(subsub)])
    effects.sort()
    if effects:
        return effects[-1][1:]


def apply_effect(sub_diff,effect):
    for key,val in sub_diff.items():
        if key[effect[0]] == effect[1]:
            sub_diff[key]-=effect[2]


def process_subdict(sub_diff):
    sub_med = np.median(np.array(list(sub_diff.values())))
    for key,val in sub_diff.items():
        sub_diff[key]-=sub_med

    steps = []
    str_eff = find_strongest_eff(sub_diff)
    if str_eff:
        while abs(str_eff[2]) > 0.01:
            steps.append(str_eff)
            apply_effect(sub_diff,str_eff)
            str_eff = find_strongest_eff(sub_diff)
    return sub_med,steps


def flatten_effdict(effects):
    effects.sort()
    if not effects:
        return
    last_eff = effects[0]
    new_list = []
    for eff in effects[1:]:
        if eff[:2] == last_eff[:2]:
            last_eff[2] = last_eff[2]+eff[2]
        else:
            new_list.append(last_eff)
            last_eff = eff
    
    new_list.append(last_eff)

    return new_list
    

def per_modifloc(dict_canon, dict_modif):
    all_effects = []
    for i in range(mer_size):
        diff_dict = {}
        for m_kmer,m_val in dict_modif.items():
            if m_kmer[i] != 'M':
                continue
            if m_kmer.count('M') > 1:
                continue
            canon = replace_edge(m_kmer.replace(mod_motif_M,mod_motif))
            diff = dict_canon[canon]-m_val
            diff_dict[m_kmer]=diff
        if diff_dict:    
            effect_all, effect_per = process_subdict(diff_dict)
            effect_per = flatten_effdict(effect_per)
            all_effects.append([i,effect_all,effect_per])
    return all_effects


def effs_to_dict(all_effects):
    eff_dict = {}
    for modpos in all_effects:
        eff_dict[str(modpos[0])+'M'] = modpos[1]
        for basepos in modpos[2]:
            eff_dict[str(modpos[0])+'M'+str(basepos[0])+str(basepos[1])] = basepos[2]
    return eff_dict


def impute_mod(eff_dict,kmer):
    delta = 0
    lid = ''
    for mi,mbase in enumerate(kmer):
        if mbase != 'M':
            continue
        else:
            if str(mi)+'M' in eff_dict:
                delta += eff_dict[str(mi)+'M']
            for oi,obase in enumerate(kmer):
                lid = str(mi)+'M'+str(oi)+obase
                if lid in eff_dict:
                    delta += eff_dict[lid]
    return delta


# def sample_from_dict(d, sample=10):
#     keys = random.sample(list(d), sample)
#     values = [d[k] for k in keys]
#     return dict(zip(keys, values))


def proper_random(max_picks,bases='ACGT'):#subsample=.4):
    # Generate all possible 4bp sequences
    env_len = mer_size-len(mod_motif)
    env_mers = [''.join(p) for p in itertools.product(bases, repeat=env_len)]
    # Ignore kmers with >1 motif for now
    env_mers = [p for p in env_mers if p.count(mod_motif) == 0]
    # Generate all possible prefix and postfix mers around the motif
    env_comb = [p for p in itertools.product(env_mers, repeat=2)]
    # Randomize order
    random.shuffle(env_comb)
    # Determine total number of kmers to sample
    # max_picks = np.power(len(bases),env_len) * (mer_size-1) * subsample
    # max_picks = len(env_mers) * (mer_size-1) * subsample
    # print(len(env_mers) * (mer_size-1),max_picks)
    kmers = set()
    contexts = set()
    # Add each combination to the kmer set
    for env in env_comb:
        env_seq = env[0]+mod_motif_M+env[1]
        contexts.add(env_seq)
        add_mers = [env_seq[i:i+mer_size] for i in range(mer_size-1)]
        kmers.update(add_mers)
        # Break if we have more than specified
        if len(kmers) >= max_picks:
            # Technically we may break after a few more kmers more than intended
            # but proper fixing would require checking and skipping sequences
            # until we find a perfect fit
            break
    return kmers, contexts


def impute_table(dict_canon,dict_modif,eff_dict,keepers={}):
    imputes = []
    for mod_kmer in dict_modif:
        canon = replace_edge(mod_kmer.replace(mod_motif_M,mod_motif))
        if mod_kmer in keepers:
            imputes.append([mod_kmer.index('M'),dict_canon[canon],dict_modif[mod_kmer],keepers[mod_kmer],keepers[mod_kmer]-dict_canon[canon]])
        else:
            imputes.append([mod_kmer.index('M'),dict_canon[canon],dict_modif[mod_kmer],dict_canon[canon]-impute_mod(eff_dict,mod_kmer),impute_mod(eff_dict,mod_kmer)])
    return np.array(imputes)


def plot_imputes(imputes,prefix='',saveas=None):
    def plot_modpos(x,y):
        for i in range(5):#np.sort(np.unique(imputes[:,0]).astype(int))[::-1]:
            i = 4-i
            pos_imputes = imputes[imputes[:,0]==i]
            plt.scatter(pos_imputes[:,x],pos_imputes[:,y],c=4-pos_imputes[:,0],alpha=1,s=1, vmin=0, vmax=9,label='.'*i+'CG'+'.'*(4-i))
        plt.legend()
    # Force same range for all plots
        plt.xlim(55,125)
        plt.ylim(55,125)
    
    plt.figure(figsize=(4,4))
    plt.set_cmap('tab10')
    plt.plot([np.min(imputes[:,1:3]),np.max(imputes[:,1:3])],[np.min(imputes[:,1:3]),np.max(imputes[:,1:3])],alpha=1)
    plot_modpos(1,2)
    # plt.title(prefix+'Input model')
    plt.title(prefix)
    plt.xlabel('Canonical value')
    plt.ylabel('Modified value')
    plt.tight_layout()
    if saveas:
        plt.savefig(saveas+'_cm.pdf')
    plt.figure(figsize=(4,4))
    plt.plot([np.min(imputes[:,1:4:2]),np.max(imputes[:,1:4:2])],[np.min(imputes[:,1:4:2]),np.max(imputes[:,1:4:2])])
    plot_modpos(1,3)
    # plt.title(prefix+'Imputed model')
    plt.title(prefix)
    plt.xlabel('Canonical value')
    plt.ylabel('Imputed modified value')
    plt.tight_layout()
    if saveas:
        plt.savefig(saveas+'_ci.pdf')
    plt.figure(figsize=(4,4))
    plt.plot([np.min(imputes[:,2:4]),np.max(imputes[:,2:4])],[np.min(imputes[:,2:4]),np.max(imputes[:,2:4])])
    plot_modpos(2,3)
    # plt.title(prefix+'Input mod vs imputed mod')
    plt.title(prefix)
    plt.xlabel('True modified value')
    plt.ylabel('Imputed modified value')
    plt.tight_layout()
    if saveas:
        plt.savefig(saveas+'_mi.pdf')


def encode_table(dict_modif,bases='ACGTM'):
    ohar = np.zeros([len(dict_modif),mer_size*len(bases)+4])
    for i,(kmer,val) in enumerate(dict_modif.items()):
        
        for j,base in enumerate(kmer):
            ohar[i,j*len(bases)+bases.index(base)] = 1
        ohar[i,-4] = kmer.find('M')
        ohar[i,-3] = dict_canon[replace_edge(kmer.replace(mod_motif_M,mod_motif))]
        ohar[i,-2] = val
        ohar[i,-1] = val - dict_canon[replace_edge(kmer.replace(mod_motif_M,mod_motif))]
        # print(i,kmer,val,ohar[i])

    # print(ohar)
    return ohar


def my_shifter(dict_canon,train):
    # all_effects = per_modifloc(dict_canon, test)
    # eff_dict = effs_to_dict(all_effects)

    sub_effects = per_modifloc(dict_canon, train)
    eff_dict = effs_to_dict(sub_effects)

    # plot_imputes(train,eff_dict)

    # plot_imputes(dict_canon,test,eff_dict)

    # plt.show()
    return(eff_dict)


def generate_table(dict_canon,eff_dict,dict_canon_sd,keepers={}):
    model_table = []
    bases = ['A','C','G','T','M']
    bases.sort()
    for p in itertools.product(bases, repeat=mer_size):
        kmer = ''.join(p)
        canon = replace_edge(kmer.replace(mod_motif_M,mod_motif))
        exp_val, exp_sd, info = 0,0,'E'

        if kmer in keepers:
            # Don't impute
            exp_val, exp_sd, info = keepers[kmer], dict_canon_sd[kmer], 'Kept'
        elif 'M' in canon:
            # kmer still has an M somewhere so it wasn't really a CG position
            canon = canon.replace("M","C")
            exp_val, exp_sd, info = dict_canon[canon],dict_canon_sd[canon],'LoneM'
        elif kmer.count(mod_motif_M) == 0:
            exp_val, exp_sd, info = dict_canon[canon],dict_canon_sd[canon],'Canon'
        else:
            new_val = dict_canon[canon]-impute_mod(eff_dict,kmer)
            exp_val, exp_sd, info = new_val,dict_canon_sd[canon],'ModFix'
        model_table.append([kmer, exp_val, exp_sd, info])
    return model_table


def write_model(new_table,out_path):
    with open(out_path,'w') as out_file:
        out_file.write('#Re-quant adjusted model based on pre-existing HMM model\n')
        # out_file.write('#training_subset\t'+str(len(kmer_train))+'\n')
        # out_file.write('#dampened_prediction\t'+str(damper)+'\n')
        out_file.write('#model_name\tr9.4_450bps.'+mod_alphabet+'.6mer.requant.model\n')
        out_file.write('#kit\tr9.4_450bps\n')
        out_file.write('#strand\ttemplate\n')
        out_file.write('#alphabet\t'+mod_alphabet+'\n')
        for kmer_line in new_table:
            out_file.write('\t'.join([str(x)[:7] for x in kmer_line])+'\n')


if __name__== "__main__":
    dict_canon, dict_modif, dict_main_sd = load_polishmodel(sys.argv[1])

    eff_dict = my_shifter(dict_canon,dict_modif)
    imputes = impute_table(dict_canon,dict_modif,eff_dict)
    plot_imputes(imputes,saveas=sys.argv[2]+'_full')
    new_table = generate_table(dict_canon,eff_dict,dict_main_sd)
    write_model(new_table,sys.argv[2])

    # # # Subsample training data semi realistically
    for subsample in [5,10,25,50,75,90,95]:
        keys,contexts = proper_random(len(dict_modif)*(0.01*subsample))
        values = [dict_modif[k] for k in keys]
        print(len(values))
        sub_dict_modif = dict(zip(keys, values))

        sub_eff_dict = my_shifter(dict_canon,sub_dict_modif)
        sub_imputes = impute_table(dict_canon,dict_modif,sub_eff_dict)
        plot_imputes(sub_imputes,prefix='Trained on: '+str(subsample)+'%',saveas=sys.argv[2]+'_sub'+str(subsample))


    # plt.show()
