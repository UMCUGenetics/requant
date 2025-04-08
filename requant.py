### ReQuant implementation ###
### Docstrings were generated using co-pilot ###

import sys
import argparse
import itertools
import numpy as np

mer_size = 6
mod_motif, mod_motif_M, mod_alphabet = 'CG', 'MG', 'cpg'


def replace_edge(kmer):
    """
    Disabled/deprecated, meant to catch half motif at end of kmer e.g. xxxxxM(G)
    """
    return kmer


def load_polishmodel(model_path):
    """
    Load the partially trained model from the given file path.

    Args:
        model_path (str): The file path of the model.

    Returns:
        tuple: A tuple containing three dictionaries:
            - dict_main: A dictionary of canonical kmers and their values.
            - dict_alt: A dictionary of modified kmers and their values.
            - dict_main_sd: A dictionary of canonical kmers and their standard deviations.
    """
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
                canon = replace_edge(kmer.replace(mod_motif_M, mod_motif))
                # kmer still has an M somewhere so it wasn't really a proper modification motif
                if 'M' in canon:
                    continue

                # Looks like a properly trained modified kmer, keep it
                dict_alt[kmer] = float(splitline[1])

            # This is a canonical kmer, other dict
            else:
                dict_main[kmer] = float(splitline[1])
            dict_main_sd[kmer] = float(splitline[2])

    return dict_main, dict_alt, dict_main_sd


def find_strongest_eff(sub_diff):
    """
    Find the strongest effect of each base by position in the kmer.

    Args:
        sub_diff (dict): A dictionary of differences between canonical and modified kmers.

    Returns:
        list: A sorted list containing the strongest effect of each base by position in the kmer.
            Each element in the list is a list containing the following information:
            - The index of the base.
            - The base itself.
            - The median difference.
    """
    effects = []
    for base in 'ACGT':
        for index in range(mer_size):
            subsub = [val for kmer, val in sub_diff.items() if kmer[index] == base]
            if subsub:
                effects.append([np.abs(np.median(subsub)), index, base, np.median(subsub)])
    effects.sort()
    if effects:
        return effects[-1][1:]


def apply_effect(sub_diff, effect):
    """
    Subtract the effect from the differences that are still left.

    Args:
        sub_diff (dict): A dictionary of unresolved differences between canonical and modified kmers.
        effect (list): A list containing the effect to be applied.
            The effect is represented as a list containing the following information:
            - The index of the base.
            - The base itself.
            - The effect value.
    """
    for key, val in sub_diff.items():
        if key[effect[0]] == effect[1]:
            sub_diff[key] -= effect[2]


def process_subdict(sub_diff):
    """
    Process the differences between canonical and modified kmers.

    Args:
        sub_diff (dict): A dictionary of differences between canonical and modified kmers.

    Returns:
        tuple: A tuple containing two elements:
            - The median difference of the sub-dictionary.
            - A list of effects per base.
    """
    sub_med = np.median(np.array(list(sub_diff.values())))
    for key, val in sub_diff.items():
        sub_diff[key] -= sub_med

    steps = []
    str_eff = find_strongest_eff(sub_diff)
    if str_eff:
        while abs(str_eff[2]) > 0.01:
            steps.append(str_eff)
            apply_effect(sub_diff, str_eff)
            str_eff = find_strongest_eff(sub_diff)
    return sub_med, steps


def flatten_effdict(effects):
    """
    Flatten the effects dictionary.

    Args:
        effects (list): A list of effects.

    Returns:
        list: A flattened list of effects.
    """
    effects.sort()
    if not effects:
        return
    last_eff = effects[0]
    new_list = []
    for eff in effects[1:]:
        if eff[:2] == last_eff[:2]:
            last_eff[2] = last_eff[2] + eff[2]
        else:
            new_list.append(last_eff)
            last_eff = eff

    new_list.append(last_eff)

    return new_list


def per_modifloc(dict_canon, dict_modif):
    """
    Determine the effects of each base by their position in the kmer.

    Args:
        dict_canon (dict): A dictionary of canonical kmers and their values.
        dict_modif (dict): A dictionary of modified kmers and their values.

    Returns:
        list: A list of effects of each base by their position in the kmer.
            Each element in the list is a list containing the following information:
            - The position of the motif in the k-mer.
            - The median effect of the motif at that position on all matching kmers.
            - The effect of bases at other positions in the k-mer.
    """
    all_effects = []
    untrained = 0
    for i in range(mer_size):
        diff_dict = {}
        for m_kmer, m_val in dict_modif.items():
            if m_kmer[i] != 'M':
                continue
            if m_kmer.count('M') > 1:
                continue
            canon = replace_edge(m_kmer.replace(mod_motif_M, mod_motif))
            diff = dict_canon[canon] - m_val

            # If mod and canon are exactly the same, mod wasn't actually trained
            if diff == 0:
                untrained += 1
                continue
            diff_dict[m_kmer] = diff
        if diff_dict:
            effect_all, effect_per = process_subdict(diff_dict)
            effect_per = flatten_effdict(effect_per)
            all_effects.append([i, effect_all, effect_per])
    print('Mod dict size:', len(dict_canon), 'Untrained:', untrained)
    return all_effects


def effs_to_dict(all_effects):
    """
    Convert the effects list to a dictionary.

    Args:
        all_effects (list): A list of effects.

    Returns:
        dict: A dictionary of effects.
    """
    eff_dict = {}
    for modpos in all_effects:
        eff_dict[str(modpos[0]) + 'M'] = modpos[1]
        for basepos in modpos[2]:
            eff_dict[str(modpos[0]) + 'M' + str(basepos[0]) + str(basepos[1])] = basepos[2]
    return eff_dict


def impute_mod(eff_dict, kmer):
    """
    Impute the modification effect for the given kmer.

    Args:
        eff_dict (dict): A dictionary of effects.
        kmer (str): The kmer to impute the modification effect for.

    Returns:
        float: The imputed modification effect.
    """
    delta = 0
    lid = ''
    # For every position of M in the k-mer (so multiple M's are handled too)
    for mi, mbase in enumerate(kmer):
        if mbase != 'M':
            # Not a modified base so nothing to do here
            continue
        else:
            if str(mi) + 'M' in eff_dict:
                # Add effect of having a modified base at index mi
                delta += eff_dict[str(mi) + 'M']
            for oi, obase in enumerate(kmer):
                # Add effect of every other base in kmer
                lid = str(mi) + 'M' + str(oi) + obase
                if lid in eff_dict:
                    delta += eff_dict[lid]
    return delta


def my_shifter(dict_canon, train):
    """
    Determine the effects on each base by position in the kmer.

    Args:
        dict_canon (dict): A dictionary of canonical kmers and their values.
        train (dict): A dictionary of modified kmers and their values.

    Returns:
        dict: A dictionary of effects.
    """
    sub_effects = per_modifloc(dict_canon, train)
    eff_dict = effs_to_dict(sub_effects)
    return eff_dict


def generate_table(dict_canon, eff_dict, dict_canon_sd, keepers={}):
    """
    Generate a table of kmers with their expected values, standard deviations, and the imputaton decision.

    Args:
        dict_canon (dict): A dictionary of canonical kmers and their values.
        eff_dict (dict): A dictionary of effects.
        dict_canon_sd (dict): A dictionary of canonical kmers and their standard deviations.
        keepers (dict): A dictionary of kmers to be kept without imputation.

    Returns:
        list: A list of kmers with their expected values, standard deviations, and the imputaton decision.
    """
    model_table = []
    bases = ['A', 'C', 'G', 'T', 'M']
    bases.sort()
    for p in itertools.product(bases, repeat=mer_size):
        kmer = ''.join(p)
        canon = replace_edge(kmer.replace(mod_motif_M, mod_motif))
        exp_val, exp_sd, info = 0, 0, 'E'

        if kmer in keepers:
            # Don't impute
            exp_val, exp_sd, info = keepers[kmer], dict_canon_sd[kmer], 'Kept'
        elif 'M' in canon:
            # kmer still has an M somewhere so it wasn't really a CG position
            canon = canon.replace("M", "C")
            exp_val, exp_sd, info = dict_canon[canon], dict_canon_sd[canon], 'LoneM'
        elif kmer.count(mod_motif_M) == 0:
            exp_val, exp_sd, info = dict_canon[canon], dict_canon_sd[canon], 'Canon'
        else:
            new_val = dict_canon[canon] - impute_mod(eff_dict, kmer)
            exp_val, exp_sd, info = new_val, dict_canon_sd[canon], 'ModFix'
        model_table.append([kmer, exp_val, exp_sd, info])
    return model_table


def write_model(new_table, out_path):
    """
    Write the model table to a file.

    Args:
        new_table (list): A list of kmers with their expected values, standard deviations, and information.
        out_path (str): The file path to write the model table to.
    """
    with open(out_path, 'w') as out_file:
        out_file.write('#Re-quant adjusted model based on pre-existing HMM model\n')
        out_file.write('#model_name\tr9.4_450bps.' + mod_alphabet + '.' + str(mer_size) + 'mer.requant.model\n')
        out_file.write('#kit\tr9.4_450bps\n')
        out_file.write('#strand\ttemplate\n')
        out_file.write('#alphabet\t' + mod_alphabet + '\n')
        for kmer_line in new_table:
            out_file.write(kmer_line[0] + '\t' + '\t'.join([str(x)[:7] for x in kmer_line[1:]]) + '\n')


def write_effs(eff_dict, out_path):
    """
    Write the effects dictionary to a file.

    Args:
        eff_dict (dict): A dictionary of effects.
        out_path (str): The file path to write the effects dictionary to.
    """
    sorted_eff_dict = sorted(eff_dict.items(), key=lambda x: abs(x[1]), reverse=True)
    with open(out_path, 'w') as out_file:
        for x in sorted_eff_dict:
            out_file.write(x[0] + '\t' + str(x[1]) + '\n')


def run_requant(f_np_model, f_imputed, mod, kmer_size=6):
    """
    Run the requantization process.

    Args:
        f_np_model (str): The file path of the partially trained nanopolish model.
        f_imputed (str): The base path for the output model files.
        mod (list): A list containing the modification definition [canonical, modified, alphabet].
        kmer_size (int, optional): The size of the kmer. Defaults to 6.
    """
    # Not much point passing these around all the time
    global mod_motif, mod_motif_M, mod_alphabet, mer_size
    mer_size = kmer_size
    # Setup modification motif information
    mod_motif, mod_motif_M, mod_alphabet = mod
    # Load the partially trained model
    dict_canon, dict_modif, dict_main_sd = load_polishmodel(f_np_model)

    # Determine the effect of on each base by position in the kmer
    eff_dict = my_shifter(dict_canon, dict_modif)
    # Write the effects to a table
    write_effs(eff_dict, f_imputed + '.rules')

    # Impute the full model but retain trained kmers
    add_table = generate_table(dict_canon, eff_dict, dict_main_sd, dict_modif)
    # Write imputed model to a table
    write_model(add_table, f_imputed + '.added.model')

    # Impute the full model and overwrite trained kmers
    rep_table = generate_table(dict_canon, eff_dict, dict_main_sd)
    # Write replaced model to a table
    write_model(rep_table, f_imputed + '.replaced.model')


def main(args):
    """
    The main function that parses command line arguments and runs the requantization process.

    Args:
        args (list): A list of command line arguments.
    """
    parser = argparse.ArgumentParser(prog='PROG')
    parser.set_defaults(func=lambda args: parser.print_help())
    parser.add_argument('input_model', help='input (partially trained) nanopolish model')
    parser.add_argument('output_base', help='output model base path, will add suffixes: .added.model and .replaced.model')
    parser.add_argument('-m', nargs=3, default=['CG', 'MG', 'cpg'], help='modification definition [CG MG cpg]',
                        metavar=('can', 'mod', 'alph'))
    parser.add_argument('-k', type=int, default=6, help='k-mer size for model [6]',
                        metavar=('k-mer_size'))

    args = parser.parse_args(args)
    run_requant(args.input_model, args.output_base, args.m, args.k)


if __name__ == '__main__':
    main(sys.argv[1:])