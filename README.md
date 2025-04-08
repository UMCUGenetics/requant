# Preface
This work aims to impute missing modified k-mer values in table based models used for modfication calling on nanopore sequencing data.

This repo/branch/release contains all scripts as ran for the ReQuant manuscript:

>__ReQuant: Improved base modification calling by k-mer value imputation__
>
>Roy Straver, Carlo Vermeulen, Joe R. Verity-Legg, Marc Pagès-Gallego, Dieter G.G. Stoker, Alexander van Oudenaarden and Jeroen de Ridder

Most of this was developed on a laptop, then relevant scripts were copied to an HPC to apply them to more data. Lots of minor edits happened to make things work for particular datasets hence some scripts occur multiple times in various folders, reflecting what analysis exactly they were involved in.

For general usage, I suggest running `requant.py` and ignoring the rest. After various versions to produce plots I removed all the irrelevant code for the algorithm itself and turned it into a command-line tool. I tested `requant.py` on a few models to ensure it produces the exact same output as the models in the paper, but it can now change the modification and k-mer size using arguments rather than having to change hard-coded values. All the other files are just included for reproducability.


# Usage
Running `requant.py` requires a python environment with numpy installed. Among others as can be found in the .yml files in this repo, these versions worked at the time of writing:
python 3.9
numpy 1.21.5

The arguments for ReQuant can be requested using -h and should be self explanatory:
```python
>python requant.py -h

usage: PROG [-h] [-m can mod alph] [-k k-mer_size] input_model output_base

positional arguments:
  input_model      input (partially trained) nanopolish model or similar
  output_base      output model base path, ReQuant will add suffixes:
                   .added.model and .replaced.model

optional arguments:
  -h, --help       show this help message and exit
  -m can mod alph  modification definition [CG MG cpg]
  -k k-mer_size    k-mer size for model [6]
```


# Reproducing data/plots
Any .ipynb was ran locally on a laptop to make some plots or get some numbers from the data. These are all located in the `./laptop/` directory, together with a conda export of the environment they were ran in (`conda_env_nanopore.yml`).
 - `event_align.ipynb`: Show a comparison of raw data between modified and unmodified data as prepared with scripts in `./hpc/np_requant/event_align`.
 - `plot_output_model.ipynb`: Create scatterplots from (imputed) models in nanopolish' table-ish format, used to create plots for figure 2.
 - `remora_callstat.ipynb`: Show remora having difficulties with limited training data as produced with scripts in `./hpc/remora_tests`.
 - `requant_callstat_2.ipynb`: Show nanopolish calling results with and without ReQuant applied, data generated with scripts in `hpc/np_requant`.
 - `rules_stat.ipynb`: Show the rules that ReQuant determined for each run and compare the resulting models directly, data produced with scripts in `hpc/np_requant`.
 - `f5c_stat.ipynb`: Show plots on the R10 modification calling results as produced with scripts in `./hpc/f5c_methcall`.

The `mark_fasta.py` and `split_fasta.py` were both ran on the same laptop and the fasta files they produced were copied to the HPC for further processing. 
 - The `mark_fasta.py` script marks random occurrences of the motif until the target number of k-mers is covered by at least one context each. 
 - Remora seemed hardcoded to assume all occurrences of a motif are either affected, or all of them are not. Hence the `split_fasta.py` just splits the genome in some sections that cover roughly the target percentage of k-mers with a motif, writes a file with the sequence for each section, and we randomly picked a few of these sections to train and test with.
For reproducability the actual fasta files used for the ReQuant paper are zipped in `./hpc/refcuts.tar.xz`.


The actual version of ReQuant ran to produce (most of) the data is `./hpc/np_requant/requant_impute.py`. This and nearly all other scripts in `./hpc/` and subdirectories of it were ran in the conda environment exported as `./hpc/conda_env_poretools.yml`, except for f5c related steps, which used conda environment `./hpc/conda_env_remora.yml`.
The difference between the `./hpc/`'s `requant.py` and `requant_impute.py` is mostly the k-mer sampling and plotting functionality. These were useful for development but useless and rather slow for bigger data analyses, hence `requant_impute.py` implements the high level logic without any of this, and imports the lower level functions from `requant.py`.
Should you want to reproduce the GC motif analysis instead of CG, or change the k-mer size, you'll need to adjust both of these files to match.

For example, `requant.py` can toggle from CG to GC by uncommenting the last line from these, and change the k-mer size by altering the first line:
```
mer_size = 6
mod_motif, mod_motif_M, mod_alphabet = 'CG', 'MG', 'cpg'
# mod_motif, mod_motif_M, mod_alphabet = 'GC', 'GM', 'gpc'
```

Similarly for `requant_impute`, toggle with these lines:
```
motif, motif_mod = 'CG', 'MG'
# motif, motif_mod = 'GC', 'GM'
```

Beyond this ReQuant is indifferent to the actual motifs to work on, so feel free to experiment.

Lots of the `.sh` scripts in `/hpc/` are just (somewhat generalized) job (submission) scripts. Anything ending with `_subloop.sh` is there to submit processing of repeats for various k-mer "coverage" values by submitting the similarly named but non-`_subloop` `.sh` script. Finally, statistics for each of these were obtained using `callstat_2` scripts.

References to files are not very organized, they are just where the data happened to be. If you really want to reproduce the manuscript to the dot, you'll have to alter these paths as necessary. The `nanopolish_og` that ran in many of these scripts is just an unaltered version of Nanopolish.
