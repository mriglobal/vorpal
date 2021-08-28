﻿# **Vorpal**
<p align=center>
<img src="vorpal_graphic.png"/>
</p>
A collection of scripts that implements the [Vorpal RNA feature extraction algorithm](https://www.biorxiv.org/content/10.1101/2020.02.28.969782v1).

To ensure the functionality of these scripts, it is strongly advised to use the included `vorpal.yml`
to reconstruct the conda environment for which the dependencies are known to fuction correctly. In newer versions of pandas, numpy, and scipy the implementation of sparse dataframes has changed and will break this implementation of the Vorpal. We plan to fix this at some point but, for now, the easiest thing to do is create a specific Vorpal conda environment for the workflow.

First, clone the library:

    git clone https://www.github.com/mriglobal/vorpal.git

Assuming you have conda installed already (in case you don't, please go to https://docs.conda.io/en/latest/ for instructions on installation) create the Vorpal conda environment: `conda env create --file vorpal.yml`
Activate your environment: `conda activate vorpal`

You should be able to execute the scripts in the library now.

## Models

If you quickly want to evaluate the models skip to the [Evaluating the Models](#evaluating-the-models) section to get started.

## Overview

1. kmercountouter_sparse.py

This script will take a specified reference sequence, `-r`, and use the length of that sequence to bin the sequences specified in the `--seqs` directory. Those that are within the length criteria of `-p` will proceed to canonical K-mer counting where K size is specified in `-k`. The `-n` and `-m`  arguments allow for a stratified resampling of the input sequences to be employed before counting, where N is the number of samples drawn (with replacement) and M is a tabular text file containing "accession" and "groups" columns assigning a group label for each input record. The ouput will be a pickled sparse dataframe object of shape K-mers x Sequence IDs. This orientation is for convenience in the high-frequency filtering stage.

    usage: kmercountouter_sparse.py [-h] -r R --seqs SEQS [-n N] [-m M] [-k [K]]
                                    [-p P]
    
    Kmer Counter for RNA viruses
    
    optional arguments:
      -h, --help   show this help message and exit
      -r R         Reference genome
      --seqs SEQS  Folder for Sequences
      -n N         Number of samples for stratified resampling.
      -m M         Metadata file containing group labels for resampling.
      -k [K]       K size for counting: default 17
      -p P         Percent Variance in reference length for replicon binning

2. hammingclusters_fast.py
 
This script will take the sparse dataframe output from `kmercountout_sparse.py` and perform hierarchical clustering on the resulting K-mers. There are many user parameters to specify. First, `-q` specifies the quantile for K-mer frequency cutoff to proceed to clustering. In other words, a .95 quantile means that only the top 5% of abundant K-mer will be clustered. The abundance calculation is defined in the manuscript. This parameter is used mostly as a tool to make calculation of the linkage function tractable given the memory requirements. The most important parameter regarding motif creation is `-n`, which specifies where the resulting tree should be cut to form clusters. The parameter is specified as average number of degenerate bases desired, and that number is converted to a fraction of the length of the K-mer size. Various additional arguments are for managing memory contraints. Memory mapping of the distance matrix with --temp allows for memory to be freed for the linkage calculation. Both the K-mer frequency counting and distance matrix calculations can been chunked as well using the `-c` and `--mem` arguments respectively.

    usage: hammingclusters_fast.py [-h] -p P [-n N] [-q Q] [-c C] [--temp TEMP]
                                   [--mem MEM]
    
    Create Degenerate Primers from Kmer Clusters
    
    optional arguments:
      -h, --help   show this help message and exit
      -p P         Pickled Kmer Sparse DataFrame
      -n N         Average Number of Allowed Degenerate Bases
      -q Q         Quantile of clustered Kmers
      -c C         Number of chunks to split matrix into for count processing.
                   Default 0
      --temp TEMP  Location of temp directory for memory mapped matrix. (Default:
                   None)
      --mem MEM    Amount of memory allocated to process distance matrix chunks.
                   (In MiB)

 3. referencemapping_mp.py
 
 This script uses a [Biopython SeqUtils](https://biopython.org/docs/1.75/api/Bio.SeqUtils.html#Bio.SeqUtils.nt_search) function to map the degenerate motifs fasta file, `-k`, resulting from `hammingclusters_fast.py`  to the concatenated sequences specified by `-r`. The resulting alignments (forward strand only) are output as a series of discrete BED format files for each unique sequence specified in `-r`, to the current working directory. Chunking is also possible using the `-c` argument (another memory constraint work-around).
 

     usage: referencemapping_mp.py [-h] -r R -k K [-t T] [-c C]
    
    Mapping Degenerate Kmers to Reference Sequences
    
    optional arguments:
      -h, --help  show this help message and exit
      -r R        Concatenated References Fasta
      -k K        Degenerate Kmers Fasta
      -t T        Number of threads. Default all.
      -c C        Number of chunks to break references into and write to disk.
                  (Default=1)

## Model Fitting Routines

Scripts handling model fitting routines are currently:

 1. binary_vorpal_model.py
 2. binary_vorpal_model_ElasticNet.py
 3. binary_vorpal_model_GSS.py
 
 The former two were used in the original, explanatory modeling work. Here, we'll focus on the last one, binary_vorpal_model_GSS.py, which was used in the predictive modeling effort.

All three scripts above share similar command line signatures. A directory where the BED format files produced by referencemapping_mp.py is specified by `--beds`.  A tabular format metadata file containing minimally, "accession", "label", and in the case of the Group-Shuffle-Split script, "groups" field is specified by the `-m` argument. A truncated example of this file is provided below:

    accession	label	species	groups
    KJ156866.1	1	1335626	1_1335626
    KM027262.1	1	1335626	1_1335626
    AY545916.1	1	694009	1_694009
    DQ415907.1	1	290028	1_290028
    DQ415909.1	1	290028	1_290028
    MG923478.1	1	1335626	1_1335626
    MF598625.1	1	1335626	1_1335626
    KR011265.1	1	1335626	1_1335626
    KF923895.1	1	694003	1_694003

Arguments specific to the GSS fitting routine involve cross validation parameters. The `-s` argument specifies the fraction of total groups withheld for each validation set split. The `-n` argument specifies the number of cross validation splits that will be performed. Finally, the `-r` argument implements a stratified resampling of the input data using the "groups" field where the input data is resample the number of times provided to `-r`. The remaining arguments regard coordinate descent parameters for liblinear, as well as number of processes to use for cross validation.

    usage: binary_vorpal_model_GSS.py [-h] --beds BEDS -m M [-o O] [-s S] [-n N]
                                      [--RVDB] [-i I] [-p P] [-t T] [-r R]
    
    Takes feature-labeled .bed files and corresponding meta data as inputs to
    generate a sparse model of phenotype predictors.
    
    optional arguments:
      -h, --help   show this help message and exit
      --beds BEDS  Directory containing .bed files.
      -m M         Meta data and groups table for genomic records.
      -o O         Prefix for output files.
      -s S         Fraction size for group splits. Default: 0.10.
      -n N         Number of splits for groups splits. Default: 100
      --RVDB       Flag for RVDB fasta headers.
      -i I         Number of iterations for coordinate descent.
      -p P         Number of processors to use. Default: Max available
      -t T         Min loss tolerance for stopping. Default: .00000001
      -r R         Number of resampled rows using stratified groups. (Default is
                   no resample)

## Example Workflow
The workflow for extracting features and fitting new models is outlined here. Before you begin should have a specified reference genome from the viral clade of interest (this is mostly for segment binning purposes. The output k-mers file name will also inherit the fasta header name from this specified reference sequence.) You will need a multifasta of sequences you are interested in counting k-mers from. You will also need a tab separated file that contains at least an "accession" column containing the accession numbers of the sequences in the multi-fasta as well as a "label" column that contains a binary class label for each accession number i.e. a 0 or a 1. If you want to take advantage of stratified resampling, and Group-Shuffle-Split cross validation utilized in `binary_vorpal_model_GSS.py`, you will also have to have a "groups" column in this text file that contains discrete group labels for each respective record. As an example, the data labels file used to fit the model in /models directory is provided in the /data directory in the repository with the name `RVDB14_Coronavirus_meta_data_group_human_pathogen.tsv`. The sequences of script execution for feature extraction and model fitting are:

 1. `kmercountouter_sparse.py`
 2. `hammingclusters_fast.py` (you'll need a lot of RAM for this)
 3. `referencemapping_mp.py`
 4. `binary_vorpal_model_GSS.py`

## Evaluating Sequences with the Models

Included with the repository is a Jupyter Notebook template `vorpal_model_eval.ipynb`. This notebook contains commented instructions for taking the models deployed as part of the repository and quickly using it to evaluate new sequences. A quick workflow to evaluate the deployed model in /models is demonstrated through this example:

 1. In a desired location create a directory to contain the bed files used to assign the motifs. `mkdir 4.0_degen; cd 4.0_degen`
 2. Activate the vorpal environment `conda activate vorpal`
 3. Use `referencemapping_mp.py` to map the motif fasta `SARS_ref_15mers_0_1000resample_sparse4.0_quantile.90degenerate_primers.fasta` located in the /models directory using the `-k` argument to a  desired reference sequence multi-fasta. For the sake of this example we can use the `RVDB14_complete_coronavirus_fixed.fasta` sequences used to train the included model. This file can be found in the /data repository directory and specified to `referencemapping_mp.py` using the `-r` argument.
 4. Launch juptyerhub `juptyerhub`
 5. Make sure the *Vorpal* git repository is in a directory that is accessible by Jupyterhub and navigate to the repository directory on your file system. Launch `vorpal_model_eval.ipynb`
 6. Follow the instructions in the comment text to assign the variables to the correct file locations. For example, the persistent model object file `RVDB14_Coronavirus_meta_data_group_human_pathogen_15mer_4.0_rep3_CLF.joblib` should be assigned to the `clf` variable and the `RVDB14_Coronavirus_meta_data_group_human_pathogen_15mer_4.0_rep3_feature_array.pickle` file should be assigned to the `feature vector` variable. Training and test set data label files `RVDB14_Coronavirus_meta_data_group_human_pathogen.tsv` and `SARS2_SADS_test_set.tsv` are assigned to their respective variables `meta` and `test_meta`. The complete file path to these file locations will depend on location in your file system you have cloned the repository to. Finally, direct the `os.chdir` call in either the "training" or "test" cells to the directories that contain the bed files output in step 3.
 7. Run the cells of the notebook and observe the results. Convenient summary outputs to tabular files for further evaluation are provided in the cells labeled "misclassifieds out" and "test set probabilities out".
 8. If desired new sequences can predicted on using the "predict on new sequence" cell. Motifs should be mapped to these new sequences using the same procedure above with `referencemapping_mp.py` and directory containing the output beds specified in the `os.chdir` call in the same cell.

> Written with [StackEdit](https://stackedit.io/).






