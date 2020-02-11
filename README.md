# **Vorpal**
A collection of scripts that implements the Vorpal RNA feature extraction algorithm.

To ensure the functionality of these scripts, it is strongly advised to use the included `vorpal.yml`
to reconstruct the conda environment for which the dependencies are known to fuction correctly. In newer versions of pandas, numpy, and scipy the implementation of sparse dataframes has changed and will break this implementation of the Vorpal. We plan to fix this at some point but, for now, the easiest thing to do is create a specific Vorpal conda environment for the workflow.

First, clone the library:

    git clone https://www.github.com/mriglobal/vorpal.git

Assuming you have conda installed already (in case you don't please go to https://docs.conda.io/en/latest/ for instructions on installation) create the Vorpal conda environment: `conda env create --file vorpal.yml`
Activate your environment: `conda activate vorpal`

You should be able to execute the scripts in the library now.

## Overview

 1. kmercountouter_sparse.py

>     usage: kmercountouter_sparse.py [-h] -r R --seqs SEQS [-k [K]] [-p P]
>     
>     Kmer Counter for RNA viruses
>     
>     optional arguments:
>       -h, --help   show this help message and exit
>       -r R         Reference genome
>       --seqs SEQS  Folder for Sequences
>       -k [K]       K size for counting: default 18
>       -p P         Percent Variance in reference length for replicon binning
> 

 2. hammingclusters_fast.py

>      usage: hammingclusters_fast.py [-h] -p P [-n N] [-q Q] [-c C] [--temp TEMP]
>                                    [--mem MEM]
>     
>     Create Degenerate Primers from Kmer Clusters
>     
>     optional arguments:   -h, --help   show this help message and exit  
>     -p P         Pickled Kmer Sparse DataFrame   
>     -n N         Average Number of Allowed Degenerate Bases   
>     -q Q         Quantile of clustered Kmers   
>     -c C         Number of chunks to split matrix into for count processing.
>                    Default 0   
>     --temp TEMP  Location of temp directory for memory mapped matrix. (Default:
>                    None)   
>     --mem MEM    Amount of memory allocated to process distance matrix chunks.
>                    (In MiB)

 3. 

> Written with [StackEdit](https://stackedit.io/).


