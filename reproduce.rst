Reproducing Results in the eGRM Manuscript
========

This is a tutorial introducing how to reproduce simulation results in the eGRM manuscript, using the bin/simulate script.

The output is a folder named as [name], including two files:

-   results.p: a pickle file storing a python dict object containing results from this simulation

-   simulation.trees: a tskit tree sequence file containing the tree sequence of the study panel

Usage example:

    simulate --nrow simulate --name simulation01 --demo ooa --nrow 2 --ns 1000 0 --ns_ref 0 1000 --time_move 8000 

This will simulate a grid of 2 rows and 1 column of demes, 
the population size history of each deme follows the pre-defined out-of-Africa demography.
The two demes were splitted 8000 generations ago, from a single deme with the same demography.
The first deme contains 1000 study panel haplotypes, and the second deme contains 1000 reference panel haplotypes.
The reference panel only affects results related to the K_obs_imputed (standard GRM with imputation).


Command Line Tools
-----------------

There are two command line tools:

    trees2egrm [--input INPUT] [--output OUTPUT]
    
    trees2mtmrca [--input INPUT] [--output OUTPUT]

Where INPUT is the tree sequence file prefix (so that the full name should be "INPUT.trees"), and OUTPUT is the output file prefix.

Optional parameters:

    [--c_extension] or [--c]

This specifies whether to use the C exntension model to accelerate the algorithm.
Usually this makes it ~10 times faster.
Recommended whenever the C environment is available.

    [--skip_first_tree] or [--sft]

This option skips the first tree in the tree sequence.
This is often useful because RELATE and some other tools always output tree sequences from 0bp, even when the genotype data starts from within the chromosome.

    [--run_var] or [--var]

This option is only for trees2egrm, not trees2mtmrca.
With this option turned on, the algorithm will output the varGRM in addition to eGRM, while roughly doubling the compuation time.

    [--left LEFT] [--right RIGHT]

The leftmost and rightmost positions (in bp) between which the eGRM or mTMRCA is computed.

    [--rlim RLIM] [--alim ALIM]

This option is only for trees2egrm, not trees2mtmrca.
RLIM and ALIM are the most recent and most ancient times (in generations) between which the eGRM is computed.

The output of trees2egrm will be two (or three, if with --var option) files in numpy NPY format: 

-   OUTPUT.npy, which contains the eGRM matrix;

-   OUTPUT_mu.npy, which contains a single number of the measure of the tree sequence (i.e., the expected number of mutations on this tree sequence);

-   OUTPUT_var.npy, which contains the varGRM matrix, if the --var option is selected.

The output of trees2mtrmca will be two files in numpy NPY format: 

-   OUTPUT.npy, which contains the mTMRCA matrix;

-   OUTPUT_l.npy, which contains a single number of the number of base pairs of the tree sequence;


Python Functions
-----------------

    varGRM_C(trees)
    
    varGRM(trees)

The C and non-C versions of the eGRM algorithm. The input is a tskit TreeSequence object.
See the source code for a complete explanation of its parameters.

    mTMRCA_C(trees)
    
    mTMRCA(trees)

The C and non-C versions of the mTMRCA algorithm. The input is a tskit TreeSequence object.
See the source code for a complete explanation of its parameters.


Reproducing Results in the paper
-----------------

There is an additional commandline tool

    simulate 

which is included in the package, but not installed by default. You may manually run this script.

A complete explanation of its parameters and output files can be found at

    simulate -h


Support
-------

If you are having issues, please let us know.
Email the author: caoqifan@usc.edu

