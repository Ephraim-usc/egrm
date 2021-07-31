Reproducing Results in the eGRM Manuscript
========

This is a tutorial introducing how to reproduce simulation results in the eGRM manuscript, using the [bin/simulate] script.

The output is a folder named as [name], including two files:

-   results.p: a pickle file storing a python dict object containing statistical results from this simulation

-   simulation.trees: a tskit tree sequence file containing the tree sequence of the study panel

Usage example:

    simulate --nrow simulate --name simulation01 --demo ooa --nrow 2 --ns 1000 0 --ns_ref 0 1000 --time_move 8000 

This will simulate a grid of 2 rows and 1 column of demes, 
the population size history of each deme follows the pre-defined out-of-Africa demography.
The two demes were splitted 8000 generations ago, from a single deme with the same demography.
The first deme contains 1000 study panel haplotypes, and the second deme contains 1000 reference panel haplotypes.
The reference panel only affects results related to the K_obs_imputed (standard GRM with imputation).

To run the [bin/simulate] script, make sure you have the following commands ready:

    ~/bin/plink2    if --run_prune
    
    ~/bin/impute2    if --run_impute
    
    ~/bin/relate    if --run_relate
    
    ~/bin/tsinfer and ~/bin/tsdate    if --run_tsinfer

To reproduce results published in the manuscript, please submit a series of [simulate] jobs with appropriate parameters,
and collect results in the [results.p] file in each output directory.


Command Line Parameters
-----------------

    --name NAME

This is only required parameter. The output folder will be called NAME.

    --run_egrm

If included, the [EK] matrix and related results will be computed.

    --run_relate

If included, the [EK_relate] matrix and related results will be computed.

    --run_tsinfer

If included, the [EK_tsinfer] matrix and related results will be computed.

    --run_mtmrca

If included, the [mtmrca] matrix and related results will be computed.

    --run_all

Including this option is equivalent to including --run_egrm, --run_relate, --run_tsinfer and --run_mtmrca at the same time.

    --demo DEMO

The demography to be used. Use "--demo ooa" to let each deme follow the out-of_Africa demography.
Use "--demo constant" to let each deme to have a constant population size decided by [pop_size].
Any other DEMO values will be treated as "constant".

    --nrow NROW

Number or rows in the grid of demes.

    --ncol NCOL

Number or columns in the grid of demes.

    --migration_rate MIGRATION_RATE

Migration rate between neighboring demes.

    --pop_size POP_SIZE

Population size of each deme. Only takes effect when "--demo constant"

    --N, --n, --ns

These three pameters specify the number of sample haplotypes in each deme. Only one of these three should be provided. 
--N specifices the total number of study haplotypes.
--n specifices the (same) number of study haplotypes in each deme.
--ns specifices a list containing numbers of study haplotypes in each deme. List should be the same length as [nrow]*[ncol].

    --ns_ref

Specifies a list containing numbers of reference haplotypes in each deme.




Support
-------

If you are having issues, please let us know.
Email the author: caoqifan@usc.edu

