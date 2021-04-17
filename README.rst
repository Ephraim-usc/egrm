Expected Genetic Relationship Matrix
========

Expected Genetic Relationship Matrix (EGRM) is the expected value of the Genetic Relationship Matrix (GRM) on unknown SNPs 
given the complete genealogical tree of a sample of individuals, derived under the coalescent theory framework.


Installation
------------

Install from PyPI (not available yet):

    pip install egrm

Or download the package and install from local:

    git clone https://github.com/Ephraim-usc/egrm.git
    
    pip install ./egrm


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
(1) OUTPUT.npy, which contains the mTMRCA matrix;
(2) OUTPUT_l.npy, which contains a single number of the number of base pairs of the tree sequence;



Python Functions
-----------------

    simulate( ... )

This will run the whole simulation process (of haplotypic data, of phenotypes and of observed SNPs) and 
store all output data in a nested dictionary structure.
It has many optional parameters to setup the simulation conditions.

    getEK(tree)

This will compute the expected GRM based on the tskit.Tree object.
The output is a numpy matrix.

    getEK_trees(trees)

This will compute the expected GRM based on the tskit.TreeSequence object.
It is just the weighted mean (by genomic interval length) of all outputs of getEK on each tree.
The output is a numpy matrix.

Support
-------

If you are having issues, please let us know.
Email the author: caoqifan@usc.edu

