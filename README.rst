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

Where INPUT is the tree sequence file prefix (so that the full name should be "INPUT.trees").
And OUTPUT is the output file prefix.
Optional parameters:

    [--c_extension] or [--c]

This specifies whether to use the C exntension model to accelerate the algorithm.
Usually this makes it ~10 times faster.
Recommended whenever the C environment is available.



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

