Expected Genetic Relationship Matrix
========

Expected Genetic Relationship Matrix (EGRM) is the expected value of the Genetic Relationship Matrix (GRM) on unknown SNPs 
given the complete genealogical tree of a sample of individuals, derived under the coalescent theory framework.


Installation
------------

Install from PyPI:

    pip install egrm

Or download the package and install from local:

    pip install ./egrm


Command Line Tools
-----------------

There are two command line tools:

    trees2egrm [--input INPUT]

Where INPUT is the tree sequence filename which must end with '.trees'.
It outputs the CRM matrix in a file with the same name but '.egrm' suffix.

    workflow [--out OUT] [--name NAME] [--gcta] [--relate]

Where OUT is the output directory, and NAME is the prefix of all output files associated with this simulation.
Results are written into OUT/NAME, Running log info is written into OUT/NAME.log.
This workflow simulates the haplotypic and phenotypic data, and runs phenotype imputation based on several relationship matrices.
Include --gcta to run the GCTA analysis.
Include --relate to run relate tree reconstruction (which is required by Km_relate).
Use 

    workflow -h

to see more parameters for the simulation.



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

