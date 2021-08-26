Expected Genetic Relationship Matrix
========

Expected Genetic Relationship Matrix (EGRM) is the expected value of the Genetic Relationship Matrix (GRM) on unknown SNPs 
given the complete genealogical tree of a sample of individuals, derived under the coalescent theory framework.

This method is described in the following [paper](https://www.biorxiv.org/content/10.1101/2021.08.18.456747v1.abstract):
Fan, Caoqi, Nicholas Mancuso, and Charleston WK Chiang. "A genealogical estimate of genetic relationships." bioRxiv (2021).
Please cite our paper if you use our method.


Installation
------------

Install from PyPI (not available yet):

    pip install egrm

Or download the package and install from local:

    git clone https://github.com/Ephraim-usc/egrm.git
    
    pip install ./egrm


Command Line Tools
------------------
The software to compute the eGRM from tskit tree sequence output is handled by `trees2egrm`. Its usage is given by,

    usage: trees2egrm [-h] [--output OUTPUT] [--c-extension] [--skip-first-tree] [--run-var] [--genetic-map GENETIC_MAP] [--left LEFT]
                      [--right RIGHT] [--rlim RLIM] [--alim ALIM] [--verbose] [--haploid] [--output-format {gcta,numpy}]
                      input

    Construct eGRM matrix from tree sequence data

    positional arguments:
      input                 path to ts-kit tree sequence file

    optional arguments:
      -h, --help            show this help message and exit
      --output OUTPUT, --o OUTPUT
                            output file prefix (default: egrm)
      --c-extension, --c    acceleration by C extension (default: False)
      --skip-first-tree, --sft
                            discard the first tree in the tree sequence (default: False)
      --run-var, --var      compute varGRM in addition to eGRM (default: False)
      --genetic-map GENETIC_MAP, --map GENETIC_MAP
                            map file fullname (default: None)
      --left LEFT, --l LEFT
                            leftmost genomic position to be included (default: 0)
      --right RIGHT, --r RIGHT
                            rightmost genomic position to be included (default: inf)
      --rlim RLIM           most recent time limit (default: 0)
      --alim ALIM           most ancient time limit (default: inf)
      --verbose             verbose logging. Includes debug info. (default: False)
      --haploid             output eGRM over haploids. Default is diploid/genotype eGRM. (default: False)
      --output-format {gcta,numpy}, --f {gcta,numpy}
                            output format of eGRM (default: gcta)

Where `input` is the tree sequence file prefix (so that the full name should be "INPUT.trees"), and `OUTPUT` is the output file prefix.

Optional parameters:

    [--c_extension] or [--c]

This specifies whether to use the C exntension model to accelerate the algorithm.
Usually this makes it ~10 times faster.
Recommended whenever the C environment is available.

    [--skip_first_tree] or [--sft]

This option skips the first tree in the tree sequence.
This is often useful because RELATE and some other tools always output tree sequences from 0bp, even when the genotype data starts from within the chromosome.

    [--run_var] or [--var]

With this option turned on, the algorithm will output the varGRM in addition to eGRM, while roughly doubling the compuation time.

    [--genetic_map] or [--map]

A (comma/space/tab separated) three-column file with first column specifying the physical position in bp and the third column specifying the genetic position in cM. The second column is not used. The first line will always be ignored as the header.

    [--left LEFT] [--right RIGHT]

The leftmost and rightmost positions (in bp) between which the eGRM is computed.

    [--rlim RLIM] [--alim ALIM]

RLIM and ALIM are the most recent and most ancient times (in generations) between which the eGRM is computed.

If output-format is set to `gcta`, the eGRM will be output into GCTA format (.grm.bin, .grm.N.bin and .grm.id files):

-   OUTPUT.grm.bin, which contains the eGRM matrix;

-   OUTPUT.grm.N.bin, which contains all the same number of the measure of the tree sequence (i.e., the expected number of mutations on this tree sequence);

-   OUTPUT.grm.id, which contains dummy ids for the samples.

-   OUTPUT_var(.grm.bin/.grm.N.bin/.grm.id), for the varGRM matrix, if the --var option is selected.

If If output-format is set to `numpy`, he output will be two (or three, if using --var option) files in numpy NPY format: 

-   OUTPUT.npy, which contains the eGRM matrix;

-   OUTPUT_mu.npy, which contains a single number of the measure of the tree sequence (i.e., the expected number of mutations on this tree sequence);

-   OUTPUT_var.npy, which contains the varGRM matrix, if the --var option is selected.


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

    manuscript/simulate 

which is included in the package, but not installed by default. You may manually run this script.

A complete explanation of its parameters and output files can be found at

    manuscript/reproduce.rst


Support
-------

If you are having issues, please let us know.
Email the author: caoqifan@usc.edu

