Coalescent Relationship Matrix
========

Coalescent Relationship Matrix (CRM) is the expected Genetic Relationship Matrix (GRM) on unknown SNPs 
given the complete genealogical tree, derived under the coalescent theory framework.


Installation
------------

Install from PyPI:

    pip install crm

Or download the package and install from local:

    pip install ./crm


Command Line Tool
-----------------

There are two command line tools:

    trees2crm [--input INPUT]

Where INPUT is the tree sequence filename which must end with '.trees'.
It outputs the CRM matrix in a file with the same name but '.crm' suffix.

    workflow [--out OUT] [--name NAME] [--gcta] [--relate]

Where OUT is the output directory, and NAME is the prefix of all output files associated with this simulation.
Include --gcta to run the GCTA analysis.
Include --relate to run relate tree reconstruction (which is required by Km_relate).
Use 

    workflow -h

to see more parameters for the simulation.


Support
-------

If you are having issues, please let us know.
Email the author: caoqifan@usc.edu

