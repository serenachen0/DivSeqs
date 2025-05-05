# DivSeqs
To construct an initial enzyme library of maximally diverse homologs, we provide a python script that applies HHfilter (from Söding Lab) to screen an MSA based on a defined minimum coverage with the query sequence and an increasing maximum pairwise sequence identity cutoff value.

Note: install HH-suite3 before running the python script. Follow the instructions at https://github.com/soedinglab/hh-suite for installation.

File required: an MSA file in A3M format

**To run the python script:** 

python3 hhfilter-dynamic.py \<hhfilter\> \<inpmsa\> --cov (optional) --idmin (optional) --idmax (optional) --outdir (optional) --pid (optional)

- hhfilter: path to the hhfilter executable
- inpmsa: path to the input MSA file in A3M format
- cov: minimum sequence coverage with the query sequence (the first sequence in the input MSA file)
- idmin: minimum pairwise sequence identity cutoff in dynamic filtering
- idmax: maximum pairwise sequence identity cutoff in dynamic filtering
- outdir: path to the output directory
- pid: protein ID to be used for output file name and plot

For help: python3 hhfilter-dynamic.py -h

**To run the python script using the provided example:**

Sort the homologs in the MSA using a minimum of 50% coverage with the query (first) sequence and a maximum pairwise sequence identity ranging from 25% to 95%
   
   python3 hhfilter-dynamic.py \<hhfilter\> ./example/NylC_inpmsa.a3m --cov 50 --idmin 25 --idmax 95 --outdir ./example --pid NylC

Outputs of the python script include

1. an .a3m file with sorted homologs based on increasing pairwise sequence identity cutoff at the defined minimum sequence coverage

2. a .pdf file showing a histogram of the homologs sorted by pairwise sequence identity

Outputs for the example can be found in example


