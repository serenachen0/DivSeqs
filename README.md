# DivSeqs
To construct an initial enzyme library of maximally diverse homologs, we provide a python script that applies HHfilter (from Söding Lab) to screen a MSA based on a defined minimum coverage with the query sequence and an increasing maximum pairwise sequence identity cutoff value.

Note: intall HH-suite3 before running the python script. Follow the instructions at https://github.com/soedinglab/hh-suite for installation.

File required: a MSA file in A3M format

To run the python script: 

python3 hhfilter-dynamic.py \<hhfilter executable\> \<input MSA\> \<minimum sequence coverage with query\> \<minimum pairwise sequence identity cutoff in dynamic filtering\> \<output directory\> \<protein ID\>

For help: python3 hhfilter-dynamic.py -h

To run the python script using the provided examples:

1. NylC— sort the homologs in the MSA using a minimum of 50% coverage with the query (first) sequence and a maximum pairwise sequence identity ranging from 25% to 95%
   
   python3 hhfilter-dynamic.py \<hhfilter executable\> ./examples/inputs/NylC_inpmsa.a3m 50 25 ./examples/outputs NylC

2. PghP— sort the homologs in the MSA using a minimum of 50% coverage with the query (first) sequence and a maximum pairwise sequence identity ranging from 0% to 95%
   
   python3 hhfilter-dynamic.py \<hhfilter executable\> ./examples/inputs/PghP_inpmsa.a3m 50 0 ./examples/outputs PghP

Outputs of the python script include

1. an .a3m file with sorted homologs based on increasing pairwise sequence identity cutoff at the defined minimum sequence coverage

2. a .pdf file showing a histogram of the homologs sorted by pairwise sequence identity

Outputs for the two examples can be found in examples/outputs/


