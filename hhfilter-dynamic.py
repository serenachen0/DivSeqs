#!/usr/bin/env python3
"""
(1) Dynamic filtering and selection of diverse homologs using different pairwise sequence identity or sequence coverage cutoffs
(2) Generate a plot of number of unique homologs in each cutoff bin
"""

import argparse, os
import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt

def create_parser():
    parser = argparse.ArgumentParser(
        description="Dynamic filtering and selection of diverse homologs from a MSA in A3M format using HHfilter"
    )
    parser.add_argument(
        "hhfilter",
        type=str,
        help="Path to the hhfilter executable"
    )
    parser.add_argument(
        "inpmsa",
        type=str,
        help="Path to the input MSA file in A3M format"
    )
    parser.add_argument(
        "cov",
        type=int,
        help="Minimum coverage with the first (query) sequence in the input MSA file",
    )
    parser.add_argument(
        "idmin",
        type=int,
        help="Minimum pairwise sequence identity cutoff in dynamic filtering",
    )
    parser.add_argument(
        "outdir",
        type=str,
        help="Path to the output file directory"
    )
    parser.add_argument(
        "pid",
        type=str,
        help="Protein ID to be used for output file name and plot"
    )

    return parser

def run(args):
    exe = args.hhfilter
    msa = args.inpmsa
    cov = args.cov
    id_min = args.idmin
    outdir = args.outdir
    pid = args.pid
    
    id_stride = 5
    
    #output files
    outfl_a3m = outdir + "/" + pid + "_sorted_seqs.a3m" #store the homologs based on increasing pairwise seq iden at the defined min sequence coverage
    outfl_pdf = outdir + "/" + pid + "_sorted_hist.pdf" #show the histogram of the homologs sorted by pairwise seq iden

    fh_out = open(outfl_a3m, "w+")
    
    seqs = 0
    iden = id_min
    seqs_list = []
    unique_list = []
    for iden in range(id_min, 100, id_stride):
        #run HHfilter
        tmp = outdir + "/tmp.a3m"
        os.system(exe + ' -id ' + str(iden) + ' -cov ' + str(cov) + ' -i ' + msa + ' -o ' + tmp)
        
        #filter out existing/duplicate sequences
        unique_seq_count = 0
        fasta_sequences = SeqIO.parse(open(tmp),'fasta')
        for i, fasta in enumerate(fasta_sequences):
            id_line, sequence = fasta.description, str(fasta.seq)
 
            if not sequence in seqs_list:
                fh_out.write(">" + id_line + "\n" + sequence + "\n")
                seqs_list.append(sequence)
                seqs += 1
                unique_seq_count += 1

        unique_list.append(unique_seq_count)

    print ("Total number of unique sequences = {}".format(seqs))
    fh_out.close()
    os.remove(tmp)

    #plot number of unique homologs in each pairwise seq id cutoff
    plt.rcParams.update({'font.size': 16})
    ind = np.arange(id_min, 100, id_stride) # the x locations for the groups
    width = 3.0 # the width of the bars
    p = plt.bar(ind, unique_list, width, color="#1E90FF", alpha=1, edgecolor="#1E90FF")
    for i, x in enumerate(ind):
        plt.text(x, unique_list[i], "%d"%unique_list[i], ha="center", va="bottom", fontsize=12)

    plt.title(pid + '\n Minimum coverage with the query = ' + str(cov) + '%')
    plt.xlabel('Pairwise sequence identity (%)')
    plt.ylabel('Unique homolog count')
    plt.xticks(np.arange(id_min, 100, 10))
    plt.savefig(outfl_pdf, dpi=600, bbox_inches='tight', format="pdf")
    plt.clf()
    plt.close()
    
def main():
    parser = create_parser()
    args = parser.parse_args()
    run(args)

if __name__ == "__main__":
    main()
