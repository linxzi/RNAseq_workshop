#!/usr/bin/env python

import os, sys

from sequence_to_kmer_list import *
from fastq_file_to_sequence_list import *


## method: count_kmers(kmer_list)
##
##  Counts the frequency of each kmer in the given list of kmers
##
##  input parameters:
##
##  kmer_list : list of kmers (type: list)
##               ie.  ["GATC", "TCGA", "GATC", ...]
##
##
##  returns kmer_counts_dict : dict containing ( kmer : count )
##                    ie.  {  "GATC" : 2,
##                            "TCGA" : 1,
##                             ...       }


def count_kmers(kmer_list):

    kmer_count_dict = dict()

    ##################
    ## Step 2:
    ## begin your code
    
    for item in kmer_list:
        for kmer in item:
            if kmer in kmer_count_dict:
                kmer_count_dict[kmer] += 1
            else:
                kmer_count_dict[kmer] = 1


    ## end your code
    ################

    return kmer_count_dict


def main():

    progname = sys.argv[0]

    usage = "\n\n\tusage: {} filename.fastq kmer_length num_top_kmers_show\n\n\n".format(
        progname
    )

    if len(sys.argv) < 4:
        sys.stderr.write(usage)
        sys.exit(1)

    # capture command-line arguments
    fastq_filename = sys.argv[1]
    kmer_length = int(sys.argv[2])
    num_top_kmers_show = int(sys.argv[3])

    seq_list = seq_list_from_fastq_file(fastq_filename)

    all_kmers = list()

    #######################
    ## Step 1:
    ## begin your code, populate 'all_kmers' list with the
    ## collection of kmers from all sequences

    for seq in seq_list:
        all_kmers.append(sequence_to_kmer_list(seq,kmer_length))
        #print(all_kmers)


    ## end your code
    #######################

    kmer_count_dict = count_kmers(
        all_kmers
    )  # see step 2 above. You implement this. :-)

    unique_kmers = list(kmer_count_dict.keys())

    #########################
    ## Step 3: sort unique_kmers by abundance descendingly
    ## (Note, you can run and test without first implementing Step 3)
    ## begin your code       hint: see the built-in 'sorted' method documentation

    kmer_count_dict = dict(sorted(kmer_count_dict.items(), key=lambda item: item[1],reverse=True))
    unique_kmers = list(kmer_count_dict.keys())

    ## end your code

    #########################
    ## calculate entropy for each kmer and add to dict
    import math
    entropys = {}
    #n = 0
    for kmer in kmer_count_dict.keys():
        nucleotides = {'A':0, 'C':0, 'G':0, 'T':0}
        for nucleotide in kmer:
            if nucleotide in nucleotides:
                nucleotides[nucleotide] += 1
        entroA = nucleotides['A']/len(kmer)
        #print(nucleotides['A'])
        #print(entroA)
        entroC = nucleotides['C']/len(kmer)
        #print(nucleotides['C'])
        #print(entroC)
        entroG = nucleotides['G']/len(kmer)
        #print(nucleotides['G'])
        #print(entroG)
        entroT = nucleotides['T']/len(kmer)
        #print(nucleotides['T'])
        #print(entroT)
        lists = [entroA, entroC, entroG, entroT]
        i = 0
        while i < len(lists):
            #print(lists[i])
            if lists[i] != 0:
                lists[i] = abs(lists[i] * math.log2(lists[i]))
            else:
                lists[i] = 0
    
            #print("after log")
            #print(lists[i])
            i = i + 1

        #print(lists)
        entropys[kmer] = lists[0] + lists[1] + lists[2] + lists[3]
        #n = n + 1
        #if n == 2:
        #    break
    
    ## end entropy code


    ## printing the num top kmers to show
    top_kmers_show = unique_kmers[0:num_top_kmers_show]

    for kmer in top_kmers_show:
        print("{}: {}".format(kmer, kmer_count_dict[kmer]))
        print("{}: {}".format(kmer, entropys[kmer]))

    sys.exit(0)  # always good practice to indicate worked ok!


if __name__ == "__main__":
    main()
