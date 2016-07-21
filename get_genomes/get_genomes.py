#!/usr/bin/env python
"""
version 0.4
minor revisions:
 - added options for individual file output
"""

from Bio import Entrez
import sys
import os
import time
## dd/mm/yyyy format
datetimetag = time.strftime("%Y%m%d%I%M%S")
#
Entrez.email = "AlfredTheDaring@gmail.com"


def usage():
    print('\nA script to fetch nucleotide sequences from NCBI when given a file containing NCBI accession numbers\n')
    print('Usage: python gen_genomes.py text_file_with_accessions.txt /output/directory/ grouped_or_single')
    print('("grouped" will return single fasta file; "single" returns fastas individually)')
    sys.exit()


def get_args():
    if len(sys.argv) != 4 or (sys.argv[3] != 'single' and sys.argv[3] != 'grouped'):
        usage()
    else:
        return sys.argv[1], sys.argv[2], sys.argv[3]


def parse_accession_list(pathtofile):
    accessions = []
    for line in open(pathtofile):
        li = line.strip()
        if not li.startswith("#"):
            accessions.append(line.rstrip())
    return(accessions)


def fetch_and_write_seqs(accessions, destination, outtype):
    if (outtype == 'grouped'):
        for i in accessions:
            print("fetching %s" % i)
            out_handle = open(str(destination.strip()+datetimetag+"get_genomes_result.fasta"), "a")
            sequence_handle = Entrez.efetch(db="nucleotide", id=i, rettype="fasta", retmode="text")
            for line in sequence_handle:
                out_handle.write(line)
            out_handle.close()
            sequence_handle.close()
    elif (outtype == 'single'):
        for i in accessions:
            print("fetching %s" % i)
            out_handle = open(str(destination.strip() + i + ".fasta"), "w")
            sequence_handle = Entrez.efetch(db="nucleotide", id=i, rettype="fasta", retmode="text")
            for line in sequence_handle:
                out_handle.write(line)
            out_handle.close()
            sequence_handle.close()
    else:
        raise ValueError("Check your third argument! It must be either 'single' or 'grouped'")
    return(out_handle.name)

##########################################################################

if __name__ == '__main__':
    inputlist, outputdirectory, outtype = get_args()
    if not os.path.isdir(outputdirectory):
        print("creating %s" % outputdirectory)
        os.mkdir(outputdirectory)
    accession_list = parse_accession_list(inputlist)
    outputpath = fetch_and_write_seqs(accession_list, outputdirectory, outtype)
