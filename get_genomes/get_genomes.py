#!/usr/bin/env python
"""
version 0.3
minor revisions:
 - tried to make it testable
""" 

from Bio import Entrez
import sys
import os
import time
## dd/mm/yyyy format
datetimetag=time.strftime("%Y%m%d%I%M%S")
#
Entrez.email = "AlfredTheDaring@gmail.com"

def usage():
    print('\nA script to fetch nucleotide sequences from NCBI when given a\
           file containing NCBI accession numbers\n')
    print('Usage: python gen_genomes.py text_file_with_accessions.txt\
           /output/directory/')
    sys.exit()


def get_args():
	if len(sys.argv) != 3:
            usage()
	else:
            return sys.argv[1], sys.argv[2]


def parse_accession_list(pathtofile):
    accessions=[]
    for line in open(pathtofile):
        li=line.strip()
        if not li.startswith("#"):
            accessions.append(line.rstrip())
    return(accessions)


def fetch_and_write_seqs(accessions, destination):
    for i in accessions:
        print("fetching %s" %i)
        out_handle=open(str(destination.strip()+datetimetag+"get_genomes_result.fasta"), "a")
        sequence_handle= Entrez.efetch(db="nucleotide", id=i, rettype="fasta", retmode="text")
        for line in  sequence_handle:
            out_handle.write(line)
        out_handle.close()
        sequence_handle.close()
    return(out_handle.name)

##########################################################################

if __name__ == '__main__':
    inputlist, outputdirectory = get_args()
    if not os.path.isdir(outputdirectory):
        print("creating %s" % outputdirectory)
        os.mkdir(outputdirectory)
    accession_list = parse_accession_list(inputlist)
    outputpath = fetch_and_write_seqs(accession_list, outputdirectory)
