#!/usr/bin/env python
"""
version 0.5
minor revisions:
 - added genbank option
 - reversed single vs grouped naming convention just to confuse people
 - probably ruined whatever good thing this had going by overcomplicating it.
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
    print('Usage: python gen_genomes.py text_file_with_accessions.txt /output/directory/ multiple_or_single fasta_or_gb')
    print('("single" will return single file; "multiple" returns  individually)')
    print('(if using "gb", will return genbank results as well)')
    sys.exit()


def get_args():
    if len(sys.argv) != 5 or (sys.argv[3] != 'multiple' and sys.argv[3] != 'single') or (sys.argv[4] != 'fasta' and sys.argv[4] != 'gb'):
        usage()
    else:
        return sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]


def parse_accession_list(pathtofile):
    accessions = []
    for line in open(pathtofile):
        li = line.strip()
        if not li.startswith("#"):
            accessions.append(line.rstrip())
    return(accessions)


def fetch_and_write_seqs(accessions, destination, outtype):
    if (outtype == 'single'):
        for i in accessions:
            print("fetching %s" % i)
            out_handle = open(str(destination.strip()+datetimetag+"get_genomes_result.fasta"), "a")
            sequence_handle = Entrez.efetch(db="nucleotide", id=i, rettype="fasta", retmode="text")
            for line in sequence_handle:
                out_handle.write(line)
            out_handle.close()
            sequence_handle.close()
    elif (outtype == 'multiple'):
        for i in accessions:
            print("fetching %s" % i)
            out_handle = open(str(destination.strip() + i + ".fasta"), "w")
            sequence_handle = Entrez.efetch(db="nucleotide", id=i, rettype="fasta", retmode="text")
            for line in sequence_handle:
                out_handle.write(line)
            out_handle.close()
            sequence_handle.close()
    else:
        raise ValueError("Check your third argument! It must be either 'multiple' or 'single'")
    return(out_handle.name)


def fetch_and_write_gb(accessions, destination, outtype):
    if (outtype == 'single'):
        for i in accessions:
            print("fetching %s" % i)
            out_handle = open(str(destination.strip()+datetimetag+"get_genomes_result.gb"), "a")
            sequence_handle = Entrez.efetch(db="nucleotide", id=i, rettype="gbwithparts",
                                            retmode="text")
            for line in sequence_handle:
                out_handle.write(line)
            out_handle.close()
            sequence_handle.close()
    elif (outtype == 'multiple'):
        for i in accessions:
            print("fetching %s" % i)
            out_handle = open(str(destination.strip() + i + ".gb"), "w")
            sequence_handle = Entrez.efetch(db="nucleotide", id=i, rettype="gbwithparts",
                                            retmode="text")
            for line in sequence_handle:
                out_handle.write(line)
            out_handle.close()
            sequence_handle.close()
    else:
        raise ValueError("Check your third argument! It must be either 'multiple' or 'single'")
    return(out_handle.name)


##########################################################################

if __name__ == '__main__':
    inputlist, outputdirectory, outtype, outfmt = get_args()
    if not os.path.isdir(outputdirectory):
        print("creating %s" % outputdirectory)
        os.mkdir(outputdirectory)
    accession_list = parse_accession_list(inputlist)
    outputpath = fetch_and_write_seqs(accession_list, outputdirectory, outtype)
    if outfmt == "gb":
            outputpath2 = fetch_and_write_gb(accession_list, outputdirectory, outtype)
