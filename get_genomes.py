#!/usr/bin/env python
"""
version 0.2
minor revisions:
 - added functionality to remove commented lines (WHOLE LINES, that is)
 
 quick and dirty script to fetch genomes from NCBI when given a file cintaining NCBI accession numbers
 usage:
    python gen_genomes.py text_file_with_accessions.txt /output/directory/ 

"""

from Bio import Entrez #, SeqIO
import sys
#import time
Entrez.email = "nickp60@gmail.com"

#filelist=open(sys.argv[1])
outputdirectory=sys.argv[2]
#accessions=filelist.readlines()
accessions=[]
for line in open(sys.argv[1]):
    li=line.strip()
    if not li.startswith("#"):
        accessions.append(line.rstrip())
for i in accessions:
    print("fetching %s" %i)
    out_handle=open(str(outputdirectory+i.strip()+".fasta"), "w")
    sequence_handle= Entrez.efetch(db="nucleotide", id=i, rettype="fasta", retmode="text")
    for line in  sequence_handle:
        out_handle.write(line)
    out_handle.close()
    sequence_handle.close()
