#!/usr/bin/env python

# quick and dirty script to fetch genes from NCBI when given a file cintaining NCBI accession numbers, blast them against a local database (from makeblastdb), and concatenate the results into a single csv  
# usage python snagnblast.py text_file_with_accessions.txt /BLAST/directory/ /output/directory/

import os
#import sys
import re
import datetime
import subprocess
import argparse
from Bio import SeqIO,Entrez
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
import numpy as np
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
DEBUG=True


#%% 
#define inputs
remake_blast_db=False
nuc_flag=False
if DEBUG:
    genelist = os.path.expanduser("~/GitHub/FB/Ecoli_comparative_genomics/data/test_virgenes.txt")
    blastdb = os.path.expanduser("~/BLAST/env_Coli")
    output = os.path.expanduser("~/GitHub/FB/GitHub/FB/Ecoli_comparative_genomics/")  
    score_min = 70
else:
    parser = argparse.ArgumentParser(description="heres where the gelp message goes")
    parser.add_argument("genelist", help="file")
    parser.add_argument("blastdb", help="target am genome to")
    parser.add_argument("-o","--output", help="fasta returned")
    parser.add_argument("-s","--score_min", help="T if you waabases")

    args = parser.parse_args()    

    genelist = args.genelist
    blastdb = args.blastdb
    output   = args.output
    score_min = args.score_min
#%%  Grab sequences
Entrez.email = "nickp60@gmail.com"

max_results=200


#out_handle=open(str(output.strip()+".fasta"), "w")
genes=open(genelist, "r")
accessions=genes.readlines()
sequence_handle= Entrez.efetch(db="nucleotide", id=accessions, rettype="fasta")
seqs=SeqIO.parse(sequence_handle, "fasta")
fasta_output= open("sequences.fa", "w") 
SeqIO.write(seqs, fasta_output, "fasta")
fasta_output.close()

#%% Blast query against target
 #path to blast database for target 
# build commandline call
blast_cline = NcbiblastnCommandline(query=fasta_output.name, 
                                    db= blastdb,  evalue=10,
                                    outfmt=7, out=str(os.path.dirname(blastdb)+os.path.sep+"results.tab"))
#%% last I checked, I couldnt figure out how to add this through the python blast api                                    
add_params=" -num_threads 4 -max_target_seqs 2000 -task dc-megablast"
#print(str(blast_cline)+add_params)
blast_command=str(str(blast_cline)+add_params)
print("Running BLAST search...")
subprocess.Popen(blast_command, stdout=subprocess.PIPE,  shell=True).stdout.read()

#%%
colnames=["query id", "subject id", "% identity", "alignment length", "mismatches", "gap opens", "q. start", "q. end", "s. start", "s. end", "evalue", "bit score"]
csv_results=pd.read_csv(open(str(os.path.dirname(blastdb)+os.path.sep+"results.tab")), comment="#", sep="\t" , names=colnames  )
