#!/usr/bin/env python
"""
version 0.3
Minor version changes:
 -now merges accession file with output to transfer over any other metadata and writes back out
 - text input should ignore hash-commented lines
 
 
 Quick and dirty script to fetch genes from NCBI when given a file cintaining NCBI accession numbers, 
 blast them against a local database (from makeblastdb), and write out the results as a csv.

USAGE:
 $ python snagnblast.py accessions.txt_or_accessions.csv /BLAST/directory/ /output/directory/
"""
import os
#import sys
#import re
#import datetime
import subprocess
import argparse
from Bio import SeqIO,Entrez
#from Bio.SeqRecord import SeqRecord
#from Bio.Seq import Seq
import pandas as pd
#import numpy as np
#from Bio.Alphabet import IUPAC
#from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
DEBUG=True


#%% 
#define inputs
remake_blast_db=False
nuc_flag=False
if DEBUG:
    genelist = os.path.expanduser("~/GitHub/FB/Ecoli_comparative_genomics/data/test_virgenes.tsv")
    blastdb = os.path.expanduser("~/BLAST/env_Coli")
    output = os.path.expanduser("~/GitHub/FB/Ecoli_comparative_genomics/")  
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
#%%  open accessions file, determine type, and parse
Entrez.email = "nickp60@gmail.com"

genes=open(genelist, "r")
if genes.name.endswith("txt"):
    genelist_type="txt"
    #accessions=genes.readlines()
    accessions=[]
    for line in genes:
        li=line.strip()
        if not li.startswith("#"):
            accessions.append(line.rstrip())
elif genes.name.endswith( "tsv"):  #if the input is tabular, accesions must be in the first column
    genelist_type="delim"
    n=("accession","name","phenotype",	"function","genome",	"note","source")
    genedf= pd.read_csv(genes, sep="\t", names=n, index_col=False)
    accessions= genedf.iloc[1:, 0].tolist()
    accessions = [x for x in accessions if str(x) != 'nan']
elif genes.name.endswith("csv"):
    genelist_type="delim"
    n=("accession",	"name","phenotype",	"function",	"genome",	"note",	"source")
    genedf= pd.read_csv(genes, sep=",", names=n)
    accessions= genedf.iloc[1:, 0].tolist()
    accessions = [x for x in accessions if str(x) != 'nan']
else:
    print("REading error")
#%% Grab sequences from NCBI, write out resulting fasta file
sequence_handle= Entrez.efetch(db="nucleotide", id=accessions, rettype="fasta")
seqs=SeqIO.parse(sequence_handle, "fasta")
with open(str(output+"sequences.fa"), "w") as fasta_output: 
    SeqIO.write(seqs, fasta_output, "fasta")


#%% Blast query against target

# build commandline call
output_path_tab=str(output+"dcmegablast_results.tab")
output_path_csv=str(output+"dcmegablast_results.csv")
blast_cline = NcbiblastnCommandline(query=fasta_output.name, 
                                    db= blastdb,  evalue=10,
                                    outfmt=7, out=output_path_tab)
#%% last I checked, I couldnt figure out how to add this through the python blast api                                    
add_params=" -num_threads 4 -max_target_seqs 2000 -task dc-megablast"
#print(str(blast_cline)+add_params)
blast_command=str(str(blast_cline)+add_params)
print("Running BLAST search...")
subprocess.Popen(blast_command, stdout=subprocess.PIPE,  shell=True).stdout.read()

#%% parse output
colnames=["query_id", "subject_id", "identity_perc", "alignment_length", "mismatches", "gap_opens", "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score"]
csv_results=pd.read_csv(open(output_path_tab), comment="#", sep="\t" , names=colnames)
#This regex will probably break things rather badly before too long...
# it looks for capital letter and numbers, dot, number, ie SHH11555JJ8.99 
csv_results["accession"]=csv_results.query_id.str.extract('(?P<accession>[A-Z \d]*\.\d*)') 
#%% write out results with new headers or with new headers and merged metadat from accessions.tab
if genelist_type=="delim":
    results_annotated=pd.merge(csv_results, genedf, how="left",  on="accession" )
    results_annotated.to_csv(open(output_path_csv, "w") )
else:
    csv_results.to_csv(open(output_path_csv, "w") )
