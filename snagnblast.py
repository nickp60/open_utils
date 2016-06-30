#!/usr/bin/env python
"""
version 0.5
Minor version changes:
 - fixed off by one
 -redid output naming
 
 Quick and dirty script to fetch genes from NCBI when given a file cintaining NCBI accession numbers, 
 blast them against a local database (from makeblastdb), and write out the results as a csv.

USAGE:
 $ python snagnblast.py accessions.txt_or_accessions.csv /BLAST/directory/ /output/directory/
"""
import os
#import sys
#import re
import datetime
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
from Bio.Blast.Applications import NcbitblastxCommandline
DEBUG=False


#%% 
#define inputs
remake_blast_db=False
nuc_flag=False
if DEBUG:
    genelist = os.path.expanduser("~/GitHub/FB/Ecoli_comparative_genomics/data/test_virgenes.csv")
    blastdb = os.path.expanduser("~/BLAST/env_Coli")
    output = os.path.expanduser("~/GitHub/FB/Ecoli_comparative_genomics/results/")  
    score_min = 70
    blasttype="tblastx"
else:
    parser = argparse.ArgumentParser(description="This script takes a list of gene accessions from either a text file or a csv, grabs the sequencs from NCBI, and proceeds to use either blastn or tblastx to detect the presence of the genes in a custom database")
    parser.add_argument("genelist", help="file containing gene accessions.  if delimited, use the headers in the example file as a template")
    parser.add_argument("blastdb", help="blastdb of interest")
    parser.add_argument("-o","--output", help="directory in which to place the output files")
    parser.add_argument("-s","--score_min", help="not currently used; will be used to determinine a scoring threshold")
    parser.add_argument("-t","--blast_type", help="blastn or tblastx")

    args = parser.parse_args()    

    genelist = args.genelist
    blastdb = args.blastdb
    blasttype=args.blast_type
    output   = args.output
    score_min = args.score_min
date=str(datetime.datetime.now().strftime('%Y%m%d'))
#%%  open accessions file, determine type, and parse
Entrez.email = "nickp60@gmail.com"
print("reading in gene list")
genes=open(genelist, "r")
if genes.name.endswith("txt"):
    genelist_type="txt"
    print("gene list is a text file")
    #accessions=genes.readlines()
    accessions=[]
    for line in genes:
        li=line.strip()
        if not li.startswith("#"):
            accessions.append(line.rstrip())
elif genes.name.endswith( "tsv"):  #if the input is tabular, accesions must be in the first column
    genelist_type="delim"
    print("gene list is a tab-delimited file")
    n=("accession","name","phenotype",	"function","genome",	"note","source")
    genedf= pd.read_csv(genes, sep="\t", names=n, index_col=False)
    accessions= genedf.iloc[0:, 0].tolist()
    accessions = [x for x in accessions if str(x) != 'nan']
elif genes.name.endswith("csv"):
    genelist_type="delim"
    print("gene list is a comma-deliminated file")
    n=("accession",	"name","phenotype",	"function",	"genome",	"note",	"source")
    genedf= pd.read_csv(genes, sep=",")
    accessions= genedf.iloc[0:, 0].tolist()
    accessions = [x for x in accessions if str(x) != 'nan']
else:
    print("Reading error")
#%% Grab sequences from NCBI, write out resulting fasta file
print("\n\nFetching %i accessions from NCBI" %len(accessions))
sequence_handle= Entrez.efetch(db="nucleotide", id=accessions, rettype="fasta")
seqs=SeqIO.parse(sequence_handle, "fasta")
with open(str(os.path.join(output, date)+"_sequences.fa"), "w") as fasta_output: 
    SeqIO.write(seqs, fasta_output, "fasta")
#%%
sequences_fasta=open(str(os.path.join(output, date)+"_sequences.fa"), "r")    
entrez_results = list(SeqIO.parse(sequences_fasta, "fasta"))
print("returned %i accessions from ncbi" %len(entrez_results))
if(len(accessions)!= len(entrez_results)):
    print("Warning! not all accessions were found!")
sequences_fasta.close()
#%% Blast query against target
def run_blastn():
    # build commandline call
    output_path_tab=str(os.path.join(output, date)+"_dcmegablast_results.tab")
    blast_cline = NcbiblastnCommandline(query=fasta_output.name, 
                                        db= blastdb,  evalue=10,
                                        outfmt=7, out=output_path_tab)
    add_params=" -num_threads 4 -max_target_seqs 2000 -task dc-megablast"
    blast_command=str(str(blast_cline)+add_params)
    print("Running blastn search...")
#    subprocess.Popen(blast_command, stdout=subprocess.PIPE,  shell=True).stdout.read()
    subprocess.call(blast_command, shell=True)
    return(output_path_tab)

def run_tblastx():
    # build commandline call
    output_path_tab=str(os.path.join(output, date)+"_tblastx_results.tab")
    blast_cline = NcbitblastxCommandline(query=fasta_output.name, 
                                        db= blastdb,  evalue=10,
                                        outfmt=7, out=output_path_tab)
    add_params=" -num_threads 4 -max_target_seqs 2000  -query_gencode 11 -db_gencode 11"
    blast_command=str(str(blast_cline)+add_params)
    print("Running tblastx search...")
#    subprocess.Popen(blast_command, stdout=subprocess.PIPE,  shell=True).stdout.read()
    subprocess.call(blast_command, shell=True)
    return(output_path_tab)

#%% Execute
if blasttype=="blastn":
    output_path_tab = run_blastn()
elif blasttype=="tblastx":
    output_path_tab = run_tblastx()
else:
    print("you need to use either blastn or tblastx, sorry!")
    

#%% parse output
print("cleaning up the csv output")
colnames=["query_id", "subject_id", "identity_perc", "alignment_length", "mismatches", "gap_opens", "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score"]
csv_results=pd.read_csv(open(output_path_tab), comment="#", sep="\t" , names=colnames)
#This regex will probably break things rather badly before too long...
# it looks for capital letter and numbers, dot, number, ie SHH11555JJ8.99 
csv_results["accession"]=csv_results.query_id.str.extract('(?P<accession>[A-Z _\d]*\.\d*)') 
#%% write out results with new headers or with new headers and merged metadat from accessions.tab
output_path_csv = str(os.path.splitext(output_path_tab)[0]+".csv")
if genelist_type=="delim":
    results_annotated=pd.merge(csv_results, genedf, how="left",  on="accession" )
    results_annotated.to_csv(open(output_path_csv, "w") )
else:
    csv_results.to_csv(open(output_path_csv, "w") )
