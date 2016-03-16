#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 24 13:39:09 2015
######################## Genome Re-Annotation Tool ##########################
#####################         By Nick Waters 2015        ######################
##################               Brinsmade Lab                #################
@author: Nick Waters

#                             Version 4.3
#                               20160315
#
#   About:  this script is essentially a wrapper for NCBI BLAST+, the 
#           standalone version. Give it a genome(.gb), and a fasta file to blast against
#   Major revisions for 4.0:
#   -incorporated reciprocal blas
#   -argparse
#   -better debug tools
#   -fixed date using strtime
#   -fixed the hack to recast malformed dataframe (index, record iterations, etc)
#
#   Minor version revisions:
#   -fixed error bug looking for "i" instead of "locus"
#
# 
#   Requires: -installation of BioPython and NCBI Blast+ standalone suite
#             -make a directory called "BLAST" in your home folder
#             
#   Questions? nickp60@gmail.com
#             
#   USAGE: python reannotate.py input.gb target.fasta
#
#   OUTPUT: this will return a csv with the genome reannotations for used with 
#           downstream R applications.  Additionally, it will output the BLAST
#           XML for use with anything else.
#           Additionally, it will spit out a protein (or DNA) fasta with your 
#           genome
#   Known Bugs:  Might crash and burn if attempted on a PC.
                But who knows? I tried to use the os.path stuff for handling 
                path names...
              
"""
import os
#import sys
import re
import datetime
import subprocess
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
import numpy as np
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline
DEBUG=False
#%% 
#define inputs
remake_blast_db=False
nuc_flag=True
if DEBUG:
    input_genome = os.path.expanduser("~/GitHub/BlastDBs/uams1.gb")
    input_target_fasta = os.path.expanduser("~/GitHub/BlastDBs/CP000253_8325.fasta")
    nuc_flag=True
    remake_blast_db=True
    pattern='(.*gene=)(.*?)](.*)'
else:
    parser = argparse.ArgumentParser(description='Reciprocal blast for crude genome reannotation. \
    Requires a .gb file as the query and a protein fasta as the target for reannotation')
    parser.add_argument("query_db", help="path to .gb or .gbk Genbank genome file")
    parser.add_argument("target_db", help="target amino acid fasta file to compare your query genome to")
    parser.add_argument("-n","--nucleotide", help="T if you need a nucleotide fasta returned", type=bool)
    parser.add_argument("-r","--remake", help="T if you need a nucleotide fasta returned", type=bool)
    parser.add_argument("-g","--grep_pattern", default='(.*gene=)(.*?)](.*)', help="grep pattern for isolating the locus_tag in the target fasta", type=str)
    args = parser.parse_args()    
    input_genome = args.query_db
    input_target_fasta = args.target_db
    nuc_flag= args.nucleotide
    pattern=args.grep_pattern
    remake_blast_db= args.remake

#%%
#print(str("input_target_fasta is "+ input_target_fasta)) # sanity check 
gb=os.path.basename(input_genome).partition(".")[0]
#print(str("gb is "+ gb)) # sanity check 
input_target=os.path.basename(input_target_fasta).partition(".")[0]
#print(str("input_target is "+ input_target)) # sanity check 
bdb_name= str(input_target+"_db") #name for blast database to be created
#print(str("bdb_name is "+ bdb_name)) # sanity check 

now=datetime.datetime.now()
home=os.path.expanduser("~")
subdirname=os.path.join( home,"%s_reannotate_output" %str(now.strftime("%Y%m%d")))
if os.path.exists(subdirname): pass 
else:
    print("Making directory %s..." % subdirname)
    os.makedirs(subdirname)
input_handle = open(input_genome,"r") #  input database
aaoutput_handle = open(os.path.join(subdirname, gb+"aa.fasta"),"w") #  #  Amino Acid output
aaoutput_path=os.path.join(subdirname, gb+"aa.fasta")  #  Amino Acid output path
# for constructing a blastDB
output_genfasta_handle = open(os.path.join(subdirname, gb+"_genomic.fasta"), "w")
#test_path= "../testFasta.fasta"
#test_handle=open(test_path, "w")
blast_path= os.path.join(subdirname, gb+"_vs_"+bdb_name+".xml")   #  Blast output path xml
blast_path_recip= os.path.join(subdirname, bdb_name+"_vs_"+gb+".xml")   #  Blast output path xml
blastcsv_path=os.path.join(subdirname,gb+"_vs_"+bdb_name+".csv")#  Blast output path csv
finaldbcsv_path=os.path.join(subdirname,gb+"_with_"+bdb_name+"anno.csv") #  Final output 
    
#%%     make dataframe from genbank file   takes ~25" 
#       ADJUST TO CLEAN UP QUALIFIERS IF NEEDED
    #  See the bioPython manual to add more fields here if desired
genbankdf=pd.DataFrame()  ## genbank List Long
print(str("Reading .gb file..."))
#with input_handle as input_genome_handle:
#    record = next(SeqIO.parse(input_genome_handle, "genbank"))

#TODO: make this less awful
genbankdf=pd.DataFrame()  ## genbank List Long
for record in SeqIO.parse(input_handle, "genbank"):
    for index, feature in enumerate(record.features):
        if feature.type != "CDS":
            continue
        try:
            if feature.qualifiers.has_key("old_locus_tag"):
                genbankdf.loc[index, "old_locus_tag"]=str(feature.qualifiers["old_locus_tag"]).replace("\'", "").replace("[", "",).replace("]", "")
            else:
                genbankdf.loc[index, "old_locus_tag"]="none"
            genbankdf.loc[index, "locus_tag"]=str(feature.qualifiers["locus_tag"]).replace("\'", "").replace("[", "",).replace("]", "")
            genbankdf.loc[index, "product"]      =str(feature.qualifiers["product"]).replace("\'", "").replace("[", "",).replace("]", "")
            genbankdf.loc[index, "protein_id"]   =str(feature.qualifiers["protein_id"]).replace("\'", "").replace("[", "",).replace("]", "") 
            genbankdf.loc[index, "db_xref"]      =str(feature.qualifiers["db_xref"]).replace("\'", "").replace("[", "",).replace("]", "")
            genbankdf.loc[index, "seq"]          =str(feature.location.extract(record).seq)
            genbankdf.loc[index, "translation"]  =str(feature.qualifiers["translation"]).replace("\'", "").replace("[", "",).replace("]", "")        
            genbankdf.loc[index, "type"]="CDS"
        except KeyError:
            genbankdf.loc[index, "type"] = "pseudo" 
genbankdf.reset_index(level=0, inplace=True)
# set up recipient structures
#print(genbankdf.loc[genbankdf['locus_tag'] == 'QV15_00005'])   ### just a test      


if DEBUG:
    testdf=genbankdf[6:15]    ####Use this to debug creation of SeqRecord dictionary
##Set up recipient structures
nseqList=[]    
aaseqList=[]
 #%%  make that dataframe into a nucleotide fasta (ie, convert genbank to fasta)

def prepare_for_blastn(x):    
    for content in x.itertuples():
        a=Seq(str(content[7]), IUPAC.unambiguous_dna)
        b=SeqRecord(seq=a, id= str(content[1]),
                    description=str(content[3]+
                    "|"+str(content[2])+
                    '|'+str(content[4])+
                    '|'+str(content[5])))
        nseqList.append(b)

if nuc_flag=="T":
    noutput_handle=open(os.path.join(subdirname, gb+"n.fasta"),"w") #  Nucleotide output
    noutput_path = os.path.join(subdirname, gb+"n.fasta")   #  Nucleotide output
    print(str("Writing out nucleotide .fasta to\n"+ noutput_path+"..."))
    prepare_for_blastn(genbankdf)
    SeqIO.write(nseqList, noutput_handle, "fasta")
    noutput_handle.close()
else: pass
#%%  make that dataframe into a protein fasta  
print(str("Writing out protein .fasta to \n"+ aaoutput_path+"..."))
def prepare_for_blastaa(x):
    for content in x.itertuples():
        for locus_tag in content[3]:
            if locus_tag in aaseqList:
                print("warning!  Duplicate %s entry" %locus_tag)
                next
        a=Seq(str(content[8]), IUPAC.protein)
        b=SeqRecord(seq=a, id= str(content[1]),
                    description=str(content[3]+
                    "|"+str(content[2])+
                    '|'+str(content[4])+
                    '|'+str(content[5])))
        aaseqList.append(b)
prepare_for_blastaa(genbankdf)
SeqIO.write(aaseqList, aaoutput_handle, "fasta")

#  close open handles
input_handle.close()
aaoutput_handle.close()
#%%  Make a database from Fasta file using subprocess to pipe to shell
# if you have issues, print(make_db_command) and run the result in terminal
make_db_command_target= str("makeblastdb -in "+input_target_fasta+
    " -input_type fasta -dbtype prot -title "+
    input_target+
    " -parse_seqids -out "+  os.path.join(home, "BLAST", bdb_name))
# for reciprocal blasting 
make_db_command_query= str("makeblastdb -in "+aaoutput_handle.name+
    " -input_type fasta -dbtype prot -title "+
    gb+
    " -parse_seqids -out "+  os.path.join(home, "BLAST",  str(gb+"_db")))

#%%#debug command
#make_db_command='makeblastdb -in ~/testdb.fasta -input_type fasta -dbtype prot -title testdb -parse_seqids -out ~/BLAST/testdb'
# helpful thing:
#http://stackoverflow.com/questions/24340877/why-does-this-bash-call-from-python-not-work
#if os.path.isfile(os.path.dirname(str('~/BLAST/'+bdb_name+ '.psq'))):

#target
test_path_target= str( os.path.join(home, "BLAST", bdb_name+'.psq')) #
if os.path.exists(test_path_target) and remake_blast_db is False:
    print "BLAST Database Already Exists"
else:
    print(str("Creating BLAST database for "+bdb_name+"..."))
    subprocess.Popen(make_db_command_target, stdout=subprocess.PIPE,  shell=True).stdout.read()
#query
test_path_query= str( os.path.join(home, "BLAST", gb+"_db.psq")) #
if os.path.exists(test_path_query)and remake_blast_db is False:
    print "BLAST Database Already Exists"
else:
    print(str("Creating BLAST database for "+gb+"..."))
    subprocess.Popen(make_db_command_query, stdout=subprocess.PIPE,  shell=True).stdout.read()
#%% Add  to home directory to BLASTDB
#TODO test if already there
#test_if_in_BLASTDB='''echo "$BLASTDB" | grep -e '$HOME/BLAST' '''
#test_if_in_BLASTDB=str('''echo "$BLASTDB" | grep -e "var" ''')

#subprocess.Popen(test_if_in_BLASTDB, stdout=subprocess.PIPE,  shell=True).stdout.read()



#%blastdbset_command='export BLASTDB=$BLASTDB:$HOME/BLAST'#subprocess.Popen(blastdbset_command, stdout=subprocess.PIPE,  shell=True).stdout.read()
#%% TEST if worked
#subprocess.Popen("export", stdout=subprocess.PIPE, shell=True).stdout.read()

#%% Blast query against target
 #path to blast database for target 
bdb_path=os.path.join(os.path.abspath("BLAST"), bdb_name)
bdb_path_recip=os.path.join(os.path.abspath("BLAST"), gb)
# build commandline call
blast_cline = NcbiblastpCommandline(query=aaoutput_path, 
                                    db= bdb_path,  evalue=0.001,
                                    outfmt=5, out=blast_path)
blast_cline_recip = NcbiblastpCommandline(query=input_target_fasta, 
                                    db= str(bdb_path_recip+"_db"),  evalue=0.001,
                                    outfmt=5, out=blast_path_recip)
# last I checked, I couldnt figure out how to add this through the python blast api                                    
add_params=" -num_threads 2 -max_target_seqs 1 -max_hsps_per_subject 1"
#print(str(blast_cline)+add_params)
blast_command=str(str(blast_cline)+add_params)
blast_command_recip=str(str(blast_cline_recip)+add_params)
print("Running BLAST search...")
subprocess.Popen(blast_command, stdout=subprocess.PIPE,  shell=True).stdout.read()
print(str("Writing out XML to "+blast_path+"..."))
print("Running reciprocal BLAST search...")
subprocess.Popen(blast_command_recip, stdout=subprocess.PIPE,  shell=True).stdout.read()
print(str("Writing out XML to "+blast_path_recip+"..."))
#%%Parse ( from http://stackoverflow.com/questions/1684470/biopython-extracting-sequence-ids-from-a-blast-output-file)          
print("Parsing XML")
# for iteration in query against fasta, return hit and stat
record = NCBIXML.parse(open(blast_path))
record_recip = NCBIXML.parse(open(blast_path_recip))
def xml_to_df(record):
    df=pd.DataFrame(columns=["locus_tag", "match_in_target","bit_score","gaps","ortho"])
    counter=1    
    for rec in record:
        locus=str(rec.query.partition(" ")[2].partition("|")[0])
        if locus in df.locus_tag.values:
            print("Caution, possible duplicate %s;  locus_tag must be unique" %locus)
            line_index=counter
            df.loc[line_index, "locus_tag"]=locus
            counter=counter+1            
            continue
        else:
           for alignment in rec.alignments:
               hsp = alignment.hsps[0] #zero allows selection of first HSP only, as max_hsps_per_subject seems to fail 
               line_index=counter
               df.loc[line_index, "locus_tag"]=locus
               df.loc[line_index, "match_in_target"]=alignment.hit_def
               df.loc[line_index, "bit_score"]      =hsp.bits
               df.loc[line_index, "gaps"]           =hsp.gaps
               counter=counter+1            
    return(df)
querydf =  xml_to_df(record)
targetdf = xml_to_df(record_recip)
#%%

d2 = pd.DataFrame(targetdf.match_in_target.str.split("|").tolist(), columns=genbankdf.columns[2:6],index=np.arange(targetdf.shape[0])+1 )
#%%targetdf.match_in_target=d2.locus_tag
targetdf=targetdf.sort_values(by="match_in_target")
targetdf.reset_index(level=0, inplace=True)
try:
    querydf=querydf.sort_values(by="locus_tag")
except AttributeError:
    print("pandas is outdated; please upgrade\n")

#%%
#TODO: use bit score and gaps to create homologue filtering schema


#TODO pretest for uniqueness?
len(targetdf.match_in_target)!=len(set(targetdf.match_in_target))
#%% judge reciprocity
#    for j in targetdf.loc[1:len(targetdf)].index:
duplicates=list()
non_reciprocal_loci=list()
for i in querydf.locus_tag: #for each locus tag in query,
    index_in_querydf=querydf[querydf['locus_tag']==i].index.tolist() #get the index in the query df
    if len(index_in_querydf) is not 1:# address issues where something is clearly broken (0), or multiple cases found, also indicating something is broken.  If either, break
        if len(index_in_querydf) == 0:
           print(i+" not found in query?")
           break
        if len(index_in_querydf) >= 2:
           print("multiple instances of "+i+"found, using first; handling of multiple hsps is unavailible currently")
           continue
    querydf_vs_targetdf_match=str(querydf["match_in_target"][index_in_querydf].item()) #find match query vs target (forward search)
    index_in_targetdf=targetdf[targetdf["locus_tag"]==querydf_vs_targetdf_match].index.tolist()  #find index for match in target vs query df
    if len(index_in_targetdf)==0: # if locus tag is absent from target df (ie, non reciprocal result, where best matches arent eachother)
        print("non-reciprocal matching for %s" %querydf["locus_tag"][index_in_querydf].item())
        querydf["ortho"][index_in_querydf]="not_reciprocal"
        non_reciprocal_loci.append(i)
        continue
    if len(index_in_targetdf)>1: # if locus tag is present more than once in target df other
        print("multiple hits for %s" %querydf["locus_tag"][index_in_querydf].item())
    targetdf_vs_querydf_match=str(targetdf["locus_tag"][index_in_targetdf[0]]) #try first
    if querydf_vs_targetdf_match==targetdf_vs_querydf_match:
        querydf["ortho"][index_in_querydf]="reciprocal"
        continue
    elif len(index_in_targetdf)>1:
        print("checking multiple reciprocal hits for %s" %i)
        for i in range(1,len(index_in_targetdf)):
            targetdf_vs_querydf_match=str(targetdf["locus_tag"][index_in_targetdf[i]] )#try first
            if querydf_vs_targetdf_match==targetdf_vs_querydf_match:
                querydf["ortho"][index_in_querydf]="reciprocal"
            else:
                print("womp")
    else:
        print("womp")


#%%
#TODO make sure this handles all cases well

new=pd.merge(genbankdf, querydf, on='locus_tag')#, suffixes=['_left', '_right'])
#%% separate match_in_target with grep

#TODO can this be sped up?
if len(pattern)>0:
    genelist=list()
    for line in new.match_in_target:
        target_locus_tag=re.sub(pattern,r"\2", str(line))
        genelist.append(target_locus_tag)
    new['target_locus_tag']=genelist

#%% Write new database to output path
print("Writing results as CSV")
new.to_csv(finaldbcsv_path)
#TODO write out log file with params, etc
