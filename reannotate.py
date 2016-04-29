#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 24 13:39:09 2015
######################## Genome Re-Annotation Tool ##########################
#####################         By Nick Waters 2015        ######################
##################               Brinsmade Lab                #################
@author: Nick Waters

#                             Version 4.7
#                               20160413
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
#   - created option to rename genomic fasta headers
#TODO
#   make renamer work with subprocess; had issues with using awk from within python
# 
#   Requires: -installation of BioPython and NCBI Blast+ standalone suite
#             -make a directory called "BLAST" in your home folder
#             
#   Questions? nickp60@gmail.com
#             
#   USAGE: python reannotate.py input.gb target.fasta new_grep_pattern*
# *optional
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
nuc_flag=False
if DEBUG:
    input_genome = os.path.expanduser("~/GitHub/BlastDBs/genbank_genomes/uams1.gb")
    input_target_fasta = os.path.expanduser("~/GitHub/BlastDBs/fasta_genomes/MRSA252.txt")
    nuc_flag=True
    remake_blast_db=True
    pattern='(.*gene=)(.*?)](.*)'
    rename_fa=True
else:
    parser = argparse.ArgumentParser(description='Reciprocal blast for crude genome reannotation. \
    Requires a .gb file as the query and a protein fasta as the target for reannotation')
    parser.add_argument("query_db", help="path to .gb or .gbk Genbank genome file")
    parser.add_argument("target_db", help="target amino acid fasta file to compare your query genome to")
    parser.add_argument("-n","--nucleotide", help="T if you need a nucleotide fasta returned")
    parser.add_argument("-r","--remake", help="T if you want to remake the blast databases")
    parser.add_argument("-g","--grep_pattern", default='(.*gene=)(.*?)](.*)', help="grep pattern for isolating the locus_tag in the target fasta", type=str)
    parser.add_argument("-f","--rename_fa", help="T  if you would like to try to simplify the names in the resulting genomic fasta file")
    args = parser.parse_args()    
    input_genome = args.query_db
    input_target_fasta = args.target_db
    if str(args.nucleotide).lower() =="t" or str(args.nucleotide).lower() =="true" :
        nuc_flag=True
    else:
        nuc_flag=False
    if str(args.rename_fa).lower() =="t" or str(args.nucleotide).lower() =="true" :
        rename_fa=True
    else:
        rename_fa=False
    pattern=args.grep_pattern
    if str(args.remake).lower() =="t" or str(args.remake).lower() =="true" :
        remake_blast_db=True
    else:
        remake_blast_db=False


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
aaoutput_handle = open(os.path.join(subdirname, gb+".fa"),"w") #  #  Amino Acid output
aaoutput_path=os.path.join(subdirname, gb+".fa")  #  Amino Acid output path
# for constructing a blastDB
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
qualifiers=["old_locus_tag","locus_tag","protein_id", "db_xref", "translation", "gene"]
index=0
for record in SeqIO.parse(input_handle, "genbank"):
    for feature in record.features:
        index=index+1
        if feature.type != "CDS":
            continue
        try:
            for qualifier in qualifiers:
                try:
                    genbankdf.loc[index,qualifier]=str(feature.qualifiers[qualifier]).replace("\'", "").replace("[", "",).replace("]", "")
                except KeyError:
                     genbankdf.loc[index,qualifier]='none'
                genbankdf.loc[index, "seq"] =str(feature.location.extract(record).seq)
                genbankdf.loc[index, "type"]="CDS"
        except KeyError:
            genbankdf.loc[index, "type"] = "pseudo" 

#%%
#this should handle casses (like MW2) where there is no locus_tag designation
#ie, this should rarely be called 
if ((genbankdf.locus_tag=='none').all() & (genbankdf.gene=="none").all()):
    invalid_genbank=KeyError("Genbank file must have at either unique locus_tag column or gene column")
    raise invalid_genbank
elif ((genbankdf.locus_tag=='none').all() & (genbankdf.gene!="none").all()):
    genbankdf["locus_tag"]=genbankdf["gene"]
else:
    pass
if (genbankdf.old_locus_tag=='none').all():
    genbankdf["old_locus_tag"]=genbankdf["locus_tag"]

if (genbankdf.locus_tag=='none').any():
    print("Warning! removing %i entries lacking locus_tag.." %len(genbankdf[genbankdf.locus_tag=='none']))
    genbankdf=genbankdf[genbankdf.locus_tag!="none"]

genbankdf.reset_index(level=0, inplace=True)
#  we split with pipes "|" later on;  this rmoves any from locus_tags , bd_xref's, and old_locus_tags
try:
    genbankdf.db_xref=genbankdf.db_xref.str.replace("|","_")
    genbankdf.locus_tag=genbankdf.locus_tag.str.replace("|","_")
    genbankdf.old_locus_tag=genbankdf.old_locus_tag.str.replace("|","_")
except KeyError:
    pass
#print(genbankdf.loc[genbankdf['locus_tag'] == 'QV15_00005'])   ### just a test      
#%%#%%#  make a fasta for genomic sequence(s), clean up header in resulting file
output_genfasta_handle = open(os.path.join(subdirname, gb+"_genomic.fasta"), "w")
SeqIO.convert(input_handle.name, "genbank", output_genfasta_handle, "fasta")
output_genfasta_handle.close()
#%%

if rename_fa:
    print("After running this script, run the following output as a command, starting with 'awk' and ending with _renamed.fasta:")
    import string
    input_handle = open(input_genome,"r") #  input database
    names=[]
    for record in SeqIO.parse(input_handle, "genbank"):
        print(record.id)
        names.append(record.id)
    if len(names)>1:
        starts_at=names[0][-1]
        base_name=names[0].partition(".")[0][0:-1] # removes number, period, and version number
        length_of_names=len(names)
        version_of_genbank=re.sub("(.*)\.(.*)$",r"\2", names[0])
        if not string.ascii_letters.find(starts_at)==-1:
            name_not_digit= KeyError("not actally a key error, but at this point, fastas can only be renamed if the sequence ends in a digit :(")
            raise name_not_digit
        if not string.ascii_letters.find(starts_at)==-1:
            version_not_digit= KeyError("not actally a key error, but at this point, fastas can only be renamed if the version ends in a digit :(")
            raise version_not_digit
        if not int(starts_at)==1:
            version_not_digit= KeyError("not actally a key error, but at this point, fastas can only be renamed if the version starts with 1 :(")
            raise version_not_digit
        renaming_script=str(r"""awk '/^>/{print ">"""+base_name+"""" ++i "."""+version_of_genbank+r""" "; next}{print}' < """+output_genfasta_handle.name +" > "+str(output_genfasta_handle.name.partition(".fasta")[0]+"_renamed.fasta"))
        print(renaming_script)
    elif len(names)==1:
        name= names[0]
        renaming_script=str(r"""awk '/^>/{print ">"""+name+""""; next}{print}' < """+output_genfasta_handle.name +" > "+str(output_genfasta_handle.name.partition(".fasta")[0]+"_renamed.fasta"))
        print(renaming_script)
    else:
        print("something funny happened when preparing the renaming script")
print("\n\n")

# TODO
#    proc=subprocess.Popen(renaming_script, shell=True)
#    proc.wait()
#    print proc.returncode
#%%
# set up recipient structures

if DEBUG:
    testdf=genbankdf[6:15]    ####Use this to debug creation of SeqRecord dictionary
##Set up recipient structures
nseqList=[]    
aaseqList=[]
 #%%  make that dataframe into a nucleotide fasta (ie, convert genbank to fasta)

def prepare_for_blastn(x):    
    for content in x.itertuples():
        a=Seq(str(content.seq), IUPAC.unambiguous_dna)
        b=SeqRecord(seq=a, id= str(content.Index), #id
                    description=str(content.locus_tag+ 
                    "|"+str(content.old_locus_tag)+
                    '|'+str(content.db_xref)))
        nseqList.append(b)
if nuc_flag==True:
    noutput_handle=open(os.path.join(subdirname, gb+".fn"),"w") #  Nucleotide output
    noutput_path = os.path.join(subdirname, gb+".fn")   #  Nucleotide output
    print(str("Writing out nucleotide .fn to\n"+ noutput_path+"..."))
    prepare_for_blastn(genbankdf)
    SeqIO.write(nseqList, noutput_handle, "fasta")
    noutput_handle.close()
else: pass
#%%  make that dataframe into a protein fasta  
print(str("Writing out protein .fn to \n"+ aaoutput_path+"..."))
def prepare_for_blastaa(x):
    for content in x.itertuples():
        for locus_tag in content.locus_tag:
            if locus_tag in aaseqList:
                print("warning!  Duplicate %s entry" %locus_tag)
                next
        a=Seq(str(content.translation), IUPAC.protein)
        b=SeqRecord(seq=a, id= str(content.Index), #id
                    description=str(content.locus_tag+ 
                    "|"+str(content.old_locus_tag)+
                    '|'+str(content.db_xref)))
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
                                    db= bdb_path,  evalue=0.01,
                                    outfmt=5, out=blast_path)
blast_cline_recip = NcbiblastpCommandline(query=input_target_fasta, 
                                    db= str(bdb_path_recip+"_db"),  evalue=0.01,
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
#create dataframe with gene information
d2 = pd.DataFrame(targetdf.match_in_target.str.split("|").tolist(), 
                  columns=["locus_tag","old_locus_tag", "db_xref"],index=np.arange(targetdf.shape[0])+1 )
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
print("there were %i not found in query" %(abs(len(querydf)-len(genbankdf))))
new=pd.merge(genbankdf, querydf,how="left", on='locus_tag')#, suffixes=['_left', '_right'])
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
