#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 24 13:39:09 2015
######################## Database Re-Annotation Tool ##########################
#####################         By Nick Waters 2015         #####################
##################               Brinsmade Lab                #################
@author: Nick Waters

#                             Version 3.2
#                               20151217
#
#   About:  this script is essentially a wrapper for NCBI BLAST+, the 
#           standalone version. Give it a genome(.gb), and a fasta file to blast against
#   Minor version revisions:  previos version had to be run with the nulceotide ouput flag being true. Oops...
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
#           Additionally, it will spit out a protein (or DNA) fasta wuth your 
#           genome
#   Known Bugs:  Might crash and burn if attempted on a PC.
                But who knows? I tried to use the os.path stuff for handling 
                path names...
              
"""
import os
import sys
import datetime
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
from Bio.Alphabet import IUPAC
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastpCommandline
#%%  Get sys args from command line (^see usage^)
input_genome = sys.argv[1]
#print(str("input_genome is "+ input_genome))   # sanity check 
input_target_fasta = sys.argv[2]
try:nuc_flag=sys.argv[3] 
except: nuc_flag='F'
#%%
#print(str("input_target_fasta is "+ input_target_fasta)) # sanity check 
#kegg_org="sar"
#gb= input_genome.split("/")[len(input_genome.split("/"))-1].partition(".")[0]
gb=os.path.basename(input_genome).partition(".")[0]
#input_genome.partition(".")[0]  # removes ".gb"
#print(str("gb is "+ gb)) # sanity check 
input_target=os.path.basename(input_target_fasta).partition(".")[0]
#print(str("input_target is "+ input_target)) # sanity check 
bdb_name= str(input_target+"_db") #name for blast database to be created
#print(str("bdb_name is "+ bdb_name)) # sanity check 

#%% debug set; adjust as needed to twea
'''
input_genome = os.path.expanduser("~/GitHub/Py/uams1.gb")
input_target_fasta = os.path.expanduser("~/Downloads/subtilis168.txt")
nuc_flag='T'
'''
#%%
now=datetime.datetime.now()
home=os.path.expanduser("~")
subdirname=os.path.join( home,"%s%s%s_reannotate_output" %(now.year,now.month, now.day))
if os.path.exists(subdirname): pass 
else:
    print("Making directory %s..." % subdirname)
    os.makedirs(subdirname)
input_handle = open(input_genome,"r") #  input database
aaoutput_handle = open(os.path.join(subdirname, gb+"aa.fasta"),"w") #  #  Amino Acid output
aaoutput_path=os.path.join(subdirname, gb+"aa.fasta")  #  Amino Acid output path
#test_path= "../testFasta.fasta"
#test_handle=open(test_path, "w")
blast_path= os.path.join(subdirname, gb+"_vs_"+bdb_name+".xml")   #  Blast output path xml
blastcsv_path=os.path.join(subdirname,gb+"_vs_"+bdb_name+".csv")#  Blast output path csv
finaldbcsv_path=os.path.join(subdirname,gb+"_with_"+bdb_name+"anno.csv") #  Final output 
    
#%%     make dataframe from genbank file   takes ~25" 
#       ADJUST TO CLEAN UP QUALIFIERS IF NEEDED
    #  See the bioPython manual to add more fields here if desired
gll=pd.DataFrame()  ## genbank List Long
print(str("Reading .gb file..."))
for rec in SeqIO.parse(input_handle, "genbank"):
    if rec.features:
        for feature in rec.features:
            if feature.type == "CDS":
                for qual in feature.qualifiers["locus_tag"]:
                    try:
                        gll[qual]=(str(feature.qualifiers["old_locus_tag"]).replace("\'", "").replace("[", "",).replace("]", ""),
                        str(feature.qualifiers["product"]).replace("\'", "").replace("[", "",).replace("]", ""),
                        str(feature.qualifiers["protein_id"]).replace("\'", "").replace("[", "",).replace("]", ""), 
                        str(feature.qualifiers["db_xref"]).replace("\'", "").replace("[", "",).replace("]", ""), 
                        str(feature.location.extract(rec).seq),
                        str(feature.qualifiers["translation"]).replace("\'", "").replace("[", "",).replace("]", "") )              
                    except KeyError:
                        gll[qual] = "pseudo" 
#%%      shake it out from wide to tall data, sanity check, and 
#       set up recipient structures
gl=gll.transpose()  # genbank List
gl.reset_index(level=0, inplace=True)  # make numeric index
gl.columns=('locus_tag', 'old_locus_tag','product', 'protein_id', 'db_xref','sequence', "aaseq") # name dem cols
    
#print(gl.loc[gl['locus_tag'] == 'QV15_00005'])   ### just a test      


#testdf=gl[6:15]    ####Use this to debug creation of SeqRecord dictionary
##Set up recipient structures
nseqDict={}
nseqList=[]    
aaseqDict={}
aaseqList=[]
 #%%  make that dataframe into a nucleotide fasta (ie, convert genbank to fasta)

def preparefor_blastn(x):    
    for content in x.itertuples():
        a=Seq(str(content[5]), IUPAC.unambiguous_dna)
#### clean up, clean up, all the brackets everywhere    
        b=SeqRecord(seq=a, id= content[1],
                    description=str(content[2]).replace("\'", "").replace("[", "",).replace("]", "")+
                    '|'+str(content[3]).replace("\'", "").replace("[", "",).replace("]", ""))
        nseqList.append(b)
#%%
if nuc_flag=="T":
    noutput_handle=open(os.path.join(subdirname, gb+"n.fasta"),"w") #  Nucleotide output
    noutput_path = os.path.join(subdirname, gb+"n.fasta")   #  Nucleotide output
    print(str("Writing out nucleotide .fasta to\n"+ noutput_path+"..."))
    preparefor_blastn(gl)
    SeqIO.write(nseqList, noutput_handle, "fasta")
    noutput_handle.close()
else: pass
''
#%%  make that dataframe into a protein fasta  
print(str("Writing out protein .fasta to \n"+ aaoutput_path+"..."))
def preparefor_blastaa(x):
    for index, row in x.iterrows():
        for locus_tag in row.locus_tag:
            if locus_tag in aaseqList:
                print("warning!  Duplicate %s entry" %locus_tag)
        a=Seq(str(row['aaseq']), IUPAC.protein)
        b=SeqRecord(seq=a, id= row['locus_tag'],
                    description=row['protein_id']+
                    "|"+row['old_locus_tag']+
                    '|'+row['product'])
        aaseqDict[index]=b ### just for funsies, in case you want a dictionary
        aaseqList.append(b)
preparefor_blastaa(gl)
SeqIO.write(aaseqList, aaoutput_handle, "fasta")    #####USE WITH CAUTION- BIG FILE

 #%%  close open handles
input_handle.close()
aaoutput_handle.close()

#%%  Make a database from Fasta file using subprocess to pipe to shell
# if you have issues, print(make_db_command) and run the result in terminal
make_db_command= str("makeblastdb -in "+input_target_fasta+
    " -input_type fasta -dbtype prot -title "+
    input_target+
    " -parse_seqids -out "+  os.path.join(home, "BLAST", bdb_name))
#debug command
#make_db_command='makeblastdb -in ~/testdb.fasta -input_type fasta -dbtype prot -title testdb -parse_seqids -out ~/BLAST/testdb'
# helpful thing:
#http://stackoverflow.com/questions/24340877/why-does-this-bash-call-from-python-not-work
#if os.path.isfile(os.path.dirname(str('~/BLAST/'+bdb_name+ '.psq'))):
test_path= str( os.path.join(home, "BLAST", bdb_name+'.psq')) #
if os.path.exists(test_path):
    print "BLAST Database Already Exists"
else:
    print(str("Creating BLAST database for "+bdb_name+"..."))
    subprocess.Popen(make_db_command, stdout=subprocess.PIPE,  shell=True).stdout.read()

     
#%%  BBBBBBLLLLLLLAAAAAAASSSSSSTTTTT!!!!!!!!!!!!!!!!!!
#Debug command
#blastp -out test.xml -outfmt 5 -query ./GitHub/BlastDBs/qv15_13255.fasta -db ~/BLAST/testdb -evalue 0.001 -num_threads 2 -max_target_seqs 1
bdb_path=os.path.join(os.path.abspath("BLAST"), bdb_name)
blast_cline = NcbiblastpCommandline(query=aaoutput_path, 
                                    db= bdb_path,  evalue=0.001,
                                    outfmt=5, out=blast_path)
add_params=" -num_threads 2 -max_target_seqs 1"
#print(str("\t ++++++copy and paste the 2 lines below into your terminal to run BLAST: ++++++").upper())
#print("cd ~/GitHub/BlastDBs/")
#print(str(blast_cline)+add_params)
blast_command=str(str(blast_cline)+add_params)
print("Running BLAST search...")
subprocess.Popen(blast_command, stdout=subprocess.PIPE,  shell=True).stdout.read()
print(str("Writing out XML to "+blast_path+"..."))

#stdout, stderr = blast_cline()


#%%  make a recipient dictionary ("not found" will be overwritten if match exists)
newdb={}
for i in gl.locus_tag:
    newdb[i]="not found"
#%%Parse ( from http://stackoverflow.com/questions/1684470/biopython-extracting-sequence-ids-from-a-blast-output-file)          
print("Parsing XML")
for record in NCBIXML.parse(open(blast_path)): #print record
    querykey=record.query.partition(" ")[0]
    for align in record.alignments :
        if querykey in newdb:
            newdb[querykey]= (align.title)

#%% make dataframe, match, amd merge over the "qv_anno" column
newdf=pd.DataFrame(newdb.items(), columns=['locus_tag', str(bdb_name+'_anno')]).sort(columns='locus_tag', axis=0, ascending=True)
locus_tag=gl
new=pd.merge(newdf, locus_tag, on='locus_tag')#, suffixes=['_left', '_right'])

#%% Write new database to output path
print("Writing results as CSV")
new.to_csv(finaldbcsv_path)

