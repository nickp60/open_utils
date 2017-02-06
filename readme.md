# open_utils
### by Nicholas Waters

This is the repo for my tools with realtively broad applicability.  See individual subdirectories for detailed readme's for the individual applications.  Assuming I have gotten around to writing one. If not, check the header of the script for usage, etc.

NOTE: 20160629 I changed the directory structure slightly.  Sorry about that; I wanted to have individual directories for sample data, etc

### gbparseR
Use to parse genbank files into fasta, ffa, fna, faa, gff, and a handy csv

### illumina2metaphlan
convert 16s metagenomic output from illumina to a format that mimicks metaphlan for use with the rest of biobakery

### reannotate
given a genbank genome of interest and a .faa for reference, this performs reciprocal blasting to "reannotate" and provide orthology information

### snagnblast
given a list of gene accesions, this retrieves them and blasts against a local BLAST database using either blastn or tblastx.

### get_genomes
_DEPRECIATED_ - added to pyutilsnrw repo which is available to install via pypi. given a list of genome accessions, this retrieves them from NCBI and writes out a fna file

### dnaGrenade
generates pseudoreads from a fasta file of loci of interest

### clermontPCR
This is a simple search tool to check for the presence/absence of genes used as part of the 2013 Clermont Phylotyping scheme.
