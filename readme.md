# open_utils

This is the repo for my tools with realtively broad applicability.  See individual subdirectories for detailed readme's for the individual applications.  Assuming I have gotten around to writing one. If not, check the header of the script or run with `-h` for usage, etc.

## Installation
`open_utils` can be installed from pip or by running `python setup.py install` from the cloned repository.  That will install entrypoints for the following programs:

### vcfortless
Given a genbank annotation file and a vcf file, generates a nice flatfile to tell Ts/Tv, which SNPs are intergenic, etc

### ExtractRegion
A relatively fast way to extract regions of interest from a fasta file.  Nice and pipe-able.

### snagnblast
given a list of gene accesions, this retrieves them and blasts against a local BLAST database using either blastn or tblastx.

## Other tools
The following other tools are available in the repository, but are not installed as part of the package by default.

### gbparseR
Use to parse genbank files into fasta, ffa, fna, faa, gff, and a handy csv

### illumina2metaphlan
convert 16s metagenomic output from illumina to a format that mimicks metaphlan for use with the rest of biobakery

### reannotate
given a genbank genome of interest and a .faa for reference, this performs reciprocal blasting to "reannotate" and provide orthology information

### get_genomes
_DEPRECIATED_ - added to pyutilsnrw repo which is available to install via pypi. given a list of genome accessions, this retrieves them from NCBI and writes out a fna file

### dnaGrenade
generates pseudoreads from a fasta file of loci of interest

### clermontPCR (renamed to EzClermont and moved to own repo)
This is a simple search tool to check for the presence/absence of genes used as part of the 2013 Clermont Phylotyping scheme.
