# open_utils
### by Nicholas Waters

This is the repo for my tools with realtively broad applicability.  See individual subdirectories for detailed radme's for the individual applications.  Assuming I have gotten around to writing one. If not, check the header of the script for usage, etc.

NOTE: 20160629 I changed the directory structure slightly.  Sorry about that; I wanted to have individual directories for sample data, etc

### gbparseR
Use to parse genbank files into fasta, ffa, fna, faa, gff, and a handy csv

### illumina2metaphlan
convert 16s metagenomic output from illumina to a format that mimicks metaphlan for use with the rest of biobakery

### reannotate
given a genbank genome of interest and a .faa for reference, this performs reciprocal blasting to "reannote" and provide othology information

### snagnblast
given a list of gene accesions, this retrieves them and blasts against a local BLAST database using either blastn or tblastx. Snagnblast multi is for use with multiple processors, but still needs a bit of a tuneup to pool jobs instead of queueing.

### get_genomes
given a list of genome accessions, this retrieves them from NCBI and writes out a fna file

### dnaGrenade
generates pseudoreads from a fasta file of loci of interest

