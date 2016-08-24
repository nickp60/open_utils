# gapped_alignment_vis.R

I am rather intersted in situations resulting in semi-global (free end, gapped) alignments, and have written this for quick comparisons.

## Usage
Requires biostrings
```
Rscript gapped_alignment_vis.R query.fasta subject.fasta
```
## Details
This does a few things:

* if you have several subject files (ie, a contigs multifasta), Each subject will be searched against
* matches are found for both forward and reverse-complimented query; the one with the higher score is outputed
* for convenience, if the query is longer than the subject, they will be switched such that the longer of the two will be displayed on the top in the alignment
* score threshold is calculated as 1/2 of the length of the pattern.

## Output
This just prints the alignment to stdout.  so usedul, right?

## Example output
```
nicholas@nicholinux[nicholas] Rscript ~/GitHub/open_utils/gapped_allignment_vis/gapped_alignment_vis.R ~/GitHub/FB/Ecoli_comparative_genomics/scripts/riboSeed_pipeline/batch_coli/20160815_ribosomal_RNA_16S226621_228162_riboSnag.fasta ~/riboSeed_sandbox/seed2/results/SPAdes_20160815_ribosomal_RNA_16S226621_228162_riboSnag/contigs.fasta
[1] "NODE_1_length_473_cov_8.52941"
[1] "using the reverse compliment"
[1] 282.039
CTAATCCCATCTGGGCACATCCGATGGCAAGAGGCCCGAAGGTCCCCCTCTTTGGTCTTGCGACGTTATGCGGTATTAGC
                                                                                
--------------------------------------------------------------------------------

TACCGTTTCCAGTAGTTATCCCCCTCCATCAGGCAGTTTCCCAGACATTACTCACCCGTCCGCCACTCGTCAGCGAAACA
                                                                                
--------------------------------------------------------------------------------

GCAAGCTGTTTCCTGTTACCGTTCGACTTGCATGTGTTTAGGCCTGCCGCCAGCGTTCAATCTGAGCCATGATCAAACTC
                                     |||||||||||||||||||||||||||||||||||||||||||
-------------------------------------TTAGGCCTGCCGCCAGCGTTCAATCTGAGCCATGATCAAACTC

TTCAATTTAAAAGTTTGATGCTCAAAGAATTAAACTTCGTAATGAATTACGTGTTCACTCTTGAGACTTGGTATTCATTT
||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
TTCAATTTAAAAGTTTGATGCTCAAAGAATTAAACTTCGTAATGAATTACGTGTTCACTCTTGAGACTTGGTATTCATTT

TTCGTCTTGCGACGTTAAGAATCCGTATCTTCGAGTGCCCACACAGATTGTCTGATAAATTGTTAAAGAGCAGTGCCGCT
||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||   |  
TTCGTCTTGCGACGTTAAGAATCCGTATCTTCGAGTGCCCACACAGATTGTCTGATAAATTGTTAAAGAGCAGTTGTGAC

TCGCTTTT--TCTCAGCGGCGCGGGGTGT-GCATAATACGCTTTCCCGCTACAGAGTCAAGCATTTTTTAC-----TTTT
 ||  |||   ||||  | |||| ||||  | ||| ||||||||||   | ||||||||| | |   || |     ||||
GCGGCTTTCAGCTCACTGTCGCGAGGTGGCGTATATTACGCTTTCCTCTTTCAGAGTCAACCCTGAATTTCAGGATTTTT

CT------------------------------------------------------------------------------
||                                                                              
CTTCTTCAACCGAACCGACTGTTTGTGTGAAGTGATTCACATCCGCCGTGTCGATGGAGGCGCATTATAGGGAGTCGGCT


```