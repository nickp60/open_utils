#!/bin/bash
mkdir test_data
cd test_data
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/019/385/GCA_000019385.1_ASM1938v1/GCA_000019385.1_ASM1938v1_genomic.gbff.gz -O ATCC8739.gbk.gz
cp ~/GitHub/soil-persistent-ecoli/Chapter-comparative-genomics/data/2020-01-clonality-core/core.vcf ./core.vcf
