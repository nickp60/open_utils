#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 11:44:26 2016
version 0.5
Minor version changes:
 - works for multiple contigs

TODO:

USAGE:
 $ python extract_region.py genome.fasta 3:111 > extractedregion.fasta
 $ python extract_region.py genome.fasta chrom1:3:111 > extractedregionfromcontig.fasta

"""
import argparse
import sys
import re
DEBUG = False
#%%
#define inputs
if DEBUG:
    input_seqs_path = "multientry.fasta"
    coords = "'NODE_2_length_777_cov_4.02:5:342'"
else:
    parser = argparse.ArgumentParser(description="Given a set of coords (and an optional chromosome or contig), extract a region of the genome.  It goes to standard out, which can easily be redirected to a file with '>')
    parser.add_argument("fasta_file", help="input genome")
    parser.add_argument("coords", help="start:end or chromosome:start:end")
    args = parser.parse_args()
    input_seqs_path = args.fasta_file
    coords = args.coords


def parse_coords(coords):
    coords = coords.replace("'","")
    if len(coords.split(":")) == 2:
        start = int(coords.split(":")[0]) - 1  # account here for 0index
        end = int(coords.split(":")[1])
        chrom = None
    elif len(coords.split(":")) == 3:
        chrom = coords.split(":")[0]
        start = int(coords.split(":")[1]) - 1  # account here for 0index
        end = int(coords.split(":")[2])
    else:
        raise ValueError("must have either chrom:start:end or start:end")
    return(chrom, start, end)


def get_chrom_of_interest(chrom, lines):
    startindex = None
    endindex = None
    maxlines = len(lines)
    while startindex is None:
        for index, line in enumerate(lines):
            if line[0] == ">":
                try:
                    if len(re.search(chrom, line).group(0)) > 0:
                        print("extracting from %s" % chrom)
                        startindex = index
                except AttributeError:
                    continue
            if index == (maxlines - 1) and startindex is None:
                raise ValueError("%s not found in any of the fasta headers!" % chrom)
    while endindex is None:
        for nline in range(startindex+1, maxlines):
            if lines[nline][0] == ">" or int(nline) == int(maxlines-1):
                endindex = nline
    if startindex is None or endindex is None:
        raise ValueError("%s was not found in any fasta headers!" % chrom)
    return(startindex, endindex)

#%%


def main():
    chrom, start, end = parse_coords(coords)
    with open(input_seqs_path, "r") as file_handle:
        lines = file_handle.readlines()
    if lines[0][0] != ">":
        raise ValueError("not a valid fasta!")
    if chrom is not None:
        chrom_start, chrom_end = get_chrom_of_interest(chrom, lines)
    else:
        chrom_start, chrom_end = 0, len(lines)

    seq = ""
    header = str(lines[0 + chrom_start].strip().split(" ")[0] + "_" +
                 str(int(start+1)) + ":" + str(end))
    for line in range(1, len(lines[chrom_start:chrom_end])):
        if len(seq) > end:  # only read in up to what you need
            break
        if lines[line][0] == ">":
            raise ValueError("only for use with single fastas")
        seq = str(seq + lines[line].strip())
    linewidth = 80  # standard
    subset = seq[start:end]
#    print(header)
    newstring = ""
    for i in range(0, len(subset)):
        if i >= linewidth and (i % linewidth == 0):
            newstring = str(newstring+"\n"+subset[i])
        else:
            newstring = str(newstring+subset[i])
#    print(newstring)
    sys.stdout.write(header)
    sys.stdout.write("\n")
    sys.stdout.write(newstring)
    sys.stdout.write("\n")  # to get rid of '%' at end shown in terminal

if __name__ == '__main__':
    main()
