#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 11:44:26 2016
version 0.4
Minor version changes:
 - fixed bug line 47: if lines[line] == ">" needed to be if lines[line][0] == ">":


TODO:

USAGE:
 $ python extract_region.py genome.fasta 3:111 > extractedregion.fasta
"""
import argparse
import sys
DEBUG = False
#%%
#define inputs
if DEBUG:
    input_seqs_path = "grab.fasta"
    coords = "5:342"
else:
    parser = argparse.ArgumentParser(description="Given a set of coords, extract a region of the genome.")
    parser.add_argument("fasta_file", help="input genome")
    parser.add_argument("coords", help="start:end")
    args = parser.parse_args()
    input_seqs_path = args.fasta_file
    coords = args.coords


def main():
    start = int(coords.split(":")[0]) - 1  # account here for 0index
    end = int(coords.split(":")[1])
    with open(input_seqs_path, "r") as file_handle:
        lines = file_handle.readlines()
    if lines[0][0] != ">":
        raise ValueError("not a valid fasta!")
    seq = ""
    # print out adjusted coords
    header = str(lines[0].split(" ")[0] + "_" + str(int(start+1)) + ":" + str(end))
    for line in range(1, len(lines)):
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
