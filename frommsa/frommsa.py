#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
import argparse
import sys
import os
from Bio.Seq import Seq
from Bio import AlignIO

def get_args():
    parser = argparse.ArgumentParser(
        description="Given a MSA, convert to phylip etc.")
    parser.add_argument('msa', help="path to msa")
    parser.add_argument(
        '-f', '--outfmt',
        default="phylip",
        help="output format, anything handled by biopython")
    parser.add_argument(
        '-o', '--outfile',
        help="output file; ignore to use stdout")
    return(parser.parse_args())


if __name__ == "__main__":
    args = get_args()
    outname=args.msa.replace('.msa', '.phy')
    with open(args.msa, "r") as inf:
        alignments = AlignIO.parse(inf, 'fasta')
        if args.outfile is None:
            args.outfile = sys.stdout
            AlignIO.write(alignments, sys.stdout, args.outfmt)
        else:
            with open(args.outfile, "w") as outf:
                AlignIO.write(alignments, outf, args.outfmt)
