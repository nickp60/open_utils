#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 11:44:26 2016
version 0.1
Minor version changes:
 - reads from stdin

TODO:

USAGE:
 $ extractFromMSA.py genome.fasta  > extractedregionfromcontig.fasta

"""
import argparse
import sys
import re
import logging
import os
from Bio import SeqIO


def get_args():
    parser = argparse.ArgumentParser(
        description="subset a list of chromosomes or sequences from a " +
        "multifasta, be it a MSA or whatever. It goes to " +
        "standard out, which can easily be redirected to a file with '>'")
    parser.add_argument("fasta_file",
                        help="input files (multifasta or aligned multifasta")
    parser.add_argument("-e", "--end_char",
                        help="if, for example, getting Lys3 in headerline  " +
                        "'something@Lys3_c12.fasta', end_char is '_'; " +
                        "This char terminates search;",
                        default="_")
    parser.add_argument("names", nargs="+")
    args = parser.parse_args()
    return(args)


def extract_record(path, chrom, end_char):
    rec_of_interest = None
    recs = []
    hits = 0
    with open(path, "r") as fasta:
        for record in SeqIO.parse(fasta, "fasta"):
            if str(chrom.strip() + end_char) in record.id:
                rec_of_interest = record
                recs.append(record.id)
                hits = hits + 1
            else:
                pass
    if hits == 0:
        raise ValueError("No hits found for %s", chrom)
    elif hits > 1:
        raise ValueError("Multiple hits found for %s: \n%s", chrom,
                         "\n".join(recs))
    else:
        return(rec_of_interest)


def main():
    logger = logging.getLogger('extract_region')
    logging.basicConfig(stream=sys.stderr, level=logging.DEBUG,
                        format='%(name)s (%(levelname)s): %(message)s')
    logger.setLevel(logging.DEBUG)
    logger.debug("All settings used:")
    args = get_args()
    for k, v in sorted(vars(args).items()):
        logger.debug("{0}: {1}".format(k, v))

    for i in args.names:
        logger.debug("searching for %s", i)
        try:
            record = extract_record(
                path=args.fasta_file, chrom=i, end_char=args.end_char)
        except Exception as e:
            logger.error(str(type(e)) + ": " + str(e))
            sys.exit(1)
        SeqIO.write(record, sys.stdout, "fasta")


if __name__ == '__main__':
    main()
