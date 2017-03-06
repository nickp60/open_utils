#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 18 11:44:26 2016
version 0.6
Minor version changes:
 - reads from stdin

TODO:

USAGE:
 $ python extract_region.py genome.fasta 3:111 > extractedregion.fasta
 $ python extract_region.py genome.fasta chrom1:3:111 > extractedregionfromcontig.fasta
 $ echo "mySuperGene@chrom1:3:111" | xargs extract_region.py genome.fasta  > extractedregionfromcontig.fasta

"""
import argparse
import sys
import logging
import os

def get_args():
    parser = argparse.ArgumentParser(
        description="given a contigs file from spades, this renames each " +
        "this renames the header line to include the file's basename")
    parser.add_argument("-i", "--infile", dest="infile", help="input genome", required=True)
    parser.add_argument("-n", "--name", dest="name", help="name", required=True)
    parser.add_argument("-o", "--outfile", help=" output file",
                        dest="outfile",
                        default=os.path.join(os.getcwd(),
                                             "spadesRenamed.fa"))
    args = parser.parse_args()
    return(args)


def renameContigs(infile, name, outfile):
    counter=1
    # name = os.path.splitext(os.path.basename(infile))
    with open(outfile, "a") as f:
        for line in open(infile):
            if '>' in line:
                f.write('>' + name + "_c" + str(counter) + " " +
                        line[1:len(line)])
                counter = counter + 1
            else:
                f.write(line)


def main():
    logger = logging.getLogger('filterBamToSam')
    logging.basicConfig(stream=sys.stderr, level=logging.DEBUG,
                        format='%(name)s (%(levelname)s): %(message)s')
    logger.setLevel(logging.DEBUG)
    logger.debug("All settings used:")
    args = get_args()
    for k, v in sorted(vars(args).items()):
        logger.debug("{0}: {1}".format(k, v))
    renameContigs(args.infile, args.name, args.outfile)


if __name__ == '__main__':
    main()
