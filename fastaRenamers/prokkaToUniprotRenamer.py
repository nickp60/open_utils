#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
version 0.1
"""
import argparse
import sys
import logging
import os


def get_args():
    parser = argparse.ArgumentParser(
        description="given a contigs file from spades, this renames each " +
        "this renames the header line to include the file's basename")
    parser.add_argument("-i", "--infile", dest="infile",
                        help="input .faa from prokka", required=True)
    parser.add_argument("-c", "--clobber", help="overwrite existing outfile",
                        dest="clobber",
                        default=False,
                        action="store_true")
    parser.add_argument("-o", "--outfile", help=" output file",
                        dest="outfile",
                        default=os.path.join(os.getcwd(),
                                             "prokka_Uniprot_renamed.fa"))
    args = parser.parse_args()
    return(args)


def renameContigs(infile, outfile):
    counter = 1
    inname = os.path.splitext(os.path.basename(infile))[0][0:50]
    with open(outfile, "a") as f:
        for line in open(infile):
            if '>' in line:
                locus_tag = line[1:len(line)].split()[0]
                f.write('>sp|' + locus_tag + "|" + line[1:len(line)].strip() + " " +
                        "OS=metagenome_" + inname + " PE=1 SV=1")
                counter = counter + 1
            else:
                f.write(line)


def main():
    logger = logging.getLogger('prokkaToUniprotRenamer')
    logging.basicConfig(stream=sys.stderr, level=logging.DEBUG,
                        format='%(name)s (%(levelname)s): %(message)s')
    logger.setLevel(logging.DEBUG)
    logger.debug("All settings used:")
    args = get_args()
    for k, v in sorted(vars(args).items()):
        logger.debug("{0}: {1}".format(k, v))
    if os.path.isfile(args.outfile):
        if args.clobber:
            os.unlink(args.outfile)
        else:
            logger.error("Output file exists! exiting")
            sys.exit(1)
    renameContigs(args.infile, args.outfile)
    logger.info("Finished")


if __name__ == '__main__':
    main()
