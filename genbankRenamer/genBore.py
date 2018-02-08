#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
import argparse
import sys
import logging
import os

def get_args():
    parser = argparse.ArgumentParser(
        description="given a genbank with long LOCUS names, rename each " +
        "to include the file's basename")
    parser.add_argument("-i", "--infile",
                        dest="infile",
                        help="input genbank", required=True)
    parser.add_argument("-n", "--name",
                        dest="name", help="name", required=False)
    parser.add_argument("-c", "--clobber",
                        help="overwrite existing outfile",
                        dest="clobber",
                        default=False,
                        action="store_true")
    parser.add_argument("-o", "--outfile",
                        help=" output file; otherwise, uses stdout",
                        dest="outfile")
    args = parser.parse_args()
    return(args)


def renameContigs(infile, name, outfile):
    """fix genbank LOCUS when the contig had a long name
    LOCUS       NODE_1_length_426502_cov_39.297426502 bp   DNA linear
    00:06      LOCUS
    06:12      spaces
    12:??      Locus name
    ??:??      space
    ??:40      Length of sequence, right-justified
    40:44      space, bp, space
    44:47      Blank, ss-, ds-, ms-
    47:54      Blank, DNA, RNA, tRNA, mRNA, uRNA, snRNA, cDNA
    54:55      space
    55:63      Blank (implies linear), linear or circular
    63:64      space
    64:67      The division code (e.g. BCT, VRL, INV)
    67:68      space
    68:79      Date, in the form dd-MMM-yyyy (e.g., 15-MAR-1991)
     """
    counter = 1
    if name is None:
        name = os.path.splitext(os.path.basename(infile))[0]
    if len(name) > 15:
        name = name[0:15]
    # handle noth file output and stdout
    f = open(target, 'w') if outfile else sys.stdout
    for line in open(infile):
        if line.startswith("LOCUS"):
            if "length" in line:
                # must be 10 chars
                len_maybe = line.split("length")[1].split("_")[1].strip()
            padded_counter = str(counter).zfill(5)
            new_name = str(name + "_c" + padded_counter)
            new_header = str('LOCUS ' +            # 6
                            " " * 6 +              # 6
                             new_name.ljust(18) +  # 28
                             len_maybe.rjust(10) + # 10  40
                             " bp " +              # 4
                             "   " +               # 3
                             " DNA   " +           # 7   54
                             " " * 10 +            # 10  64
                             "BCT"                 # 3   67
                             " 01-APR-1492\n")     # 12  79
            counter = counter + 1
            f.write(new_header)
        else:
            f.write(line)
    if f is not sys.stdout:
        f.close()


def main():
    logger = logging.getLogger('getBore')
    logging.basicConfig(stream=sys.stderr, level=logging.DEBUG,
                        format='%(name)s (%(levelname)s): %(message)s')
    logger.setLevel(logging.DEBUG)
    logger.debug("All settings used:")
    args = get_args()
    for k, v in sorted(vars(args).items()):
        logger.debug("{0}: {1}".format(k, v))
    if args.outfile is not None:
        if os.path.isfile(args.outfile):
            if args.clobber:
                os.unlink(args.outfile)
            else:
                logger.error("Ouput file exists! exiting")
                sys.exit(1)
    renameContigs(args.infile, args.name, args.outfile)
    logger.info("Finished")


if __name__ == '__main__':
    main()
