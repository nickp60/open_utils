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
        description="given a contigs file from spades, this renames each " +
        "this renames the header line to include the file's basename")
    parser.add_argument("-i", "--infile", dest="infile", help="input genome", required=True)
    parser.add_argument("-n", "--name", dest="name", help="name", required=True)
    parser.add_argument("-c", "--clobber", help="overwrite existing outfile",
                        dest="clobber",
                        default=False,
                        action="store_true")
    parser.add_argument("-o", "--outfile", help=" output file",
                        dest="outfile",
                        default=os.path.join(os.getcwd(),
                                             "spadesRenamed.fa"))
    args = parser.parse_args()
    return(args)


def renameContigs(infile, name, outfile):
    counter = 1
    # name = os.path.splitext(os.path.basename(infile))
    with open(outfile, "a") as f:
        for line in open(infile):
            if '>' in line:
                f.write('>' + name + "_c" + str(counter) + " " + name + " " +
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
