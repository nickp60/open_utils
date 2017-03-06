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
import re
import logging
import os
from Bio.Seq import Seq

def get_args():
    parser = argparse.ArgumentParser(
        description="Given a set of coords (and an optional chromosome " +
        "or contig), extract a region of the genome.  It goes to " +
        "standard out, which can easily be redirected to a file with '>'")
    parser.add_argument('coords', nargs='?', type=str,
                        default=sys.stdin,
                        help="start:end or chromosome:start:end; " +
                        "can be read from standard in.  If reading from " +
                        "stdin, can take an additional 'name' attribute: " +
                        "name@chromosome:start:end;")
    parser.add_argument("fasta_file", help="input genome")
    parser.add_argument("-n", "--name", help="An optional name to be added " +
                        "to the fasta header on output",
                        default="")
    args = parser.parse_args()
    return(args)


def parse_coords(coords, from_stdin=False):
    coords = coords.replace("'", "")
    if from_stdin and "@" in coords:
        name, coords = coords.split("@")
        # name = name + "_"
    else:
        name = ''
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
    return(name, chrom, start, end)


def get_chrom_of_interest(chrom, lines, warn_multiple=False, logger=None):
    """ By defualt, this will get the first occurance of the chromorome
    of interest.  So, a search for "chrom13" will be returned even if you
    wanted "chrom133".  the "warn multiple flad
    """
    assert logger is not None, "must use logging"
    startindex = None
    endindex = None
    maxlines = len(lines)
    hits = 0
    # while startindex is None:
    for index, line in enumerate(lines):
        if line[0] == ">":
            try:
                if len(re.search(chrom, line).group(0)) > 0:
                    logger.info("Hit for %s found: %s" % (chrom, line.strip()))
                    startindex = index
                    hits = hits + 1
            except AttributeError:
                continue
        if startindex is not None:
            for index2, line2 in enumerate(lines[index + 1: ]):
                if line2[0] == ">":
                    try:
                        if len(re.search(chrom, line2).group(0)) > 0:
                            logger.info("Hit for %s found: %s" % (chrom, line2.strip()))
                            hits = hits + 1
                    except AttributeError:
                        continue
            if hits != 1:  # this is depreciated
                raise ValueError("Hits for %s were found in mutiple headers!" % chrom)
            else:
                break
    if index == (maxlines - 1) and startindex is None:
        raise ValueError("%s not found in any headers!" % chrom)
    # while endindex is None:
    for nline in range(startindex + 1, maxlines):
        if lines[nline][0] == ">" or int(nline) == int(maxlines - 1):
            endindex = nline
    return(startindex, endindex)

#%%


def main():
    logger = logging.getLogger('extract_region')
    logging.basicConfig(stream=sys.stderr, level=logging.DEBUG,
                        format='%(name)s (%(levelname)s): %(message)s')
    logger.setLevel(logging.DEBUG)
    logger.debug("All settings used:")
    args = get_args()
    for k, v in sorted(vars(args).items()):
        logger.debug("{0}: {1}".format(k, v))

    COORDS_FROM_STDIN = False
    if type(args.coords) == file:
        # with open(args.coords) as c:
        coord_string = args.coords.readlines()[0].strip()
        COORDS_FROM_STDIN = True
    else:
        coord_string = args.coords

    name, chrom, start, end = parse_coords(coord_string,
                                           from_stdin=COORDS_FROM_STDIN)
    logger.info(name)
    logger.info(chrom)
    logger.info(start)
    logger.info(end)
    if name == "" and args.name != "":
        name = args.name
    with open(args.fasta_file, "r") as file_handle:
        lines = file_handle.readlines()
    if lines[0][0] != ">":
        raise ValueError("not a valid fasta!")
    if chrom is not None:
        try:
            chrom_start, chrom_end = get_chrom_of_interest(
                chrom, lines, logger=logger)
        except Exception as e:
            logger.error(str(type(e)) + ": " + str(e))
            sys.exit(1)
    else:
        chrom_start, chrom_end = 0, len(lines)

    seq = ""
    header = str(lines[0 + chrom_start].strip().split(" ")[0] + "_" +
                 name + "_" + str(int(start + 1)) + ":" + str(end))
    for line in range(1, len(lines[chrom_start:chrom_end])):
        if len(seq) > end:  # only read in up to what you need
            break
        # if lines[line][0] == ">":
        #     raise ValueError("only for use with single fastas")
        seq = str(seq + lines[line].strip())
    linewidth = 80  # standard
    if name.endswith("_RC"):
        logger.debug("getting RC")
        try:
            subset = str(Seq(seq[start:end]).reverse_complement())
        except Exception as e:
            logger.error("could not get reverse complement with BioPython module")
            subset = seq[start: end]
    else:
        subset = seq[start:end]
    newstring = ""
    for i in range(0, len(subset)):
        if i >= linewidth and (i % linewidth == 0):
            newstring = str(newstring + "\n" + subset[i])
        else:
            newstring = str(newstring + subset[i])
    sys.stdout.write(header)
    sys.stdout.write("\n")
    sys.stdout.write(newstring)
    sys.stdout.write("\n")  # to get rid of '%' at end shown in terminal

if __name__ == '__main__':
    main()
