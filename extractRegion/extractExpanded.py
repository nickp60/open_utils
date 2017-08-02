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
import io
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import SimpleFastaParser


from extractRegion import parse_coords, extract_record, \
    check_missing_or_duplicated

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
    parser.add_argument("-l", "--list", dest="reglist",
                        help="region list; " +
                        "if you have a lot of regions, this file " +
                        "will be used instead of the coords arg")
    parser.add_argument("-e", "--extra", dest="extra",
                        help=" Get a region surrounding the region given," +
                        "by this many bases "
                        "default: %(default)s", type=int, default=0)
    parser.add_argument("-q", "--quick", dest="quick",
                        help="uses the experimental SimpleFastaParser;",
                        action="store_false", default=True)

    parser.add_argument("-n", "--name", help="An optional name to be added " +
                        "to the fasta header on output",
                        default="")
    args = parser.parse_args()
    return(args)


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
    if sys.stdin.isatty() and args.reglist is None:
        # with open(args.coords) as c:
        coord_string_list = [args.coords.readlines()[0].strip()]
        COORDS_FROM_STDIN = True
    else:
        if args.reglist is not None:
            logger.info("reading coords from file")
            coord_string_list = open(args.reglist).readlines()
            COORDS_FROM_STDIN = True
        else:
            coord_string_list = [args.coords]
    #hold all them names and stuff
    nameChromStartEnd = []
    for coords in [x.strip() for x in coord_string_list if x.strip() ]:
        name, chrom, start, end = parse_coords(coords,
                                               from_stdin=COORDS_FROM_STDIN,
                                               logger=logger)
        nameChromStartEnd.append([name, chrom, start, end])
    # cache record list
    rec_ids = []
    with open(args.fasta_file, "r") as f:
        for rec in SeqIO.parse(f, 'fasta'):
            rec_ids.append(str(rec.id + " "))
    check_missing_or_duplicated(fastaids=rec_ids,
                                extids=[x[1] for x in nameChromStartEnd])
    nregions = len(nameChromStartEnd)
    for idx, ncse in enumerate(nameChromStartEnd):
        logger.debug("processing entry %i of %i:", idx, nregions)
        logger.debug(ncse)
        if ncse[0] == "" and args.name != "":
            ncse[0] = args.name
        logger.debug("extracting sequence named %s", ncse[1])
        if ncse[1] is not None:
            if args.quick:
                try:
                    record = qextract_record(path=args.fasta_file,
                                             chrom=ncse[1].strip())
                except Exception as e:
                    logger.error(str(type(e)) + ": " + str(e))
                    sys.exit(1)
            else:
                try:
                    record = extract_record(path=args.fasta_file,
                                            chrom=ncse[1].strip())
                except Exception as e:
                    logger.error(str(type(e)) + ": " + str(e))
                    sys.exit(1)
        else:
            with open(args.fasta_file) as f:
                records = list(SeqIO.parse(f, "fasta"))
            if len(records) != 1:
                logger.error("fasta contains multiple records! must " +
                             "provide a sequence ID")
                sys.exit(1)
            else:
                record = records[0]

        ## get the region
        # seqid = "{0}_{1}{2}:{3}".format(record.id, name, start, end)
        seqid = "{0}@{1}:{2}:{3}".format(
            ncse[0], ncse[1],
            ncse[2] + args.extra,
            ncse[3] + args.extra)
        if ncse[0].endswith("-RC_"):
            logger.debug("getting RC")
            try:
                subrec = SeqRecord(record.seq[ncse[2]: ncse[3]].reverse_complement(),
                                   id=seqid, description="")
            except Exception as e:
                logger.error(str(type(e)) + ": " + str(e))
                sys.exit(1)
        else:
            logger.debug("getting sequences")
            try:
                subrec = SeqRecord(record.seq[ncse[2]: ncse[3]],
                                   id=seqid, description="")
            except Exception as e:
                logger.error(str(type(e)) + ": " + str(e))
                sys.exit(1)
        SeqIO.write(subrec, sys.stdout, "fasta")


if __name__ == '__main__':
    # help(SimpleFastaParser)
    # sys.exit(0)
    main()
