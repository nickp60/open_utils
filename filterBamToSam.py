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
import pysam

def get_args():
    parser = argparse.ArgumentParser(
        description="given a bam file, this will index, sort, and filter " +
        "out any reads with an alignment score less than --min_AS_score " +
        "(100 by default)")
    parser.add_argument("bam_file", help="input genome")
    parser.add_argument("-a", "--min_AS_score", dest="AS_min",
                        help="region list; " +
                        "if you have a lot of regions, this file " +
                        "will be used instead of the coords arg",
                        default=100)

    parser.add_argument("-o", "--outfile", help=" output file",
                        dest="outfile",
                        default=os.path.join(os.getcwd(),
                                             "filteredOutput.sam"))
    args = parser.parse_args()
    return(args)


def filter_bam_AS(inbam, outsam, score, logger=None):
    """ This is needed because bwa cannot filter based n alignment score
    for paired reads.
    https://sourceforge.net/p/bio-bwa/mailman/message/31968535/
    Given a  bam file from bwa (has "AS" tags), write out
    reads with AS higher than --score to outsam
    read count from https://www.biostars.org/p/1890/
    """
    notag = 0
    written = 0
    unwritten = 0  # cant read my mind
    sortf = os.path.join(os.path.dirname(inbam),
                        os.path.splitext(inbam)[0] + "_sorted.bam")
    print(sortf)
    pysam.sort(inbam, "-o", sortf)
    pysam.index(sortf)
    bam = pysam.AlignmentFile(sortf, "rb")
    osam = pysam.Samfile(outsam, 'wh', template=bam)
    for read in bam.fetch():
        if read.has_tag('AS'):
            if read.get_tag('AS') >= score:
                osam.write(read)
                written = written + 1
            else:
                unwritten = unwritten + 1
                pass
        else:
            notag = notag + 1
            pass
    bam.close()
    logger.debug("Reads before filtering: %i", written + unwritten)
    logger.debug("Reads after filtering: %i", written)
    if notag != 0:
        logger.debug("Reads lacking alignment score: %i", notag)


def main():
    logger = logging.getLogger('filterBamToSam')
    logging.basicConfig(stream=sys.stderr, level=logging.DEBUG,
                        format='%(name)s (%(levelname)s): %(message)s')
    logger.setLevel(logging.DEBUG)
    logger.debug("All settings used:")
    args = get_args()
    for k, v in sorted(vars(args).items()):
        logger.debug("{0}: {1}".format(k, v))
        if args.bam_file.endswith(".sam"):
            logger.error("This acts on BAM files; " +
                         "it looks like you've supplied a SAM File")
            sys.exit(1)
    if os.path.exists(args.outfile):
        logger.error("Output file already exists.  Exiting...")
        sys.exit(1)
    # try:
    filter_bam_AS(inbam=os.path.abspath(os.path.expanduser(args.bam_file)),
                  outsam=args.outfile,
                  score=args.AS_min,
                  logger=logger)
    # except pysam.utils.SamtoolsError as e:
    #     logger.error(e)
    #     sys.exit(1)


if __name__ == '__main__':
    main()
