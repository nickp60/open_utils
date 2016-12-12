#!/usr/bin/env python3
#-*- coding: utf-8 -*-
import re
import argparse
import sys
import Bio
import itertools
import logging

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


class PcrHit(object):
    newid = itertools.count()

    def __init__(self, index=None,
                 template_orientation=None,
                 template_id=None,
                 true_hit=False,
                 F_start=None,
                 R_start=None,
                 F_end=None,
                 R_end=None,
                 F_hit=None,
                 R_hit=None,
                 partial=False):
        # int: unique identifier for cluster
        self.index = next(PcrHit.newid)
        self.template_orientation = template_orientation
        self.template_id = template_id
        self.true_hit = true_hit
        self.F_start = F_start
        self.F_end = F_end
        self.R_start = R_start
        self.R_end = R_end
        self.F_hit = F_hit
        self.R_hit = R_hit
        self.partial = partial
        self.parse_hits()

    def parse_hits(self):
        if self.F_start is not None and self.F_end is not None:
            self.F_hit = True
        if self.R_start is not None and self.R_end is not None:
            self.R_hit = True
        if self.F_hit and self.R_hit:
            self.true_hit = True
        elif (self.F_hit or self.R_hit) and self.partial:
            self.true_hit = True
        elif (self.F_hit or self.R_hit):
            self.true_hit = False


def get_args():  # pragma: no cover
    parser = argparse.ArgumentParser(
        description="run a 'PCR' to get clermont types",
        add_help=False)  # to allow for custom help
    parser.add_argument("contigs", action="store",
                        help="FASTA formatted genome or set of contigs")

    # # taking a hint from http://stackoverflow.com/questions/24180527
    # requiredNamed = parser.add_argument_group('required named arguments')
    # requiredNamed.add_argument("-F", "--fastq1", dest='fastq1', action="store",
    #                            help="forward fastq reads, can be compressed",
    #                            type=str, default="", required=True)
    # # had to make this faux "optional" parse so that the named required ones
    # # above get listed first
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("-p", "--partial", dest='partial',
                          action="store_true",
                          help="If scanning contigs, breaks between " +
                          "contigs could potentially contain your " +
                          "sequence of interest.  if --partial, partial " +
                          "matches that could be intereupted by contig " +
                          "breaks are reported",
                          default=False)
    optional.add_argument("-c", "--no_control", dest='no_control',
                          action="store_true",
                          help="ignore failure of control PCR",
                          default=False)
    # had to make this explicitly to call it a faux optional arg
    optional.add_argument("-h", "--help",
                          action="help", default=argparse.SUPPRESS,
                          help="Displays this help message")
    args = parser.parse_args()
    return args


def get_matches(seq_list, fwd, rev, expected_size, partial=False, strand="+", logger=None):
    """given a seqence list  and regex compilations of your primers
    return the matches
    """
    assert logger is not None, "must use logger!"
    assert strand in ["-", "+"], "strand must be either + or -"
    if strand == "-":
        seqs = [x.reverse_complement() for x in seq_list]
    else:
        seqs = seq_list
    matches = []
    for i in seqs:
        coords_F = None
        coords_R = None
        try:
            coords_F = fwd.search(str(i.seq)).span()
            logger.info("F match!")
        except:
            pass
        try:
            coords_R = rev.search(str(i.seq)).span()
            logger.info("R match!")
        except:
            pass

        if coords_F is not None and coords_R is not None:
            logger.info("Match found on %s (%s)" % (i.id, strand))
            matches.append(PcrHit(
                template_orientation=strand,
                template_id=i.id,
                F_start=coords_F[0],
                R_start=coords_R[0],
                F_end=coords_F[1],
                R_end=coords_R[1],
                partial=partial
            ))
        elif coords_F is not None and coords_R is None:
            if len(i.seq[coords_F[0]:]) > expected_size:
                logger.info("Possible match on %s (%s)" % (i.id, strand))
                matches.append(PcrHit(
                    template_orientation=strand,
                    template_id=i.id,
                    F_start=coords_F[0],
                    R_start=None,
                    F_end=coords_F[1],
                    R_end=None,
                    partial=partial
                ))
            else:
                pass
        elif coords_R is not None and coords_F is None:
            if not coords_R[0] < expected_size:
                logger.info("Possible match on %s (%s)" % (i.id, strand))
                matches.append(PcrHit(
                    template_orientation=strand,
                    template_id=i.id,
                    R_start=coords_R[0],
                    F_start=None,
                    R_end=coords_R[1],
                    F_end=None,
                    partial=partial
                ))

            else:
                pass
    return(matches)


def main():
    args = get_args()
    chuA_1b = "ATGGTACCGGACGAACCAAC"
    chuA_2  = "TGCCGCCAGTACCAAAGACA"
    yjaA_1b = "CAAACGTGAAGTGTCAGGAG"
    yjaA_2b = "AATGCGTTCCTCAACCTGTG"
    TspE4_C2 = "CACTATTCGTAAGGTCATCC"
    TspE4C2_2b = "AGTTTATCGCTGCGGGTCGC"
    AceK_f = "AACGCTATTCGCCAGCTTGC"
    ArpA1_r = "TCTCCCCATACCGTACGCTA"
    ArpAgpE_f = "GATTCCATCTTGTCAAAATATGCC"
    ArpAgpE_r = "GAAAAGAAAAAGAATTCCCAAGAG"
    trpBA_f = "CGGCGATAAAGACATCTTCAC"
    trpBA_r = "GCAACGCGGCCTGGCGGAAG"

    quad_primers = {"chu": [chuA_1b, chuA_2, 288],
                    "yjaA": [yjaA_1b, yjaA_2b, 211],
                    "TspE4": [TspE4_C2, TspE4C2_2b, 152],
                    "AceK": [AceK_f, ArpA1_r, 400]}

    controls = [trpBA_f, trpBA_r]
    logger = logging.getLogger('root')
    with open(args.contigs, 'r') as fasta:
        seqs = list(SeqIO.parse(fasta, 'fasta'))
    # Start with the Control primers before trying to
    fwd = re.compile(controls[0], re.IGNORECASE)
    rev = re.compile(str(SeqRecord(Seq(controls[1]).reverse_complement()).seq),
                     re.IGNORECASE)
    forward_control_matches = get_matches(seq_list=seqs, fwd=fwd, rev=rev,
                                          partial=args.partial,
                                          expected_size=489,
                                          strand='+',
                                          logger=logger)
    reverse_control_matches = get_matches(seq_list=seqs, fwd=fwd, rev=rev,
                                          partial=args.partial,
                                          expected_size=489,
                                          strand='-',
                                          logger=logger)

    if len(forward_control_matches) and len(reverse_control_matches) == 0:
        if not args.no_control:
            logger.info("No matches found for control PCR.  Exiting")
            sys.exit(1)
        else:
            logger.info("No matches found for control PCR, but continuing analysis")
    # run Clermont Typing
    for key, val in quad_primers.items():
        fwd = re.compile(val[0], re.IGNORECASE)
        rev = re.compile(str(SeqRecord(Seq(val[1]).reverse_complement()).seq),
                         re.IGNORECASE)
        fwd_matches = get_matches(seq_list=seqs, fwd=fwd, rev=rev,
                                  partial=args.partial,
                                  expected_size=val[2],
                                  strand='+',
                                  logger=logger)

        rev_matches = get_matches(seq_list=seqs, fwd=fwd, rev=rev,
                                  partial=args.partial,
                                  expected_size=val[2],
                                  strand='-',
                                  logger=logger)

        if len(fwd_matches) or len(rev_matches) != 0:
            print("%s: +" % key)
        else:
            print("%s: -" % key)


if __name__ == "__main__":
    main()
