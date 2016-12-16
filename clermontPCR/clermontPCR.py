#!/usr/bin/env python3
#-*- coding: utf-8 -*-
import re
import argparse
import sys
import unittest
import itertools

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
        # int: unique identifier for match
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


def get_matches(seq_list, fwd_primer, rev_primer, expected_size,
                partial=False, strand="+"):
    """given a seqence list  and regex compilations of your primers
    return the matches
    """
    # assert logger is not None, "must use logger!"
    assert strand in ["-", "+"], "strand must be either + or -"
    assert isinstance(seq_list[0], SeqRecord), "must submit list of SeqRecords"
    if strand == "+":
        fwd = re.compile(fwd_primer, re.IGNORECASE)
        rev = re.compile(str(SeqRecord(Seq(rev_primer).reverse_complement()).seq),
                         re.IGNORECASE)
    else:
        fwd = re.compile(rev_primer, re.IGNORECASE)
        rev = re.compile(str(SeqRecord(Seq(fwd_primer).reverse_complement()).seq),
                         re.IGNORECASE)
    matches = []
    for i in seq_list:
        coords_F = None
        coords_R = None
        try:
            coords_F = fwd.search(str(i.seq)).span()
            print("F match!")
        except:
            pass
        try:
            coords_R = rev.search(str(i.seq)).span()
            print("R match!")
        except:
            pass

        if coords_F is not None and coords_R is not None:
            print("Match found on %s (%s)" % (i.id, strand))
            matches.append(PcrHit(
                template_orientation=strand,
                template_id=i.id,
                F_start=coords_F[0],
                R_start=coords_R[0],
                F_end=coords_F[1],
                R_end=coords_R[1],
                partial=partial
            ))
        elif coords_F is not None:
            if len(i.seq[coords_F[0]:]) > expected_size:
                print("Possible match on %s (%s)" % (i.id, strand))
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
        elif coords_R is not None:
            if not coords_R[0] < expected_size:
                print("Possible match on %s (%s)" % (i.id, strand))
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
            # print("No hits on %s" % i.id)
            pass
    return(matches)


def interpret_hits(arpA, chu, TspE4, yjaA):
    if arpA:
        if chu:
            if TspE4 and yjaA:
                result = "U"
            elif not TspE4 and not yjaA:
                result = "E/D"
            elif TspE4:
                result = "E/D"
            else:
                assert yjaA, "error interpretting results!"
                result = "E/Cryptic"
        else:
            if TspE4 and yjaA:
                result = "U"
            elif not TspE4 and not yjaA:
                result = "A"
            elif TspE4:
                result = "B1"
            else:
                assert yjaA, "error interpretting results"
                result = "A/C"
    else:
        if chu:
            if yjaA or TspE4:
                result = "B2"
            else:
                result = "F"
        else:
            if yjaA:
                result = "cryptic"
            else:
                result = "U"
    return(result)


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
    # lets not be hasty
    # ArpAgpE_f = "GATTCCATCTTGTCAAAATATGCC"
    # ArpAgpE_r = "GAAAAGAAAAAGAATTCCCAAGAG"
    trpBA_f = "CGGCGATAAAGACATCTTCAC"
    trpBA_r = "GCAACGCGGCCTGGCGGAAG"

    quad_primers = {"chu": [chuA_1b, chuA_2, 288],
                    "yjaA": [yjaA_1b, yjaA_2b, 211],
                    "TspE4": [TspE4_C2, TspE4C2_2b, 152],
                    "arpA": [AceK_f, ArpA1_r, 400]}

    controls = [trpBA_f, trpBA_r]
    print("Reading in sequence(s)")
    with open(args.contigs, 'r') as fasta:
        seqs = list(SeqIO.parse(fasta, 'fasta'))

    # Start with the Control primers before trying anythong else
    print("Running Control PCR")
    forward_control_matches = get_matches(seq_list=seqs,
                                          fwd_primer=controls[0],
                                          rev_primer=controls[1],
                                          partial=args.partial,
                                          expected_size=489,
                                          strand='+')
    reverse_control_matches = get_matches(seq_list=seqs,
                                          fwd_primer=controls[0],
                                          rev_primer=controls[1],
                                          partial=args.partial,
                                          expected_size=489,
                                          strand='-')

    if (
            len(forward_control_matches) == 0 and
            len(reverse_control_matches) == 0):
        if not args.no_control:
            print("No matches found for control PCR.  Exiting")
            sys.exit(1)
        else:
            print("No matches found for control PCR, but continuing analysis")
    else:
        pass
    # run Clermont Typing
    print("Running Quadriplex PCR")
    profile = ""
    for key, val in sorted(quad_primers.items()):
        print("Scanning %s" % key)
        # fwd = re.compile(val[0], re.IGNORECASE)
        # rev = re.compile(str(SeqRecord(Seq(val[1]).reverse_complement()).seq),
        #                  re.IGNORECASE)
        fwd_matches = get_matches(seq_list=seqs,
                                  fwd_primer=val[0],
                                  rev_primer=val[1],
                                  partial=args.partial,
                                  expected_size=val[2],
                                  strand='+')

        rev_matches = get_matches(seq_list=seqs,
                                  fwd_primer=val[0],
                                  rev_primer=val[1],
                                  partial=args.partial,
                                  expected_size=val[2],
                                  strand='-')
        if len(fwd_matches) != 0 or len(rev_matches) != 0:
            # print("%s: +" % key)
            profile = "{0}\n{1}: +".format(profile, key)
            val.append(True)
        else:
            # print("%s: -" % key)
            profile = "{0}\n{1}: -".format(profile, key)
            val.append(False)
    print("\n-------- Results -------")
    print(profile)
    print("--------   --    -------")
    Clermont_type = interpret_hits(arpA=quad_primers['arpA'][3],
                                   chu=quad_primers['chu'][3],
                                   TspE4=quad_primers['TspE4'][3],
                                   yjaA=quad_primers['yjaA'][3])
    print("Clermont type: %s" % Clermont_type)
    print("------------------------\n")


class clermontTestCase(unittest.TestCase):
    """
    """
    def test_interpret(self):

        ref = ['A', 'A/C', 'B1', 'A/C', 'E/D', 'E/D', 'E/Cryptic', 'E/D',
               'E/D', 'F', 'B2', 'B2', 'B2', 'cryptic', 'E/Cryptic', 'U']
        test = []
        # if True:
        # A's
        test.append(interpret_hits(arpA=True, chu=False, yjaA=False, TspE4=False))
        test.append(interpret_hits(arpA=True, chu=False, yjaA=True, TspE4=False))
        # B1
        test.append(interpret_hits(arpA=True, chu=False, yjaA=False, TspE4=True))
        # C
        test.append(interpret_hits(arpA=True, chu=False, yjaA=True, TspE4=False))
        # E
        test.append(interpret_hits(arpA=True, chu=True, yjaA=False, TspE4=False))
        test.append(interpret_hits(arpA=True, chu=True, yjaA=False, TspE4=True))
        test.append(interpret_hits(arpA=True, chu=True, yjaA=True, TspE4=False))
        # D
        test.append(interpret_hits(arpA=True, chu=True, yjaA=False, TspE4=False))
        test.append(interpret_hits(arpA=True, chu=True, yjaA=False, TspE4=True))
        # F
        test.append(interpret_hits(arpA=False, chu=True, yjaA=False, TspE4=False))
        # B2
        test.append(interpret_hits(arpA=False, chu=True, yjaA=True, TspE4=False))
        test.append(interpret_hits(arpA=False, chu=True, yjaA=False, TspE4=True))
        test.append(interpret_hits(arpA=False, chu=True, yjaA=True, TspE4=True))
        # cryptic
        test.append(interpret_hits(arpA=False, chu=False, yjaA=True, TspE4=False))
        test.append(interpret_hits(arpA=True, chu=True, yjaA=True, TspE4=False))
        # unknown
        test.append(interpret_hits(arpA=False, chu=False, yjaA=False, TspE4=False))
        print(ref)
        print(test)
        self.assertEqual(ref, test)

    def test_get_matches(self):
        test_seq = [SeqRecord(
            Seq("GCACAGTCGATCAAAATTTTTGCAGTCGACTGGACTGACTGTCGGATCTCAGTCAT"))]
        test_fwd = "GCACAG"
        test_rev = "TGACTG"
        self.assertEqual(test_fwd.search(str(test_seq[0].seq)).span(), (0, 6))
        forward_control_matches = get_matches(
            seq_list=test_seq,
            fwd_primer=test_fwd,
            rev=test_rev,
            partial=True,
            expected_size=60,
            strand='+')
        self.assertEqual(forward_control_matches[0].R_end, 55)


if __name__ == "__main__":
    main()
