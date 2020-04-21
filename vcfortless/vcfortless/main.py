#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
import argparse
import vcf
import sys
import gzip
import os
from collections import namedtuple
from operator import attrgetter
# import io
from Bio import SeqIO
from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord
# from Bio.SeqIO.FastaIO import SimpleFastaParse

# BI: whether a allele is bialleleic 1 or not 0
# CG: change: transversion TV or TS
# SYN: 0 synonymous, 1 non-synonymous
# EFF: Effect, one of Premature Stop Codon, putative PROmoter disruption inTERgenic, inTRAgenic
Result = namedtuple('Result', 'CHROM  POS REF ALT BI CHG SYN EFF')
DEBUG=False

def get_args():
    parser = argparse.ArgumentParser(
        description="Given a VCF and a genbank file, writes out a report")
    parser.add_argument('vcf', help="path to vcf")
    parser.add_argument("gbk", help="path to genbank")
    parser.add_argument("-t", "--trans_table",
                        help="translation table; default 11",
                        type=int, default=11)
    parser.add_argument("-f", "--feature",
                        help="which feature to consider: gene or cds",
                        choices=["gene", "cds"],
                        default="gene")
    parser.add_argument("-b", "--binding_width",
                        help="width to mark whether snps might disrupt RBS",
                        type=int, default=15)
    args = parser.parse_args()
    return(args)

def tvts(ref, alt):
    valid = ["A", "T", "C", "G"]
    for nuc in [ref, alt]:
        if nuc not in valid:
            sys.stderr.write("Non-standard nucleotide: %s\n" % nuc)
            return ("x")
    refs = {
        "A": {"A": "-",
              "T": "tv",
              "C": "tv",
              "G": "ts"},
        "T": {"A": "tv",
              "T": "-",
              "C": "ts",
              "G": "tv"},
        "C": {"A": "tv",
              "T": "ts",
              "C": "-",
              "G": "tv"},
        "G": {"A": "ts",
              "T": "tv",
              "C": "tv",
              "G": "-"}
    }
    result = refs[ref.upper()][alt.upper()]
    return (result)


def subs_nuc(refseq, start, end, pos, alt):
    #gene ------
    # ----;----;
    #        *
    # snp: 8
    # start = 5
    # start = 4 in python, but biopython uses 0-based, but vcf used 1-based
    # end = 10
    # region seq[5 - 1: 8 - 1 ] + SNP + seq[ SNP : end ]
    thisseq =  refseq[start : pos ] + alt +refseq[pos + 1  : end ]
    # print(refseq[start: end])
    # print(thisseq)
    assert len(thisseq) == len(refseq[start: end]), "bad reconstruction of  reference; ref length is %i and reconstructed length is %i" %(len(refseq[start: end]), len(thisseq))
    return thisseq

def test_subs_nuc_psc():
    refseq = "ATGCCCAAATTTTACTAG"
    mutseq = "ATGCCCAAATTTTAGTAG"
    newseq =  subs_nuc(refseq, start=0, end=18, pos=14, alt="G")
    assert   mutseq == newseq, "error in sub_nuc function"

def test_subs_nuc_norm():
    refseq = "ATGCCCAAATTTTACTAG"
    mutseq = "ATGCCCAAATTTTATTAG"
    assert mutseq == subs_nuc(refseq, start=0, end=18, pos=14, alt="T"), "error in sub_nuc function"

def test_pmc():
    refseq = "ATGCCCAAATTTTACTAG"
    mutseq = "ATGCCCAAATTTTAGTAG"
    thisp = Seq(mutseq).translate(table=11, to_stop=True)
    refp = Seq(refseq).translate(table=11, to_stop=True)
    print(thisp)
    print(refp)
    assert len(refp) != len(thisp), "error detecting premature stop codon!"

def process_region(args, vcf_data, chrom, start, end, rec, strand, is_locus=False):
    if is_locus:
        assert rec is not None, "must provide rec for loci"
        assert strand is not None, "must provide feature for loci"
        nucseq  = rec.seq[start: end]
        if strand == 1:
            nucseqp = nucseq.translate(table=args.trans_table, to_stop=True)
        else:
            nucseqp = nucseq.reverse_complement().translate(table=args.trans_table, to_stop=True)
    these_vcfs = vcf_data[chrom][start: end]
    ignored = 0
    for pos, ref, altlist, PROCESS in these_vcfs:
        if len(ref) != 1:
            ignored = ignored + 1
            continue
        if not PROCESS:
            continue
        bialleleic = False
        if len(altlist) > 1:
            biallelic = True

        for alt in altlist:
            if len(alt) > 1:
                ignored = ignored + 1
                continue
            thiststv = tvts(ref, str(alt))
            if is_locus:
                try:
                    thisseq = subs_nuc(rec.seq, start, end, pos, str(alt))
                except AssertionError:
                    sys.stderr.write("start: %i; end %i; pos: %i ; alt: %s\n" %(start, end, pos, str(alt)))
                    sys.exit(1)
                
                assert len(thisseq) == len(nucseq), "bad reconstruction of  reference"
                if strand == 1:
                    thisseqp = thisseq.translate(table=args.trans_table, to_stop=True)
                else:
                    thisseqp = thisseq.reverse_complement().translate(table=args.trans_table, to_stop=True)
                if DEBUG:
                    print(nucseq)
                    print(thisseq)

                    print(rec.seq[start : pos ] + ref)
                    print(nucseqp)
                    print(thisseqp)
                SYN = 1
                EFF = "TRA"
                if nucseqp != thisseqp:
                    SYN = 0
                    if len(thisseqp) != len(nucseqp):
                        EFF = "PSC"
                # back to 1-indexed
                thisres = Result(chrom, pos+1, ref, alt, bialleleic, thiststv, SYN, EFF)
                sys.stdout.write("%s\t%i\t%s\t%s\t%i\t%s\t%i\t%s\n" % thisres)
            else:
                # process intergenic region
                EFF = "TER"
                if (
                        (pos - start) < args.binding_width or
                        (end - pos) < args.binding_width
                ):
                    EFF = "PRO"
                thisres = Result(chrom, pos, ref, alt, 0, thiststv, 0, EFF)
                sys.stdout.write("%s\t%i\t%s\t%s\t%i\t%s\t%i\t%s\n" % thisres)
    return ignored


def main(args=None):
    """
    """
    if args is None:
        args = get_args()

    gbk_open_fun = open
    vcf_open_fun = open
    if os.path.splitext(args.gbk)[-1] in ['.gz', '.gzip']:
        gbk_open_fun = gzip.open
    if os.path.splitext(args.vcf)[-1] in ['.gz', '.gzip']:
        vcf_open_fun = gzip.open
    vcf_reader = vcf.Reader(vcf_open_fun(args.vcf, 'r'))
    found_one = False
    chroms = []


    sys.stderr.write("Getting IDs from Genbank\n")
    vcf_data = {}
    with gbk_open_fun(args.gbk, "r") as ingbk:
        for rec in SeqIO.parse(ingbk, "genbank"):
            # print(rec.features[2])
            # sys.exit()
            vcf_data[rec.id.split(".")[0]] = []

    sys.stderr.write("Reading in vcf\n")
    # we do this weird counter thing so that we have an entry for each position in the genome
    # its a dumb idea until you have to deal with subsets of this list, in which case trading off the ram for the speed

    # 
    prev_pos = 0  # we keep track of previous position se we know when to reset the counter for new contigs
    for i, v in enumerate(vcf_reader):
        # if (i % 1000) == 0:
        #     sys.stderr.write(str(i) + " ")
        # here wer set to counter 
        if v.POS < prev_pos or i == 0:
            counter = 1
        if v.POS > 200000000:
            sys.stderr.write("Warning: long sequence detected, only processing the first 20Mb")
            break
        # this pads out for the non-snp regions
        while counter != v.POS and counter < v.POS:
            vcf_data[v.CHROM].append([counter, "-", "-", False])
            counter = counter + 1
        assert counter == v.POS, "error syncing counters;\n -chrom: %s \n-position: %i \n -previous: %i \n -counter: %i" % (v.CHROM, v.POS, prev_pos, counter)
        # make 0-indexed
        vcf_data[v.CHROM].append([v.POS-1, v.REF, v.ALT, True])
        counter = counter + 1
        prev_pos = v.POS
    last_gene_end = 0
    # first process all the coding sequences, then hit the remaining intergenic loci
    ignored_positons = 0
    with gbk_open_fun(args.gbk, "r") as ingbk:
        for rec in SeqIO.parse(ingbk, "genbank"):
            thischrom = rec.id.split(".")[0]
            sys.stderr.write("Processing %s\n" % thischrom)
            for feat in rec.features:
                #if feat.type not in ["source"]:
                if feat.type == args.feature:
                    # process coding region
                    ig = process_region(
                        args, vcf_data,
                        chrom=thischrom,
                        start=feat.location.start,
                        end=feat.location.end,
                        rec=rec,
                        strand=feat.strand,
                        is_locus=True)
                    ignored_positons = ignored_positons + ig
                    ig = process_region(
                        args, vcf_data,
                        chrom=thischrom,
                        start=last_gene_end,
                        end=feat.location.start,
                        rec=rec,
                        strand=feat.strand,
                        is_locus=False)
                    ignored_positons = ignored_positons + ig
                    # keep track of where the last gene ended
                    last_gene_end = feat.location.end
    if ignored_positons != 0:
        sys.stderr.write("ignored %d complex entries\n" %ignored_positons)

if __name__ == '__main__':
    main()
