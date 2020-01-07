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
# from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord
# from Bio.SeqIO.FastaIO import SimpleFastaParse

# BI: whether a allele is bialleleic 1 or not 0
# CG: change: transversion TV or TS
# SYN: 0 synonymous, 1 non-synonymous
# EFF: Effect, one of Premature Stop Codon, putative PROmoter disruption inTERgenic, inTRAgenic
Result = namedtuple('Result', 'CHROM  POS REF ALT BI CHG SYN EFF')


def get_args():
    parser = argparse.ArgumentParser(
        description="Given a VCF and a genbank file, writes out a report")
    parser.add_argument('vcf', help="path to vcf")
    parser.add_argument("gbk", help="path to genbank")
    parser.add_argument("-t", "--trans_table",
                        help="translation table; default 11",
                        type=int, default=11)
    parser.add_argument("-b", "--binding_width",
                        help="width to mark whether snps might disrupt RBS",
                        type=int, default=15)
    args = parser.parse_args()
    return(args)

def tvts(ref, alt):
    valid = ["A", "T", "C", "G"]
    for nuc in [ref, alt]:
        if nuc not in valid:
            sys.stderr.write("Non-standard nucleotide: %s" % nuc)
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
    for pos, ref, altlist, PROCESS in these_vcfs:
        if not PROCESS:
            continue
        bialleleic = False
        if len(altlist) > 1:
            biallelic = True

        for alt in altlist:
            thiststv = tvts(ref, str(alt))
            if is_locus:
                #gene ------
                # ----;----;
                #        *
                # snp: 8
                # start = 5
                # start = 4 in python, buy biopython uses trad numbering
                # end = 10
                # region seq[5 - 1: 8 - 1 ] + SNP + seq[ SNP : end ]
                thisseq =  rec.seq[start - 1: pos - 1] + ref +rec.seq[pos : end]
                if strand == 1:
                    thisseqp = thisseq.translate(table=args.trans_table, to_stop=True)
                else:
                    thisseqp = thisseq.reverse_complement().translate(table=args.trans_table, to_stop=True)
                SYN = 1
                EFF = "TER"
                if nucseqp != thisseqp:
                    SYN = 0
                    if len(thisseqp) != len(nucseqp):
                        EFF = "PSC"
                thisres = Result(chrom, pos, ref, alt, bialleleic, thiststv, SYN, EFF)
                sys.stdout.write("%s\t%i\t%s\t%s\t%i\t%s\t%i\t%s\n" % thisres)
            else:
                # process intergenic region
                EFF = "TRA"
                if (
                        (pos - start) < args.binding_width or
                        (end - pos) < args.binding_width
                ):
                    EFF = "PRO"
                thisres = Result(chrom, pos, "-", "-", 0, "-", 0, EFF)
                sys.stdout.write("%s\t%i\t%s\t%s\t%i\t%s\t%i\t%s\n" % thisres)



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
    counter = 1
    for i, v in enumerate(vcf_reader):
        # if (i % 1000) == 0:
        #     sys.stderr.write(str(i) + " ")
        if v.POS > 100000000:
            break
        # if v.CHROM not in vcf_data.keys():
        #     vcf_data[v.CHROM] = []
        while counter != v.POS and counter < v.POS:
            vcf_data[v.CHROM].append([counter, "-", "-", False])
            counter = counter + 1
        assert counter == v.POS, "error syncing counters"
        vcf_data[v.CHROM].append([v.POS, v.REF, v.ALT, True])
        counter = counter + 1

    last_gene_end = 0
    # first process all the coding sequences, then hit the remaining intergenic loci
    with gbk_open_fun(args.gbk, "r") as ingbk:
        for rec in SeqIO.parse(ingbk, "genbank"):
            thischrom = rec.id.split(".")[0]
            sys.stderr.write("Processing %s\n" % thischrom)
            for feat in rec.features:
                if feat.type not in ["source"]:
                    # process coding region
                    process_region(args, vcf_data,
                                   chrom=thischrom,
                                   start=feat.location.start,
                                   end=feat.location.end,
                                   rec=rec,
                                   strand=feat.strand,
                                   is_locus=True)
                    # process previous intergenic region
                    process_region(args, vcf_data,
                                   chrom=thischrom,
                                   start=last_gene_end,
                                   end=feat.location.start,
                                   rec=rec,
                                   strand=feat.strand,
                                   is_locus=False)
                    # keep track of where the last gene ended
                    last_gene_end = feat.location.end

if __name__ == '__main__':
    main()
