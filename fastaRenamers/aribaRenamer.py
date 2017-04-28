#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
version 0.1
"""
import argparse
import sys
import logging
import os
import shutil
import subprocess

from Bio import SeqIO
from Bio import Entrez


def get_args():
    parser = argparse.ArgumentParser(
        description="given a protein fasta file (from PHASTER, writes out a " +
        "sanitized nuc fasta and a tsv file with te annotations")
    parser.add_argument("-i", "--infile", dest="infile",
                        help="input protein .fa", required=True)
    # parser.add_argument("-c", "--clobber", help="overwrite existing outfile",
    #                     dest="clobber",
    #                     default=False,
    #                     action="store_true")
    parser.add_argument("-s", "--start_from", action="store", dest="start_from",
                        help="if you get interupted, restart from this index",
                        default=0, type=int)
    parser.add_argument("-e", "--email", action="store", dest="email",
                        help="email address",
                        default='joe_smith@mail.gov', type=str)
    parser.add_argument("-o", "--outpre",
                        help=" output file prefix (no extension)",
                        dest="outpre",
                        default=os.path.join(os.getcwd(),
                                             "prokka_Uniprot_renamed.fa"))
    args = parser.parse_args()
    return(args)


def fetch_and_write_seqs(accessions, destination, region, db='nucleotide',
                         outfmt='fasta', concat=False, logger=None):
    """
    now should be able to detect and handle ftp calls to NCBI only
    takes an accession list
    if first item listed is an ftp address, it will run FTP calls for all the
    rest of the items. if not, the files will be fetched from entrex and if
    concat, smooshed into one file. Returns 0 if succesful
    """
    if len(accessions) > 1 and region:
        logger.info("region will be ignored when using multiple accessions!")
    if outfmt == 'fasta':
        rettype = 'fasta'
    elif outfmt == "gb":
        rettype = "gbwithparts"
    else:
        logger.error("only supports fasta and gb")
        sys.exit(1)
    i = accessions[0]
    if region:
        outf = str(destination.strip() +
                   i + ":" + region + "." + outfmt)
        if len(region.strip().split(":")) != 2:
            logger.error("Region must be two colon-sparated integers!")
            sys.exit(1)
        reg_start, reg_end = region.strip().split(":")
        logger.debug("fetching %s from %s to %s as %s",
                     i, reg_start, reg_end, outfmt)
        out_handle = open(outf, "w")
        sequence_handle = Entrez.efetch(
            db=db, id=i, rettype=rettype,
            retmode="text",
            seq_start=reg_start,
            seq_stop=reg_end)
    else:
        outf = str(destination.strip() +
                   i + "." + outfmt)
        logger.debug("fetching %s as %s", i, outfmt)
        out_handle = open(outf, "w")
        sequence_handle = Entrez.efetch(
            db=db, id=i, rettype=rettype, retmode="text")
    for line in sequence_handle:
        out_handle.write(line)
    out_handle.close()
    sequence_handle.close()
    return outf


def getAccessionAndCoords(rec, logger):
    """ use efetch to get the nucleotide sequence from a protein genbank.
relies on having curl and grep, so its pretty janky.
    """
    identifier = rec.id.split("|")[3]
    logger.debug("Finding  coding sequence associated with %s", identifier)
    cmd = 'curl -L "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id={0}&rettype=gp" | grep coded_by'.format(identifier)
    # {0}".format(shutil.which(backtranseq))
    try:
        res = subprocess.run(cmd,
                             shell=sys.platform != "win32",
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             check=True)
    except:
        logger.warning("Error trying to get coords for %s using id %s" % (
            rec, identifier))
        logger.warning("with the command:")
        logger.warning(cmd)
        return None, None, None
    outline = res.stdout.decode("utf-8").replace(" ", "")
    # /coded_by="complement(NC_004337.2:1806396..1807301)"
    outline = outline.split('"')[1]
    accession, outline = outline.split(':')
    accession = accession.replace("complement(", "")
    accession = accession.replace("join(", "")
    outline = outline.replace(')', "")
    outline = outline.replace('..', ":")
    return(identifier, accession, outline)


def main():
    logger = logging.getLogger('aribaRenamer')
    logging.basicConfig(stream=sys.stderr, level=logging.DEBUG,
                        format='%(name)s (%(levelname)s): %(message)s')
    logger.setLevel(logging.DEBUG)
    logger.debug("All settings used:")
    args = get_args()
    for k, v in sorted(vars(args).items()):
        logger.debug("{0}: {1}".format(k, v))
    outdir = args.outpre + "_results"
    name = os.path.basename(outdir)
    outfasta = os.path.join(outdir, name + ".fa")
    outtsv = os.path.join(outdir, name + ".tsv")
    outtemp = os.path.join(outdir, name + "_temp")
    logger.info("temp folder: %s", outtemp)
    Entrez.email = args.email
    if os.path.isdir(outdir):
        logger.error("Output file exists! exiting")
        sys.exit(1)
    os.makedirs(outdir)
    os.makedirs(outtemp)
    #  start processing the records
    entries = 0
    records = SeqIO.parse(args.infile, "fasta")
    for rec in records:
        entries = entries + 1
    tsv_lines = []
    fastas = []
    records = SeqIO.parse(args.infile, "fasta")
    for idx, i in enumerate(records):
        if idx < args.start_from:
            continue
        logger.debug("processing record %i/%i", idx, entries)
        # if idx < 5:
        #     break
        protacc, acc, region = getAccessionAndCoords(i, logger=logger)
        if not protacc:
            continue
        try:
            output_fasta = fetch_and_write_seqs(
                accessions=[acc],
                destination=os.path.join(outtemp, ""),
                region=region, db='nucleotide',
                outfmt='fasta', concat=False)
        except Exception as e:
            logger.warning("something went wrong while fetching seqs:")
            logger.warning(e)
            continue

        fastas.append(output_fasta)
        # again, see the spec on the ariba website
        tsv_lines.append("{0}\t{1}\t{2}\t{3}\n".format(
            # because some phages come from the same genome, we have to append the index to keep it unique
            protacc + "_" + str(idx),  # unique identifier
            1,  # 0 is noncoding, 1 is cds
            0,  # 0 is presence/absense, 1 is variant
            i.id + "_" + acc + ":" + region))  # description

    # write out tsv
    assert len(fastas) == len(tsv_lines), \
        "Unequal number of fastas and accessions. Did one fail?"
    with open(outtsv, "a") as outtsv:
        for line in tsv_lines:
            outtsv.write(line)
    with open(outfasta, "a") as outfasta:
        for idx, fasta in enumerate(fastas):
            with open(fasta, "r") as fhandle:
                seq = list(SeqIO.parse(fhandle, "fasta"))
                if len(seq) != 1:
                    logger.warning("multiple entries returned in error, " +
                                   "only writing out first")
                simplified_rec = seq[0]
                simplified_rec.description = ""
                simplified_rec.name = None
                simplified_rec.id = tsv_lines[idx].split("\t")[0]
                SeqIO.write(simplified_rec, outfasta, "fasta")
    logger.info("Finished")


if __name__ == '__main__':
    main()
