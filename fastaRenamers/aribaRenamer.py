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
import glob
import multiprocessing
import unittest

from Bio import SeqIO
from Bio import Entrez


def get_args():
    parser = argparse.ArgumentParser(
        description="given a protein fasta file (from PHASTER, writes out a " +
        "sanitized nuc fasta and a tsv file with te annotations")
    parser.add_argument("-i", "--infile", dest="infile",
                        help="input protein .fa", required=True)
    parser.add_argument("-c", "--cores", help="for multiprocessing. This " +
                        "influences how many ncbi calls get sent, so 2 " +
                        "is probably a good number.",
                        dest="cores",
                        type=int, default=1)
    parser.add_argument("-s", "--start_from", action="store",
                        dest="start_from",
                        help="if you get interupted, restart from this index",
                        default=0, type=int)
    parser.add_argument("-j", "--just_aggregate", action="store_true",
                        dest="just_aggregate",
                        help="if you get interupted, after fetching all the " +
                        "accession, rerun with this to do the name " +
                        "conversion and aggregations for all the files " +
                        " ending in *_status",
                        default=False)
    parser.add_argument("-e", "--email", action="store", dest="email",
                        help="email address",
                        required=True, default='', type=str)
    parser.add_argument("-o", "--outdir",
                        help=" output direcrory",
                        dest="outdir",
                        required=True)
    parser.add_argument("-u", "--unittest",
                        help=" run unittest and exit",
                        dest="unittest",
                        action="store_true")
    args = parser.parse_args()
    return(args)


def getAccessionAndCoords(rec, idx, nseqs, outdir):
    """ use efetch to get the nucleotide sequence from a protein genbank.
relies on having curl and grep, so its pretty janky.
    return codes :
    0 all good
    1 error with curl
    2 error parsing
    """
    identifier = rec.id.split("|")[3]
    resfile = os.path.join(outdir, str(idx) + "_status")
    print("Getting Accessions for %s, item %d of %d" % (identifier, idx + 1, nseqs))
    cmd = 'curl -L "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id={0}&rettype=gp" | grep coded_by'.format(identifier)
    try:
        res = subprocess.run(cmd,
                             shell=sys.platform != "win32",
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             check=True)
    except:
        with open(resfile, "w") as resout:
            resout.write("{0}\t{1}\t{2}".format(
                identifier, "NULL", "could not get accession"))
        return 1
    outline = res.stdout.decode("utf-8").replace(" ", "")
    # /coded_by="complement(NC_004337.2:1806396..1807301)"
    outline = outline.split('"')[1]
    try:
        accession, outline = outline.split(':')
    except:
        with open(resfile, "w") as resout:
            resout.write("{0}\t{1}\t{2}".format(
                identifier,
                "NULL",
                str("could not split : for  %s" % outline)))
        return 2
    accession = accession.replace("complement(", "")
    accession = accession.replace("join(", "")
    outline = outline.replace(')', "")
    outline = outline.replace('..', ":")
    with open(resfile, "w") as resout:
        resout.write("{0}\t{1}\t{2}".format(
            identifier, accession, outline))
    return 0


def fetch_and_write_fasta(accessions, destination, region, db='nucleotide',
                          concat=False):
    """
    now should be able to detect and handle ftp calls to NCBI only
    takes an accession list
    if first item listed is an ftp address, it will run FTP calls for all the
    rest of the items. if not, the files will be fetched from entrex and if
    concat, smooshed into one file. Returns 0 if succesful
    return codes:
    0) all is well
    1) input is messed up
    2) error fetching sequence
    3) error on output
    """
    rettype = 'fasta'
    i = accessions
    if region:
        outf = str(destination.strip() +
                   i + ":" + region + ".fasta")
        if len(region.strip().split(":")) != 2:
            return 1
        reg_start, reg_end = region.strip().split(":")
        try:
            out_handle = open(outf, "w")
        except Exception as e:
            print(e)
            return 3
        try:
            sequence_handle = Entrez.efetch(
                db=db, id=i, rettype=rettype,
                retmode="text",
                seq_start=reg_start,
                seq_stop=reg_end)
        except Exception as e:
            print(e)
            out_handle.close()
            return 2
    else:
        outf = str(destination.strip() +
                   i + ".fasta")
        try:
            out_handle = open(outf, "w")
        except Exception as e:
            print(e)
            return 3
        try:
            sequence_handle = Entrez.efetch(
                db=db, id=i, rettype=rettype, retmode="text")
        except Exception as e:
            print(e)
            out_handle.close()
            return 2
    try:
        for line in sequence_handle:
            out_handle.write(line)
        out_handle.close()
    except Exception as e:
        print(e)
        sequence_handle.close()
        out_handle.close()
        return 3
    sequence_handle.close()
    return 0


# import time


# def RateLimited(maxPerSecond):
#     minInterval = 1.0 / float(maxPerSecond)

#     def decorate(func):
#         lastTimeCalled = [0.0]

#         def rateLimitedFunction(*args, **kargs):
#             elapsed = time.clock() - lastTimeCalled[0]
#             leftToWait = minInterval - elapsed
#             if leftToWait > 0:
#                 time.sleep(leftToWait)
#             ret = func(*args, **kargs)
#             lastTimeCalled[0] = time.clock()
#             return ret
#         return rateLimitedFunction
#     return decorate


def main(args):
    logger = logging.getLogger('aribaRenamer')
    logging.basicConfig(stream=sys.stderr, level=logging.DEBUG,
                        format='%(name)s (%(levelname)s): %(message)s')
    logger.setLevel(logging.DEBUG)
    logger.debug("All settings used:")
    for k, v in sorted(vars(args).items()):
        logger.debug("{0}: {1}".format(k, v))
    # outdir = args.outpre + "_results"
    output_root = os.path.abspath(os.path.expanduser(args.outdir))
    bigoutfasta = os.path.join(output_root, "results.fa")
    outtsv = os.path.join(output_root, "results.tsv")
    outerr = os.path.join(output_root, "errors.tsv")
    outtemp_get_accessions = os.path.join(output_root, "temp_get_accessions")
    outtemp = os.path.join(output_root, "temp")
    logger.info("temp folder: %s", outtemp)
    Entrez.email = args.email
    if not args.just_aggregate:
        if os.path.isdir(output_root):
            logger.error("Output file exists! exiting")
            sys.exit(1)
        os.makedirs(output_root)
        os.makedirs(outtemp)
        os.makedirs(outtemp_get_accessions)
        #  start processing the records
    entries = 0
    records = SeqIO.parse(args.infile, "fasta")
    for rec in records:
        entries = entries + 1
    failed_records = []
    if not args.just_aggregate:
        logger.info("Getting the nucleotide accessions associated with " +
                    "the protein IDs.  This can take a while.  " +
                    "The best way to check progress is by the number " +
                    "of output files in the temp dir")
        records = SeqIO.parse(args.infile, "fasta")
        # run all the getAccessionCoords
        pool = multiprocessing.Pool(processes=args.cores)
        results = [
            pool.apply_async(getAccessionAndCoords,
                             (rec,),
                             {"outdir": outtemp_get_accessions,
                              "nseqs": entries,
                              "idx": idx})
            for idx, rec in enumerate(records) if idx > args.start_from]
        pool.close()
        pool.join()
        logger.info("Sum of return codes (should be 0):")
        logger.info(sum([r.get() for r in results]))
    # get all the results files
    acc_results_files = glob.glob(os.path.join(outtemp_get_accessions, "") +
                                  "*_status")
    acc_results = {}
    for acc_res in acc_results_files:
        with open(acc_res, "r") as f:
            lines = f.read().split('\n')
            if len(lines) != 1:
                logger.error("how did we get a multiline result file for %s?",
                             acc_res)
            protacc, acc, region = lines[0].strip().split("\t")
            if "NULL" in acc:
                logger.warning("Error getting coordinates for accession %s:",
                               protacc)
                logger.warning(region)
                failed_records.append(
                    [protacc, "error fetching or parsing protein accession"])
                continue
            acc_results[protacc] = [acc, region]
    logger.info("Fetching the nucleotide sequences.  This can take a while." +
                "  The best way to check progress is to count the " +
                "output files")
    # run all the getAccessionCoords
    pool2 = multiprocessing.Pool(processes=args.cores)
    results_fetch = [
        pool2.apply_async(fetch_and_write_fasta,
                          (rec,),
                          {"destination": os.path.join(outtemp, ""),
                           "region": reg,
                           "db": 'nucleotide',
                           "concat": False})
        for rec, reg in acc_results.values()]
    pool2.close()
    pool2.join()
    logger.info("Sum of return codes from fetching " +
                "nucleotide sequences (should be 0):")
    logger.info(sum([r.get() for r in results_fetch]))

    fetched_files = glob.glob(str(os.path.join(outtemp, "") + "*fasta"))
    logger.info("Writing out final results ")
    idx = 0
    idxs = len(acc_results)
    for protacc, recreg in acc_results.items():
        # check for output fasta
        idx = idx + 1
        logger.info("Processing %s / %s:%s, item %d of %d",
                    protacc, recreg[0], recreg[1], idx, idxs)
        target_fasta = "{0}{1}:{2}.fasta".format(
            os.path.join(outtemp, ""),
            recreg[0], recreg[1])
        if target_fasta not in fetched_files:
            logger.warning("No data could be fetched for protein %s / nuc %s",
                           protacc, recreg[0])
            failed_records.append(
                [protacc,
                 "error getting nucleotide sequences from {0}:{1}".format(
                     recreg[0], recreg[1])])
            continue
        tsv_line = "{0}\t{1}\t{2}\t{3}\n".format(
            # because some phages come from the same genome,
            # we have to append the index to keep it unique
            protacc + "_" + str(idx),  # unique identifier
            1,  # 0 is noncoding, 1 is cds
            0,  # 0 is presence/absense, 1 is variant
            recreg[0] + ":" + recreg[1])  # description
            # i.id + "_" + recreg[0] + ":" + recreg[1])  # description
        try:
            with open(bigoutfasta, "a") as ofasta:
                with open(target_fasta, "r") as fhandle:
                    seq = list(SeqIO.parse(fhandle, "fasta"))
                    if len(seq) != 1:
                        logger.warning("multiple entries returned in %s",
                                       target_fasta)
                        logger.warning("only writing out first")
                    simplified_rec = seq[0]
                    simplified_rec.description = ""
                    simplified_rec.name = None
                    simplified_rec.id = tsv_line.split("\t")[0]
                    SeqIO.write(simplified_rec, ofasta, "fasta")
        except Exception as e:
            logger.error("error with writing out %s with the proper header",
                         target_fasta)
            logger.error(e)
            failed_records.append(tsv_line)
            continue
        with open(outtsv, "a") as otsv:
            otsv.write(tsv_line)

    with open(outerr, "a") as oerr:
        for errec in failed_records:
            oerr.write("{0}\t{1}".format(errec[0], errec[1]))
    logger.info("Finished!")
    logger.info("Check the errors.tsv files to figure out what happened " +
                "with any accessions skipped")


if __name__ == '__main__':
    args = get_args()

    if args.unittest:
        unittest.main()
    else:
        main(args)


class aribaRenamer(unittest.TestCase):
    """ test for ariba
    """
    def setUp(self):
        self.testfile = os.path.join(os.path.dirname(__file__),
                                     "ariba_test.db")
        self.testdir = os.path.dirname(__file__)
        with open(self.testfile, "r") as t:
            self.reclist = list(SeqIO.parse(t, "fasta"))
        self.acc_dict = {'NP_708230.1': ['NC_004337.2', '2476842:2477336'],
                         'NP_707631.1': ['NC_004337.2', '1806396:1807301'],
                         'NP_752289.1': ['NC_004431.1', '326209:330324'],
                         'NP_752209.1': ['NC_004431.1', '253783:254469']}

    def test_getAccessions(self):
        """ Test all three return codes for dataset
        """
        # return code 0

        code0 = getAccessionAndCoords(rec=self.reclist[0],
                                      idx=1, nseqs=6, outdir=self.testdir)
        self.assertEqual(0, code0)
        # return code 1
        code1 = getAccessionAndCoords(rec=self.reclist[1],
                                      idx=1, nseqs=6, outdir=self.testdir)
        self.assertEqual(1, code1)
        # return code 2
        code2 = getAccessionAndCoords(rec=self.reclist[2],
                                      idx=1, nseqs=6, outdir=self.testdir)
        self.assertEqual(2, code2)
        os.unlink(os.path.join(self.testdir, "1_status"))

    def test_fetch(self):
        Entrez.email = "nickp60@hotmail.com",
        fcode0 = fetch_and_write_fasta(
            accessions=self.acc_dict["NP_708230.1"][0],
            destination=self.testdir,
            region=self.acc_dict["NP_708230.1"][1], db='nucleotide',
            concat=False)
        self.assertEqual(0, fcode0)

        fcode1 = fetch_and_write_fasta(
            accessions=self.acc_dict["NP_708230.1"][0],
            destination=self.testdir,
            region="1100:2200:3366",
            db='nucleotide',
            concat=False)
        self.assertEqual(1, fcode1)

        fcode2 = fetch_and_write_fasta(
            accessions=self.acc_dict["NP_708230.1"][0],
            destination=self.testdir,
            region=self.acc_dict["NP_708230.1"][1], db='porcupine',
            concat=False)
        self.assertEqual(2, fcode2)

        fcode3 = fetch_and_write_fasta(
            accessions=self.acc_dict["NP_708230.1"][0],
            destination=os.path.join(self.testdir, "porcupine_recipes", ""),
            region=self.acc_dict["NP_708230.1"][1], db='nucleotide',
            concat=False)
        self.assertEqual(3, fcode3)
