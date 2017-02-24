#!/usr/bin/env python3
"""
"""
import os
import sys
import shutil
import datetime
import subprocess
import argparse
import multiprocessing
import logging
from operator import itemgetter
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio import AlignIO


def get_args(DEBUG=False):
    parser = argparse.ArgumentParser(
        description="This does some simple blasting to get region of interest")
    parser.add_argument("big_fasta",
                        help="")
    parser.add_argument("-o", "--output", dest='output',
                        help="directory in which to place the output files",
                        default=os.path.join(os.getcwd(), "orthoMSA"))
    parser.add_argument("-s", "--score_min",
                        help="not currently used; will be used to determinine \
                        an optional scoring threshold")
    # parser.add_argument("--prank_exe", dest="prank_exe",
    #                     action="store", default="prank",
    #                     help="Path to PRANK executable; " +
    #                     "default: %(default)s")
    parser.add_argument("--mafft_exe", dest="mafft_exe",
                        action="store", default="mafft",
                        help="Path to MAFFT executable; " +
                        "default: %(default)s")
    parser.add_argument("-t", "--blast_type",
                        help="blastn or tblastx", default="tblastx")
    args = parser.parse_args()
    return(args)


def prepare_mafft_cmd(outdir, combined_fastas, mafft_exe,
                      add_args="", outfile_name="best_MSA",
                      logger=None):
    """returns command line for constructing MSA with
    mafft and the path to results file
    """
    assert logger is not None, "Must use logger"
    if not os.path.exists(outdir):
        raise FileNotFoundError("output directory not found!")
    mafft_cmd = "{0} {1} {2} > {3}".format(
        mafft_exe, add_args, combined_fastas,
        os.path.join(outdir, outfile_name))
    logger.debug("MAFFT command: \n %s", mafft_cmd)
    return (mafft_cmd, os.path.join(outdir, outfile_name))


def make_msa(msa_tool, unaligned_seqs, outname, prank_exe, mafft_exe,
             args, outdir, logger=None):
    """returns msa cmd and results path
    """
    assert logger is not None, "Must use logger"
    if shutil.which(mafft_exe) is not None:
        msa_cmd, results_path = prepare_mafft_cmd(
            outdir=outdir,
            outfile_name=outname,
            combined_fastas=unaligned_seqs,
            mafft_exe=shutil.which(mafft_exe),
            add_args=args,
            logger=logger)
    else:
        raise ValueError("Construction of MSA skipped because " +
                         "%s is not a valid executable!", mafft_exe)
    return(msa_cmd, results_path)


def process_msa_list(msa_list, name=""):
    df_list = []
    for msa in msa_list:
        print(msa)
        alignment = AlignIO.read(msa, "fasta")
        names = [x.id for x in alignment]
        align_df = pd.DataFrame(np.array(alignment, dtype=str),
                                # rows=range(1, len(alignment)),
                                index=[x.id for x in alignment])
        # >gene=RpoS_NP_417221.1-RC@Lys36_c1 :96444:97433
        genomes = [x.split("@")[1].split("_")[0] for x in names]
        genes = [x.split("@")[0].split("_")[0] for x in names]
        # print(genomes)
        # print(genes)
        assert all([genes[0] == x for x in genes[1:]]), "Alignment must be of the same gene!"
        # print("passed assertion")
        transp = align_df.T
        transp.columns = genomes
        transp["gene"] = genes[0]
        transp["gene_index"] = range(0, len(transp))
        # print(transp.head())
        df_list.append(transp)
    print(df_list[0].head())
    for x in df_list:
        print(x.shape)
        dups = x.columns[x.columns.duplicated()]
        x.drop(dups, 1, inplace=True)
        print(x.shape)
        print("dropped dup col: %s" % "\n".join(dups))
    df_big = pd.concat(df_list)
    print(df_big.shape)
    print(df_big.dtypes)
    print(df_big.head())
    for col in df_big.columns:
        df_big[col] = df_big[col].astype(str)
        print(df_big.head())
        df_big.to_csv(os.path.join(
            args.output, "combined_" + name + "_msas.tab"), sep="\t")


def main(args):
    logger = logging.getLogger('prepareOrthoMSA')
    logging.basicConfig(stream=sys.stderr, level=logging.DEBUG,
                        format='%(name)s (%(levelname)s): %(message)s')
    logger.setLevel(logging.DEBUG)
    logger.debug("All settings used:")
    for k, v in sorted(vars(args).items()):
        logger.debug("{0}: {1}".format(k, v))
    if not os.path.isdir(args.output):
        logger.info("creating output directory %s", args.output)
        os.makedirs(args.output)
    else:
        logger.error("Output Directory already exists!")
        sys.exit(1)
    date = str(datetime.datetime.now().strftime('%Y%m%d'))
    output_root = os.path.abspath(os.path.expanduser(args.output))
    with open(args.big_fasta, 'r') as fasta:
        records = list(SeqIO.parse(fasta, 'fasta'))

    genome_entries = {}
    gene_entries = {}
    for rec in records:
        # genomes = [x.split("@")[1].split("_")[0] for x in names]
        # genes = [x.split("@")[0] for x in names]
        genome = rec.id.split("@")[1].split("_")[0]
        gene = rec.id.split("@")[0].replace("gene=", "").split("_")[0]
        if genome in genome_entries:
            genome_entries[genome].append((gene, rec))
        else:
            genome_entries[genome] = [(gene, rec)]
        if gene in gene_entries:
            gene_entries[gene].append((genome, rec))
        else:
            gene_entries[gene] = [(genome, rec)]
    for key, items in genome_entries.items():
        outpath = os.path.join(args.output, key + "_regions.fasta")
        with open(outpath, "w") as outp:
            for item in items:
            # for item in sorted(items):
                # print(item[1].seq.translate())
                # SeqIO.write(SeqRecord(item[1].seq.translate(),
                #                       id=item[1].id), outp, "fasta")
                SeqIO.write(item[1], outp, "fasta")
    msas = []
    for key, items in gene_entries.items():
        outpath = os.path.join(output_root, key + "_regions.fasta")
        with open(outpath, "w") as outp:
            for item in items:
                SeqIO.write(item[1], outp, "fasta")
        msa_cmd, results_path = make_msa(msa_tool="mafft",
                                         unaligned_seqs=outpath,
                                         prank_exe='',
                                         args='',
                                         outname=os.path.join(
                                             output_root,
                                             key + "_alligned_regions.fasta"),
                                         mafft_exe=args.mafft_exe,
                                         outdir=output_root,
                                         logger=logger)

        logger.info("Running %s for MSA of %s", "mafft", key)
        subprocess.run(msa_cmd,
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)
        msas.append(results_path)
    ######
    prot_msas = []
    for key, items in gene_entries.items():
        outpath_prot = os.path.join(output_root, key + "_regions.faa")
        with open(outpath_prot, "w") as outprot:
            for item in items:
                SeqIO.write(SeqRecord(item[1].seq.translate(),
                                      id=item[1].id), outprot, "fasta")

                # SeqIO.write(item[1], outp, "fasta")
        msa_cmd_prot, results_path_prot = make_msa(msa_tool="mafft",
                                                   unaligned_seqs=outpath_prot,
                                                   prank_exe='',
                                                   args='',
                                                   outname=os.path.join(
                                                       output_root,
                                                       key + "_alligned_regions.faa"),
                                                   mafft_exe=args.mafft_exe,
                                                   outdir=output_root,
                                                   logger=logger)

        logger.info("Running %s for MSA of %s", "mafft", key)
        subprocess.run(msa_cmd_prot,
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)
        prot_msas.append(results_path_prot)
    #####
    process_msa_list(msa_list=msas, name="nuc")
    process_msa_list(msa_list=prot_msas, name="prot")


if __name__ == '__main__':
    assert ((sys.version_info[0] == 3) and (sys.version_info[1] >= 5)), \
        "Must use python3.5 or higher!"
    args = get_args()
    main(args=args)
