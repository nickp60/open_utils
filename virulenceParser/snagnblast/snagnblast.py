#!/usr/bin/env python
"""
branch of snagnblast version 0.5
version 0.9
Minor version changes:
- revampted to use multiprocessing.pool, added logging

 script to fetch genes from NCBI when given a file cintaining NCBI accession numbers,
 blast them against a local database (from makeblastdb), and write out the results
 as a csv.

USAGE:
 $ python snagnblast.py accessions.txt_or_accessions.csv /BLAST/directory/ /output/directory/
"""
import os
import sys
import datetime
import subprocess
import argparse
import multiprocessing
import logging
import pandas as pd
from Bio import SeqIO, Entrez
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbitblastxCommandline


def get_args(DEBUG=False):
    if DEBUG:
        args = argparse.Namespace(
            genelist=os.path.expanduser(
                "~/GitHub/FB/Ecoli_comparative_genomics/data/test_virgenes2.csv"),
            blastdb=os.path.expanduser("~/BLAST/env_Coli"),
            output=os.path.expanduser("~/GitHub/FB/Ecoli_comparative_genomics/results/1/"),
            score_min=70,
            blasttype="tblastx",
            cores=2,
            query_bin=20)
    else:
        parser = argparse.ArgumentParser(description="This script takes a list of gene accessions \
        from either a text file or a csv, grabs the \
        sequencs from NCBI, and proceeds to use either \
        blastn or tblastx to detect the presence of the \
        genes in a custom database")
        parser.add_argument("genelist",
                            help="file containing gene accessions. if delimited, use \
                            the headers in the example file as a template")
        parser.add_argument("blastdb", help="blastdb of interest")
        parser.add_argument("-o", "--output", dest='output',
                            help="directory in which to place the output files",
                            default=os.path.join(os.getcwd(), "snagnblast_results"))
        parser.add_argument("-s", "--score_min",
                            help="not currently used; will be used to determinine \
        an optional scoring threshold")
        parser.add_argument("-t", "--blast_type", help="blastn or tblastx", default="tblastx")
        parser.add_argument("-c", "--cores", help="number of cores to use (at your own risk)",
                            type=int, default=2)
        # parser.add_argument("-b", "--query_bin", help="number of sequences to put per blast job",
        # type=int)

        args = parser.parse_args()
    return(args)


def get_accessions(genelist):
    print("reading in gene list")
    genes = open(genelist, "r")
    if genes.name.endswith("txt"):
        raise ValueError("only csv's taken as input now, sorry!")
    elif genes.name.endswith("tsv"):
    # if the input is tabular, accesions must be in the first column
        raise ValueError("only csv's taken as input now, sorry!")
    elif genes.name.endswith("csv"):
        print("gene list is a comma-deliminated file")
        genedf = pd.read_csv(genes, sep=",")
        accessions = genedf.iloc[0:, 0].tolist()
        accessions = [x for x in accessions if str(x) != 'nan']
    else:
        raise ValueError("Reading error: must be either a csv, tab, or txt file")
    return(accessions)


def run_entrez(accessions, output, date, logger):
    """ Grab sequences from NCBI, write out resulting fasta file
    """
    logger.info("\n\nFetching %i accessions from NCBI" % len(accessions))
    Entrez.email = "alfredTheDaring@gmail.com"
    try:
        sequence_handle = Entrez.efetch(db="nucleotide",
                                        id=accessions, rettype="fasta")
    except Exception as e:
        logger.error("Error fetching accessions from Entrez")
        raise ValueError(e)
    seqs = SeqIO.parse(sequence_handle, "fasta")
    # logger.info("writing results to %s" % str(os.path.join(output, date) +
    #                                           "_sequences.fa"))
    outfile = str(os.path.join(
        output, date) + "_sequences.fa")
    with open(outfile, "w") as fasta_output:
        count = SeqIO.write(seqs, fasta_output, "fasta")
    logger.info("wrote %i records to %s", count, outfile)
    return(outfile)


def parse_multifasta(accessions, entrez_file, output, date, logger):
    logger.info("Splitting multifasta")
    sequences_fasta = open(entrez_file, "r")
    entrez_results = list(SeqIO.parse(sequences_fasta, "fasta"))
    if(len(accessions) != len(entrez_results)):
        logger.warning("Warning! not all accessions were found!")
    if not os.path.isdir(os.path.join(output, "accession_fastas")):
        os.makedirs(os.path.join(output, "accession_fastas"))
    for rec in entrez_results:
        dest = os.path.join(output, str("accession_fastas" +
                                        os.path.sep +
                                        rec.id + ".fasta"))
        logger.debug("writing %s to %s", rec.id, dest)
        with open(dest, 'w') as out:
            SeqIO.write(rec, out, "fasta")

    sequences_fasta.close()
    return(os.path.join(output, "accession_fastas"))


def make_blast_cmds(filename_list, blast_type, output, blastdb, date):
    """given a file, make a blast cmd, and return path to output csv
    """
    blast_cmds = []
    blast_outputs = []
    for f in filename_list:
        output_path_tab = str(
            os.path.join(output, date) + "_dcmegablast_results_" +
            os.path.basename(f) + ".tab")
        if blast_type == 'blastn':
            blast_cline = NcbiblastnCommandline(query=f,
                                                db=blastdb, evalue=10,
                                                outfmt=6, out=output_path_tab)
            add_params = str(" -num_threads 1 -max_target_seqs " +
                             "2000 -task dc-megablast")
        elif blast_type == 'tblastx':
            blast_cline = NcbitblastxCommandline(query=f,
                                                 db=blastdb, evalue=10,
                                                 outfmt=6, out=output_path_tab)
            add_params = str(" -num_threads 1 -max_target_seqs 2000 " +
                             "-query_gencode 11 -db_gencode 11")
        else:
            raise ValueError("must use either blastn or tblastx")
        blast_command = str(str(blast_cline) + add_params)
        blast_cmds.append(blast_command)
        blast_outputs.append(output_path_tab)
    return(blast_cmds, blast_outputs)


def get_complete_paths_of_files(directory):
    filenames = []
    shortnames = [i for i in os.listdir(directory) if not os.path.isdir(os.path.join(directory, i))]
    for i in shortnames:
        filenames.append(os.path.join(directory, i))
    return(filenames)


def merge_outfiles(filelist, outfile_name):
    """
    """
    # only grab .tab files, ie, the blast output
    filelist = [i for i in filelist if i.split(".")[-1:] == ['tab']]
    if len(filelist) == 1:
        print("only one file found! no merging needed")
        return(filelist)
    else:
        print("merging all the blast results to %s" % outfile_name)
        nfiles = len(filelist)
        fout = open(outfile_name, "a")
        # first file:
        for line in open(filelist[0]):
            fout.write(line)
        #  now the rest:
        for num in range(1, nfiles):
            f = open(filelist[num])
            for line in f:
                fout.write(line)
            f.close()  # not really needed
        fout.close()
    return(outfile_name)


def cleanup_output_to_csv(infile, accession_pattern='(?P<accession>[A-Z _\d]*\.\d*)'):
    #% parse output
    print("cleaning up the csv output")
    colnames = ["query_id", "subject_id", "identity_perc", "alignment_length",
                "mismatches", "gap_opens", "q_start", "q_end", "s_start",
                "s_end", "evalue", "bit_score"]
    csv_results = pd.read_csv(
        open(infile), comment="#", sep="\t", names=colnames)
    #This regex will probably break things rather badly before too long...
    # it looks for capital letter and numbers, dot, number, ie SHH11555JJ8.99
    csv_results["accession"] = \
        csv_results.query_id.str.extract(accession_pattern)
    # write out results with new headers or with
    # new headers and merged metadat from accessions.tab
    genes = open(args.genelist, "r")
    genedf = pd.read_csv(genes, sep=",")
    output_path_csv = str(os.path.splitext(infile)[0] + ".csv")
    results_annotated = pd.merge(csv_results, genedf,
                                 how="left", on="accession")
    results_annotated.to_csv(open(output_path_csv, "w"))
    print("wrote final csv to %s" % output_path_csv)
#%%


def main():
    args = get_args(DEBUG=False)
    logger = logging.getLogger('snagnblast')
    logger.setLevel(logging.DEBUG)
    logger.debug("All settings used:")
    for k, v in sorted(vars(args).items()):
        logger.debug("{0}: {1}".format(k, v))
    if not os.path.isdir(args.output):
        logger.info("creating output directory %s", args.output)
        os.makedirs(args.output)
    else:
        logging.info("Using existing output_directory")
    date = str(datetime.datetime.now().strftime('%Y%m%d'))
    logger.info("reading in genelist: %s", args.genelist)
    genes = get_accessions(args.genelist)  # parse the input file
    logger.debug("Gene accessions found:")
    logger.debug(" ".join([x for x in genes]))
    try:
        entrez_outfile = run_entrez(genes, output=args.output, date=date,
                                    logger=logger)
    except ValueError as e:
        logger.error(e)
        sys.exit(1)

    accs_dir = parse_multifasta(genes, output=args.output,
                                entrez_file=entrez_outfile,
                                date=date, logger=logger)
    files = get_complete_paths_of_files(accs_dir)

    commands, paths_to_outputs = make_blast_cmds(
        filename_list=files, blast_type=args.blast_type,
        output=args.output, blastdb=args.blastdb, date=date)
    ###
    logger.info("running %s commands", args.blast_type)
    pool = multiprocessing.Pool(processes=args.cores)
    logger.debug("Running the following commands in parallel " +
                 "(this could take a while):")
    logger.debug("\n" + "\n".join([x for x in commands]))
    results = [
        pool.apply_async(subprocess.run,
                         (cmd,),
                         {"shell": sys.platform != "win32",
                          "stdout": subprocess.PIPE,
                          "stderr": subprocess.PIPE,
                          "check": True})
        for cmd in commands]
    pool.close()
    pool.join()
    reslist = []
    reslist.append([r.get() for r in results])

    merged_fasta = merge_outfiles(
        paths_to_outputs,
        os.path.join(args.output,
                     str(date + "_results_merged.csv")))
    cleanup_output_to_csv(merged_fasta)
    print("<insert summary here>. have a great day!")


if __name__ == '__main__':
    assert ((sys.version_info[0] == 3) and (sys.version_info[1] >= 5)), \
        "Must use python3.5 or higher!"
    main()
