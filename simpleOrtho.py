#!/usr/bin/env python
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
import pandas as pd
from Bio.Blast.Applications import NcbiblastxCommandline


def get_args(DEBUG=False):
    parser = argparse.ArgumentParser(
        description="This does some simple blasting to get region of interest")
    parser.add_argument("db_aa",
                        help="file containing gene accessions. if delimited," +
                        " use the headers in the example file as a template")
    parser.add_argument("genomes_dir", help="dir with and only iwth genomes")
    parser.add_argument("-o", "--output", dest='output',
                        help="directory in which to place the output files",
                        default=os.path.join(os.getcwd(), "simpleOrtho"))
    parser.add_argument("-s", "--score_min",
                        help="not currently used; will be used to determinine \
                        an optional scoring threshold")
    parser.add_argument("-t", "--blast_type",
                        help="blastn or tblastx", default="tblastx")
    args = parser.parse_args()
    return(args)


def setup_blast_db(input_file, input_type="fasta", dbtype="prot",
                   title="blastdb", out="blastdb",
                   makeblastdb_exe='', logger=None):
    """
    This runs make blast db with the given parameters
    requires logging, os, subprocess, shutil
    """
    if makeblastdb_exe == '':
        makeblastdb_exe = shutil.which("makeblastdb")
        if logger:
            logger.info("makeblastdb executable: %s", makeblastdb_exe)
    makedbcmd = str("{0} -in {1} -input_type {2} -dbtype {3} " +
                    "-out {4}").format(makeblastdb_exe,
                                       input_file,
                                       input_type,
                                       dbtype,
                                       out)
    if logger:
        logger.info("Making blast db: {0}".format(makedbcmd))
    try:
        subprocess.run(makedbcmd, shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE, check=True)
        logging.debug("BLAST database '{0}' created here: {1}".format(
            title, out))
        return 0
    except:
        if logger:
            logging.error("Something bad happened when trying to make " +
                          "a blast database")
        sys.exit(1)


def make_blast_cmds(filename_list, blast_type, output, blastdb, date):
    """given a file, make a blast cmd, and return path to output csv
    """
    blast_cmds = []
    blast_outputs = []
    for f in filename_list:
        output_path_tab = str(
            os.path.join(output, date) + "_simpleOrtho_results_" +
            os.path.basename(f) + ".tab")
        if blast_type == 'blastx':
            blast_cline = NcbiblastxCommandline(query=f,
                                                db=blastdb, evalue=10,
                                                outfmt=6, out=output_path_tab)
            add_params = str(" -num_threads 1 ")
        else:
            raise ValueError("must use blastx")
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
        # print("only one file found! no merging needed")
        return(filelist)
    else:
        # print("merging all the blast results to %s" % outfile_name)
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


def filter_BLAST_csv(csv_results, logger=None):
    """ results from pd.read_csv with default BLAST output 6 columns
    """
    assert logger is not None, "must use a logger"
    logger.debug("shape of blast results")
    logger.debug(csv_results.shape)
    idx = csv_results.groupby(['subject_id'])['identity_perc'].transform(max) == csv_results['identity_perc']
    best_hits_df = csv_results[idx]
    logger.debug("Shape of best hits:")
    logger.debug(best_hits_df.shape)
    new_df = best_hits_df.query('identity_perc >= 97 and evalue < .0001')
    return new_df


def cleanup_output_to_csv(infile, accession_pattern='(?P<accession>[A-Z _\d]*\.\d*)', logger=None):
    #% parse output
    assert logger is not None, "must use a logger"
    logger.debug("cleaning up the csv output")
    colnames = ["query_id", "subject_id", "identity_perc", "alignment_length",
                "mismatches", "gap_opens", "q_start", "q_end", "s_start",
                "s_end", "evalue", "bit_score"]
    raw_csv_results = pd.read_csv(
        open(infile), comment="#", sep="\t", names=colnames)
    # csv_results['genome'], csv_results['contig'] = zip(*csv_results['subject_id'].map(lambda x: x.split('_')))
    raw_csv_results['genome'] = raw_csv_results.subject_id.str.split('_').str.get(0)

    csv_results = filter_BLAST_csv(raw_csv_results, logger=logger)
    filtered_path = os.path.splitext(infile)[0] + "_filtered" + os.path.splitext(infile)[1]
    csv_results.to_csv(filtered_path, sep='\t')
    for index, row in csv_results.iterrows():
        if row['q_start'] > row['q_end']:
            logger.debug("hit is on the (-) strand")
            line = "{0}-RC@{1} :{2}:{3}".format(
                row['subject_id'],
                row['query_id'],
                row['q_end'],
                row['q_start'])
        else:
            line = "{0}@{1} :{2}:{3}".format(
                row['subject_id'],
                row['query_id'],
                row['q_start'],
                row['q_end'])
        sys.stdout.write(line + "\n")
    return(filtered_path)


def main(args):
    logger = logging.getLogger('simpleOrtho')
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
    logger.info("Creating protein BLAST database")
    db_dir = os.path.join(args.output,
                          os.path.splitext(os.path.basename(args.db_aa))[0])
    os.makedirs(db_dir)
    db = os.path.join(db_dir,
                      os.path.splitext(os.path.basename(args.db_aa))[0])

    setup_blast_db(input_file=args.db_aa, input_type="fasta", dbtype="prot",
                   out=db, logger=logger)

    genomes = get_complete_paths_of_files(args.genomes_dir)

    commands, paths_to_outputs = make_blast_cmds(
        filename_list=genomes, blast_type="blastx",
        output=args.output, blastdb=db, date=date)
    ###
    logger.info("running %s commands", args.blast_type)
    pool = multiprocessing.Pool()
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

    filtered_outputs = []
    for i in paths_to_outputs:
        filtered_path = cleanup_output_to_csv(i, logger=logger)
        filtered_outputs.append(filtered_path)


if __name__ == '__main__':
    assert ((sys.version_info[0] == 3) and (sys.version_info[1] >= 5)), \
        "Must use python3.5 or higher!"
    args = get_args()
    main(args=args)
