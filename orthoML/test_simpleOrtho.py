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
import pandas as pd
from Bio.Blast.Applications import NcbiblastxCommandline
from Bio.Blast.Applications import NcbitblastnCommandline


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
    # parser.add_argument("-t", "--blast_type",
    #                     help="blastn or tblastx", default="tblastx")
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


def make_prot_nuc_recip_blast_cmds(
        query_list, date,
        output, subject_file=None, logger=None):
    """given a file, make a blast cmd, and return path to output csv
    Only works is query_list is nucleotide and subject_file is protein
    """
    assert logger is not None, "must use logging"
    logger.info("Creating protein BLAST database")
    db_dir = os.path.join(output,
                          os.path.splitext(os.path.basename(subject_file))[0])
    os.makedirs(db_dir)
    protdb = os.path.join(db_dir,
                          os.path.splitext(os.path.basename(subject_file))[0])

    setup_blast_db(input_file=args.db_aa, input_type="fasta", dbtype="prot",
                   out=protdb, logger=logger)
    blast_cmds = []
    blast_outputs = []
    recip_blast_outputs = []
    for f in query_list:
        # run forward, nuc aganst prot, blast
        output_path_tab = str(
            os.path.join(output, date) + "_simpleOrtho_results_" +
            os.path.basename(f) + "_vs_protdb.tab")
        blast_cline = NcbiblastxCommandline(query=f,
                                            db=protdb, evalue=.001,
                                            outfmt=6, out=output_path_tab)
        add_params = str(" -num_threads 1 -num_alignments 20")
        blast_command = str(str(blast_cline) + add_params)
        blast_cmds.append(blast_command)
        blast_outputs.append(output_path_tab)
        # run reverse, prot against nuc, blast
        recip_output_path_tab = str(
            os.path.join(output, date) + "_simpleOrtho_results_" +
            "prot_vs_" + os.path.basename(f) + ".tab")
        recip_blast_cline = NcbitblastnCommandline(query=subject_file,
                                                   subject=f,
                                                   evalue=.001,
                                                   outfmt=6, out=recip_output_path_tab)
        recip_blast_command = str(str(recip_blast_cline) + add_params)
        blast_cmds.append(recip_blast_command)
        recip_blast_outputs.append(recip_output_path_tab)

    return(blast_cmds, blast_outputs, recip_blast_outputs)


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


def BLAST_tab_to_df(path):
    colnames = ["query_id", "subject_id", "identity_perc", "alignment_length",
                "mismatches", "gap_opens", "q_start", "q_end", "s_start",
                "s_end", "evalue", "bit_score"]
    raw_csv_results = pd.read_csv(
        open(path), comment="#", sep="\t", names=colnames)
    return raw_csv_results


def filter_recip_BLAST_df(df1, df2, logger=None):
    """ results from pd.read_csv with default BLAST output 6 columns
    df1 must be genomes against genes, and df2 must be genes against genomes,
    because we have to split the names so all all the contigs are recognized
    as coming from one genome.  returns a df
    """
    assert logger is not None, "must use a logger"
    logger.debug("shape of blast results")
    logger.debug("shape of recip blast results")
    df1['genome'] = df1.query_id.str.split('_').str.get(0)
    df2['genome'] = df2.subject_id.str.split('_').str.get(0)
    logger.debug(df1.shape)
    logger.debug(df2.shape)
    # recip structure
    filtered = pd.DataFrame(columns=df1.columns)
    unq_subject = df1.subject_id.unique()
    unq_query = df1.genome.unique()
    recip_hits = []
    nonrecip_hits = []
    for gene in unq_subject:
        for genome in unq_query:
            logger.debug("Checking %s in %s for reciprocity" % (gene, genome))
            tempdf1 = df1.loc[(df1["subject_id"] == gene) &
                              (df1["genome"] == genome), ]
            tempdf2 = df2.loc[(df2["query_id"] == gene) &
                              (df2["genome"] == genome), ]
            if tempdf1.empty or tempdf2.empty:
                logger.info("skipping %s in %s", gene, genome)
            else:
                subset1 = tempdf1.loc[(tempdf1["identity_perc"] > 90) &
                                      (tempdf1["bit_score"] == tempdf1["bit_score"].max())]
                subset2 = tempdf2.loc[(tempdf2["identity_perc"] > 90) &
                                      (tempdf2["bit_score"] == tempdf2["bit_score"].max())]
                logger.debug("grouped df shape: ")
                logger.debug(tempdf1.shape)
                logger.debug("grouped df2 shape: %s")
                logger.debug(tempdf2.shape)
                if subset1.empty or subset2.empty:
                    logger.info("No reciprocol hits for %s in %s", gene, genome)
                    nonrecip_hits.append([gene, genome])
                else:
                    # logger.debug(tempdf1)
                    # logger.debug("tempdf2")
                    # logger.debug(tempdf2)
                    # logger.debug("subset1")
                    # logger.debug(subset1)
                    # logger.debug("subset2")
                    # logger.debug(subset2)
                    if subset1.iloc[0]["query_id"] == subset2.iloc[0]["subject_id"]:
                        recip_hits.append([gene, genome])
                        filtered = filtered.append(subset1)
                        logger.info("Reciprocol hits for %s in %s!", gene, genome)
                    else:
                        nonrecip_hits.append([gene, genome])
                        logger.info("No reciprocol hits for %s in %s", gene, genome)

            # logger.debug(subset.shape)
    logger.debug("Non-reciprocal genes:")
    logger.debug(nonrecip_hits)
    logger.debug("Reciprocal genes:")
    logger.debug(recip_hits)
    logger.debug("filtered shape:")
    logger.debug(filtered.shape)
    return(filtered)

    # idx = csv_results.groupby(
    #     ['subject_id'])['identity_perc'].transform(max) == csv_results['identity_perc']
    # best_hits_df = csv_results[idx]
    # logger.debug("Shape of best hits:")
    # logger.debug(best_hits_df.shape)
    # new_df = best_hits_df.query('identity_perc >= 95 and evalue < .0001')
    # # new_df = best_hits_df.query('identity_perc >= 95')
    # return new_df


def write_pipe_extract_cmds(df, outfile, logger=None):
    #% parse output
    assert logger is not None, "must use a logger"
    logger.debug("cleaning up the csv output")
    with open(outfile, "a") as outf:
        for index, row in df.iterrows():
            if row['q_start'] > row['q_end']:
                logger.debug("hit is on the (-) strand")
                line = "{0}-RC@{1} :{2}:{3}".format(
                    row['subject_id'],
                    row['query_id'],
                    int(row['q_end']),
                    int(row['q_start']))
            else:
                line = "{0}@{1} :{2}:{3}".format(
                    row['subject_id'],
                    row['query_id'],
                    int(row['q_start']),
                    int(row['q_end']))
            sys.stdout.write(line + "\n")
            outf.write(line + "\n")


def main(args):
    logger = logging.getLogger('simpleOrtho')
    logging.basicConfig(stream=sys.stderr, level=logging.DEBUG,
                        format='%(name)s (%(levelname)s): %(message)s')
    logger.setLevel(logging.DEBUG)
    logger.debug("All settings used:")
    output_root = os.path.abspath(os.path.expanduser(args.output))
    for k, v in sorted(vars(args).items()):
        logger.debug("{0}: {1}".format(k, v))
    if not os.path.isdir(output_root):
        logger.info("creating output directory %s", output_root)
        os.makedirs(output_root)
    else:
        logger.error("Output Directory already exists!")
        sys.exit(1)
    date = str(datetime.datetime.now().strftime('%Y%m%d'))

    genomes = get_complete_paths_of_files(args.genomes_dir)

    commands, paths_to_outputs, paths_to_recip_outputs = \
        make_prot_nuc_recip_blast_cmds(
            query_list=genomes,
            subject_file=args.db_aa,
            output=output_root, date=date,
            logger=logger)
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
    merged_tab = os.path.join(output_root,
                              "merged_results.tab")
    recip_merged_tab = os.path.join(output_root,
                                    "recip_merged_results.tab")
    merge_outfiles(filelist=paths_to_outputs,
                   outfile_name=merged_tab)
    merge_outfiles(filelist=paths_to_recip_outputs,
                   outfile_name=recip_merged_tab)
    resultsdf = BLAST_tab_to_df(merged_tab)
    recip_resultsdf = BLAST_tab_to_df(recip_merged_tab)
    filtered_hits = filter_recip_BLAST_df(df1=resultsdf, df2=recip_resultsdf, logger=logger)
    write_pipe_extract_cmds(outfile=os.path.join(output_root, "simpleOrtho_regions.txt"),
                            df=filtered_hits, logger=logger)

if __name__ == '__main__':
    assert ((sys.version_info[0] == 3) and (sys.version_info[1] >= 5)), \
        "Must use python3.5 or higher!"
    args = get_args()
    main(args=args)
