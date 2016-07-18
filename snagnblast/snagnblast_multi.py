#!/usr/bin/env python
"""
branch of snagnblast version 0.5
version 0.7
Minor version changes:
 - pep8 kinda? I think?


 script to fetch genes from NCBI when given a file cintaining NCBI accession numbers,
 blast them against a local database (from makeblastdb), and write out the results
 as a csv.

USAGE:
 $ python snagnblast.py accessions.txt_or_accessions.csv /BLAST/directory/ /output/directory/
"""
import os
#import sys
import datetime
import subprocess
import argparse
from Bio import SeqIO, Entrez
from multiprocessing import Process, Queue, current_process, cpu_count
import pandas as pd
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbitblastxCommandline
import md5
DEBUG = False
#%%
#define inputs
if DEBUG:
    genelist = os.path.expanduser("~/GitHub/FB/Ecoli_comparative_genomics/data/test_virgenes2.csv")
    blastdb = os.path.expanduser("~/BLAST/env_Coli")
    output = os.path.expanduser("~/GitHub/FB/Ecoli_comparative_genomics/results/1/")
    score_min = 70
    blasttype = "tblastx"
    cores = 2
    query_bin = 20
else:
    parser = argparse.ArgumentParser(description="This script takes a list of gene accessions \
                                                 from either a text file or a csv, grabs the \
                                                 sequencs from NCBI, and proceeds to use either \
                                                 blastn or tblastx to detect the presence of the \
                                                 genes in a custom database")
    parser.add_argument("genelist", help="file containing gene accessions. if delimited, use \
                                            the headers in the example file as a template")
    parser.add_argument("blastdb", help="blastdb of interest")
    parser.add_argument("-o", "--output", help="directory in which to place the output files")
    parser.add_argument("-s", "--score_min", help="not currently used; will be used to determinine \
                                                   an optional scoring threshold")
    parser.add_argument("-t", "--blast_type", help="blastn or tblastx")
    parser.add_argument("-c", "--cores", help="number of cores to use (at your own risk)",
                        type=int)
    parser.add_argument("-b", "--query_bin", help="number of sequences to put per blast job",
                        type=int)

    args = parser.parse_args()

    genelist = args.genelist
    blastdb = args.blastdb
    blasttype = args.blast_type
    output = args.output
    score_min = args.score_min
    cores = args.cores
    query_bin = args.query_bin
date = str(datetime.datetime.now().strftime('%Y%m%d'))
if cores > cpu_count():
    raise ValueError("number of cores specified exceeds available cores!")
if not os.path.isdir(output):
    print("creating %s" % output)
    os.mkdir(output)
#%%  open accessions file, determine type, and parse


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


def run_entrez(accessions):
    """ Grab sequences from NCBI, write out resulting fasta file
    """
    print("\n\nFetching %i accessions from NCBI" % len(accessions))
    Entrez.email = "alfredTheDaring@gmail.com"
    sequence_handle = Entrez.efetch(db="nucleotide", id=accessions, rettype="fasta")
    seqs = SeqIO.parse(sequence_handle, "fasta")
    print("writing results to %s" % str(os.path.join(output, date) + "_sequences.fa"))
    with open(str(os.path.join(output, date)+"_sequences.fa"), "w") as fasta_output:
        SeqIO.write(seqs, fasta_output, "fasta")
    return(1)


def parse_multifasta(accessions):
    sequences_fasta = open(str(os.path.join(output, date)+"_sequences.fa"), "r")
    entrez_results = list(SeqIO.parse(sequences_fasta, "fasta"))
    nseqs = len(entrez_results)
    print("returned %i accessions from NCBI" % nseqs)
    if(len(accessions) != len(entrez_results)):
        print("Warning! not all accessions were found!")
    sequences_fasta.close()
    return(nseqs)


def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = iterator.next()
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch


def split_fasta(nseqs, cores):
    """ build set: bin by query_bin, devide up into number of cores
    """
    all_seqs = SeqIO.parse(open(str(os.path.join(output, date)+"_sequences.fa"), "r"), "fasta")
    set_size = int(round(nseqs/cores)+1)  # add 1 to avoid round-down errors
    nfiles = cores  # at least 1 file per core
    while set_size > query_bin:
        nfiles = nfiles * 2  # reduce set by half, double the number of files
        # if too large, just cut in half.  best method? probably not
        set_size = int(round(set_size / 2))
        temp_output_prefix = os.path.join(output, str(date + "_snagnblast_temp_dir"), "")
    if not os.path.isdir(temp_output_prefix):
        os.mkdir(temp_output_prefix)
    for i, batch in enumerate(batch_iterator(all_seqs, set_size)):
        filename = os.path.join(temp_output_prefix, "sequences_%i.fna" % (i + 1))
        handle = open(filename, "w")
        count = SeqIO.write(batch, handle, "fasta")
        handle.close()
        print("Wrote %i records to %s" % (count, filename))
    return(os.path.join(os.path.split(filename)[0], ""))


def run_blastn(filename):
    # build commandline call
    output_path_tab = str(os.path.join(output, date) + "_dcmegablast_results_" +
                          str(md5.new(str(datetime.datetime.now())).hexdigest()) + ".tab")
    blast_cline = NcbiblastnCommandline(query=filename,
                                        db=blastdb,  evalue=10,
                                        outfmt=6, out=output_path_tab)
    add_params = " -num_threads 1 -max_target_seqs 2000 -task dc-megablast"
    blast_command = str(str(blast_cline) + add_params)
    print("Running blastn search for %s..." % filename)
    subprocess.call(blast_command, shell=True)
    return(output_path_tab)


def run_tblastx(filename):
    # build commandline call
    output_path_tab = str(os.path.join(output, date) + "_tblastx_results_" +
                          str(md5.new(str(datetime.datetime.now())).hexdigest()) + ".tab")
    blast_cline = NcbitblastxCommandline(query=filename,
                                         db=blastdb, evalue=10,
                                         outfmt=6, out=output_path_tab)
    add_params = " -num_threads 1 -max_target_seqs 2000 -query_gencode 11 -db_gencode 11"
    blast_command = str(str(blast_cline)+add_params)
    print("Running tblastx search for %s ..." % filename)
    subprocess.call(blast_command, shell=True)
    return(output_path_tab)


def multiblast_worker(blast_method, work_queue, done_queue):
    ''' this is going to be how blast is execulted for multiprocessing, stolen laregely
        from http://toastdriven.com/blog/2008/nov/11/brief-introduction -multiprocessing/
    '''
    try:
        if blast_method == "blastn":
            for fasta in iter(work_queue.get, 'STOP'):
                outfile = run_blastn(fasta)
                done_queue.put(outfile)
        elif blast_method == "tblastx":
            for fasta in iter(work_queue.get, 'STOP'):
                outfile = run_tblastx(fasta)
                done_queue.put(outfile)
        else:
            raise ValueError("you need to use either blastn or tblastx, sorry!")
    except Exception, e:
        done_queue.put("%s failed on %s with: %s" % (current_process().name, file, e.message))
    return True


def get_complete_paths_of_files(directory):
    filenames = []
    shortnames = [i for i in os.listdir(directory) if not os.path.isdir(os.path.join(directory, i))]
    for i in shortnames:
        filenames.append(os.path.join(directory, i))
    return(filenames)


def artstop():
    """artificial stoppiing function
    """
    if 27 == 27:
        raise ValueError("thats all for now folks!")


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
    colnames = ["query_id", "subject_id", "identity_perc", "alignment_length", "mismatches",
                "gap_opens", "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score"]
    csv_results = pd.read_csv(open(infile), comment="#", sep="\t", names=colnames)
    #This regex will probably break things rather badly before too long...
    # it looks for capital letter and numbers, dot, number, ie SHH11555JJ8.99
    csv_results["accession"] = csv_results.query_id.str.extract(accession_pattern)
    # write out results with new headers or with new headers and merged metadat from accessions.tab
    genes = open(genelist, "r")
    genedf = pd.read_csv(genes, sep=",")
    output_path_csv = str(os.path.splitext(infile)[0]+".csv")
    results_annotated = pd.merge(csv_results, genedf, how="left",  on="accession")
    results_annotated.to_csv(open(output_path_csv, "w"))
    print("wrote final csv to %s" % output_path_csv)
#%%


def main():
    genes = get_accessions(genelist)  # parse the input file
    entrez_return_code = run_entrez(genes)  # get data from entrez , write out
    print(entrez_return_code)
    no_seqs = parse_multifasta(genes)  # counts sequences in fasta
    # split multifasta into smaller files according to nmber of cores, etc
    temp_directory = split_fasta(no_seqs, cores)
    files = get_complete_paths_of_files(temp_directory)
    workers = cores
    work_queue = Queue()
    done_queue = Queue()
    processes = []

    print("queued files:")
    for fasta in files:   # put the split up fastas in  the queue
        print("\t %s" % fasta)
        work_queue.put(fasta)

    for w in xrange(workers):  # for each core
        # do some blast magic
        p = Process(target=multiblast_worker, args=(blasttype, work_queue, done_queue))
        p.start()  # run the process
        processes.append(p)  # log completed process
        work_queue.put('STOP')  # at the end of working queue, put stop as the last entry

    for p in processes:  # for all processes, join when finished
        p.join()

    done_queue.put('STOP')  # at the end of done queue, put stop as the last entry

    for fasta in iter(done_queue.get, 'STOP'):  # keep the user posted on whats crackin (craic-in?)
        print ("%s completed" % fasta)
    output_dir_files = get_complete_paths_of_files(output)
    merged_fasta = merge_outfiles(output_dir_files,
                                  os.path.join(output, str(date+"_results_merged.csv")))
    cleanup_output_to_csv(merged_fasta)
    print("<insert summary here>. have a great day!")


if __name__ == '__main__':
    main()
